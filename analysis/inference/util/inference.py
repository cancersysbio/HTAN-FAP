# Use pyabc env
import os
os.environ['OPENBLAS_NUM_THREADS'] = '5'
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import sys
import warnings
import pyabc
from scipy.stats import gaussian_kde
from statsmodels.distributions.empirical_distribution import ECDF
from sklearn.metrics import mean_squared_error
import pickle
import shutil
from scipy.special import kl_div

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# MUST LOAD AFTER
import h5py

if os.path.exists('/home/rschenck/oak_dir/cluster/outputs/inferences/inferences_set4/util'):
	sys.path.append('/home/rschenck/oak_dir/cluster/outputs/inferences/inferences_set4/util')
else:
	sys.path.append('/home/users/rschenck/oak/cluster/outputs/inferences/inferences_set4/util')
from util.IO import read_yaml_to_dict, get_parameter_output, run_model_job, get_parameter_output_inference, h5_tree, read_h5, process_data
from util.model import model

class inference_Model:

	def __init__(self, name = "test", summary_stats={}, max_t=100):
		self.name = name
		
		self.summary_stats = summary_stats
		self.max_t = max_t

		self.prior = None
		self.transition = None

		self.outfile = None

		self.history = None

		self.model_used = None
	
	def count_distance(self, x, y):
		if x['num_clonal_subclonal_total'][0] != -1:
			d0 = x['num_clonal_subclonal_total']
			d = y['num_clonal_subclonal_total']
			
			mse = mean_squared_error(d0, d)
			rmse = np.sqrt(mse)
		else:
			rmse = 100000000

		return rmse

	def probs_distance(self, x, y):
		if x['probs'][0] != -1:
			d0 = x['probs']
			d = y['probs']

			# Sum over the absolute difference at each idx
			# This is the L1 norm
			# l1 = np.sum(np.abs(d0 - d))

			mse = mean_squared_error(d0, d)
			rmse = np.sqrt(mse)
		else:
			rmse = 100000000
		return rmse
	
	def density_distance(self, x, y):
		if x['density'][0] != -1:
			d0 = x['density']
			d = y['density']

			# Sum over the absolute difference at each idx
			# This is the L1 norm
			# l1 = np.sum(np.abs(d0 - d))
			mse = mean_squared_error(d0, d)
			rmse = np.sqrt(mse)
		else:
			rmse = 100000000
			
		return rmse

	def kldiv_distance(self, x, y):
		d = sum(kl_div(x['density'], y['density']))
		return(d)
	
	def run_inference(self, samnam, npoints=100, max_nr_populations=3, epsilon_diff=1., parent_output="./", outputdb="inference.db", cores=100, continue_abc=False):
		#### Setup priors
		discrete_domain_founders = np.arange(1, 20) # Number of seeds domain
		discrete_domain_polyp_mut = np.arange(10,80)
		discrete_init_time = np.arange(1, (self.max_t-1) * 2)

		self.prior = pyabc.Distribution(
			s_coef=pyabc.random_variables.LowerBoundDecorator( pyabc.RV("halfnorm", 0., 0.25), 0. ),
			# num_seeds=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('poisson', 4 ), 1 ),
			num_seeds=pyabc.RV('rv_discrete', values=(discrete_domain_founders, [1 / len(discrete_domain_founders)] * len(discrete_domain_founders))),
			# polyp_init_time=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('poisson', 2 ), 1),
			polyp_init_time=pyabc.RV('rv_discrete', values=(discrete_init_time, [1 / len(discrete_init_time)] * len(discrete_init_time))),
			s_coef_polyp=pyabc.random_variables.LowerBoundDecorator( pyabc.RV("halfnorm", 0.,  0.25), 0. ),
			mut_rate_polyp=pyabc.RV('rv_discrete', values=(discrete_domain_polyp_mut, [1 / len(discrete_domain_polyp_mut)] * len(discrete_domain_polyp_mut))),
			# mut_rate_polyp=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('poisson', 36 ), 1),
			polyp_birth_rate=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('norm', 0.034*3. , 0.01), 0.021) # Lower bound on this must be > death
		)

		#### transition kernels
		self.transition = pyabc.AggregatedTransition(
			mapping={
				's_coef': pyabc.MultivariateNormalTransition(),
				# 'num_seeds': pyabc.MultivariateNormalTransition(),
				'num_seeds': pyabc.DiscreteJumpTransition(
				    domain=discrete_domain_founders, p_stay=0.8),
				# 'polyp_init_time': pyabc.MultivariateNormalTransition(),
				'polyp_init_time': pyabc.DiscreteJumpTransition(
				    domain=discrete_init_time, p_stay=0.8),
				's_coef_polyp': pyabc.MultivariateNormalTransition(),
				# 'mut_rate_polyp': pyabc.MultivariateNormalTransition(),
				'mut_rate_polyp': pyabc.DiscreteJumpTransition(domain=discrete_domain_polyp_mut, p_stay=0.8),
				'polyp_birth_rate': pyabc.MultivariateNormalTransition()
			}
		)

		# Gets dictionary of parameters from template
		template = read_yaml_to_dict("/oak/stanford/groups/ccurtis2/users/rschenck/cluster/util/conf.yaml")

		# # Setup the model with the parameters that are constant
		# To handle in class: outfile and seed
		const = {
			'outfile' : parent_output,
			'initSize' : 100000,  # Initial size
			'birth_rate' : float( 0.034 ), # Normal fission rate
			'death_rate' : float( 0.01 ), # Initial death rate
			'mut_rate' : int( 36  ), # Initial mutation rate
			'adv_mut_rate' : float( 0.05  ),
			'adv_mut_rate_polyp' : float( 0.05  ),
			'polyp_death_rate' : 0.02,
			'max_t' : self.max_t, # Age of patient
			'final_pop_size' : 1000000, # This is maximum of the population
			'runID': samnam,
			'verbose': 1
		}

		# Update with constants
		for key, value in const.items():
			template[key] = value

		# Init model class
		m = model(template)
		self.model_used = m

		# dist = pyabc.AggregatedDistance([self.count_distance, self.probs_distance, self.density_distance], weights= np.ones(3)/3 )
		dist = self.kldiv_distance

		# Run ABC
		multi_sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=cores)
		abc = pyabc.ABCSMC(m.run_model,
					self.prior,
					dist,
					transitions=self.transition,
					population_size=npoints,
					sampler = multi_sampler
					)

		if continue_abc==False:
			# Setup the database
			abc.new(f"sqlite:///{outputdb}", self.summary_stats)
			self.outfile = outputdb

			self.history = abc.run(max_nr_populations=max_nr_populations, min_eps_diff=epsilon_diff)
		else:
			abc.load(f"sqlite:///{outputdb}", 1)

			self.outfile = outputdb

			self.history = abc.run(max_nr_populations=max_nr_populations, min_eps_diff=epsilon_diff)

		return self.history, self.prior, self.transition