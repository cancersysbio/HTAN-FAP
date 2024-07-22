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
import argparse
from scipy.stats import gaussian_kde
from statsmodels.distributions.empirical_distribution import ECDF
from sklearn.metrics import mean_squared_error
import pickle
import shutil

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# MUST LOAD AFTER
import h5py

if os.path.exists('/home/rschenck/oak_dir/cluster/outputs/inferences/inferences_set4/util'):
	sys.path.append('/home/rschenck/oak_dir/cluster/outputs/inferences/inferences_set4/util')
else:
	sys.path.append('/home/users/rschenck/oak/cluster/outputs/inferences/inferences_set4/util')
from util.IO import read_yaml_to_dict, get_parameter_output, run_model_job, get_parameter_output_inference, h5_tree, read_h5, process_data
from util.inference import inference_Model

warnings.filterwarnings("ignore")

# def calculate_densities(arr, x_range = np.linspace(0,1,21)):
# 	# Define the range
# 	x_range = np.linspace(0, 1, 50)

# 	# Calculate the KDE and density
# 	kde = gaussian_kde(arr)
# 	density = kde(x_range)

# 	return density

# def calculate_ecdf(arr, points = np.linspace(0,1,21)):
# 	ecdf = ECDF(arr)
	
# 	# Eval ECDF at points
# 	ecdf_values = ecdf(points)
	
# 	return ecdf_values

# def get_sample_summary_data(muts):
# 	# Take only autosomes
# 	# Return num total mutations, num subclonal mutations, num clonal mutations, array of mutation frequencies as vaf and ccfs
# 	muts = muts[(muts['group'] == 'Polyp') & (muts['mut_type'] == 'Somatic')].reset_index(drop=True)
# 	muts['purity_clonal'] = muts['clonality'].to_list()
# 	c = muts.groupby('clonality')['purity_clonal'].size().reset_index(drop=False)
	
# 	num_clonal_subclonal_total = np.zeros(3)
# 	clonal = c[c['clonality']=='CLONAL']['purity_clonal'].to_list()
# 	subclonal = c[c['clonality']=='SUBCLONAL']['purity_clonal'].to_list()
# 	if len(clonal)>0:
# 		num_clonal_subclonal_total[0] = clonal[0]
# 	if len(subclonal)>0:
# 		num_clonal_subclonal_total[1] = subclonal[0]

# 	num_clonal_subclonal_total[2] = muts.shape[0]

# 	ccfs = muts['ppVAF'].to_numpy()
# 	vafs = muts['vaf'].to_numpy()

# 	densities = calculate_densities(ccfs)

# 	probs = calculate_ecdf(ccfs)

# 	return(num_clonal_subclonal_total, ccfs, vafs, probs, densities)

# class model:

# 	def __init__(self, parameters, julia='/home/rschenck/.juliaup/bin/julia', model='/home/rschenck/oak_dir/cluster/model/run_model.jl'):
# 		self.parameters = parameters
# 		self.julia = julia
# 		self.model = model

# 	def run_model(self, params):
# 		# init_s, n, init_t, pol_s, pol_mu, pol_b = params['s_coef'], params['num_seeds'], params['polyp_init_time'], params['s_coef_polyp'], params['mut_rate_polyp'], params['polyp_birth_rate']

# 		# Set these parameters from inference run
# 		# params is provided from pyabc
# 		parameters_this_kernel = self.parameters.copy()
# 		# uniq_fn = f"{int(np.random.randint(0,10000))}_{int(params['num_seeds'])}_{int(params['polyp_init_time'])}_{int(params['mut_rate_polyp'])}_{round(params['polyp_birth_rate'],3)}"
# 		uniq_fn = ""
# 		for key, value in params.items():
# 			parameters_this_kernel[key] = value
# 			uniq_fn += str(np.round(value, 4)).replace('.','')
# 		uniq_fn += str(np.random.randint(0,10000))

# 		# Set max_t based on the init time of the kernel
# 		# tick is now in 6 month intervals!!!!!!!
# 		parameters_this_kernel['max_t'] = (int(parameters_this_kernel['max_t'])*2)
		
# 		# Set outputdir for this SINGLE model run
# 		thisoutfile = f"{'/'.join(parameters_this_kernel['outfile'].rstrip('/').split('/'))}/infer_tmp_{uniq_fn}"
# 		parameters_this_kernel['outfile'] = thisoutfile
# 		parameters_this_kernel['randomSeed'] = np.random.randint(43, 1000000) # Anything but 42
# 		# parameters_this_kernel['verbose'] = 1
# 		print(parameters_this_kernel['outfile'])

# 		# Create directory if it doesn't exist
# 		if os.path.exists(parameters_this_kernel['outfile']) == False:
# 			os.mkdir(parameters_this_kernel['outfile'])

# 		modellog = f"{parameters_this_kernel['outfile']}/run.log"
		
# 		# To run the model
# 		exitcode = run_model_job(parameters_this_kernel, julia=self.julia, model=self.model, model_log=modellog)
		
# 		# Process output data
# 		# To load the model results
# 		h5 = parameters_this_kernel['outfile'] + '/out.hdf5'
# 		muts = parameters_this_kernel['outfile'] + '/muts.tsv'

# 		out = read_h5(h5)
# 		# print(h5_tree(out))
# 		muts = pd.read_csv(muts, sep='\t')

# 		# Process data
# 		# Label "Germline" and "Somatic"
# 		# Convert time for inutero to months
# 		# Convert VAF to ppVAF
# 		popdf, muts = process_data(out, muts, parameters_this_kernel['outfile'] + "/", min_vaf=0.01)
# 		num_clonal_subclonal_total, ccfs, vafs, probs, densities = get_sample_summary_data(muts)

# 		# Put sum stats in a dictionary
# 		these_sum_stats = {'num_clonal_subclonal_total':num_clonal_subclonal_total, 'probs':probs, 'density':densities}

# 		# Remove the output file
# 		if os.path.exists(parameters_this_kernel['outfile']):
# 			shutil.rmtree(parameters_this_kernel['outfile'])

# 		return these_sum_stats

# class inference_Model:

# 	def __init__(self, name = "test", summary_stats={}, max_t=100):
# 		self.name = name
		
# 		self.summary_stats = summary_stats
# 		self.max_t = max_t

# 		self.prior = None
# 		self.transition = None

# 		self.outfile = None

# 		self.history = None

# 		self.model_used = None
	
# 	def count_distance(self, x, y):
# 		d0 = x['num_clonal_subclonal_total']
# 		d = y['num_clonal_subclonal_total']
		
# 		mse = mean_squared_error(d0, d)
# 		rmse = np.sqrt(mse)

# 		return rmse

# 	def probs_distance(self, x, y):
# 		d0 = x['probs']
# 		d = y['probs']

# 		# Sum over the absolute difference at each idx
# 		# This is the L1 norm
# 		# l1 = np.sum(np.abs(d0 - d))

# 		mse = mean_squared_error(d0, d)
# 		rmse = np.sqrt(mse)
# 		return rmse
	
# 	def density_distance(self, x, y):
# 		d0 = x['density']
# 		d = y['density']

# 		# Sum over the absolute difference at each idx
# 		# This is the L1 norm
# 		# l1 = np.sum(np.abs(d0 - d))
# 		mse = mean_squared_error(d0, d)
# 		rmse = np.sqrt(mse)
# 		return rmse
	
# 	def run_inference(self, samnam, npoints=100, max_nr_populations=3, epsilon_diff=0.01, parent_output="./", outputdb="inference.db"):
# 		#### Setup priors
# 		discrete_domain_founders = np.arange(1, 200) # Number of seeds domain
# 		discrete_domain_polyp_mut = np.arange(10,80)
# 		discrete_init_time = np.arange(1, (self.max_t-1) * 2)

# 		self.prior = pyabc.Distribution(
# 			s_coef=pyabc.random_variables.LowerBoundDecorator( pyabc.RV("halfnorm", 0.15), 0. ),
# 			# num_seeds=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('poisson', 4 ), 1 ),
# 			num_seeds=pyabc.RV('rv_discrete', values=(discrete_domain_founders, [1 / len(discrete_domain_founders)] * len(discrete_domain_founders))),
# 			# polyp_init_time=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('poisson', 2 ), 1),
# 			polyp_init_time=pyabc.RV('rv_discrete', values=(discrete_init_time, [1 / len(discrete_init_time)] * len(discrete_init_time))),
# 			s_coef_polyp=pyabc.random_variables.LowerBoundDecorator( pyabc.RV("halfnorm", 0.15), 0. ),
# 			mut_rate_polyp=pyabc.RV('rv_discrete', values=(discrete_domain_polyp_mut, [1 / len(discrete_domain_polyp_mut)] * len(discrete_domain_polyp_mut))),
# 			# mut_rate_polyp=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('poisson', 36 ), 1),
# 			polyp_birth_rate=pyabc.random_variables.LowerBoundDecorator( pyabc.RV('norm', 0.034*10. , 0.01), 0.052) # Lower bound on this must be > death
# 		)

# 		#### transition kernels
# 		self.transition = pyabc.AggregatedTransition(
# 			mapping={
# 				's_coef': pyabc.MultivariateNormalTransition(),
# 				# 'num_seeds': pyabc.MultivariateNormalTransition(),
# 				'num_seeds': pyabc.DiscreteJumpTransition(
# 				    domain=discrete_domain_founders, p_stay=0.6),
# 				# 'polyp_init_time': pyabc.MultivariateNormalTransition(),
# 				'polyp_init_time': pyabc.DiscreteJumpTransition(
# 				    domain=discrete_domain_founders, p_stay=0.6),
# 				's_coef_polyp': pyabc.MultivariateNormalTransition(),
# 				# 'mut_rate_polyp': pyabc.MultivariateNormalTransition(),
# 				'mut_rate_polyp': pyabc.DiscreteJumpTransition(domain=discrete_init_time, p_stay=0.6),
# 				'polyp_birth_rate': pyabc.MultivariateNormalTransition()
# 			}
# 		)

# 		# Gets dictionary of parameters from template
# 		template = read_yaml_to_dict("/home/rschenck/oak_dir/cluster/util/conf.yaml")

# 		# # Setup the model with the parameters that are constant
# 		# To handle in class: outfile and seed
# 		const = {
# 			'outfile' : parent_output,
# 			'initSize' : 1000000,  # Initial size
# 			'birth_rate' : float( 0.034 ), # Normal fission rate
# 			'death_rate' : float( 0.01 ), # Initial death rate
# 			'mut_rate' : int( 36  ), # Initial mutation rate
# 			'adv_mut_rate' : float( 0.02  ),
# 			'adv_mut_rate_polyp' : float( 0.02  ),
# 			'polyp_death_rate' : 0.05,
# 			'max_t' : self.max_t*2., # Age of patient
# 			'final_pop_size' : 1000000, # This is maximum of the population
# 			'runID': samnam,
# 			'verbose': 1
# 		}

# 		# Update with constants
# 		for key, value in const.items():
# 			template[key] = value

# 		# Init model class
# 		m = model(template)
# 		self.model_used = m

# 		dist = pyabc.AggregatedDistance([self.count_distance, self.probs_distance, self.density_distance], weights= np.ones(3)/3 )

# 		# Run ABC
# 		multi_sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=100)
# 		abc = pyabc.ABCSMC(m.run_model,
# 					self.prior,
# 					dist,
# 					transitions=self.transition,
# 					population_size=npoints,
# 					sampler = multi_sampler
# 					)
		
# 		# Setup the database
# 		abc.new(f"sqlite:///{outputdb}", self.summary_stats)
# 		self.outfile = outputdb

# 		self.history = abc.run(max_nr_populations=max_nr_populations, min_eps_diff=epsilon_diff)

# 		return self.history, self.prior, self.transition

def main():
	parser = argparse.ArgumentParser(description='Run inference.')
	parser.add_argument('--npoints', default=10, dest='npoints', type=int,
						help='Number of simulations per iteration (default:250)')
	parser.add_argument('--epsilon_diff', default=1, dest='epsilon_diff', type=float,
						help='Difference in between runs (default:0.000000001)')
	parser.add_argument('--max_nr_pop', default=1, dest='maxiters', type=int,
						help='Maximum number of inference iterations (default:10)')
	parser.add_argument('--sample_name', default='A002C203', dest='sample_name', type=str,
						help='Maximum number of inference iterations (default:10)')
	parser.add_argument('--patient_id', default='A002', dest='patient', type=str,
						help='Output path for results')
	parser.add_argument('--wxs', default='wgs', dest='wxs', type=str,)
	parser.add_argument('--outputdir', default='/home/rschenck/oak_dir/cluster/outputs/inferences/inferences_set4/', dest='outputdir', type=str,
						help='Output path for results, for a group of inferences this should be the top level directory.')
	parser.add_argument('--data_sum_stats', default='/home/rschenck/oak_dir/cluster/data/summary_stats/wgs/A002C203.hdf5', dest='d_sum_stats', type=str, help='Path to data summary statistics')
	parser.add_argument('--continue_inf', default=False, type=bool, dest='continue_inf', help="Continue inference already started.")
	
	# Execute the parse_args() method
	args = parser.parse_args()

	# pat_var = args.pat_var
	npoints = int(args.npoints)
	epsilon_diff = float(args.epsilon_diff)
	maxiters = int(args.maxiters)
	sample_to_fit = args.sample_name
	wxs = args.wxs
	outputdb = args.outputdir + wxs + '/' + sample_to_fit + '/' + '%s.db'%(sample_to_fit)
	outputpkl = args.outputdir + wxs + '/' + sample_to_fit + '/' + '%s.pkl'%(sample_to_fit)
	patient = args.patient
	data_sum_stats = args.d_sum_stats
	continue_inf = args.continue_inf

	parent_output = args.outputdir + wxs + '/' + sample_to_fit + '/'

	if os.path.exists(parent_output.replace(sample_to_fit, '')) == False:
		os.mkdir(parent_output.replace(sample_to_fit, ''))
	if os.path.exists(parent_output) == False:
		os.mkdir(parent_output)

	# template_path = '/Users/ryanschenck/Dropbox/Projects/FAP_Evo/cluster/util/conf.yaml'
	# output_path = '/Users/ryanschenck/Dropbox/Projects/FAP_Evo/cluster/outputs/'
	# julia = '/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia'
	# model = '/Users/ryanschenck/Dropbox/Projects/FAP_Evo/cluster/model/run_model.jl'

	# Get summary stats
	f = read_h5(data_sum_stats)
	sam = os.path.basename(data_sum_stats).replace('.hdf5','')
	ppvaf = f[f'{sam}/ccfs'][0:]
	colectomy_age = f[f'{sam}/colectomy_age'][0:]
	density = f[f'{sam}/density'][0:]
	num_clonal_subclonal_total = f[f'{sam}/num_clonal_subclonal_total'][0:]
	probs = f[f'{sam}/probs'][0:]
	stage = f[f'{sam}/stage']
	vafs = f[f'{sam}/vafs'][0:]

	summary_stats = {
		'num_clonal_subclonal_total':num_clonal_subclonal_total,
		'probs':probs,
		'density':density
	}

	# Initialize the inference class
	inference = inference_Model(sam, summary_stats=summary_stats, max_t=colectomy_age)

	# Run the inference
	history, prior, transition = inference.run_inference(sam, npoints=npoints, max_nr_populations=maxiters, epsilon_diff=epsilon_diff, 
													  outputdb=outputdb, parent_output=parent_output, continue_abc=continue_inf, cores=20)

	runID = history.id
	print(runID)

	# Iterate over args and print to log
	outputfile = open(parent_output + '/log.runID_%s.txt'%(history.id), 'w')

	# outputfile.write( 'Stage: %s\n' % (stage) )

	# Close log file
	for arg in vars(args):
		print ( arg, getattr(args, arg) )
		outputfile.write( "%s: %s\n"%(arg, getattr(args, arg)) )
	outputfile.write( 'History ID: %s\n' % (history.id) )
	outputfile.close()

	pickle_out_data = {'Sample':sample_to_fit, 'Run_ID': history.id, 'prior_funcs':prior, 'transitions':transition}
	outputpkl = parent_output + '/' + '%s.runID_%s.pkl'%(sample_to_fit, history.id)
	pickle.dump(pickle_out_data, open(outputpkl, 'wb'))

if __name__=="__main__":
	main()