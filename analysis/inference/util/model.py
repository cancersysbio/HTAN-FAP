import os
os.environ['OPENBLAS_NUM_THREADS'] = '5'

# Use pyabc env
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import sys
import warnings
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

warnings.filterwarnings("ignore")

def calculate_densities(arr, x_range = np.linspace(0,1,21)):
	# Define the range
	x_range = np.linspace(0, 1, 50)

	# Calculate the KDE and density
	kde = gaussian_kde(arr)
	density = kde(x_range)

	return density

def calculate_ecdf(arr, points = np.linspace(0,1,21)):
	ecdf = ECDF(arr)
	
	# Eval ECDF at points
	ecdf_values = ecdf(points)
	
	return ecdf_values

def get_sample_summary_data(muts):
	# Take only autosomes
	# Return num total mutations, num subclonal mutations, num clonal mutations, array of mutation frequencies as vaf and ccfs
	muts = muts[(muts['group'] == 'Polyp') & (muts['mut_type'] == 'Somatic')].reset_index(drop=True)
	muts['purity_clonal'] = muts['clonality'].to_list()
	c = muts.groupby('clonality')['purity_clonal'].size().reset_index(drop=False)
	
	num_clonal_subclonal_total = np.zeros(3)
	clonal = c[c['clonality']=='CLONAL']['purity_clonal'].to_list()
	subclonal = c[c['clonality']=='SUBCLONAL']['purity_clonal'].to_list()
	if len(clonal)>0:
		num_clonal_subclonal_total[0] = clonal[0]
	if len(subclonal)>0:
		num_clonal_subclonal_total[1] = subclonal[0]

	num_clonal_subclonal_total[2] = muts.shape[0]

	ccfs = muts['ppVAF'].to_numpy()
	vafs = muts['vaf'].to_numpy()

	densities = calculate_densities(ccfs)

	probs = calculate_ecdf(ccfs)

	return(num_clonal_subclonal_total, ccfs, vafs, probs, densities)

class model:

	def __init__(self, parameters, julia='/home/rschenck/.juliaup/bin/julia', model='/home/rschenck/oak_dir/cluster/model/run_model.jl'):
		if os.path.exists(julia)==False:
			mpath = os.path.dirname(os.path.abspath("__file__"))
			julia_path = mpath.split('cluster/')[0]
			julia = os.path.join(julia_path, "../tools/julia-1.9.4/bin/julia")

		if os.path.exists(model)==False:
			mpath = os.path.dirname(os.path.abspath("__file__"))
			model = mpath.split('cluster/')[0] + "cluster/model/run_model.jl"

		self.parameters = parameters
		self.julia = julia
		self.model = model

	def run_model(self, params):
		# init_s, n, init_t, pol_s, pol_mu, pol_b = params['s_coef'], params['num_seeds'], params['polyp_init_time'], params['s_coef_polyp'], params['mut_rate_polyp'], params['polyp_birth_rate']

		# Set these parameters from inference run
		# params is provided from pyabc
		parameters_this_kernel = self.parameters.copy()
		# uniq_fn = f"{int(np.random.randint(0,10000))}_{int(params['num_seeds'])}_{int(params['polyp_init_time'])}_{int(params['mut_rate_polyp'])}_{round(params['polyp_birth_rate'],3)}"
		uniq_fn = ""
		for key, value in params.items():
			parameters_this_kernel[key] = value
			uniq_fn += str(np.round(value, 4)).replace('.','')
		uniq_fn += str(np.random.randint(0,10000))

		# Set max_t based on the init time of the kernel
		# tick is now in 6 month intervals!!!!!!!
		parameters_this_kernel['max_t'] = (int(parameters_this_kernel['max_t'])*2)
		
		# Set outputdir for this SINGLE model run
		thisoutfile = f"{'/'.join(parameters_this_kernel['outfile'].rstrip('/').split('/'))}/infer_tmp_{uniq_fn}"
		parameters_this_kernel['outfile'] = thisoutfile
		parameters_this_kernel['randomSeed'] = np.random.randint(43, 1000000) # Anything but 42
		# parameters_this_kernel['verbose'] = 1
		# print(parameters_this_kernel['outfile'])

		# Create directory if it doesn't exist
		if os.path.exists(parameters_this_kernel['outfile']) == False:
			os.mkdir(parameters_this_kernel['outfile'])

		modellog = f"{parameters_this_kernel['outfile']}/run.log"
		
		try:
			# To run the model
			# print(parameters_this_kernel)
			exitcode = run_model_job(parameters_this_kernel, julia=self.julia, model=self.model, model_log=modellog)
			
			# Process output data
			# To load the model results
			h5 = parameters_this_kernel['outfile'] + '/out.hdf5'
			muts = parameters_this_kernel['outfile'] + '/muts.tsv'

			out = read_h5(h5)
			# print(h5_tree(out))
			muts = pd.read_csv(muts, sep='\t')

			# Process data
			# Label "Germline" and "Somatic"
			# Convert time for inutero to months
			# Convert VAF to ppVAF
			popdf, muts = process_data(out, muts, parameters_this_kernel['outfile'] + "/", min_vaf=0.01)
			num_clonal_subclonal_total, ccfs, vafs, probs, densities = get_sample_summary_data(muts)

			# Put sum stats in a dictionary
			these_sum_stats = {'num_clonal_subclonal_total':num_clonal_subclonal_total, 'probs':probs, 'density':densities}

			# Remove the output file
			if os.path.exists(parameters_this_kernel['outfile']):
				shutil.rmtree(parameters_this_kernel['outfile'])
				
		except Exception as e:
			print(e)
			print(parameters_this_kernel)
			# PLACEHOLDER!!!!! WHILE DETERMINING ERROR WHERE OUTPUTS DO NOT GET MADE (SUSPECT an out of memory error)
			these_sum_stats = {'num_clonal_subclonal_total':np.asarray([10000000000]), 'probs':np.asarray([10000000000]), 'density':np.asarray([10000000000])}
			
			# Remove the output file
			if os.path.exists(parameters_this_kernel['outfile']):
				shutil.rmtree(parameters_this_kernel['outfile'])

		return these_sum_stats