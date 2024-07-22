import yaml
import os
os.environ['OPENBLAS_NUM_THREADS'] = '5'
import subprocess as sp
import pandas as pd
import numpy as np

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# MUST LOAD AFTER
import h5py

def get_parameter_output(template_path, outpath='./model_outputs/logistic_growth/eval_pop_size/', runID='output_1', **kwargs):
    """
    Reads a YAML template, updates parameters with additional keyword arguments,
    and saves the updated parameters as a new YAML file.
    
    Args:
        template_path (str): Path to the template YAML file.
        outpath (str): Output path for the parameters file.
        runID (str): Run ID to name the output directory.
        makedir: Make directory and yaml on disk
        **kwargs: Additional parameters to update in the template.
        
    Returns:
        dict: Updated parameters dictionary.
    """
    # Load the template YAML into a dictionary
    parameters = read_yaml_to_dict(template_path)
    
    if parameters is None:
        print("Failed to load template YAML.")
        return None
    
    # Update parameters dictionary with any additional keyword arguments
    for key, value in kwargs.items():
        parameters[key] = value
    
    # Build the output file path
    parameters['outfile'] = os.path.join(outpath, runID, '')
    
    # Ensure the directory exists
    if not os.path.exists(parameters['outfile']):
        os.makedirs(parameters['outfile'], exist_ok=False)
    
    # Save updated parameters as a YAML file
    with open(os.path.join(parameters['outfile'], 'params.yaml'), 'w') as yaml_file:
        yaml.dump(parameters, yaml_file, default_flow_style=False)
    
    return parameters

def get_parameter_output_inference(parameters_loaded, outpath='./model_outputs/logistic_growth/eval_pop_size/', runID='output_1', **kwargs):
    # Load the template YAML into a dictionary
    parameters = parameters_loaded
    
    if parameters is None:
        print("Failed to load template YAML.")
        return None
    
    # Update parameters dictionary with any additional keyword arguments
    for key, value in kwargs.items():
        parameters[key] = value
    
    # Build the output file path
    parameters['outfile'] = os.path.join(outpath, runID, '')
    return parameters

def h5_tree(val, pre=''):
    """
    Pretty print HDF5 file structure.
    """
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                print(pre + '└── ' + key)
                h5_tree(val, pre+'    ')
            else:
                try:
                    print(pre + '└── ' + key + ' (%d)' % len(val))
                except TypeError:
                    print(pre + '└── ' + key + ' (scalar)')
        else:
            if type(val) == h5py._hl.group.Group:
                print(pre + '├── ' + key)
                h5_tree(val, pre+'│   ')
            else:
                try:
                    print(pre + '├── ' + key + ' (%d)' % len(val))
                except TypeError:
                    print(pre + '├── ' + key + ' (scalar)')

def read_h5(filename):
    f = h5py.File(filename, 'r')
    return f

def read_yaml_to_dict(filepath):
    """
    Reads a YAML file and returns its contents as a dictionary.
    
    Args:
        filepath (str): Path to the YAML file.
        
    Returns:
        dict: Contents of the YAML file as a dictionary.
    """
    try:
        with open(filepath, 'r') as file:
            config_dict = yaml.safe_load(file)
        return config_dict
    except Exception as e:
        print(f"Error reading the YAML file: {e}")
        return None
    
def run_model_job(parameters, julia='julia', model='run_model.jl', model_log='tmp.log'):
    """
    Runs the Julia model with the specified parameters.
    
    Args:
        parameters (dict): Dictionary of parameters to pass to the Julia model.
    """
    args = [
        julia,
        model,
        '--runtype', str(parameters['runtype']),
        '--initSize', str(parameters['initSize']),
        '--birth_rate', str(parameters['birth_rate']),
        '--death_rate', str(parameters['death_rate']),
        '--mut_rate', str(parameters['mut_rate']),
        '--adv_mut_rate', str(parameters['adv_mut_rate']),
        '--s_coef', str(parameters['s_coef']),
        '--num_seeds', str(parameters['num_seeds']),
        '--final_pop_size', str(parameters['final_pop_size']),
        '--s_coef_polyp', str(parameters['s_coef_polyp']),
        '--adv_mut_rate_polyp', str(parameters['adv_mut_rate_polyp']),
        '--polyp_birth_rate', str(parameters['polyp_birth_rate']),
        '--polyp_death_rate', str(parameters['polyp_death_rate']),
        '--polyp_init_time', str(parameters['polyp_init_time']),
        '--mut_rate_polyp', str(parameters['mut_rate_polyp']),
        '--outfile', str(parameters['outfile']),
        '--seed', str(parameters['seed']),
        '--max_t', str(parameters['max_t']),
        '--verbose', str(parameters['verbose']),
        '--mean_depth', str(parameters['mean_depth']),
        '--sd_depth', str(parameters['sd_depth']),
        '--sample_size', str(parameters['sample_size'])
    ]

    cmd = [str(x) for x in args]
    # Call the Julia script with the specified parameters
    # sp.call([str(x) for x in args])
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    out, err = proc.communicate()
    
    with open(model_log, 'w') as outlog:
        outlog.write("STDOUT:\n")
        outlog.write(out.decode())
        outlog.write("\n\nSTDERR:\n")
        outlog.write(err.decode())

    return(proc.returncode)

def process_data(out, muts, outdir, min_vaf=0.01, germline_cutoff=-1.*1/12., save_files = True):
    init = out['pops/inutero'][0:]
    fission = out['pops/fission'][0:]
    polyp = out['pops/polypfission'][0:]

    init = pd.DataFrame({'Pop': init, 't': [i/12. for i in range(-len(init), 0)],
                        'Group':['Dev' for i in range(0, len(init))]})
    init['t'] = init['t']+1/12. # Sets final timepoint at 0
    fission = pd.DataFrame({'Pop': fission, 't': [i/2. for i in range(0, len(fission))],
                            'Group':['Normal' for i in range(0, len(fission))]})
    polyp = pd.DataFrame({'Pop': polyp, 't': [i/2. for i in range(0, len(polyp))],
                        'Group':['Polyp' for i in range(0, len(polyp))]})
    polyp['t'] = polyp['t'] + (fission['t'].max())

    df = pd.concat([init, fission, polyp], axis=0).reset_index(drop=True)

    all_drivers = np.concatenate([out['driver_muts/polypfission'][0:], out['driver_muts/fission'][0:], out['driver_muts/inutero'][0:]])
    muts['is_driver'] = False
    # Relabel for driver mutations
    muts.loc[muts['mutation'].isin(all_drivers), 'is_driver'] = True
    # Adjust time for inutero
    # Monthly and 6 months adjustment
    muts['t'] = muts.apply(lambda x: (x['t']-10.)/12. if x['group']=='inutero' else x['t']/2., axis=1)

    # Filter on 'sequenced' vaf
    muts = muts[muts['vaf']>=min_vaf].reset_index(drop=True)
    muts['ppVAF'] = muts['vaf'] * 2.

    muts['clonality'] = muts['ppVAF'].apply(lambda x: 'SUBCLONAL' if x < 0.85 else 'CLONAL')

    # Label germline and somatic
    # There are a few different ways of doing this
    # VAF = 0.5, clonal in inutero, or time dependency
    # This approach uses a 1 month cutoff and removes mutations with 0 ID
    germline = muts[(muts['t'] < germline_cutoff) & (muts['group']=='inutero') & (muts['is_driver']==False)]['mutation'].to_list()
    # germline = muts[(muts['clonality'] == 'CLONAL') & (muts['group']=='inutero') & (muts['is_driver']==False)]['mutation'].to_list()
    # germline = muts[(muts['raw_vaf'] == 0.5) & (muts['group']=='inutero') & (muts['is_driver']==False)]['mutation'].to_list()

    muts['mut_type'] = muts.apply(lambda x: 'Germline' if x['mutation'] in germline else 'Somatic', axis=1).reset_index(drop=True)

    # Remove 0 mutation IDs
    muts = muts[muts['mutation'] != 0].reset_index(drop=True)
    
    if save_files:
        df.to_csv(outdir + 'post_pops.tsv', sep='\t', index=False)
        muts.to_csv(outdir + 'post_muts.tsv', sep='\t', index=False)

    return df, muts