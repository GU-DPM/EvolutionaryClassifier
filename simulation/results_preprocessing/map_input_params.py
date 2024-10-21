import os
import pandas as pd
import math
import concurrent.futures

def get_param_df(filename):
  param_column_names = ['Parameter_ID',
                        'S_pop',
                        'R1_pop',
                        'R2_pop',
                        'R12_pop',
                        'g0',
                        'S_cell_sensitivity_D2',
                        'S_cell_sensitivity_D1',
                        'R1_cell_sensitivity_D2',
                        'R1_cell_sensitivity_D1',
                        'R2_cell_sensitivity_D2',
                        'R2_cell_sensitivity_D1',
                        'R12_cell_sensitivity_D2',
                        'R12_cell_sensitivity_D1',
                        'R12_transition_to_R12',
                        'R12_transition_to_R2',
                        'R12_transition_to_R1',
                        'R12_transition_to_S',
                        'R2_transition_to_R12',
                        'R2_transition_to_R2',
                        'R2_transition_to_R1',
                        'R2_transition_to_S',
                        'R1_transition_to_R12',
                        'R1_transition_to_R2',
                        'R1_transition_to_R1',
                        'R1_transition_to_S',
                        'S_transition_to_R12',
                        'S_transition_to_R2',
                        'S_transition_to_R1',
                        'S_transition_to_S',
                        ]
  param_df = pd.read_csv(filename,header=None)
  param_df.columns = param_column_names
  return param_df

def round_to_significant_digits(value, digits):
    if value == 0:
        return 0
    else:
        return round(value, digits - int(math.floor(math.log10(abs(value)))) - 1)

def gen_input_params(param_df):

  input_params_dict = {'Parameter_ID' : param_df['Parameter_ID'],
                       'initial_R1_percent' : [x/5e9 for x in param_df['R1_pop'] ],
                       'initial_R2_percent' : [x/5e9 for x in param_df['R2_pop'] ],
                       'S_sensitivity_D1_to_g0' : param_df['S_cell_sensitivity_D1']/param_df['g0'],
                       'S_sensitivity_D2_to_S_sensitivity_D1' : param_df['S_cell_sensitivity_D2']/param_df['S_cell_sensitivity_D1'],
                       'R1_sensitivity_D1_to_S_sensitivity_D1' : param_df['R1_cell_sensitivity_D1']/param_df['S_cell_sensitivity_D1'],
                       'R2_sensitivity_D2_to_S_sensitivity_D2' : param_df['R2_cell_sensitivity_D2']/param_df['S_cell_sensitivity_D2'],
                       }
  return pd.DataFrame(input_params_dict)

import math
def round_to_significant_digits(value, digits):
    if value == 0:
        return 0
    else:
        return round(value, digits - int(math.floor(math.log10(abs(value)))) - 1)

def check_input_params(input_param_df):
  # check if params are valid inputs
  valid_input_values = {'initial_R1_percent' : [0, 1e-09, 1e-07, 1e-05, 0.001, 0.1, 0.9],
                        'initial_R2_percent' : [0, 1e-09, 1e-07, 1e-05, 0.001, 0.1, 0.9],
                        'g0' : [0.001, 0.0026, 0.007, 0.0184, 0.0487, 0.129, 0.34],
                        'S_sensitivity_D1_to_g0' : [5.6e-4, 0.0054, 0.0517, 0.496, 4.77, 45.8, 440.0],
                        'S_sensitivity_D2_to_S_sensitivity_D1' : [0.0004, 0.0015, 0.0054, 0.02, 0.0737, 0.271, 1.0],
                        'R1_sensitivity_D1_to_S_sensitivity_D1' : [0, 1e-05, 9.56e-05, 0.000915, 0.0087, 0.0837, 0.8],
                        'R2_sensitivity_D2_to_S_sensitivity_D2' : [0, 1e-05, 9.56e-05, 0.000915, 0.0087, 0.0837, 0.8],
                        'S_transition_R1' : [1e-11, 2.15e-10, 4.64e-09, 1e-07, 2.15e-06, 4.64e-05, 0.001],
                        'S_transition_R2' : [1e-11, 2.15e-10, 4.64e-09, 1e-07, 2.15e-06, 4.64e-05, 0.001],
                        }
  
  valid_parameters = {'Parameter_ID' : input_param_df['Parameter_ID'] }
  for input_param in valid_input_values.keys():
    if input_param in input_param_df.columns:
      input_param_df[input_param] = input_param_df[input_param].apply(lambda x: round_to_significant_digits(x,3))
      valid_parameters[input_param] = input_param_df[input_param].apply(lambda x: x in valid_input_values[input_param]) 
  invalid_parameters_df = pd.DataFrame(valid_parameters)
  invalid_parameters_df = invalid_parameters_df[invalid_parameters_df.applymap(lambda x: x is False).any(axis=1)]
  return invalid_parameters_df

def check_output_params():
  #check transitions and rates
  return

def merge_input_output_params(output_param_df, input_param_df):
  return pd.merge(input_param_df, output_param_df, on='Parameter_ID', how='inner')

def map_parameters(sim_run_id, results_dir, mapped_dir):
  print(str(sim_run_id))
  param_file = os.path.join(sim_results_dir,'param_ALLDRUG_' + str(sim_run_id) + '.csv')

  param_df = get_param_df(param_file)
  input_param_df = gen_input_params(param_df)
  invalid_params_df = check_input_params(input_param_df)

  param_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_simParamOutput.csv'), header=True, index=False)
  input_param_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_inputParam.csv'), header=True, index=False)
  if len(invalid_params_df) > 0:
    invalid_params_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_invalidParams.csv'), header=True, index=False)

sim_results_dir = '~/sim_trial_results/'
output_dir = os.path.join(sim_results_dir,'preprocessed_sim_trial_results')
os.makedirs(output_dir, exist_ok=True)

param_files = os.listdir(sim_results_dir)
run_id_list = list(set([ os.path.splitext(os.path.basename(x))[0].split("_")[2] for x in param_files ]))

#map_parameters(run_id_list[7],sim_results_dir, output_dir)
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Submit tasks to the executor
    futures = [executor.submit(map_parameters, run_id, sim_results_dir, output_dir) for run_id in run_id_list]