import os
import pandas as pd
import math
import numpy as np
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
  if os.path.exists(filename):
    param_df = pd.read_csv(filename,header=None)
    param_df.columns = param_column_names
    return param_df
  else:
    print('no param_df ' + filename)

def round_to_significant_digits(value, digits):
    if value == 0:
        return 0
    else:
        return round(value, digits - int(math.floor(math.log10(abs(value)))) - 1)

def gen_input_params(param_df):

  input_params_dict = {'Parameter_ID' : param_df['Parameter_ID'],
                       'initial_R1_percent' : [x/5e9 for x in param_df['R1_pop'] ],
                       'initial_R2_percent' : [x/5e9 for x in param_df['R2_pop'] ],
                       'g0' : param_df['g0'],
                       'S_sensitivity_D1_to_g0' : param_df['S_cell_sensitivity_D1']/param_df['g0'],
                       'S_sensitivity_D2_to_S_sensitivity_D1' : param_df['S_cell_sensitivity_D2']/param_df['S_cell_sensitivity_D1'],
                       'R1_sensitivity_D1_to_S_sensitivity_D1' : param_df['R1_cell_sensitivity_D1']/param_df['S_cell_sensitivity_D1'],
                       'R2_sensitivity_D2_to_S_sensitivity_D2' : param_df['R2_cell_sensitivity_D2']/param_df['S_cell_sensitivity_D2'],
                       'S_transition_to_R1': param_df['S_transition_to_R1'],
                       'S_transition_to_R2': param_df['S_transition_to_R2'],
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
                        'S_transition_to_R1' : [1e-11, 2.15e-10, 4.64e-09, 1e-07, 2.15e-06, 4.64e-05, 0.001],
                        'S_transition_to_R2' : [1e-11, 2.15e-10, 4.64e-09, 1e-07, 2.15e-06, 4.64e-05, 0.001],
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

  if not os.path.exists(os.path.join(output_dir,str(sim_run_id) + '_simParamOutput.csv')):
    param_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_simParamOutput.csv'), header=True, index=False)  
  if not os.path.exists(os.path.join(output_dir,str(sim_run_id) + '_inputParam.csv')):
    input_param_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_inputParam.csv'), header=True, index=False)
  if not os.path.exists(os.path.join(output_dir,str(sim_run_id) + '_invalidParams.csv')):
    invalid_params_df = check_input_params(input_param_df)
    if len(invalid_params_df) > 0:
      invalid_params_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_invalidParams.csv'), header=True, index=False)

def get_dosage_df(filename):
  param_column_names = ['Parameter_ID', 'Strategy_index',
                        'Drug2_0','Drug1_0',
                        'Drug2_45','Drug1_45',
                        'Drug2_90','Drug1_90',
                        'Drug2_135','Drug1_135',
                        'Drug2_180','Drug1_180',
                        'Drug2_225','Drug1_225',
                        'Drug2_270','Drug1_270',
                        'Drug2_315','Drug1_315',
                        'Drug2_360','Drug1_360',
                        'Drug2_405','Drug1_405',
                        'Drug2_450','Drug1_450',
                        'Drug2_495','Drug1_495',
                        'Drug2_540','Drug1_540',
                        'Drug2_585','Drug1_585',
                        'Drug2_630','Drug1_630',
                        'Drug2_675','Drug1_675',
                        'Drug2_720','Drug1_720',
                        'Drug2_765','Drug1_765',
                        'Drug2_810','Drug1_810',
                        'Drug2_855','Drug1_855',
                        'Drug2_900','Drug1_900',
                        'Drug2_945','Drug1_945',
                        'Drug2_990','Drug1_990',
                        'Drug2_1035','Drug1_1035',
                        'Drug2_1080','Drug1_1080',
                        'Drug2_1125','Drug1_1125',
                        'Drug2_1170','Drug1_1170',
                        'Drug2_1215','Drug1_1215',
                        'Drug2_1260','Drug1_1260',
                        'Drug2_1305','Drug1_1305',
                        'Drug2_1350','Drug1_1350',
                        'Drug2_1395','Drug1_1395',
                        'Drug2_1440','Drug1_1440',
                        'Drug2_1485','Drug1_1485',
                        'Drug2_1530','Drug1_1530',
                        'Drug2_1575','Drug1_1575',
                        'Drug2_1620','Drug1_1620',
                        'Drug2_1665','Drug1_1665',
                        'Drug2_1710','Drug1_1710',
                        'Drug2_1755','Drug1_1755',                        
                        ]
                
  if os.path.exists(filename):
    dosage_df = pd.read_csv(filename,header=None)
    dosage_df.columns = param_column_names
    # only need drug 1 since it defines drug 2 dosage
    melted_dosage_df = dosage_df.melt(id_vars=['Parameter_ID', 'Strategy_index'], value_vars=[col for col in dosage_df.columns if 'Drug1' in col], var_name='timepoint', value_name='Drug1_dosage')
    melted_dosage_df['timepoint'] = melted_dosage_df['timepoint'].str.extract(r'_(\d+)').astype(int)

    # add t = 1800
    new_rows = melted_dosage_df[['Parameter_ID', 'Strategy_index']].drop_duplicates()
    new_rows['timepoint'] = 1800
    new_rows['Drug1_dosage'] = -1

    melted_dosage_df = pd.concat([melted_dosage_df, new_rows], ignore_index=True)
    melted_dosage_df = melted_dosage_df.sort_values(by=['Parameter_ID', 'Strategy_index', 'timepoint']).reset_index(drop=True)

    strategy_map = { 0 : 'CPM', 1 : 'DPM', 2 : 'DPMtrial'}
    melted_dosage_df['Strategy_index'] = melted_dosage_df['Strategy_index'].replace(strategy_map)
    new_col_name = {'Strategy_index': 'Strategy_name'}
    melted_dosage_df = melted_dosage_df.rename(columns = new_col_name)


    return melted_dosage_df
  else:
    print('no dosage_df ' + filename)


def get_pop_df(filename, initial_pop_df):
  param_column_names = ['Parameter_ID', 'Strategy_index', 
                        'Spop_45', 'R1pop_45', 'R2pop_45', 'R12pop_45',
                        'Spop_90', 'R1pop_90', 'R2pop_90', 'R12pop_90',
                        'Spop_135', 'R1pop_135', 'R2pop_135', 'R12pop_135',
                        'Spop_180', 'R1pop_180', 'R2pop_180', 'R12pop_180',
                        'Spop_225', 'R1pop_225', 'R2pop_225', 'R12pop_225',
                        'Spop_270', 'R1pop_270', 'R2pop_270', 'R12pop_270',
                        'Spop_315', 'R1pop_315', 'R2pop_315', 'R12pop_315',
                        'Spop_360', 'R1pop_360', 'R2pop_360', 'R12pop_360',
                        'Spop_405', 'R1pop_405', 'R2pop_405', 'R12pop_405',
                        'Spop_450', 'R1pop_450', 'R2pop_450', 'R12pop_450',
                        'Spop_495', 'R1pop_495', 'R2pop_495', 'R12pop_495',
                        'Spop_540', 'R1pop_540', 'R2pop_540', 'R12pop_540',
                        'Spop_585', 'R1pop_585', 'R2pop_585', 'R12pop_585',
                        'Spop_630', 'R1pop_630', 'R2pop_630', 'R12pop_630',
                        'Spop_675', 'R1pop_675', 'R2pop_675', 'R12pop_675',
                        'Spop_720', 'R1pop_720', 'R2pop_720', 'R12pop_720',
                        'Spop_765', 'R1pop_765', 'R2pop_765', 'R12pop_765',
                        'Spop_810', 'R1pop_810', 'R2pop_810', 'R12pop_810',
                        'Spop_855', 'R1pop_855', 'R2pop_855', 'R12pop_855',
                        'Spop_900', 'R1pop_900', 'R2pop_900', 'R12pop_900',
                        'Spop_945', 'R1pop_945', 'R2pop_945', 'R12pop_945',
                        'Spop_990', 'R1pop_990', 'R2pop_990', 'R12pop_990',
                        'Spop_1035', 'R1pop_1035', 'R2pop_1035', 'R12pop_1035',
                        'Spop_1080', 'R1pop_1080', 'R2pop_1080', 'R12pop_1080',
                        'Spop_1125', 'R1pop_1125', 'R2pop_1125', 'R12pop_1125',
                        'Spop_1170', 'R1pop_1170', 'R2pop_1170', 'R12pop_1170',
                        'Spop_1215', 'R1pop_1215', 'R2pop_1215', 'R12pop_1215',
                        'Spop_1260', 'R1pop_1260', 'R2pop_1260', 'R12pop_1260',
                        'Spop_1305', 'R1pop_1305', 'R2pop_1305', 'R12pop_1305',
                        'Spop_1350', 'R1pop_1350', 'R2pop_1350', 'R12pop_1350',
                        'Spop_1395', 'R1pop_1395', 'R2pop_1395', 'R12pop_1395',
                        'Spop_1440', 'R1pop_1440', 'R2pop_1440', 'R12pop_1440',
                        'Spop_1485', 'R1pop_1485', 'R2pop_1485', 'R12pop_1485',
                        'Spop_1530', 'R1pop_1530', 'R2pop_1530', 'R12pop_1530',
                        'Spop_1575', 'R1pop_1575', 'R2pop_1575', 'R12pop_1575',
                        'Spop_1620', 'R1pop_1620', 'R2pop_1620', 'R12pop_1620',
                        'Spop_1665', 'R1pop_1665', 'R2pop_1665', 'R12pop_1665',
                        'Spop_1710', 'R1pop_1710', 'R2pop_1710', 'R12pop_1710',
                        'Spop_1755', 'R1pop_1755', 'R2pop_1755', 'R12pop_1755',
                        'Spop_1800', 'R1pop_1800', 'R2pop_1800', 'R12pop_1800']

                  
  if os.path.exists(filename):
    pop_df = pd.read_csv(filename,header=None)
    pop_df.columns = param_column_names

    # add t = 0
    pop_df['Spop_0'] = initial_pop_df['S_pop']
    pop_df['R1pop_0'] = initial_pop_df['R1_pop']
    pop_df['R2pop_0'] = initial_pop_df['R2_pop']
    pop_df['R12pop_0'] = initial_pop_df['R12_pop']

    melted_Spop_df = pop_df.melt(id_vars=['Parameter_ID', 'Strategy_index'], value_vars=[col for col in pop_df.columns if 'Spop' in col], var_name='timepoint', value_name='Spop')
    melted_Spop_df['timepoint'] = melted_Spop_df['timepoint'].str.extract(r'_(\d+)').astype(int)
    melted_R1pop_df = pop_df.melt(id_vars=['Parameter_ID', 'Strategy_index'], value_vars=[col for col in pop_df.columns if 'R1pop' in col], var_name='timepoint', value_name='R1pop')
    melted_R1pop_df['timepoint'] = melted_R1pop_df['timepoint'].str.extract(r'_(\d+)').astype(int)
    melted_R2pop_df = pop_df.melt(id_vars=['Parameter_ID', 'Strategy_index'], value_vars=[col for col in pop_df.columns if 'R2pop' in col], var_name='timepoint', value_name='R2pop')
    melted_R2pop_df['timepoint'] = melted_R2pop_df['timepoint'].str.extract(r'_(\d+)').astype(int)
    melted_R12pop_df = pop_df.melt(id_vars=['Parameter_ID', 'Strategy_index'], value_vars=[col for col in pop_df.columns if 'R12pop' in col], var_name='timepoint', value_name='R12pop')
    melted_R12pop_df['timepoint'] = melted_R12pop_df['timepoint'].str.extract(r'_(\d+)').astype(int)

    merged_pop_df = pd.merge(melted_Spop_df,melted_R1pop_df, on=['Parameter_ID', 'Strategy_index', 'timepoint'])
    merged_pop_df = pd.merge(merged_pop_df,melted_R2pop_df, on=['Parameter_ID', 'Strategy_index', 'timepoint'])
    merged_pop_df = pd.merge(merged_pop_df,melted_R12pop_df, on=['Parameter_ID', 'Strategy_index', 'timepoint'])


    #melted_dosage_df = pd.concat([melted_dosage_df, new_rows], ignore_index=True)
    merged_pop_df = merged_pop_df.sort_values(by=['Parameter_ID', 'Strategy_index', 'timepoint']).reset_index(drop=True)

    strategy_map = { 0 : 'CPM', 1 : 'DPM', 2 : 'DPMtrial'}
    merged_pop_df['Strategy_index'] = merged_pop_df['Strategy_index'].replace(strategy_map)
    new_col_name = {'Strategy_index': 'Strategy_name'}
    merged_pop_df = merged_pop_df.rename(columns = new_col_name)

    return merged_pop_df

    
  else:
    print('no pop_df ' + filename)


def map_trajectories(sim_run_id, sim_results_dir, output_dir):
  if not os.path.exists(os.path.join(output_dir,str(sim_run_id) + '_simTrajectories.csv')):
    param_file = os.path.join(sim_results_dir,'param_ALLDRUG_' + str(sim_run_id) + '.csv')
    dosage_file = os.path.join(sim_results_dir,'dosage_ALLDRUG_' + str(sim_run_id) + '.csv')
    population_file = os.path.join(sim_results_dir,'pop_ALLDRUG_' + str(sim_run_id) + '.csv')

    param_df = get_param_df(param_file)
    dosage_df = get_dosage_df(dosage_file)
    pop_df = get_pop_df(population_file, param_df[['S_pop','R1_pop','R2_pop','R12_pop']])
    traj_df = pd.merge(dosage_df, pop_df, on =['Parameter_ID', 'Strategy_name', 'timepoint'] )
    traj_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_simTrajectories.csv'), header = True, index = False)
  
def get_stopt_df(filename):
  param_column_names = ['Parameter_ID', 'Survival_CPM', 'Survival_DPM', 'Survival_DPMtrial']
  if os.path.exists(filename):
    stopt_df = pd.read_csv(filename,header=None)
    stopt_df.columns = param_column_names
  return stopt_df

def collect_EC_and_survival(sim_run_id, sim_results_dir, output_dir):
  if not os.path.exists(os.path.join(output_dir,str(sim_run_id) + '_ECsurvival.csv')):
    stopt_file = os.path.join(sim_results_dir,'stopt_ALLDRUG_' + str(sim_run_id) + '.csv')
    stopt_df = get_stopt_df(stopt_file)
    run_trajectory_df = pd.read_csv(os.path.join(output_dir,str(sim_run_id) + '_simTrajectories.csv'))
    param_id_list = []
    category_list = []
    for param_id in run_trajectory_df['Parameter_ID'].unique():
      param_id_list.append(param_id)
      param_id_check = run_trajectory_df['Parameter_ID'] == param_id
      time_0_check = run_trajectory_df['timepoint'] == 0
      time_45_check = run_trajectory_df['timepoint'] == 45

      param_trajectory_0_df = run_trajectory_df[param_id_check & time_0_check]
      param_trajectory_45_df = run_trajectory_df[param_id_check & time_45_check]
      #print(param_trajectory_df)
      CPM_drug_0 = param_trajectory_0_df[param_trajectory_0_df['Strategy_name'] == 'CPM']['Drug1_dosage'].iloc[0]
      DPM_drug_0 = param_trajectory_0_df[param_trajectory_0_df['Strategy_name'] == 'DPM']['Drug1_dosage'].iloc[0]
      CPM_drug_45 = param_trajectory_45_df[param_trajectory_45_df['Strategy_name'] == 'CPM']['Drug1_dosage'].iloc[0]
      DPM_drug_45 = param_trajectory_45_df[param_trajectory_45_df['Strategy_name'] == 'DPM']['Drug1_dosage'].iloc[0]
      
      if CPM_drug_0 == DPM_drug_0 and CPM_drug_45 == DPM_drug_45:
        category_list.append('both_same')
      elif CPM_drug_0 == DPM_drug_0 and CPM_drug_45 != DPM_drug_45:
        category_list.append('first_same_only')
      elif CPM_drug_0 != DPM_drug_0 and CPM_drug_45 == DPM_drug_45:
        category_list.append('second_same_only')
      else:
        category_list.append('both diff')

    category_df = pd.DataFrame({'Parameter_ID' : param_id_list,
                                'EC_category' : category_list})

    survival_df = pd.merge(stopt_df,category_df, on = ['Parameter_ID'])
    survival_df['DPMtrail_days_improvement'] = survival_df['Survival_DPMtrial'] - survival_df['Survival_CPM']
    survival_df['DPMtrail_percent_improvement'] = (survival_df['Survival_DPMtrial'] - survival_df['Survival_CPM'])/survival_df['Survival_DPMtrial']
    survival_df['DPM_days_improvement'] = survival_df['Survival_DPM'] - survival_df['Survival_CPM']
    survival_df['DPM_percent_improvement'] = (survival_df['Survival_DPM'] - survival_df['Survival_CPM'])/survival_df['Survival_DPM']
    # survival on cpm
    # months conversion = 356/12 = 30.4
    survival_df['CPM_surv_under_6_month'] = survival_df['Survival_CPM']/30.4 <= 6
    survival_df['CPM_surv_6_to_24_month'] = np.array(list(survival_df['Survival_CPM']/30.4 > 6)) & np.array(list(survival_df['Survival_CPM']/30.4 < 24))
    survival_df['CPM_surv_over_24_month'] = survival_df['Survival_CPM']/30.4 >= 24

    # improvement on DPM
    survival_df['DPMtrial_improve_under_6_month'] = survival_df['DPMtrail_days_improvement']/30.4 <= 6
    survival_df['DPMtrial_improve_6_to_24_month'] = np.array(list(survival_df['DPMtrail_days_improvement']/30.4 > 6)) & np.array(list(survival_df['DPMtrail_days_improvement']/30.4 < 24))
    survival_df['DPMtrial_improve_over_24_month'] = survival_df['DPMtrail_days_improvement']/30.4 >= 24

    survival_df['DPM_improve_under_6_month'] = survival_df['DPM_days_improvement']/30.4 <= 6
    survival_df['DPM_improve_6_to_24_month'] = np.array(list(survival_df['DPM_days_improvement']/30.4 > 6)) & np.array(list(survival_df['DPM_days_improvement']/30.4 < 24))
    survival_df['DPM_improve_over_24_month'] = survival_df['DPM_days_improvement']/30.4 >= 24

    survival_df['DPMtrial_improve_25_to_50_percent'] = np.array(list(survival_df['DPMtrail_percent_improvement'] >= 0.25)) & np.array(list(survival_df['DPMtrail_percent_improvement'] <= 0.5))
    survival_df['DPMtrial_improve_50_to_75_percent'] = np.array(list(survival_df['DPMtrail_percent_improvement'] > 0.5)) & np.array(list(survival_df['DPMtrail_percent_improvement'] < 0.75))
    survival_df['DPMtrial_improve_over_75_percent'] = survival_df['DPMtrail_percent_improvement'] >= 0.75

    survival_df['DPM_improve_25_to_50_percent'] = np.array(list(survival_df['DPM_percent_improvement'] >= 0.25)) & np.array(list(survival_df['DPM_percent_improvement'] <= 0.5))
    survival_df['DPM_improve_50_to_75_percent'] = np.array(list(survival_df['DPM_percent_improvement'] > 0.5)) & np.array(list(survival_df['DPM_percent_improvement'] < 0.75))
    survival_df['DPM_improve_over_75_percent'] = survival_df['DPM_percent_improvement'] >= 0.75

    survival_df.to_csv(os.path.join(output_dir,str(sim_run_id) + '_ECsurvival.csv'), header = True, index = False)
    
def process_sim_output(sim_run_id, sim_results_dir, output_dir):
  map_parameters(sim_run_id, sim_results_dir, output_dir)
  map_trajectories(sim_run_id, sim_results_dir, output_dir)
  collect_EC_and_survival(sim_run_id, sim_results_dir, output_dir) 

sim_results_dir = '/home/mdm299/sim_trial_results/'
output_dir = os.path.join(sim_results_dir,'processed_results')
os.makedirs(output_dir, exist_ok=True)

param_files = os.listdir(sim_results_dir)
run_id_list = list(set([ os.path.splitext(os.path.basename(x))[0].split("_")[2] for x in param_files ]))


#map_parameters(run_id_list[7],sim_results_dir, output_dir)
with concurrent.futures.ThreadPoolExecutor() as executor:
  futures = [executor.submit(process_sim_output, run_id, sim_results_dir, output_dir) for run_id in run_id_list]
