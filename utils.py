import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass
import numpy as np
import random
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from os import path, mkdir
from collections import Counter

def get_output_basename(fname):
    basename = path.basename(fname)
    splt = path.splitext(basename)
    basename = splt[0]
    if 'pep' not in splt[1].lower():
        basename = path.splitext(basename)[0]
    return basename

def get_output_folder(outfolder, fname):
    if not outfolder:
        return path.dirname(path.realpath(fname))
    else:
        tmp_outfolder = path.join(path.dirname(path.realpath(fname)), outfolder)
        if not path.isdir(tmp_outfolder):
            mkdir(tmp_outfolder)
        return tmp_outfolder

def calc_RT(seq, RC):
    try:
        return achrom.calculate_RT(seq, RC)
    except:
        return 0
    
def is_decoy(proteins, decoy_prefix='DECOY_'):
    return all(z.startswith(decoy_prefix) for z in proteins)

def split_decoys(df):
    df['decoy1'] = df['decoy'].apply(lambda x: True if x and random.random() >= 0.5 else False)
    df['decoy2'] = df.apply(lambda x: x['decoy'] and not x['decoy1'], axis=1)
    return df

def remove_column_hit_rank(df):
    if 'hit_rank' in df.columns:
        return df[df['hit_rank'] == 1]
    else:
        print('no hit_rank column')
        return df

def parse_mods(df_raw):
    mods_counter = Counter()
    sequence, mods = df_raw['peptide'], df_raw['modifications']
    if mods:
        for mod in mods.split(','):
            mod_mass, aa_ind = mod.split('@')
            mod_mass = float(mod_mass)
            aa_ind = int(aa_ind)
            if aa_ind == 0:
                aa = 'N_term'
                mod_mass = round(mod_mass - 1.007825, 3)
            elif aa_ind == len(sequence) + 1:
                aa = 'C_term'
                mod_mass = round(mod_mass - 17.002735, 3)
            else:
                aa = sequence[aa_ind-1]
                mod_mass = round(mod_mass - mass.std_aa_mass[aa], 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] += 1
    return mods_counter

def add_mod_info(df_raw, mod):
    sequence, mods_counter = df_raw['peptide'], df_raw['mods_counter']
    mod_aa = mod.split(' at ')[1]
    if 'term' not in mod_aa and mod_aa not in sequence:
        return -1
    else:
        return mods_counter.get(mod, 0)

def prepare_mods(df):
    all_mods = set()
    for cnt in df['mods_counter'].values:
        for k in cnt.keys():
            all_mods.add(k)
    for mod in all_mods:
        df[mod] = df.apply(add_mod_info, axis=1, mod=mod)
    return df

def prepare_dataframe_xtandem(infile_path, decoy_prefix='DECOY_'):
    df1 = pepxml.DataFrame(infile_path)
    df1['length'] = df1['peptide'].apply(len)
    df1 = df1[df1['length'] >= 6]
    df1['spectrum'] = df1['spectrum'].apply(lambda x: x.split(' RTINS')[0])
    df1['RT exp'] = df1['retention_time_sec'] / 60
    df1 = df1.drop(['retention_time_sec', ], axis=1)
    df1['massdiff_int'] = df1['massdiff'].apply(lambda x: int(round(x, 0)))
    df1['massdiff_ppm'] = 1e6 * df1['massdiff'] / df1['calc_neutral_pep_mass']
    df1['decoy'] = df1['protein'].apply(is_decoy, decoy_prefix=decoy_prefix)
    df1 = split_decoys(df1)
    df1 = remove_column_hit_rank(df1)
    df1['mods_counter'] = df1.apply(parse_mods, axis=1)
    df1 = prepare_mods(df1)
    df1_f = aux.filter(df1, fdr=0.01, key='expect', is_decoy='decoy', correction=1)
    print('Default target-decoy filtering, 1%% PSM FDR: Number of target PSMs = %d' \
             % (df1_f[~df1_f['decoy']].shape[0]))
    print('Calibrating retention model...')
    retention_coefficients = achrom.get_RCs_vary_lcp(df1_f['peptide'].values, \
                                                     df1_f['RT exp'].values)
    df1_f['RT pred'] = df1_f['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
    df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
    _, _, r_value, std_value = aux.linear_regression(df1_f['RT pred'], df1_f['RT exp'])
    print('R^2 = %f , std = %f' % (r_value**2, std_value))
    df1['RT diff'] = df1['RT pred'] - df1['RT exp']
    return df1

def get_features(dataframe):
    feature_columns = dataframe.columns
    columns_to_remove = []
    for feature in feature_columns:
        if feature not in ['expect', 'hyperscore', 'calc_neutral_pep_mass', 'bscore', 'yscore', \
                            'massdiff', 'massdiff_ppm', 'nextscore', 'RT pred', 'RT diff', \
                            'sumI', 'RT exp', 'precursor_neutral_mass', \
                            'num_missed_cleavages', 'tot_num_ions', 'num_matched_ions', 'length']:
            columns_to_remove.append(feature)
    feature_columns = feature_columns.drop(columns_to_remove)
    return feature_columns

def get_X_array(df, feature_columns):
    return df.loc[:, feature_columns].values

def get_Y_array(df):
    return df.loc[:, 'decoy1'].values

def get_gbm_model(df):
    # train, test = train_test_split(df, test_size = 0.5)
#     print(type(feature_columns))
    # if type(feature_columns)==bool:
    feature_columns = get_features(df)
    # x_train = get_X_array(train, feature_columns)#train.iloc[:, 0:-1].values
    # y_train = get_Y_array(train)#train['target'].values
    x_train = get_X_array(df, feature_columns)#train.iloc[:, 0:-1].values
    y_train = get_Y_array(df)#train['target'].values
    # x_test = get_X_array(test, feature_columns)#test.iloc[:, 0:-1].values
    # y_test = get_Y_array(test)#test['target'].values
    # print(np.sum(y_train), len(y_train))
    # create dataset for lightgbm
    lgb_train = lgb.Dataset(x_train, y_train)
    # lgb_eval = lgb.Dataset(x_test, y_test, reference=lgb_train)

    # specify your configurations as a dict
    params = {
        'task': 'train',
        'boosting_type': 'gbdt',
        'objective': 'regression',
        'metric': {'l2', 'auc'},
#         'metric': 'multi_logloss',
        'num_leaves': 10,
        'learning_rate': 0.1,
        'feature_fraction': 0.95,
        'bagging_fraction': 0.95,
        'bagging_freq': 5,
        'verbose': 0,
       'min_data_in_bin': 1,
       'min_data': 1,
#        'min_sum_hessian_in_leaf': 0.001,
    }

#     print('Start training...')
    # train
    # gbm = lgb.train(params,
    #                 lgb_train,
    #                 num_boost_round=1000,
    #                 valid_sets=lgb_eval,
    #                 early_stopping_rounds=10, verbose_eval=1)
    
    
    cv_result_lgb = lgb.cv(params, 
                        lgb_train, 
                        num_boost_round=2000, 
                        nfold=9, 
                        stratified=True, 
                        early_stopping_rounds=10, 
                        verbose_eval=False, 
                        show_stdv=False)
    
    num_boost_rounds_lgb = len(cv_result_lgb['auc-mean'])
    # print(cv_result_lgb['auc-mean'])
    print('num_boost_rounds_lgb=' + str(num_boost_rounds_lgb))
    # train model
    gbm = lgb.train(params, lgb_train, num_boost_round=num_boost_rounds_lgb)
    
#     print('Save model...')
#     #save model to file
#     gbm.save_model('model.txt')

#     print('Start predicting...')
    # predict
    # y_pred = gbm.predict(X_test)
    # eval
#     print('The rmse of prediction is:', mean_squared_error(y_test, y_pred) ** 0.5)
    return gbm

def calc_PEP(df):
    feature_columns = get_features(df)
    gbm_model = get_gbm_model(df)
    x_all = get_X_array(df, feature_columns)
    df['PEP'] = gbm_model.predict(x_all, num_iteration=gbm_model.best_iteration)
    return df
