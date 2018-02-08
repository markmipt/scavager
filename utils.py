import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta
import numpy as np
import random
SEED = 42
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from os import path, mkdir
from collections import Counter, defaultdict


def keywithmaxval(d):
    #this method is much faster than using max(prots.iterkeys(), key=(lambda key: prots[key]))
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]

def get_protein_groups(df):
    pept_prots = defaultdict(set)
    prot_prots = defaultdict(set)
    prot_pepts = dict()
    for peptides, dbname in df[['peptides set', 'dbname']].values:
        prot_pepts[dbname] = peptides
        for peptide in peptides:
            pept_prots[peptide].add(dbname)
    for prots in pept_prots.values():
        for dbname in prots:
            prot_prots[dbname].update(prots)
    prot_pepts_count = dict()
    for k, v in prot_pepts.items():
        prot_pepts_count[k] = len(v)
    tostay = set()
    while pept_prots:
        bestprot = keywithmaxval(prot_pepts_count)
        tostay.add(bestprot)
        for pep in prot_pepts[bestprot]:
            for k in pept_prots[pep]:
                prot_pepts_count[k] -= 1
            del pept_prots[pep]
    df['groupleader'] = df['dbname'].apply(lambda x: x in tostay)
    df['all proteins'] = df['dbname'].apply(lambda x: ';'.join(prot_prots[x]))
    return df


def process_fasta(df, path_to_fasta):
    protsS = dict()
    for x in fasta.read(path_to_fasta):
        dbname = x[0].split(' ')[0]
        protsS[dbname] = x[1]
    df['sequence'] = df['dbname'].apply(lambda x: protsS.get(x, ''))
    return df

def get_proteins_dataframe(df1_f2, df1_peptides_f, decoy_prefix, all_decoys_2, path_to_fasta=False):
    print(path_to_fasta)
    proteins_dict = dict()
    for proteins, protein_descriptions, peptide, pep in df1_peptides_f[['protein', 'protein_descr', 'peptide', 'PEP']].values:
        for prot, prot_descr in zip(proteins, protein_descriptions):
            if prot not in proteins_dict:
                proteins_dict[prot] = dict()
                proteins_dict[prot]['dbname'] = prot
                proteins_dict[prot]['description'] = prot_descr
                proteins_dict[prot]['PSMs'] = 0
                proteins_dict[prot]['peptides set'] = set()
                proteins_dict[prot]['sequence'] = ''
                proteins_dict[prot]['NSAF'] = 0
                proteins_dict[prot]['sq'] = 0
                proteins_dict[prot]['score'] = dict()
                proteins_dict[prot]['q-value'] = 1.0
                proteins_dict[prot]['decoy'] = prot.startswith(decoy_prefix)
                proteins_dict[prot]['decoy2'] = prot in all_decoys_2
            proteins_dict[prot]['peptides set'].add(peptide)
            proteins_dict[prot]['score'][peptide] = min(proteins_dict[prot]['score'].get(peptide, 1.0), pep)

    for proteins in df1_f2[['protein']].values:
        for prot in proteins[0]:
            if prot in proteins_dict:
                proteins_dict[prot]['PSMs'] += 1
    df_proteins = pd.DataFrame.from_dict(proteins_dict, orient='index').reset_index()
    if path_to_fasta:
        df_proteins = process_fasta(df_proteins, path_to_fasta)
    df_proteins['length'] = df_proteins['sequence'].apply(len)
    df_proteins['sq'] = df_proteins.apply(calc_sq, axis=1)
    df_proteins['peptides'] = df_proteins['peptides set'].apply(len)
    df_proteins['score'] = df_proteins['score'].apply(lambda x: np.prod(list(x.values())))
    return df_proteins

def calc_sq(df_raw):
    protein = df_raw['sequence']
    peptides = df_raw['peptides set']
    if not protein:
        return 0
    psq = [False for x in protein]
    plen = len(protein)
    for pep in peptides:
        csize = len(pep)
        for j in range(plen):
            if protein[j:j+csize] == pep:
                for y in range(csize):
                    psq[j + y] = True
    return float(sum(psq)) / len(psq) * 100

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
    
def is_decoy(proteins, decoy_prefix):
    return all(z.startswith(decoy_prefix) for z in proteins)

def is_decoy_2(proteins, decoy_set):
    return all(z in decoy_set for z in proteins)

def split_decoys(df, decoy_prefix):
    all_decoys = set()
    for proteins in df[['protein']].values:
        for dbname in proteins[0]:
            if dbname.startswith(decoy_prefix):
                all_decoys.add(dbname)
    all_decoys = sorted(list(all_decoys)) # sort is done for working of random SEED
    random.seed(SEED)
    all_decoys_2 = set(random.sample(all_decoys, int(len(all_decoys)/2)))
    df['decoy2'] = df['protein'].apply(is_decoy_2, decoy_set=all_decoys_2)
    df['decoy1'] = df.apply(lambda x: x['decoy'] and not x['decoy2'], axis=1)
    return df, all_decoys_2

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
    df1, all_decoys_2 = split_decoys(df1, decoy_prefix=decoy_prefix)
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
    return df1, all_decoys_2

def get_features(dataframe):
    feature_columns = dataframe.columns
    columns_to_remove = []
    for feature in feature_columns:
        if feature not in ['expect', 'hyperscore', 'calc_neutral_pep_mass', 'bscore', 'yscore', \
                            'massdiff', 'massdiff_ppm', 'nextscore', 'RT pred', 'RT diff', \
                            'sumI', 'RT exp', 'precursor_neutral_mass', 'massdiff_int', \
                            'num_missed_cleavages', 'tot_num_ions', 'num_matched_ions', 'length']:
            if not feature.startswith('mass shift'):
                columns_to_remove.append(feature)
    feature_columns = feature_columns.drop(columns_to_remove)
    return feature_columns

def get_X_array(df, feature_columns):
    return df.loc[:, feature_columns].values

def get_Y_array(df):
    return df.loc[:, 'decoy1'].values

def get_gbm_model(df):
    feature_columns = get_features(df)
    x_train = get_X_array(df, feature_columns)
    y_train = get_Y_array(df)
    lgb_train = lgb.Dataset(x_train, y_train)

    params = {
        'task': 'train',
        'boosting_type': 'gbdt',
        'objective': 'regression',
        'metric': {'l2', 'auc'},
        'num_leaves': 10,
        'learning_rate': 0.1,
        'feature_fraction': 0.95,
        'feature_fraction_seed': SEED,
        'bagging_fraction': 0.95,
        'bagging_seed': SEED,
        'bagging_freq': 5,
        'verbose': 0,
        'min_data_in_bin': 1,
        'min_data': 1
    }
    
    cv_result_lgb = lgb.cv(params, 
                        lgb_train, 
                        num_boost_round=2000, 
                        nfold=9, 
                        stratified=True, 
                        early_stopping_rounds=10, 
                        verbose_eval=False, 
                        show_stdv=False, seed=SEED)
    
    num_boost_rounds_lgb = len(cv_result_lgb['auc-mean'])
    print('num_boost_rounds_lgb=' + str(num_boost_rounds_lgb))
    # train model
    gbm = lgb.train(params, lgb_train, num_boost_round=num_boost_rounds_lgb)
    return gbm

def calc_PEP(df):
    feature_columns = get_features(df)
    gbm_model = get_gbm_model(df)
    x_all = get_X_array(df, feature_columns)
    df['PEP'] = gbm_model.predict(x_all, num_iteration=gbm_model.best_iteration)
    pep_min = df['PEP'].min()
    df['log_score'] = np.log10(df['PEP'] - ((pep_min - 1e-15) if pep_min < 0 else 0))
    return df
