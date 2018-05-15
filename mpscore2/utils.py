from __future__ import division
import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta, mzid, parser
import numpy as np
import random
SEED = 42
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from os import path, mkdir
from collections import Counter, defaultdict
import warnings
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'

def convert_tandem_cleave_rule_to_regexp(cleavage_rule):

    def get_sense(c_term_rule, n_term_rule):
        if '{' in c_term_rule:
            return 'N'
        elif '{' in n_term_rule:
            return 'C'
        else:
            if len(c_term_rule) <= len(n_term_rule):
                return 'C'
            else:
                return 'N'

    def get_cut(cut, no_cut):
        aminoacids = set(parser.std_amino_acids)
        cut = ''.join(aminoacids & set(cut))
        if '{' in no_cut:
            no_cut = ''.join(aminoacids & set(no_cut))
            return cut, no_cut
        else:
            no_cut = ''.join(set(parser.std_amino_acids) - set(no_cut))
            return cut, no_cut

    out_rules = []
    for protease in cleavage_rule.split(','):
        protease = protease.replace('X', ''.join(parser.std_amino_acids))
        c_term_rule, n_term_rule = protease.split('|')
        sense = get_sense(c_term_rule, n_term_rule)
        if sense == 'C':
            cut, no_cut = get_cut(c_term_rule, n_term_rule)
        else:
            cut, no_cut = get_cut(n_term_rule, c_term_rule)

        if no_cut:
            if sense == 'C':
                out_rules.append('([%s](?=[^%s]))' % (cut, no_cut))
            else:
                out_rules.append('([^%s](?=[%s]))' % (no_cut, cut))
        else:
            if sense == 'C':
                out_rules.append('([%s])' % (cut, ))
            else:
                out_rules.append('(?=[%s])' % (cut, ))
    return '|'.join(out_rules)

def calc_NSAF(df):
    df['NSAF'] = df['PSMs'] / df['length']
    NSAF_sum = np.sum(df['NSAF'])
    df['NSAF'] = df['NSAF'] / NSAF_sum
    if sum(pd.notna(df['NSAF'])):
        df['LOG10_NSAF'] = np.log10(df['NSAF'])
    return df

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
    proteins_dict = dict()
    for proteins, protein_descriptions, peptide, pep in df1_peptides_f[['protein', 'protein_descr', 'peptide', 'ML score']].values:
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
    df_proteins = calc_NSAF(df_proteins)
    df_proteins['peptides'] = df_proteins['peptides set'].apply(len)
    df_proteins['score'] = df_proteins['score'].apply(lambda x: np.prod(list(x.values())))
    return df_proteins

def calc_sq(df_raw):
    protein = df_raw['sequence']
    protein = protein.replace('L', 'I')
    peptides = df_raw['peptides set']
    if not protein:
        return 0
    psq = [False for x in protein]
    plen = len(protein)
    for pep in peptides:
        pep = pep.replace('L', 'I')
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

def parse_mods_msgf(df_raw):
    mods_counter = Counter()
    sequence, mods = df_raw['peptide'], df_raw['Modification']
    if isinstance(mods, list):
        for mod in mods:
            mod_mass, aa_ind = mod['monoisotopicMassDelta'], mod['location']
            if aa_ind == 0:
                aa = 'N_term'
                mod_mass = round(mod_mass, 3)
            elif aa_ind == len(sequence) + 1:
                aa = 'C_term'
                mod_mass = round(mod_mass, 3)
            else:
                aa = sequence[aa_ind-1]
                mod_mass = round(mod_mass, 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] += 1
    return mods_counter


def add_mod_info(df_raw, mod):
    sequence, mods_counter = df_raw['peptide'], df_raw['mods_counter']
    mod_aa = mod.split(' at ')[1]
    if 'term' not in mod_aa and mod_aa not in sequence:
        return -1
    else:
        return mods_counter.get(mod, 0) >= 1

def prepare_mods(df):
    all_mods = set()
    for cnt in df['mods_counter'].values:
        for k in cnt.keys():
            all_mods.add(k)
    for mod in all_mods:
        df[mod] = df.apply(add_mod_info, axis=1, mod=mod)
    return df

def prepare_dataframe_xtandem(infile_path, decoy_prefix='DECOY_', cleavage_rule=False, allowed_peptides=False):
    if not cleavage_rule:
        cleavage_rule = parser.expasy_rules['trypsin']
    if infile_path.endswith('.pep.xml'):
        df1 = pepxml.DataFrame(infile_path)
        ftype = 'pepxml'
    else:
        df1 = mzid.DataFrame(infile_path)

    if 'Morpheus Score' in df1.columns:
        df1 = df1[df1['Morpheus Score'] != 0]
        df1['expect'] = 1 / df1['Morpheus Score']
        df1['num_missed_cleavages'] = df1['peptide'].apply(lambda x: parser.num_sites(x, rule=cleavage_rule))
        
    if 'MS-GF:EValue' in df1.columns:
        #MSGF search engine
        ftype = 'msgf'
        df1['peptide'] = df1['PeptideSequence']
        df1['num_missed_cleavages'] = df1['peptide'].apply(lambda x: parser.num_sites(x, rule=cleavage_rule))
        df1['assumed_charge'] = df1['chargeState']
        df1['spectrum'] = df1['spectrumID']
        df1['massdiff'] = (df1['experimentalMassToCharge'] - df1['calculatedMassToCharge']) * df1['assumed_charge']
        df1['calc_neutral_pep_mass'] = df1['calculatedMassToCharge'] * df1['chargeState'] - df1['chargeState'] * 1.00727649
        df1['protein'] = df1['protein description']
        df1['protein_descr'] = df1['protein description']
        df1['expect'] = df1['MS-GF:EValue']

    df1 = df1[~pd.isna(df1['peptide'])]
    if allowed_peptides:
        df1 = df1[df1['peptide'].apply(lambda x: x in allowed_peptides)]
    df1['length'] = df1['peptide'].apply(len)
    df1 = df1[df1['length'] >= 6]
    df1['spectrum'] = df1['spectrum'].apply(lambda x: x.split(' RTINS')[0])
    if 'retention_time_sec' not in df1.columns:
        if 'scan start time' in df1.columns:
            df1['RT exp'] = df1['scan start time']
            df1 = df1.drop(['scan start time', ], axis=1)
        else:
            df1['RT exp'] = 0
    else:
        df1['RT exp'] = df1['retention_time_sec'] / 60
        df1 = df1.drop(['retention_time_sec', ], axis=1)

    df1['massdiff_int'] = df1['massdiff'].apply(lambda x: int(round(x, 0)))
    df1['massdiff_ppm'] = 1e6 * (df1['massdiff'] - df1['massdiff_int'] * 1.003354)/ df1['calc_neutral_pep_mass']

    df1['decoy'] = df1['protein'].apply(is_decoy, decoy_prefix=decoy_prefix)
    df1, all_decoys_2 = split_decoys(df1, decoy_prefix=decoy_prefix)
    df1 = remove_column_hit_rank(df1)
    # try:
    if ftype == 'pepxml':
        df1['mods_counter'] = df1.apply(parse_mods, axis=1)
    elif ftype == 'msgf':
        df1['mods_counter'] = df1.apply(parse_mods_msgf, axis=1)
    df1 = prepare_mods(df1)
    # except:
    #     pass

    df1_f = aux.filter(df1, fdr=0.01, key='expect', is_decoy='decoy', correction=1, remove_decoy=True, formula=1)
    print('Default target-decoy filtering, 1%% PSM FDR: Number of target PSMs = %d' \
             % (df1_f[~df1_f['decoy']].shape[0]))
    try:
        print('Calibrating retention model...')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            retention_coefficients = achrom.get_RCs_vary_lcp(df1_f['peptide'].values, \
                                                        df1_f['RT exp'].values)
        df1_f['RT pred'] = df1_f['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
        df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
        _, _, r_value, std_value = aux.linear_regression(df1_f['RT pred'], df1_f['RT exp'])
        print('R^2 = %f , std = %f' % (r_value**2, std_value))
        df1['RT diff'] = df1['RT pred'] - df1['RT exp']
    except:
        df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, achrom.RCs_krokhin_100A_tfa))
        df1['RT diff'] = df1['RT exp']
    return df1, all_decoys_2

def get_features(dataframe):
    feature_columns = dataframe.columns
    columns_to_remove = []
    for feature in feature_columns:
        if feature not in ['expect', 'hyperscore', 'calc_neutral_pep_mass', 'bscore', 'yscore', \
                            'massdiff', 'massdiff_ppm', 'nextscore', 'RT pred', 'RT diff', \
                            'sumI', 'RT exp', 'precursor_neutral_mass', 'massdiff_int', \
                            'num_missed_cleavages', 'tot_num_ions', 'num_matched_ions', 'length', \
                            'ScoreRatio', 'Energy', 'MS2IonCurrent', 'MeanErrorTop7', 'sqMeanErrorTop7', 'StdevErrorTop7', \
                            'MS-GF:DeNovoScore', 'MS-GF:EValue', 'MS-GF:RawScore', 'MeanErrorAll', \
                            'MeanRelErrorAll', 'MeanRelErrorTop7', 'NumMatchedMainIons', 'StdevErrorAll', \
                            'StdevErrorTop7', 'StdevRelErrorAll', 'StdevRelErrorTop7', 'NTermIonCurrentRatio', \
                            'CTermIonCurrentRatio', 'ExplainedIonCurrentRatio']:
            if not feature.startswith('mass shift'):
                columns_to_remove.append(feature)
    feature_columns = feature_columns.drop(columns_to_remove)
    return feature_columns

def get_X_array(df, feature_columns):
    return df.loc[:, feature_columns].values

def get_Y_array(df):
    return df.loc[:, 'decoy1'].values

def get_cat_model(df):
    feature_columns = get_features(df)
    train, test = train_test_split(df, test_size = 0.3)
    x_train = get_X_array(train, feature_columns)
    y_train = get_Y_array(train)
    x_test = get_X_array(test, feature_columns)
    y_test = get_Y_array(test)
    # model = CatBoostClassifier(iterations=1000, learning_rate=0.05, depth=10, loss_function='Logloss', logging_level='Silent', random_seed=SEED)
    # model.fit(x_train, y_train, use_best_model=True, eval_set=(x_test, y_test))

    model = CatBoostClassifier(iterations=2500, learning_rate=0.005, depth=6, loss_function='Logloss', eval_metric='Logloss',
                               od_type='Iter', od_wait=5, random_seed=SEED, logging_level='Silent')
    model.fit(x_train, y_train, use_best_model=True, eval_set=(x_test, y_test))

    return model

def calc_PEP(df):
    feature_columns = get_features(df)
    cat_model = get_cat_model(df)
    x_all = get_X_array(df, feature_columns)
    df['ML score'] = cat_model.predict_proba(x_all)[:, 1]
    pep_min = df['ML score'].min()
    df['log_score'] = np.log10(df['ML score'] - ((pep_min - 1e-15) if pep_min < 0 else 0))
    return df

def calc_target_decoy_ratio(df):
    df = df.sort_values(by='expect', ascending=False)
    zvals = np.cumsum(df['decoy'].values) / (np.arange(1, df.shape[0] + 1, 1) - np.cumsum(df['decoy'].values))
    zvals = zvals[~np.isnan(zvals)]
    zvals = zvals[np.isfinite(zvals)]
    zvals = zvals[zvals >= 0.5]
    H1 = np.histogram(zvals, bins=np.arange(0.1, 1.1, 0.01))
    return H1[1][np.argmax(H1[0])]

def get_columns_to_output(out_type):
    if out_type == 'psm':
        return ['peptide', 'length', 'spectrum', 'q', 'ML score', 'modifications', 'assumed_charge', 'num_missed_cleavages', 'num_tol_term', 'peptide_next_aa',
         'peptide_prev_aa', 'calc_neutral_pep_mass', 'massdiff_ppm', 'massdiff_int', 'RT exp', 'RT pred', 'protein', 'protein_descr', 'decoy']
    elif out_type == 'peptide':
        return ['peptide', '#PSMs', 'length', 'spectrum', 'q', 'ML score', 'modifications', 'assumed_charge', 'num_missed_cleavages', 'num_tol_term', 'peptide_next_aa',
         'peptide_prev_aa', 'calc_neutral_pep_mass', 'massdiff_ppm', 'massdiff_int', 'RT exp', 'RT pred', 'protein', 'protein_descr', 'decoy']
    elif out_type == 'protein':
        return ['dbname','description','PSMs','peptides','NSAF','sq','score','length', 'all proteins', 'groupleader']

def calc_psms(df):
    peptides = Counter(df['peptide'])
    df['#PSMs'] = df['peptide'].apply(lambda x: peptides.get(x))
    return df