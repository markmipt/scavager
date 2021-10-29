from __future__ import division
import pandas as pd
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta, mzid, parser
import numpy as np
import random
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
import os
from collections import Counter, defaultdict
from .utils_figures import get_fdbinsize
from scipy.stats import scoreatpercentile
from sklearn.isotonic import IsotonicRegression
import logging
import warnings
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'
logger = logging.getLogger(__name__)
SEED = 42


class NoDecoyError(ValueError):
    pass


class WrongInputError(NotImplementedError):
    pass


class EmptyFileError(ValueError):
    pass


def filter_custom(df, fdr, key, is_decoy, reverse, remove_decoy, ratio, formula, correction=None, loglabel=None):
    kw = dict(key=key, is_decoy=is_decoy, reverse=reverse, full_output=True,
        remove_decoy=False, ratio=ratio, formula=formula)
    df = df.copy()
    q = aux.qvalues(df, correction=1, **kw)
    q_uncorr = aux.qvalues(df, correction=0, **kw)
    df['q'] = q['q']
    df['q_uncorrected'] = q_uncorr['q']

    if correction is not None:
        qlabel = 'q' if correction else 'q_uncorrected'
        logger.debug('Explicitly using %s for filtering.', qlabel)
    elif df['q'].min() < fdr:
        logger.debug('Successfully filtered with +1 correction (label = %s).', loglabel)
        qlabel = 'q'
    else:
        logger.info('No results for filtering with +1 correction (label = %s). Rerunning without correction...', loglabel)
        qlabel = 'q_uncorrected'
    if remove_decoy:
        df = df[~df[is_decoy]]
    return df[df[qlabel] < fdr].copy()


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


def calc_TOP3(df):
    df['TOP3'] = df['TOP3'].apply(lambda z: sum(sorted(z, reverse=True)[:3]))


def calc_NSAF(df):
    df['NSAF'] = df['PSMs'] / df['length']
    NSAF_sum = np.sum(df['NSAF'])
    df['NSAF'] = df['NSAF'] / NSAF_sum
    if sum(pd.notna(df['NSAF'])):
        df['LOG10_NSAF'] = np.log10(df['NSAF'])
    df.loc[pd.isna(df['NSAF']), 'NSAF'] = 0.0


def keywithmaxval(d):
    # this method is much faster than using max(prots.iterkeys(), key=(lambda key: prots[key]))
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def add_protein_groups(df, df_ms1_path, sort_random=True):
    # input: df - pandas dataframe which contains columns "dbname" with protein names and
    #        "peptides set" with set of peptide sequences belong to this protein and identified in MS/MS analysis.
    #        df_ms1_path - path to the table *_proteins_full_noexclusion.tsv which is output by DirectMS1 analysis.
    #        sort_random - if True, the proteins with same scores are chosen randomly. Otherwise, they are chosen
    #        in alphabetical order.
    # output: pandas dataframe with new columns "groupleader" (True for protein group leaders) and
    #         "all proteins" (list of proteins belong to protein group and separeted by ';')

    pept_prots = defaultdict(set)
    prot_prots = defaultdict(set)
    prot_pepts = dict()
    if not sort_random:
        iter_list = df.sort_values(by='dbname').reset_index(drop=True)[['peptides set', 'dbname']].values
    else:
        iter_list = df.sample(frac=1).reset_index(drop=True)[['peptides set', 'dbname']].values
    for peptides, dbname in iter_list:
        prot_pepts[dbname] = peptides
        for peptide in peptides:
            pept_prots[peptide].add(dbname)
    for prots in pept_prots.values():
        for dbname in prots:
            prot_prots[dbname].update(prots)
    prot_pepts_count = dict()
    prot_pepts_count2 = dict()

    if not df_ms1_path:
        for k, v in prot_pepts.items():
            prot_pepts_count[k] = len(v)
            prot_pepts_count2[k] = len(v)

    else:
        df_ms1 = pd.read_csv(df_ms1_path, sep='\t')
        ms1s = dict()
        for qval, prot in df_ms1[['score', 'dbname']].values:
            ms1s[prot] = float(qval)

        max_k = max(ms1s.values())
        for k, v in prot_pepts.items():
            prot_pepts_count[k] = len(v) + ms1s.get(k, 0) / max_k
            prot_pepts_count2[k] = len(v)
    tostay = set()
    while pept_prots:
        bestprot = keywithmaxval(prot_pepts_count)
        tostay.add(bestprot)
        for pep in prot_pepts[bestprot]:
            for k in pept_prots[pep]:
                prot_pepts_count[k] -= 1
                prot_pepts_count2[k] -= 1
            del pept_prots[pep]

        for k, v in list(prot_pepts_count2.items()):
            if v == 0:
                del prot_pepts_count[k]
                del prot_pepts_count2[k]
    df['groupleader'] = df['dbname'].apply(lambda x: x in tostay)
    df['all proteins'] = df['dbname'].apply(lambda x: ';'.join(prot_prots[x]))


def process_fasta(df, path_to_fasta, decoy_prefix, decoy_infix=False):
    protsS = dict()
    decoy_check_flag = False
    for x in fasta.read(path_to_fasta):
        dbname = x[0].split(' ')[0]
        if not decoy_check_flag:
            if (not decoy_infix and dbname.startswith(decoy_prefix)) or (decoy_infix and decoy_infix in dbname):
                decoy_check_flag = True
        protsS[dbname] = x[1]
    df['sequence'] = df['dbname'].apply(lambda x: protsS.get(x, protsS.get(x.split(' ')[0], '')))
    if not decoy_check_flag:
        if not decoy_infix:
            df['sequence'] = df.apply(
                lambda x: x['sequence'] if x['sequence'] else protsS.get(
                    x['dbname'].replace(decoy_prefix, ''), protsS.get(x['dbname'].split(' ')[0].replace(decoy_prefix, ''), '')),
                axis=1)
        else:
            df['sequence'] = df.apply(
                lambda x: x['sequence'] if x['sequence'] else protsS.get(
                    x['dbname'].replace(decoy_infix, ''), protsS.get(x['dbname'].split(' ')[0].replace(decoy_infix, ''), '')),
                axis=1)

    return df


def get_tag_names(columns):
    return [c for c in columns if c[:4] == 'tag_']


def get_proteins_dataframe(df1_f2, decoy_prefix, all_decoys_2, decoy_infix=False, path_to_fasta=False, pif_threshold=0):
    proteins_dict = dict()
    cols = ['protein', 'protein_descr', 'peptide', 'PEP', 'MS1Intensity', 'PIF']
    tagnames = get_tag_names(df1_f2.columns)
    tagsums = []
    for c in tagnames:
        cols.append(c)
        tagsums.append(0.)
    have_pif = 'PIF' in df1_f2.columns
    if not have_pif:
        logger.debug('PIF not found.')
        if pif_threshold > 0:
            logger.warning('PIF not found, threshold will not be applied.')
        df1_f2['PIF'] = 0.

    for tup in df1_f2[cols].values:
        # avoid star unpacking to support Python 2.7 for a little longer
        (proteins, protein_descriptions, peptide, pep, ms1_i, pif), tags = tup[:6], tup[6:]
        for prot, prot_descr in zip(proteins, protein_descriptions):
            if prot not in proteins_dict:
                proteins_dict[prot] = dict()
                proteins_dict[prot]['dbname'] = prot
                proteins_dict[prot]['description'] = prot_descr
                proteins_dict[prot]['PSMs'] = 0
                proteins_dict[prot]['peptides set'] = set()
                proteins_dict[prot]['sequence'] = ''
                proteins_dict[prot]['NSAF'] = 0
                proteins_dict[prot]['TOP3'] = []
                proteins_dict[prot]['sq'] = 0
                proteins_dict[prot]['score'] = dict()
                proteins_dict[prot]['q-value'] = 1.0
                for tag in tagnames:
                    proteins_dict[prot][tag] = 0.
                if not decoy_infix:
                    proteins_dict[prot]['decoy'] = prot.startswith(decoy_prefix)
                else:
                    proteins_dict[prot]['decoy'] = decoy_infix in prot
                proteins_dict[prot]['decoy2'] = prot in all_decoys_2
                proteins_dict[prot]['decoy1'] = proteins_dict[prot]['decoy'] and not proteins_dict[prot]['decoy2']
            proteins_dict[prot]['peptides set'].add(peptide)
            proteins_dict[prot]['TOP3'].append(ms1_i)
            proteins_dict[prot]['score'][peptide] = min(proteins_dict[prot]['score'].get(peptide, 1.0), pep)
            proteins_dict[prot]['PSMs'] += 1
            if pif > pif_threshold:
                for tag, val in zip(tagnames, tags):
                    proteins_dict[prot][tag] += val
        if not have_pif or pif > pif_threshold:
            for i, (tag, val) in enumerate(zip(tagnames, tags)):
                tagsums[i] += val

    if not have_pif:
        df1_f2.drop('PIF', axis='columns', inplace=True)
    df_proteins = pd.DataFrame.from_dict(proteins_dict, orient='index').reset_index()
    if path_to_fasta:
        df_proteins = process_fasta(df_proteins, path_to_fasta, decoy_prefix, decoy_infix)
    df_proteins['length'] = df_proteins['sequence'].apply(len)
    df_proteins['sq'] = df_proteins.apply(calc_sq, axis=1)
    df_proteins['peptides'] = df_proteins['peptides set'].apply(len)
    df_proteins['PSMs'] = df_proteins.apply(lambda x: max(x['PSMs'], x['peptides']), axis=1)
    # df_proteins.loc[:, ['PSMs', 'peptides']].max(axis=1)
    calc_NSAF(df_proteins)
    calc_TOP3(df_proteins)
    df_proteins['score'] = df_proteins['score'].apply(lambda x: np.prod(list(x.values())))
    norm = np.array(tagsums)
    df_proteins[tagnames] /= norm
    return df_proteins, norm


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


def get_output_basename(fname, suffix=''):
    basename = os.path.basename(fname)
    splt = os.path.splitext(basename)
    basename = splt[0]
    if 'pep' not in splt[1].lower():
        basename = os.path.splitext(basename)[0]
    return basename + suffix


def get_output_folder(folder, fname):
    if not folder:
        return os.path.dirname(os.path.realpath(fname))
    else:
        if not os.path.isdir(folder):
            os.mkdir(folder)
        return folder


def calc_RT(seq, RC):
    try:
        return achrom.calculate_RT(seq, RC)
    except Exception:
        return 0


def is_decoy(proteins, decoy_prefix, decoy_infix=False):
    if not decoy_infix:
        return all(z.startswith(decoy_prefix) for z in proteins)
    else:
        return all(decoy_infix in z for z in proteins)


def is_group_specific(proteins, group_prefix, group_infix, decoy_prefix, decoy_infix=None):
    if group_infix:
        return all(group_infix in z for z in proteins)
    if not decoy_infix:
        return all(z.startswith(decoy_prefix + group_prefix) or z.startswith(group_prefix) for z in proteins)
    return all(z.startswith(group_prefix) for z in proteins)


def is_decoy_2(proteins, decoy_set):
    return all(z in decoy_set for z in proteins)


def split_fasta_decoys(db, decoy_prefix, decoy_infix=None):
    decoy_dbnames = set()
    with fasta.read(db) as f:
        for protein in f:
            dbname = protein.description.split()[0]
            if (decoy_infix and decoy_infix in dbname) or dbname.startswith(decoy_prefix):
                decoy_dbnames.add(dbname)
    decoy_dbnames = sorted(decoy_dbnames)
    random.seed(SEED)
    all_decoys_2 = set(random.sample(decoy_dbnames, len(decoy_dbnames) // 2))
    logger.debug('Marking %s out of %s decoys as decoy2', len(all_decoys_2), len(decoy_dbnames))
    return all_decoys_2


def split_decoys(df):
    """Split decoys into decoy1 and decoy2 without FASTA file"""
    all_decoys = set()
    for proteins in df.loc[df.decoy, 'protein'].values:
        all_decoys.update(proteins)
    logger.debug('proteins: %s', proteins)
    all_decoys = sorted(all_decoys)  # sort is done for working of random SEED
    random.seed(SEED)
    all_decoys_2 = set(random.sample(all_decoys, int(len(all_decoys) / 2)))
    df['decoy2'] = df['protein'].apply(is_decoy_2, decoy_set=all_decoys_2)
    df['decoy1'] = df.apply(lambda x: x['decoy'] and not x['decoy2'], axis=1)
    return all_decoys_2


def remove_column_hit_rank(df):
    if 'hit_rank' in df.columns:
        return df[df['hit_rank'] == 1]
    else:
        return df


def parse_mods(df_raw):
    mods_counter = {}
    sequence, mods = df_raw['peptide'], df_raw['modifications']
    if isinstance(mods, list):
        for mod in mods:
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
                aa = sequence[aa_ind - 1]
                mod_mass = round(mod_mass - mass.std_aa_mass[aa], 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] = mods_counter.get(mod_name, 0) + 1
    return mods_counter


def parse_mods_msgf(df_raw):
    mods_counter = {}
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
                aa = sequence[aa_ind - 1]
                mod_mass = round(mod_mass, 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] = mods_counter.get(mod_name, 0) + 1
    return mods_counter


def add_mod_info(df_raw, mod):
    sequence, mods_counter = df_raw['peptide'], df_raw['mods_counter']
    mod_aa = mod.split(' at ')[1]
    if 'term' not in mod_aa and mod_aa not in sequence:
        return -1
    else:
        return (1 if mods_counter.get(mod, 0) >= 1 else 0)


def prepare_mods(df):
    all_mods = set()
    for cnt in df['mods_counter'].values:
        for k in cnt.keys():
            all_mods.add(k)
    for mod in all_mods:
        df[mod] = df.apply(add_mod_info, axis=1, mod=mod)


def prepare_dataframe(infile_path, decoy_prefix=None, decoy_infix=False, cleavage_rule=False, fdr=0.01, decoy2set=None):
    if not cleavage_rule:
        cleavage_rule = parser.expasy_rules['trypsin']
    if infile_path.lower().endswith('.pep.xml') or infile_path.lower().endswith('.pepxml'):
        df1 = pepxml.DataFrame(infile_path)
        ftype = 'pepxml'
    elif infile_path.lower().endswith('.mzid'):
        df1 = mzid.DataFrame(infile_path)
    else:
        raise WrongInputError()
    if not df1.shape[0]:
        raise EmptyFileError()

    if 'Morpheus Score' in df1.columns:
        df1 = df1[df1['Morpheus Score'] != 0]
        df1['expect'] = 1 / df1['Morpheus Score']
        df1['num_missed_cleavages'] = df1['peptide'].apply(lambda x: parser.num_sites(x, rule=cleavage_rule))

    if 'MS-GF:EValue' in df1.columns:
        # MSGF search engine
        ftype = 'msgf'
        df1['peptide'] = df1['PeptideSequence']
        df1['num_missed_cleavages'] = df1['peptide'].apply(lambda x: parser.num_sites(x, rule=cleavage_rule))
        df1['assumed_charge'] = df1['chargeState']
        df1['spectrum'] = df1['spectrumID']
        df1['massdiff'] = (df1['experimentalMassToCharge'] - df1['calculatedMassToCharge']) * df1['assumed_charge']
        df1['calc_neutral_pep_mass'] = df1['calculatedMassToCharge'] * df1['chargeState'] - df1['chargeState'] * 1.00727649
        df1['protein'] = df1['accession']
        df1['protein_descr'] = df1['protein description']
        df1['expect'] = df1['MS-GF:EValue']

    if set(df1['protein_descr'].str[0]) == {None}:
        # MSFragger
        logger.debug('Adapting MSFragger DataFrame.')
        logger.debug('Proteins before: %s', df1.loc[1, 'protein'])
        protein = df1['protein'].apply(lambda row: [x.split(None, 1) for x in row])
        df1['protein'] = protein.apply(lambda row: [x[0] for x in row])
        try:
            df1['protein_descr'] = protein.apply(lambda row: [x[1] for x in row])
        except IndexError:
            df1['protein_descr'] = protein.apply(lambda row: ['' for x in row])
        logger.debug('Proteins after: %s', df1.loc[1, 'protein'])

    # if any(None in set(df1['protein_descr'].str[0])):
    #     print('HERE')
    #     df1['protein_descr'] = df1.apply(lambda x: x['protein_descr'] if x['protein_descr'] else x['protein'], axis=1)
    df1.loc[pd.isna(df1['protein_descr']), 'protein_descr'] = df1.loc[pd.isna(df1['protein_descr']), 'protein']
    # try:
    #     df1['expect'] = 1.0 / df1['bions_score_neg'].values
    # except:
    #     pass

    df1 = df1[~pd.isna(df1['peptide'])]
    if 'MS1Intensity' not in df1:
        df1['MS1Intensity'] = 0.0
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
    df1['massdiff_ppm'] = 1e6 * (df1['massdiff'] - df1['massdiff_int'] * 1.003354) / df1['calc_neutral_pep_mass']

    df1['decoy'] = df1['protein'].apply(is_decoy, decoy_prefix=decoy_prefix, decoy_infix=decoy_infix)
    if not df1.decoy.sum():
        raise NoDecoyError()
    if decoy2set is None:
        decoy2set = split_decoys(df1)
    else:
        df1['decoy2'] = df1['protein'].apply(lambda p: all(x in decoy2set for x in p))
        df1['decoy1'] = df1['decoy'] & (~df1['decoy2'])
    df1 = remove_column_hit_rank(df1)

    if ftype == 'pepxml':
        df1['mods_counter'] = df1.apply(parse_mods, axis=1)
    elif ftype == 'msgf':
        df1['mods_counter'] = df1.apply(parse_mods_msgf, axis=1)
    prepare_mods(df1)

    pep_ratio = df1['decoy2'].sum() / df1['decoy'].sum()
    df1_f = filter_custom(df1[~df1['decoy1']], fdr=fdr, key='expect', is_decoy='decoy2',
        reverse=False, remove_decoy=False, ratio=pep_ratio, formula=1, correction=None, loglabel='PSMs default')
    num_psms_def = df1_f[~df1_f['decoy2']].shape[0]
    logger.info('Default target-decoy filtering, 1%% PSM FDR: Number of target PSMs = %d', num_psms_def)
    try:
        logger.info('Calibrating retention model...')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            retention_coefficients = achrom.get_RCs_vary_lcp(df1_f['peptide'].values, df1_f['RT exp'].values)
        df1_f['RT pred'] = df1_f['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
        df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))
        _, _, r_value, std_value = aux.linear_regression(df1_f['RT pred'], df1_f['RT exp'])
        logger.info('RT model training results: R^2 = %f , std = %f', r_value**2, std_value)
        df1['RT diff'] = df1['RT pred'] - df1['RT exp']
        logger.info('Retention model calibrated successfully.')
    except Exception:
        logger.warning('Retention times are probably missing in input file.')
        df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, achrom.RCs_krokhin_100A_tfa))
        df1['RT diff'] = df1['RT exp']
    return df1, decoy2set


_standard_features = {'calc_neutral_pep_mass', 'bscore', 'yscore',
                'massdiff', 'massdiff_ppm', 'nextscore', 'RT pred', 'RT diff',
                'sumI', 'RT exp', 'precursor_neutral_mass', 'massdiff_int',
                'num_missed_cleavages', 'tot_num_ions', 'num_matched_ions', 'length',
                'ScoreRatio', 'Energy', 'MS2IonCurrent', 'MeanErrorTop7', 'sqMeanErrorTop7', 'StdevErrorTop7',
                'MS-GF:DeNovoScore', 'MS-GF:EValue', 'MS-GF:RawScore', 'MeanErrorAll',
                'MeanRelErrorAll', 'MeanRelErrorTop7', 'NumMatchedMainIons', 'StdevErrorAll',
                'StdevErrorTop7', 'StdevRelErrorAll', 'StdevRelErrorTop7', 'NTermIonCurrentRatio',
                'CTermIonCurrentRatio', 'ExplainedIonCurrentRatio', 'fragmentMT', 'ISOWIDTHDIFF',
                'MS1Intensity', 'sumI_to_MS1Intensity', 'nextscore_std', 'IPGF', 'IPGF2', 'hyperscore', 'PIF'}


def get_features(dataframe):
    feature_columns = dataframe.columns
    columns_to_remove = []
    for feature in feature_columns:
        if feature not in _standard_features:
            if not feature.startswith('mass shift') and not feature.startswith('matched_y') and not feature.startswith('matched_b'):
                columns_to_remove.append(feature)
    feature_columns = feature_columns.drop(columns_to_remove)
    return sorted(feature_columns)


def get_X_array(df, feature_columns):
    return df.loc[:, feature_columns].values


def get_Y_array(df):
    return df.loc[:, 'decoy1'].values.astype(float)


def variant_peptides(allowed_peptides, group_prefix, group_infix):
    if allowed_peptides:
        with open(allowed_peptides) as f:
            allowed_peptides = set(pseq.strip().split()[0] for pseq in f)
    else:
        allowed_peptides = None

    return allowed_peptides, group_prefix, group_infix


def filename(outfolder, outbasename, ftype):
    type_suffix = {
        'psm_full': '_PSMs_full.tsv',
        'psm': '_PSMs.tsv',
        'pepxml': '.scavager.pep.xml',
        'peptide': '_peptides.tsv',
        'protein': '_proteins.tsv',
        'protein_group': '_protein_groups.tsv'
    }
    return os.path.join(outfolder, outbasename + type_suffix[ftype])


def get_cat_model(df, feature_columns):
    logger.info('Starting machine learning...')

    # logger.info(df.shape)
    # df['orig_spectrum'] = df['spectrum'].apply(lambda x: x.split('.')[-3])
    # df['massdiff_abs'] = df['massdiff'].abs()
    # df = df.sort_values(by='massdiff_abs')
    # df = df.drop_duplicates(subset = ['orig_spectrum', 'peptide'])
    # logger.info(df.shape)

    train, test = train_test_split(df, test_size=0.3, random_state=SEED)
    x_train = get_X_array(train, feature_columns)
    y_train = get_Y_array(train)
    x_test = get_X_array(test, feature_columns)
    y_test = get_Y_array(test)
    np.random.seed(SEED)
    # model = CatBoostClassifier(iterations=1000, learning_rate=0.05, depth=10, loss_function='Logloss', logging_level='Silent', random_seed=SEED)
    # model.fit(x_train, y_train, use_best_model=True, eval_set=(x_test, y_test))

    model = CatBoostClassifier(iterations=10000, learning_rate=0.01, depth=8, loss_function='Logloss', eval_metric='Logloss',
                               od_type='Iter', od_wait=33, random_state=SEED, logging_level='Silent')
    model.fit(x_train, y_train, use_best_model=True, eval_set=(x_test, y_test))
    best_iter = model.best_iteration_
    logger.debug('Best iteration: %d', best_iter)
    ln_rt = max(0.001, round(0.01 * best_iter / 1000, 3))
    model = CatBoostClassifier(iterations=10000, learning_rate=ln_rt, depth=8, loss_function='Logloss', eval_metric='Logloss',
                               od_type='Iter', od_wait=33, random_state=SEED, logging_level='Silent')
    model.fit(x_train, y_train, use_best_model=True, eval_set=(x_test, y_test))
    best_iter = model.best_iteration_
    logger.debug('Best iteration: %d', best_iter)
    X = get_X_array(df, feature_columns)
    y = get_Y_array(df)
    model = CatBoostClassifier(iterations=int(best_iter / 0.7), learning_rate=ln_rt, depth=8, loss_function='Logloss',
                               random_state=SEED, logging_level='Silent')
    model.fit(X, y)
    logger.debug('Feature importance:')
    for fi, fn in sorted(zip(model.feature_importances_, feature_columns), key=lambda x: x[0])[::-1]:
        logger.debug('%s: %s', fi, fn)

    logger.info('Machine learning is finished.')

    return model


def calc_PEP(df, pep_ratio=1.0, reduced=False):
    if not reduced:
        feature_columns = get_features(df)
        cat_model = get_cat_model(df, feature_columns)
        x_all = get_X_array(df, feature_columns)
        df['ML score'] = cat_model.predict_proba(x_all)[:, 1]
    else:
        df['ML score'] = df['expect']

    df0_t = df[~df['decoy']]
    df0_d = df[df['decoy']]
    df0_d = df0_d[~df0_d['decoy1']]

    binsize = min(get_fdbinsize(df0_t['ML score'].values), get_fdbinsize(df0_d['ML score'].values))
    tmp = np.concatenate([df0_t['ML score'].values, df0_d['ML score'].values])
    minv = df['ML score'].min()
    maxv = df['ML score'].max()
    lbin_s = scoreatpercentile(tmp, 1.0)
    lbin = minv
    if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
        lbin = lbin_s * 1.05
    rbin_s = scoreatpercentile(tmp, 99.0)
    rbin = maxv
    if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
        rbin = rbin_s * 1.05
    rbin += 1.5 * binsize
    logger.debug('cbins: lbin = %s, rbin = %s, binsize = %s', lbin, rbin, binsize)
    cbins = np.arange(lbin, rbin + 2 * binsize, binsize)

    H1, b1 = np.histogram(df0_d['ML score'].values, bins=cbins)
    H2, b2 = np.histogram(df0_t['ML score'].values, bins=cbins)

    H2[H2 == 0] = 1
    H1_2 = H1 * (1 + 1. / pep_ratio) / H2
    ir = IsotonicRegression(y_min=0, y_max=1.0)
    ir.fit(b1[:-1], H1_2)
    df['PEP'] = ir.predict(df['ML score'].values)

    pep_min = df['ML score'].min()
    df['log_score'] = np.log10(df['ML score'] - ((pep_min - 1e-15) if pep_min < 0 else 0))


def calc_qvals(df, ratio):
    logger.debug('Q-value calculation started...')
    df_t_1 = aux.qvalues(df[~df['decoy1']], key='ML score', is_decoy='decoy2',
        remove_decoy=False, formula=1, full_output=True, ratio=ratio, correction=1)
    df_t = aux.qvalues(df[~df['decoy1']], key='ML score', is_decoy='decoy2',
        remove_decoy=False, formula=1, full_output=True, ratio=ratio, correction=0)
    df.loc[~df['decoy1'], 'q'] = df_t_1['q']
    df.loc[~df['decoy1'], 'q_uncorrected'] = df_t['q']
    df.loc[df['decoy1'], 'q'] = None
    df.loc[df['decoy1'], 'q_uncorrected'] = None


_columns_to_output = {
    'psm_full': ['peptide', 'length', 'spectrum', 'file', 'q', 'q_uncorrected', 'ML score', 'modifications', 'assumed_charge',
         'num_missed_cleavages', 'num_tol_term', 'peptide_next_aa',
         'peptide_prev_aa', 'calc_neutral_pep_mass', 'massdiff_ppm', 'massdiff_int', 'RT exp', 'RT pred',
         'RT diff', 'protein', 'protein_descr', 'decoy', 'decoy1', 'decoy2', 'PEP',
         'MS1Intensity', 'ISOWIDTHDIFF'],
    'psm': ['peptide', 'length', 'spectrum', 'file', 'q', 'q_uncorrected','ML score', 'modifications', 'modified_peptide',
         'assumed_charge', 'num_missed_cleavages', 'num_tol_term', 'peptide_next_aa',
         'peptide_prev_aa', 'calc_neutral_pep_mass', 'massdiff_ppm', 'massdiff_int', 'RT exp', 'RT pred',
         'protein', 'protein_descr', 'decoy', 'PEP', 'MS1Intensity', 'ISOWIDTHDIFF', 'PIF'],
    'peptide': ['peptide', '#PSMs', 'length', 'spectrum', 'file', 'q', 'q_uncorrected', 'ML score', 'modifications',
         'assumed_charge', 'num_missed_cleavages', 'num_tol_term', 'peptide_next_aa',
         'peptide_prev_aa', 'calc_neutral_pep_mass', 'massdiff_ppm', 'massdiff_int', 'RT exp',
         'RT pred', 'protein', 'protein_descr', 'decoy', 'PEP', 'MS1Intensity', 'ISOWIDTHDIFF', 'PIF'],
    'protein': ['dbname', 'description', 'PSMs', 'peptides', 'NSAF', 'TOP3', 'sq', 'score', 'q', 'q_uncorrected',
         'length', 'all proteins', 'groupleader'],
    }


def get_columns_to_output(columns, out_type):
    present = set(columns)
    order = _columns_to_output[out_type]
    labels = [label for label in order if label in present]
    if out_type == 'psm_full':
        labels.extend(present.difference(order))
    if out_type != 'psm_full':
        for c in columns:
            if c[:4] == 'tag_':
                labels.append(c)
    logger.debug('Writing out %s table, q_uncorrected present: %s', out_type, 'q_uncorrected' in labels)
    return labels


def calc_psms(df, df_filt):
    peptides = Counter(df_filt['peptide'])
    df['#PSMs'] = df['peptide'].apply(lambda x: peptides.get(x, 0))
