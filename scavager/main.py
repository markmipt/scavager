from __future__ import division
import os.path
import logging
import ast

import pandas as pd
from pyteomics import auxiliary as aux
from catboost import CatBoostError
try:
    from pyteomics import pepxmltk
except ImportError:
    pepxmltk = None
from . import utils
from .utils_figures import plot_outfigures
logger = logging.getLogger(__name__)

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def process_files(args):
    """Run Scavager for multiple files (`args['file']` should be a list of file names)
    and possibly for their union.

    Parameters
    ----------

    args : dict
        A dictionary of parameters as produced from argparse in :py:func:`search.run`.

    Returns
    -------
    out : int
        Exit code. 0 for success, 1 for empty (union) result, negative for errors
        (first encountered error code is returned).
    """
    files = args['file']
    N = len(files)
    logger.info('%d file(s) to process.', N)
    cargs = args.copy()
    if args['union'] or args['quick_union']:
        if not args['database']:
            logger.error('--database is required for union calculation.')
            return -101
    if args['create_pepxml'] and pepxmltk is None:
        logger.error('pepxmltk is required for --create-pepxml. Please install it.')
        return -102
    if args['database']:
        decoy_prots_2 = utils.split_fasta_decoys(args['database'], args['prefix'], args['infix'])
    else:
        decoy_prots_2 = None
        logger.info('Database file not provided. Decoy randomization will be done per PSM file.')
    errors = 0
    retvalues = []
    if not args['quick_union']:
        for f in files:
            cargs['file'] = f
            retv = process_file(cargs, decoy2=decoy_prots_2)
            retvalues.append(retv)
            if -10 < retv < 0:
                logger.info('Stopping due to previous errors.')
                return retv
            if retv < 0:
                errors += 1
    else:
        logger.info('Skipping individual file processing.')
    if N == 1:
        return retv
    if errors >= N - 1:
        if N > 1:
            logger.info('Union will not be run because %s out of %s files were processed with errors.', errors, N)
        return retvalues
    if args['union'] or args['quick_union']:
        logger.info('Starting the union calculation...')
        psm_full_dfs = []
        for file in files:
            outfolder = utils.get_output_folder(args['output'], file)
            outbasename = utils.get_output_basename(file, args['name_suffix'])
            csvname = utils.filename(outfolder, outbasename, 'psm_full')
            try:
                df = pd.read_csv(csvname, sep='\t', converters={
                    key: ast.literal_eval for key in [
                        'protein', 'peptide_next_aa', 'peptide_prev_aa', 'num_tol_term',
                        'protein_descr', 'modifications', 'mods_counter']})
                df['file'] = os.path.basename(file)
                psm_full_dfs.append(df)
            except FileNotFoundError as e:
                logger.debug('Exception when reading %s: %s (%s)', file, e, e.args)
                logger.warning('File %s not found, skipping...', csvname)
        if not psm_full_dfs:
            logger.error('No PSMs found for union calculation.')
            retvalues.append(1)
            return retvalues
        all_psms = pd.concat(psm_full_dfs, sort=False)
        all_psms.reset_index(inplace=True, drop=True)
        logger.debug('Recovered PSMs for analysis: %s, of those: %s decoy1, %s decoy2, %s have q < %s',
            all_psms.shape, all_psms.decoy1.sum(), all_psms.decoy2.sum(),
            (all_psms['q'] < args['fdr'] / 100).sum(), args['fdr'] / 100)
        utils.prepare_mods(all_psms)
        q_label = 'q'
        if not (args['no_correction'] or args['force_correction']):
            logger.info('Using the corrected q-values for union.')
            all_psms_f2 = all_psms[(~all_psms['decoy1']) & (all_psms[q_label] < args['fdr'] / 100)]
            if not all_psms_f2.shape[0]:
                q_label = 'q_uncorrected'
                logger.info('No union results with correction. Disabling...')
        if args['no_correction']:
            q_label = 'q_uncorrected'
        logger.debug('Filtering union PSMs by %s.', q_label)
        all_psms_f2 = all_psms.loc[(~all_psms['decoy1']) & (all_psms[q_label] < args['fdr'] / 100)].copy()

        peptides, peptides_f, proteins, proteins_f, protein_groups = build_output_tables(all_psms,
            all_psms_f2, decoy_prots_2, args, 'PEP', calc_qvals=False)
        if peptides is None:
            logger.warning('No peptides identified in union.')
            retvalues.append(1)
            return retvalues

        logger.debug('Protein FDR in full table: %f%%', 100 * aux.fdr(proteins, is_decoy='decoy2'))

        write_tables(outfolder, 'union' + args['name_suffix'] + args['union_name_suffix'],
            all_psms, all_psms_f2, peptides_f, proteins_f, protein_groups)
        if args['create_pepxml']:
            pepxmltk.easy_write_pepxml(files, utils.filename(outfolder, 'union', 'pepxml'),
            set(all_psms_f2.loc[~all_psms_f2['decoy2'], 'spectrum']))

        if len(all_psms_f2[~all_psms_f2['decoy2']]) >= 3:
            plot_outfigures(all_psms, all_psms_f2[~all_psms_f2['decoy2']], peptides,
                peptides_f[~peptides_f['decoy2']],
                outfolder, 'union' + args['name_suffix'] + args['union_name_suffix'], df_proteins=proteins,
                df_proteins_f=proteins_f[~proteins_f['decoy2']],
                separate_figures=args['separate_figures'])

        logger.info('Union calculation complete.')
        retvalues.append(0)
    return retvalues


def filter_dataframe(df1, outfdr, correction, allowed_peptides, group_prefix, group_infix, group_regex, decoy_prefix, decoy_infix):
    d2 = df1['decoy2'].sum()
    d = df1['decoy'].sum()
    pep_ratio = d2 / d
    logger.debug('Peptide ratio for ML: %s (%s / %s)', pep_ratio, d2, d)
    df1_f = utils.filter_custom(df1[~df1['decoy1']], fdr=outfdr, key='expect', is_decoy='decoy2',
        reverse=False, remove_decoy=False, ratio=pep_ratio, formula=1, correction=correction, loglabel='PSMs default')
    num_psms_def = df1_f[~df1_f['decoy2']].shape[0]
    if not num_psms_def:
        logger.warning('No valid PSMs at all are found at %s%% FDR. Aborting the analysis.', 100 * outfdr)
        return df1, None
    if num_psms_def < 100:
        logger.warning('Not enough statistics for ML training (%d PSMs is less than 100).', num_psms_def)
        utils.calc_PEP(df1, pep_ratio=pep_ratio, reduced=True)
    else:
        utils.calc_PEP(df1, pep_ratio=pep_ratio)
        df1_f2 = utils.filter_custom(df1[~df1['decoy1']], fdr=outfdr, key='ML score',
            is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, formula=1, correction=correction, loglabel='PSMs')
        if df1_f2[~df1_f2['decoy2']].shape[0] < num_psms_def:
            logger.warning('Machine learning works worse than default filtering: %d vs %d PSMs.', df1_f2.shape[0], num_psms_def)
            logger.warning('Using only default search scores for machine learning...')
            utils.calc_PEP(df1, pep_ratio=pep_ratio, reduced=True)

    if allowed_peptides or group_prefix or group_infix or group_regex:
        prev_num = df1.shape[0]
        if allowed_peptides:
            df1 = df1[df1['peptide'].apply(lambda x: x in allowed_peptides)]
        else:
            logger.debug('Protein column looks like this: %s', df1['protein'].iloc[0])
            df1 = df1[df1['protein'].apply(utils.is_group_specific,
                group_prefix=group_prefix, group_infix=group_infix, group_regex=group_regex,
                decoy_prefix=decoy_prefix, decoy_infix=decoy_infix)]

        logger.info('%.1f%% of identifications were dropped during group-specific filtering.',
            (100 * float(prev_num - df1.shape[0]) / prev_num))

        if df1[df1['decoy']].shape[0] == 0:
            logger.warning('0 decoy identifications are present in the group. Please check'
            'that allowed_peptides contains decoy peptides or that decoy proteins have group_prefix/infix!')

    pep_ratio = df1['decoy2'].sum() / df1['decoy'].sum()
    logger.debug('Peptide ratio within group: %s', pep_ratio)
    df1_f2 = utils.filter_custom(df1[~df1['decoy1']], fdr=outfdr, key='ML score', is_decoy='decoy2',
        reverse=False, remove_decoy=False, ratio=pep_ratio, formula=1, correction=correction, loglabel='PSMs in group')

    return df1, df1_f2


def build_output_tables(df1, df1_f2, decoy2, args, key='ML score', calc_qvals=True):
    if args['database']:
        path_to_fasta = os.path.abspath(args['database'])
    else:
        path_to_fasta = None
    outfdr = args['fdr'] / 100
    correction = False if args['no_correction'] else (args['force_correction'] or None)
    pep_ratio = df1['decoy2'].sum() / df1['decoy'].sum()
    logger.debug('Peptide ratio for q-value calculation: %s', pep_ratio)

    if calc_qvals:
        utils.calc_qvals(df1, ratio=pep_ratio)

    if df1_f2.shape[0]:
        utils.calc_psms(df1, df1_f2)

        df_proteins, norm = utils.get_proteins_dataframe(df1_f2, decoy_prefix=args['prefix'],
            decoy_infix=args['infix'], all_decoys_2=decoy2, path_to_fasta=path_to_fasta, pif_threshold=args['pif_threshold'])
        tagnames = utils.get_tag_names(df1_f2.columns)

        logger.debug('Normalizing tag intensities by: %s', norm)
        logger.debug('Channels: %s', tagnames)

        nonnormalized = ['raw_' + n for n in tagnames]
        df1[nonnormalized] = df1[tagnames]
        df1[tagnames] /= norm
        df1_f2[tagnames] /= norm
        df1_peptides = df1.sort_values(key, ascending=True).drop_duplicates(['peptide'])
        df1_peptides_f = utils.filter_custom(df1_peptides[~df1_peptides['decoy1']], fdr=outfdr,
            key=key, is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, formula=1,
            correction=correction, loglabel='peptides')

        prot_ratio = 0.5
        df_proteins = df_proteins[~df_proteins['decoy1']]
        df_proteins_f = utils.filter_custom(df_proteins, fdr=outfdr, key='score', is_decoy='decoy2',
            reverse=False, remove_decoy=True, ratio=prot_ratio, formula=1, correction=correction, loglabel='proteins')
        utils.add_protein_groups(df_proteins_f, args['ms1'])
        df_protein_groups = df_proteins_f[df_proteins_f['groupleader']]

        logger.info('Final results at %s%% FDR level:', args['fdr'])
        logger.info('Identified PSMs: %s', df1_f2[~df1_f2['decoy2']].shape[0])
        logger.info('Identified peptides: %s', df1_peptides_f[~df1_peptides_f['decoy2']].shape[0])
        logger.info('Identified proteins: %s', df_proteins_f.shape[0])
        logger.info('Identified protein groups: %s', df_protein_groups.shape[0])
        logger.info('Processing finished.')

        return df1_peptides, df1_peptides_f, df_proteins, df_proteins_f, df_protein_groups
    else:
        logger.warning('PSMs cannot be filtered at %s%% FDR. Please increase allowed FDR.', args['fdr'])
        return (None,) * 5


def write_tables(outfolder, outbasename, df1, df1_f2, df1_peptides_f, df_proteins_f, df_protein_groups):
    df1.to_csv(utils.filename(outfolder, outbasename, 'psm_full'),
        sep='\t', index=False, columns=utils.get_columns_to_output(df1.columns, 'psm_full'))
    df1_f2[~df1_f2['decoy2']].to_csv(utils.filename(outfolder, outbasename, 'psm'),
        sep='\t', index=False, columns=utils.get_columns_to_output(df1_f2.columns, 'psm'))
    df1_peptides_f[~df1_peptides_f['decoy2']].to_csv(utils.filename(outfolder, outbasename, 'peptide'),
        sep='\t', index=False, columns=utils.get_columns_to_output(df1_peptides_f.columns, 'peptide'))

    df_proteins_f.to_csv(utils.filename(outfolder, outbasename, 'protein'), sep='\t', index=False,
            columns=utils.get_columns_to_output(df_proteins_f.columns, 'protein'))
    df_protein_groups.to_csv(utils.filename(outfolder, outbasename, 'protein_group'), sep='\t', index=False,
            columns=utils.get_columns_to_output(df_protein_groups.columns, 'protein'))
    logger.info('Output tables saved.')


def process_file(args, decoy2=None):
    """Run Scavager for a single file (`args['file']` should be a single file name).

    Parameters
    ----------

    args : dict
        A dictionary of parameters as produced from argparse in :py:func:`search.run`.
    decoy2 : set, optional
        A set of proteins labeled as "decoy2"

    Returns
    -------
    out : int
        Exit code. 0 for success, 1 for empty result, negative for errors.
    """
    fname = args['file']
    outfolder = utils.get_output_folder(args['output'], fname)
    outbasename = utils.get_output_basename(fname, args['name_suffix'])
    outfdr = args['fdr'] / 100
    correction = False if args['no_correction'] else (args['force_correction'] or None)
    logger.debug('Operating with correction = %s', correction)
    decoy_prefix = args['prefix']
    decoy_infix = args['infix']
    sf = args['separate_figures']
    logger.info('Loading file %s...', os.path.basename(fname))
    if args['enzyme']:
        cleavage_rule = utils.convert_tandem_cleave_rule_to_regexp(args['enzyme'])
    else:
        cleavage_rule = None

    try:
        df1, all_decoys_2 = utils.prepare_dataframe(fname, decoy_prefix=decoy_prefix,
            decoy_infix=decoy_infix, cleavage_rule=cleavage_rule, fdr=outfdr, decoy2set=decoy2)
    except utils.NoDecoyError:
        logger.error('No decoys were found. Please check decoy_prefix/infix parameter or your search output.')
        return -12
    except utils.WrongInputError:
        logger.error('Unsupported input file format. Use .pep.xml or .mzid files.')
        return -2
    except utils.EmptyFileError:
        logger.error('Input file %s is empty.', fname)
        return 1

    try:
        allowed_peptides, group_prefix, group_infix, group_regex = utils.variant_peptides(
            args['allowed_peptides'], args['group_prefix'], args['group_infix'], args['group_regex'])
    except ValueError as e:
        logger.error(e.args[0])
        return -3

    try:
        df1, df1_f2 = filter_dataframe(df1, outfdr, correction,
            allowed_peptides, group_prefix, group_infix, group_regex, decoy_prefix, decoy_infix)
    except CatBoostError as e:
        logger.error('There was an error in Catboost: %s', e.args)
        return -11

    if df1_f2 is None:
        return 1
    df1_peptides, df1_peptides_f, df_proteins, df_proteins_f, df_protein_groups = build_output_tables(
        df1, df1_f2, all_decoys_2, args)
    if df1_peptides is None:
        return 1

    write_tables(outfolder, outbasename, df1, df1_f2, df1_peptides_f, df_proteins_f, df_protein_groups)
    if args['create_pepxml']:
        pepxmltk.easy_write_pepxml([args['file']], utils.filename(outfolder, outbasename, 'pepxml'),
            set(df1_f2.loc[~df1_f2['decoy2'], 'spectrum']))

    if df1_f2[~df1_f2['decoy2']].shape[0] >= 3:
        plot_outfigures(df1, df1_f2[~df1_f2['decoy2']], df1_peptides, df1_peptides_f[~df1_peptides_f['decoy2']],
            outfolder, outbasename, df_proteins=df_proteins,
            df_proteins_f=df_proteins_f[~df_proteins_f['decoy2']], separate_figures=sf)
    return 0
