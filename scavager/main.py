from __future__ import division
from . import utils
from .utils_figures import plot_outfigures
from pyteomics import auxiliary as aux
import os.path
import logging
import ast
import pandas as pd

def process_files(args):
    files = args['file']
    logging.info('%d file(s) to process.', len(files))
    cargs = args.copy()
    if args['union']:
        if not args['database']:
            logging.error('--database is required with --union.')
            return

    if args['database']:
        decoy_prots_2 = utils.split_fasta_decoys(args['database'], args['prefix'], args['infix'])
    else:
        decoy_prots_2 = None
        logging.info('Database file not provided. Decoy randomization will be done per PSM file.')
    for f in files:
        cargs['file'] = f
        if process_file(cargs, decoy2=decoy_prots_2) == -1:
            logging.info('Stopping due to previous errors.')
            return
    if args['union'] and len(files) > 1:
        logging.info('Starting the union calculation...')
        psm_full_dfs = []
        for file in files:
            outfolder = utils.get_output_folder(args['output'], file)
            outbasename = utils.get_output_basename(file)
            csvname = utils.filename(outfolder, outbasename, 'psm_full')
            try:
                df = pd.read_csv(csvname, sep='\t')
                for key in ['protein', 'peptide_next_aa', 'peptide_prev_aa', 'num_tol_term', 'protein_descr', 'modifications']:
                    df[key] = df[key].apply(ast.literal_eval)
                psm_full_dfs.append(df)
            except FileNotFoundError:
                logging.warning('File %s not found, skipping...', csvname)
        all_psms = pd.concat(psm_full_dfs)
        all_psms_f2 = all_psms[(~all_psms['decoy1']) & (all_psms['q'] < args['fdr'] / 100)]

        peptides, peptides_f, proteins, proteins_f, protein_groups = build_output_tables(all_psms, all_psms_f2, decoy_prots_2, args, 'PEP', calc_qvals=False)
        if peptides is None:
            logging.warning('No peptides identified in union.')
            return 0

        logging.debug('Protein FDR in full table: %f%%', 100*aux.fdr(proteins, is_decoy='decoy2'))

        write_tables(outfolder, 'union', all_psms, all_psms_f2, peptides_f, proteins_f, protein_groups)

        plot_outfigures(all_psms, all_psms_f2[~all_psms_f2['decoy2']], peptides, peptides_f[~peptides_f['decoy2']],
            outfolder, 'union', df_proteins=proteins, df_proteins_f=proteins_f[~proteins_f['decoy2']],
            separate_figures=args['separate_figures'])

        logging.info('Union calculation complete.')


def filter_dataframe(df1, outfdr, num_psms_def, allowed_peptides, group_prefix, decoy_prefix, decoy_infix):
    pep_ratio = df1['decoy2'].sum() / df1['decoy'].sum()
    utils.calc_PEP(df1, pep_ratio=pep_ratio)
    if allowed_peptides or group_prefix:
        prev_num = df1.shape[0]
        if allowed_peptides:
            df1 = df1[df1['peptide'].apply(lambda x: x in allowed_peptides)]
        elif group_prefix:
            df1 = df1[df1['protein'].apply(utils.is_group_specific, group_prefix=group_prefix, decoy_prefix=decoy_prefix, decoy_infix=decoy_infix)]

        logging.info('%.1f%% of identifications were dropped during group-specific filtering.', (100 * float(prev_num - df1.shape[0]) / prev_num))

        if df1[df1['decoy']].shape[0] == 0:
            logging.warning('0 decoy identifications are present in the group.\n Please check\
            that allowed_peptides or group_prefix contains also decoy peptides and proteins!')

        df1_f = utils.filter_custom(df1[~df1['decoy1']], fdr=outfdr, key='expect', is_decoy='decoy2', reverse=False,
        remove_decoy=False, ratio=pep_ratio, formula=1)
        num_psms_def = df1_f[~df1_f['decoy2']].shape[0]

    df1_f2 = utils.filter_custom(df1[~df1['decoy1']], fdr=outfdr, key='ML score', is_decoy='decoy2', reverse=False,
        remove_decoy=False, ratio=pep_ratio, formula=1)

    if df1_f2[~df1_f2['decoy2']].shape[0] < num_psms_def:
        logging.warning('Machine learning works worse than default filtering: %d vs %d PSMs.', df1_f2.shape[0], num_psms_def)
        logging.warning('Using only default search scores for machine learning...')
        utils.calc_PEP(df1, pep_ratio=pep_ratio, reduced=True)
        df1_f2 = utils.filter_custom(df1[~df1['decoy1']], fdr=outfdr, key='ML score', is_decoy='decoy2', reverse=False,
        remove_decoy=False, ratio=pep_ratio, formula=1)

    return df1, df1_f2


def build_output_tables(df1, df1_f2, decoy2, args, key='ML score', calc_qvals=True):
    if args['database']:
        path_to_fasta = os.path.abspath(args['database'])
    else:
        path_to_fasta = None
    outfdr = args['fdr'] / 100
    pep_ratio = df1['decoy2'].sum() / df1['decoy'].sum()

    if calc_qvals:
        utils.calc_qvals(df1, ratio=pep_ratio)

    if df1_f2.shape[0]:
        utils.calc_psms(df1)
        df1_peptides = df1.sort_values(key, ascending=True).drop_duplicates(['peptide'])
        df1_peptides_f = utils.filter_custom(df1_peptides[~df1_peptides['decoy1']], fdr=outfdr,
            key=key, is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, formula=1)

        df_proteins = utils.get_proteins_dataframe(df1_f2, decoy_prefix=args['prefix'],
            decoy_infix=args['infix'], all_decoys_2=decoy2, path_to_fasta=path_to_fasta)
        prot_ratio = 0.5
        df_proteins = df_proteins[~df_proteins['decoy1']]
        df_proteins_f = utils.filter_custom(df_proteins, fdr=outfdr, key='score', is_decoy='decoy2',
            reverse=False, remove_decoy=True, ratio=prot_ratio, formula=1)
        utils.add_protein_groups(df_proteins_f)
        df_protein_groups = df_proteins_f[df_proteins_f['groupleader']]

        logging.info('Final results at %s%% FDR level:', args['fdr'])
        logging.info('Identified PSMs: %s', df1_f2[~df1_f2['decoy2']].shape[0])
        logging.info('Identified peptides: %s', df1_peptides_f[~df1_peptides_f['decoy2']].shape[0])
        logging.info('Identified proteins: %s', df_proteins_f.shape[0])
        logging.info('Identified protein groups: %s', df_protein_groups.shape[0])
        logging.info('The search is finished.')

        return df1_peptides, df1_peptides_f, df_proteins, df_proteins_f, df_protein_groups
    else:
        logging.error('PSMs cannot be filtered at %s%% FDR. Please increase allowed FDR.', args['fdr'])
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


def process_file(args, decoy2=None):
    fname = args['file']
    outfolder = utils.get_output_folder(args['output'], fname)
    outbasename = utils.get_output_basename(fname)
    outfdr = args['fdr'] / 100
    decoy_prefix = args['prefix']
    decoy_infix = args['infix']
    sf = args['separate_figures']
    logging.info('Loading file %s...', os.path.basename(fname))
    if args['enzyme']:
        cleavage_rule = utils.convert_tandem_cleave_rule_to_regexp(args['enzyme'])
    else:
        cleavage_rule = None

    try:
        df1, all_decoys_2, num_psms_def = utils.prepare_dataframe(fname, decoy_prefix=decoy_prefix,
            decoy_infix=decoy_infix, cleavage_rule=cleavage_rule, fdr=outfdr, decoy2set=decoy2)
    except utils.NoDecoyError:
        logging.error('No decoys were found. Please check decoy_prefix/infix parameter or your search output.')
        return -1
    except utils.WrongInputError:
        logging.error('Unsupported input file format. Use .pep.xml or .mzid files.')
        return -1

    try:
        allowed_peptides, group_prefix = utils.variant_peptides(args['allowed_peptides'], args['group_prefix'])
    except ValueError as e:
        logging.error(e.args[0])
        return -1

    df1, df1_f2 = filter_dataframe(df1, outfdr, num_psms_def, allowed_peptides, group_prefix, decoy_prefix, decoy_infix)

    df1_peptides, df1_peptides_f, df_proteins, df_proteins_f, df_protein_groups = build_output_tables(df1, df1_f2, all_decoys_2, args)
    if df1_peptides is None:
        return 0

    write_tables(outfolder, outbasename, df1, df1_f2, df1_peptides_f, df_proteins_f, df_protein_groups)

    plot_outfigures(df1, df1_f2[~df1_f2['decoy2']], df1_peptides, df1_peptides_f[~df1_peptides_f['decoy2']],
            outfolder, outbasename, df_proteins=df_proteins, df_proteins_f=df_proteins_f[~df_proteins_f['decoy2']],
            separate_figures=sf)
