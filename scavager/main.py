from __future__ import division
from . import utils
from .utils_figures import plot_outfigures
from pyteomics import auxiliary as aux
import os.path
import logging

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
        process_file(cargs, decoy2=decoy_prots_2)
    if args['union']:
        logging.info('Starting the union calculation...')
        logging.warning('Union calculation is not implemented yet.')

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
        return
    except utils.WrongInputError:
        logging.error('Unsupported input file format. Use .pep.xml or .mzid files.')
        return

    if args['allowed_peptides']:
        allowed_peptides = set([pseq.strip().split()[0] for pseq in open(args['allowed_peptides'], 'r')])
    else:
        allowed_peptides = None
    if args['group_prefix']:
        group_prefix = args['group_prefix']
    else:
        group_prefix = None
    if group_prefix and allowed_peptides:
        logging.error('Only one type of group filter can be used: --allowed-peptides or --group-prefix.')
        return

    pep_ratio = df1['decoy2'].sum() / df1['decoy'].sum()
    df1 = utils.calc_PEP(df1, pep_ratio=pep_ratio)
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

        df1_f = aux.filter(df1[~df1['decoy1']], fdr=outfdr, key='expect', is_decoy='decoy2', reverse=False,
        remove_decoy=False, ratio=pep_ratio, correction=1, formula=1)
        num_psms_def = df1_f[~df1_f['decoy2']].shape[0]
        if num_psms_def == 0:
            df1_f = aux.filter(df1[~df1['decoy1']], fdr=outfdr, key='expect', is_decoy='decoy2', reverse=False,
            remove_decoy=False, ratio=pep_ratio, correction=0, formula=1)
            num_psms_def = df1_f[~df1_f['decoy2']].shape[0]

    df1_f2 = aux.filter(df1[~df1['decoy1']], fdr=outfdr, key='ML score', is_decoy='decoy2',
        reverse=False, remove_decoy=False, ratio=pep_ratio, correction=1, formula=1)
    if df1_f2.shape[0] == 0:
        df1_f2 = aux.filter(df1[~df1['decoy1']], fdr=outfdr, key='ML score', is_decoy='decoy2',
            reverse=False, remove_decoy=False, ratio=pep_ratio, correction=0, formula=1)

    if df1_f2[~df1_f2['decoy2']].shape[0] < num_psms_def:
        logging.warning('Machine learning works worse than default filtering: %d vs %d PSMs.', df1_f2.shape[0], num_psms_def)
        logging.warning('Using only default search scores for machine learning...')
        df1 = utils.calc_PEP(df1, pep_ratio=pep_ratio, reduced=True)
        df1_f2 = aux.filter(df1[~df1['decoy1']], fdr=outfdr, key='ML score', is_decoy='decoy2',
            reverse=False, remove_decoy=False, ratio=pep_ratio, correction=1, formula=1)


    output_path_psms_full = os.path.join(outfolder, outbasename + '_PSMs_full.tsv')
    df1 = utils.calc_qvals(df1, ratio=pep_ratio)
    df1.to_csv(output_path_psms_full, sep='\t', index=False, columns=utils.get_columns_to_output(out_type='psm_full'))
    if df1_f2.shape[0] > 0:
        output_path_psms = os.path.join(outfolder, outbasename + '_PSMs.tsv')
        df1_f2[~df1_f2['decoy2']].to_csv(output_path_psms, sep='\t', index=False, columns=utils.get_columns_to_output(out_type='psm'))

        df1 = utils.calc_psms(df1)
        df1_peptides = df1.sort_values('ML score', ascending=True).drop_duplicates(['peptide'])
        df1_peptides_f = aux.filter(df1_peptides[~df1_peptides['decoy1']], fdr=outfdr,
            key='ML score', is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, correction=1, formula=1)
        if df1_peptides_f.shape[0] == 0:
            df1_peptides_f = aux.filter(df1_peptides[~df1_peptides['decoy1']], fdr=outfdr,
                key='ML score', is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, correction=0, formula=1)
        output_path_peptides = os.path.join(outfolder, outbasename + '_peptides.tsv')
        df1_peptides_f[~df1_peptides_f['decoy2']].to_csv(output_path_peptides, sep='\t', index=False,
            columns=utils.get_columns_to_output(out_type='peptide'))

        if args['database']:
            path_to_fasta = os.path.abspath(args['database'])
        else:
            path_to_fasta = None
        df_proteins = utils.get_proteins_dataframe(df1_f2, df1_peptides_f, decoy_prefix=args['prefix'],
            decoy_infix=args['infix'], all_decoys_2=all_decoys_2, path_to_fasta=path_to_fasta)
        prot_ratio = 0.5
        df_proteins = df_proteins[df_proteins.apply(lambda x: not x['decoy'] or x['decoy2'], axis=1)]
        df_proteins_f = aux.filter(df_proteins, fdr=outfdr, key='score', is_decoy='decoy2',
            reverse=False, remove_decoy=True, ratio=prot_ratio, formula=1, correction=1)
        if df_proteins_f.shape[0] == 0:
            df_proteins_f = aux.filter(df_proteins, fdr=outfdr, key='score', is_decoy='decoy2',
                reverse=False, remove_decoy=True, ratio=prot_ratio, formula=1, correction=0)
        df_proteins_f = utils.get_protein_groups(df_proteins_f)
        output_path_proteins = os.path.join(outfolder, outbasename + '_proteins.tsv')
        df_proteins_f.to_csv(output_path_proteins, sep='\t', index=False,
            columns=utils.get_columns_to_output(out_type='protein'))

        df_protein_groups = df_proteins_f[df_proteins_f['groupleader']]
        output_path_protein_groups = os.path.join(outfolder, outbasename + '_protein_groups.tsv')
        df_protein_groups.to_csv(output_path_protein_groups, sep='\t', index=False,
            columns=utils.get_columns_to_output(out_type='protein'))

        plot_outfigures(df1, df1_f2[~df1_f2['decoy2']], df1_peptides, df1_peptides_f[~df1_peptides_f['decoy2']],
            outfolder, outbasename, df_proteins=df_proteins, df_proteins_f=df_proteins_f[~df_proteins_f['decoy2']],
            separate_figures=sf)

        logging.info('Final results at %s%% FDR level:', args['fdr'])
        logging.info('Identified PSMs: %s', df1_f2[~df1_f2['decoy2']].shape[0])
        logging.info('Identified peptides: %s', df1_peptides_f[~df1_peptides_f['decoy2']].shape[0])
        logging.info('Identified proteins: %s', df_proteins_f.shape[0])
        logging.info('Identified protein groups: %s', df_protein_groups.shape[0])
        logging.info('The search is finished.')

    else:
        logging.error('PSMs cannot be filtered at %s%% FDR. Please increase allowed FDR.', args['fdr'])
