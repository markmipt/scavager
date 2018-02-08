from utils import prepare_dataframe_xtandem, calc_PEP, get_output_basename, get_output_folder, get_proteins_dataframe, get_protein_groups
from utils_figures import plot_outfigures
from pyteomics import auxiliary as aux
from os import path
from collections import defaultdict
import pandas as pd
import numpy as np

def process_file(args):
    fname = args['file']
    outfolder = get_output_folder(args['o'], fname)
    outbasename = get_output_basename(fname)
    outfdr = args['fdr'] / 100
    print('Loading file %s...' % (path.basename(fname), ))
    df1, all_decoys_2 = prepare_dataframe_xtandem(fname, decoy_prefix=args['prefix'])
    df1 = calc_PEP(df1)
    pep_ratio = np.sum(df1['decoy2'])/np.sum(df1['decoy'])

    output_path_psms_full = path.join(outfolder, outbasename + '_PSMs_full.tsv')
    df1.to_csv(output_path_psms_full, sep='\t', index=False)

    df1_f2 = aux.filter(df1[~df1['decoy1']], fdr=outfdr, key='PEP', is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, correction=1)

    output_path_psms = path.join(outfolder, outbasename + '_PSMs.tsv')
    df1_f2[~df1_f2['decoy2']].to_csv(output_path_psms, sep='\t', index=False)

    df1_peptides = df1.sort_values('PEP', ascending=True).drop_duplicates(['peptide'])
    df1_peptides_f = aux.filter(df1_peptides[~df1_peptides['decoy1']], fdr=outfdr, key='PEP', is_decoy='decoy2', reverse=False, remove_decoy=False, ratio=pep_ratio, correction=1)
    output_path_peptides = path.join(outfolder, outbasename + '_peptides.tsv')
    df1_peptides_f[~df1_peptides_f['decoy2']].to_csv(output_path_peptides, sep='\t', index=False)

    plot_outfigures(df1, df1_f2[~df1_f2['decoy2']], df1_peptides, df1_peptides_f[~df1_peptides_f['decoy2']], outfolder, outbasename)

    if args['db']:
        path_to_fasta = path.abspath(args['db'])
    else:
        path_to_fasta = args['db']
    df_proteins = get_proteins_dataframe(df1_f2, df1_peptides_f, decoy_prefix=args['prefix'], all_decoys_2=all_decoys_2, path_to_fasta=path_to_fasta)
    prot_ratio = 0.5
    df_proteins = df_proteins[df_proteins.apply(lambda x: not x['decoy'] or x['decoy2'], axis=1)]
    df_proteins = aux.filter(df_proteins, fdr=outfdr, key='score', is_decoy='decoy2', reverse=False, remove_decoy=True, ratio=prot_ratio)
    df_proteins = get_protein_groups(df_proteins)
    output_path_proteins = path.join(outfolder, outbasename + '_proteins.tsv')
    df_proteins.to_csv(output_path_proteins, sep='\t', index=False, columns = ['dbname','description','PSMs','peptides','NSAF','sq','score','length', 'all proteins', 'groupleader']) 

    df_protein_groups = df_proteins[df_proteins['groupleader']]
    output_path_protein_groups = path.join(outfolder, outbasename + '_protein_groups.tsv')
    df_protein_groups.to_csv(output_path_protein_groups, sep='\t', index=False, columns = ['dbname','description','PSMs','peptides','NSAF','sq','score','length', 'all proteins', 'groupleader']) 
