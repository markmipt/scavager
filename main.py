from utils import prepare_dataframe_xtandem, calc_PEP, get_output_basename, get_output_folder
from pyteomics import auxiliary as aux
from os import path

def process_file(args):
    fname = args['file']
    outfolder = get_output_folder(args['o'], fname)
    outbasename = get_output_basename(fname)
    print('Loading file %s...' % (path.basename(fname), ))
    df1 = prepare_dataframe_xtandem(fname, decoy_prefix='DOUBLE_')
    df1 = calc_PEP(df1)
    df1 = df1[~df1['decoy1']]
    df1_f2 = aux.filter(df1, fdr=0.01, key='PEP', is_decoy='decoy2', reverse=False, remove_decoy=True, ratio=0.5, correction=1)
    real_false = df1_f2[df1_f2['protein'].apply(lambda x: all(z.startswith('DECOY_') for z in x))].shape[0]
    real_target = df1_f2.shape[0]
    print(real_false * 100 / real_target, real_target)

    output_path_psms_full = path.join(outfolder, outbasename + '_PSMs_full.tsv')
    df1.to_csv(output_path_psms_full, sep='\t', index=False)

    output_path_psms = path.join(outfolder, outbasename + '_PSMs.tsv')
    df1_f2.to_csv(output_path_psms, sep='\t', index=False)