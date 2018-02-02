from utils import prepare_dataframe_xtandem, calc_PEP
from pyteomics import auxiliary as aux

def process_file(args):
    fname = args['file']
    df1 = prepare_dataframe_xtandem(fname, decoy_prefix='DOUBLE_')
    df1 = calc_PEP(df1)
    # df1_f2 = aux.filter(df1, fdr=0.01, key='expect', is_decoy='decoy', reverse=False, remove_decoy=True, ratio=1.0)
    df1_f2 = aux.filter(df1, fdr=0.01, key='PEP', is_decoy='decoy2', reverse=False, remove_decoy=True, ratio=0.333)
    real_false = df1_f2[df1_f2['protein'].apply(lambda x: all(z.startswith('DECOY_') for z in x))].shape[0]
    real_target = df1_f2.shape[0]
    print(real_false * 100 / real_target, real_target)