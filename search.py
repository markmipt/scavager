import argparse
import main

def run():
    parser = argparse.ArgumentParser(
        description='Search proteins using LC-MS spectra',
        epilog='''

    Example usage
    -------------
    $ search.py input.pep.xml
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', help='input pepXML file')
    parser.add_argument('-o', help='path to output folder', default=False)
    parser.add_argument('-db', help='path to fasta file. \
                        Used for sequence coverage and LFQ calculation', default=False)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-fdr', help='false discovery rate in %%', default=1.0, type=float)
    args = vars(parser.parse_args())

    main.process_file(args)
    print('The search is finished.')



if __name__ == '__main__':
    run()