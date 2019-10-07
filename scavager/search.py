from __future__ import division
import argparse
import logging
import pkg_resources
from . import main

def run():
    parser = argparse.ArgumentParser(
        description='postsearch analysis of peptides and proteins',
        epilog='''

    Example usage
    -------------
    $ scavager input.pep.xml -prefix DECOY_ -fdr 1.0
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', nargs='+', help='input pepXML file')
    decoy = parser.add_mutually_exclusive_group(required=False)
    decoy.add_argument('-p', '--prefix', help='decoy prefix', default='DECOY_')
    decoy.add_argument('-i', '--infix', help='decoy infix if database was generated by SearchGUI')
    parser.add_argument('-o', '--output', help='path to output folder')
    parser.add_argument('-x', '--create-pepxml', action='store_true',
        help='Create a pepXML file with validated PSMs (requires pepxmltk)')
    parser.add_argument('-db', '--database', help='path to fasta file. \
                        Used for sequence coverage and LFQ calculation')
    parser.add_argument('-fdr', '--fdr', help='false discovery rate in %%', default=1.0, type=float)
    parser.add_argument('-e', '--enzyme', help='Used only for msgf+ and Morpheus search engines.\
    Cleavage rule in quotes! X!Tandem style for cleavage rules. Examples:\
    "[RK]|{P}" means cleave after R and K, but not before P;\
    "[X]|[D]" means cleave before D;\
    "[RK]|{P},[M]|[X]" means mix of trypsin and cnbr', default='[RK]|{P}')
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-ap', '--allowed-peptides', help='Path to file with peptides considered '
        'in postsearch analysis. Sequences must be separated by new line. '
        'For example, it can be variant peptides and their decoys in case of proteogenomics searches'
        ' for group-specific FDR calculation')
    group.add_argument('-gp', '--group-prefix', help='Protein prefix for group-specific filtering. '
        'For example, if `mut_` prefix, peptides from mut_ or DECOY_mut_ proteins will be reported. '
        'This can useful in proteogenomic searches for group-specific FDR calculation')
    group.add_argument('-gi', '--group-infix', help='Protein infix for group-specific filtering.')
    parser.add_argument('-sf', '--separate-figures', action='store_true',
        help='save figures as separate files')
    parser.add_argument('-u', '--union', action='store_true',
        help='Produce a summary table where IDs are pooled from all files (requires -db)')
    correction = parser.add_mutually_exclusive_group(required=False)
    correction.add_argument('-c', '--force-correction', action='store_true',
        help='Force the use of "+1" correction when calculating q-values, even if it results in empty output.')
    correction.add_argument('-nc', '--no-correction', action='store_true',
        help='Disable the use of "+1" correction when calculating q-values, even if it results in highly inaccurate q-values.')
    parser.add_argument('--quick-union', action='store_true',
        help='Assume that individual files have been already processed and go straight to union calculation.')
    parser.add_argument('--debug', action='store_true', help='Enable debugging output')
    parser.add_argument('-v', '--version', action='version',
        version='%s' % (pkg_resources.require("scavager")[0], ))
    args = vars(parser.parse_args())
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=[logging.INFO, logging.DEBUG][args['debug']])
    logger = logging.getLogger(__name__)
    logger.debug('Starting with args: %s', args)
    return main.process_files(args)


if __name__ == '__main__':
    run()
