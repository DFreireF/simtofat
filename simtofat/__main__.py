import argparse
import logging as log
from datetime import datetime
import os


def main():
    scriptname = 'FAT_pySimToF' 
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, nargs='+', help='Name of the input file.')
    parser.add_argument('-l', '--lise_file', type=str, help='Name of the LISE file.')
    
    parser.add_argument('-hrm', '--harmonics', type=int, nargs='+', help='Harmonics to simulate')
    parser.add_argument('-b', '--brho', type=float, default=6.90922, help='Brho value of the reference ion beam at ESR')
    parser.add_argument('-g', '--gammat', type=float, default=1.395, help='GammaT value of ESR')
    
    parser.add_argument('-i', '--refisotope', type=str, default='72Ge', help='Isotope of study')
    parser.add_argument('-c', '--refcharge', type=float, default=32, help='Charge state of the studied isotope')
    parser.add_argument('-t', '--time', type=float, nargs='?', default=1, help='Analysis time from the begining')
    parser.add_argument('-sk', '--skip', type=float, nargs='?', default=0, help='Start of the analysis')

    parser.add_argument('-v', '--verbose',
                        help='Increase output verbosity', action='store_true')
    
    parser.add_argument('-out', '--outdir', type=str, default='.',
                                                help='output directory.')

    args = parser.parse_args()

    print(f'Running {scriptname}')
    if args.verbose: log.basicConfig(level=log.DEBUG)
    if args.outdir: outfilepath = os.path.join(args.outdir, '')

    # here we go:
    log.info(f'File {args.filename} passed for processing the information of {args.refisotope}+{args.refcharge}.')
    
    if ('txt') in args.filename[0]:
        filename_list=read_masterfile(args.filename[0])
        for filename in filename_list:do_your_stuff(filename[0], args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot, args.time, args.skip)
    else:
        for file in args.filename:
            do_your_stuff(file, args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot, args.time, args.skip)
            gApplication.Run()
    
def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.
    return [file[:-1] for file in open(master_filename).readlines()]
    
def do_your_stuff(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot, time, skip):
    date_time=datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    info_name=f'{outfilepath}{date_time}_b{brho}_g{gammat}'

if __name__ == '__main__':
    main()
