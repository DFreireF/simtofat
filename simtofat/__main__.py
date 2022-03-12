import argparse
import logging as log
from datetime import datetime
from pysimtof.importdata import *
from pysimtof.creategui import *
from iqtools import *
import os


def main():
    
    scriptname = 'FAT_pySimToF' 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('filename', type = str, nargs = '+', help = 'Name of the input file.')
    parser.add_argument('-l', '--lise_file', type = str, nargs = '?', help = 'Name of the LISE file.')
    
    parser.add_argument('-hrm', '--harmonics', type = int, nargs = '+', help = 'Harmonics to simulate.')
    parser.add_argument('-b', '--brho', type = float, default = 6.90922, help = 'Brho value of the reference ion beam at ESR.')
    parser.add_argument('-g', '--gammat', type = float, default = 1.395, help = 'GammaT value of ESR.')
    parser.add_argument('-i', '--refisotope', type = str, default = '72Ge', help = 'Isotope of study.')
    parser.add_argument('-c', '--refcharge', type = int, default = 32, help = 'Charge state of the studied isotope.')
    parser.add_argument('-d', '--ndivs', type = int, default = 4, help = 'Number of divisions in the display.')
    parser.add_argument('-o', '--dops', type = int, default = 1, help = 'Display of srf data options. 0-> constant height, else->scaled.')
    parser.add_argument('-t', '--time', type = float, nargs = '?', default = 1, help = 'Data time to analyse.')
    parser.add_argument('-sk', '--skip', type = float, nargs = '?', default = 0, help = 'Start of the analysis.')
    parser.add_argument('-bin', '--binning', type = int, nargs = '?', default = 1024, help = 'Number of frecuency bins.')
    parser.add_argument('-out', '--outdir', type = str, default = '.', help = 'output directory.')

    parser.add_argument('-v', '--verbose',
                        help = 'Increase output verbosity', action = 'store_true')

    parser.add_argument('-s', '--spdf',
                        help = 'Save canvas to pdf.', action = 'store_true')

    parser.add_argument('-r', '--sroot',
                        help = 'Save canvas to root.', action = 'store_true')
    
    parser.add_argument('-re', '--read',
                        help = 'Read frecuency and power already processed.', action = 'store_true')

    args = parser.parse_args()

    print(f'Running {scriptname}')
    if args.verbose: log.basicConfig(level = log.DEBUG)
    if args.outdir: outfilepath = os.path.join(args.outdir, '')
    
    #if args.read:
    controller(args.filename[0], args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot)
        
def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.
    return [file[:-1] for file in open(master_filename).readlines()]

def read_csv(filename):
    f = np.array([])
    p = np.array([])
    with open(filename) as f:
        cont = f.readlines()[1:]
        for l in cont:
            l = l.split('|')
            f = np.append(f, float(l[0]))
            p = np.append(p, float(l[1]))
    return (np.stack((f, p), axis = 1).reshape((len(f), 2)))

def ntcap_data(filename, nframes = 4096, lframes = 2**18, csv = False, root = False, plot = False):
    iq = get_iq_object(filename)
    print('iq got')
    iq.read_complete_file()
    print('read complete')
    xx, yy, zz = iq.get_spectrogram(lframes = 2**18, nframes = nframes)
    print('got spectr')
    xa, ya, za = get_averaged_spectrogram(xx, yy, zz, every = nframes)
    import gc
    del(xx)
    del(yy)
    del(zz)
    gc.collect()
    #for name in dir():
     #   if name == 'xx' or name == 'yy' or name == 'zz':
      #      del globals()[name]
    
    print('averaging')
    if root: write_spectrum_to_root(xa[0], za[0], filename, center = iq.center, title = filename)
    elif csv: write_spectrum_to_csv(xa[0], za[0], filename = filename, center = iq.center)
    elif plot: plot_spectrum(xa[0], za[0], cen = iq.center, filename = filename, title = filename)
    else: return  (np.stack((xa[0, :], za[0, :]), axis = 1).reshape((len(xa[0, :]), 2)))

def controller(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot):
    
    '''
    We use as experimental data a root file containing a histogram with the data
    '''

    mydata = ImportData(refisotope, refcharge, brho, gammat)
    mydata._set_secondary_args(lise_file, harmonics)
    mydata._calculate_srrf() # -> moq ; srrf
    mydata._simulated_data() # -> simulated frecs
    
    exp_data = filename   # pass directly root file               
    mycanvas = CreateGUI(refisotope, mydata.nuclei_names, ndivs, dops)
    mycanvas._view(exp_data, mydata.simulated_data_dict, filename)
    
    date_time = datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    info_name = f'{outfilepath}{date_time}_b{brho}_g{gammat}'
    if spdf: mycanvas.save_pdf(info_name)
    if sroot: mycanvas.save_root(info_name)

def controller1(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot):
    
    '''
    We read the experimental data from a csv file containing frecuency and power
    '''
    
    exp_data = read_csv(filename)
    mydata = ImportData(refisotope, refcharge, brho, gammat)
    mydata._set_secondary_args(lise_file, harmonics)
    mydata._calculate_srrf() # -> moq ; srrf
    mydata._simulated_data() # -> simulated frecs
                                                   
    mycanvas = CreateGUI(refisotope, mydata.nuclei_names, ndivs, dops)
    mycanvas._view(exp_data, mydata.simulated_data_dict, filename)
    
    date_time = datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    info_name = f'{outfilepath}{date_time}_b{brho}_g{gammat}'
    if spdf: mycanvas.save_pdf(info_name)
    if sroot: mycanvas.save_root(info_name)
    
def controller2(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot):
    
    '''
    We process the NTCAP data: ex. /lustre/ap/litv-exp/2021-05-00_E143_TwoPhotonDeday_ssanjari/NTCAP/iq/IQ_2021-05-08_19-51-35/0000050.iq.tdms
    '''
    
    exp_data = ntcap_data(filename)
    mydata = ImportData(refisotope, refcharge, brho, gammat)
    mydata._set_secondary_args(lise_file, harmonics)
    mydata._calculate_srrf() # -> moq ; srrf
    mydata._simulated_data() # -> simulated frecs
                                                   
    mycanvas = CreateGUI(refisotope, mydata.nuclei_names, ndivs, dops)
    mycanvas._view(exp_data, mydata.simulated_data_dict, filename)
    
    date_time = datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    info_name = f'{outfilepath}{date_time}_b{brho}_g{gammat}'
    if spdf: mycanvas.save_pdf(info_name)
    if sroot: mycanvas.save_root(info_name)

def controller3(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot):
    
    '''
    We process the NTCAP data: ex. /lustre/ap/litv-exp/2021-05-00_E143_TwoPhotonDeday_ssanjari/NTCAP/iq/IQ_2021-05-08_19-51-35/0000050.iq.tdms
    '''
    
    mydata = ImportData(refisotope, refcharge, brho, gammat)
    mydata._set_secondary_args(lise_file, harmonics)
    mydata._calculate_srrf() # -> moq ; srrf
    mydata._simulated_data() # -> simulated frecs

    ntcap_data(filename, root = True)
    exp_data = filename + '.root'
    mycanvas = CreateGUI(refisotope, mydata.nuclei_names, ndivs, dops)
    mycanvas._view(exp_data, mydata.simulated_data_dict, filename)
    
    date_time = datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    info_name = f'{outfilepath}{date_time}_b{brho}_g{gammat}'
    if spdf: mycanvas.save_pdf(info_name)
    if sroot: mycanvas.save_root(info_name)

if __name__ == '__main__':
    main()
