import argparse
import os
import logging as log
from datetime import datetime
from pysimtof.importdata import *
from pysimtof.creategui import *
from simtofat.model import *
from iqtools import *


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
    
    parser.add_argument('-im', '--imoq',
                        help = 'Identify moq window', action = 'store_true')
    parser.add_argument('-ir', '--iroot',
                        help = 'Identify LISE particles in root spectrum', action = 'store_true')
    parser.add_argument('-ic', '--icsv',
                        help = 'Identify LISE particles in spectrum readed from csv file', action = 'store_true')
    parser.add_argument('-in', '--intcap',
                        help = 'Identify LISE particles in a tdms file', action = 'store_true')
    parser.add_argument('-inr', '--introot',
                        help = 'Identify process tdms data, write to root, and read it and apply identification of LISE', action = 'store_true')


    args = parser.parse_args()

    print(f'Running {scriptname}')
    if args.verbose: log.basicConfig(level = log.DEBUG)
    if args.outdir: outfilepath = os.path.join(args.outdir, '')
    if args.intcap: identification_ntcap(args.filename[0], args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot, args.time, args.skip)
    if args.iroot: identification_root(args.filename[0], args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot)
    if args.icsv: identification_csv(args.filename[0], args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot)
    if args.introot: identification_ntcap_root(args.filename[0], args.lise_file, args.harmonics, args.brho, args.gammat, args.refisotope, args.refcharge, args.ndivs, args.dops, args.spdf, args.sroot)
    if args.imoq: identification_moq_window(args.moq_cen, args.moq_span, args.ref_nuclei, args.ref_charge, args.brho, args.gammat, args.lise_file, args.harmonics, args.filename[0], args.data_time, args.skip_time, args.binning)
    
    #if args.read:
    controller2
        
def read_masterfile(master_filename):
    # reads list filenames with experiment data. [:-1] to remove eol sequence.
    return [file[:-1] for file in open(master_filename).readlines()]

def ntcap_data(filename, lframes = 2**18, data_time = 1, skip_time = 2,csv = False, root = False, plot = False, read_all = False):
    iq = get_iq_object(filename)
    print('iq got')
    if not read_all:
        nframes = int(data_time * iq.fs / lframes)
        sframes = int(skip_time * iq.fs / lframes)
        iq.read(nframes = nframes, lframes = lframes, sframes = sframes)
    else:
        nframes = 4096
        iq.read_complete_file()
    print('read complete')
    xx, yy, zz = iq.get_spectrogram(nframes = nframes, lframes = lframes)
    print('spectrogram got')
    xa, ya, za = get_averaged_spectrogram(xx, yy, zz, every = nframes)
    print('average done')
    import gc
    del(xx)
    del(yy)
    del(zz)
    gc.collect()
    print('memory released')
    #for name in dir():
     #   if name == 'xx' or name == 'yy' or name == 'zz':
      #      del globals()[name]
    
    if root: write_spectrum_to_root(xa[0], za[0], filename, center = iq.center, title = filename)
    elif csv: write_spectrum_to_csv(xa[0], za[0], filename = filename, center = iq.center)
    elif plot: plot_spectrum(xa[0], za[0], cen = iq.center, filename = filename, title = filename)
    else: return  (np.stack((xa[0, :], za[0, :]), axis = 1).reshape((len(xa[0, :]), 2)))

def identification_root(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot):
    
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

def identification_ntcap(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot, data_time, skip_time):
    
    '''
    We process the NTCAP data: ex. /lustre/ap/litv-exp/2021-05-00_E143_TwoPhotonDeday_ssanjari/NTCAP/iq/IQ_2021-05-08_19-51-35/0000050.iq.tdms
    '''
    
    exp_data = ntcap_data(filename, data_time = data_time, skip_time = skip_time)
    mydata = ImportData(refisotope, refcharge, brho, gammat)
    mydata._set_secondary_args(lise_file, harmonics)
    mydata.calculate_moqs()
    mydata._calculate_srrf() # -> moq ; srrf
    mydata._simulated_data() # -> simulated frecs
                                                   
    mycanvas = CreateGUI(refisotope, mydata.nuclei_names, ndivs, dops)
    mycanvas._view(exp_data, mydata.simulated_data_dict, filename)
    
    date_time = datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    info_name = f'{outfilepath}{date_time}_b{brho}_g{gammat}'
    if spdf: mycanvas.save_pdf(info_name)
    if sroot: mycanvas.save_root(info_name)

def identification_ntcap_root(filename, lise_file, harmonics, brho, gammat, refisotope, refcharge, ndivs, dops, spdf, sroot):
    
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

def identification_moq_window(moq_cen, moq_span, ref_ion, ref_charge, brho, gammat, lise, harmonics, filename, data_time, skip_time, binning):
    
    moqs_in_window = get_all_in_moq_window(moq_cen, moq_span)
    mydata = ImportData(ref_ion, alphap)
    mydata._set_secondary_args(lise, harmonics)
    mydata._set_tertiary_args(filename, data_time, skip_time, binning)
    mydata._exp_data()
    mydata._calculate_srrf(moqs_in_window)
    mydata._simulated_data(particles = True)
    mycanvas = CreateGUI(ref_nuclei, mydata.nuclei_names, 1, 0)
    mycanvas._view(mydata.exp_data, mydata.simulated_data_dict, filename)
    '''
    Example 1: moq_cen = 2.054, moq_span = 0.002, ref_charge = 32, ref_nuclei = '72Ge', brho = 6.937117, gammat = 1.395, lise = 'e.lpp', harmonics = [200]
    filename = '/lustre/ap/litv-exp/2021-07-03_E143_TwoPhotonDecay_ssanjari/analyzers/410/410MHz-2021.07.02.17.58.30.204.tiq', data_time = 10, skip_time = 2.5, binning = 2048
    Example 2: brho = 6.980882, gammat = 1.395, harmonics = [120], #filename = '/lustre/ap/litv-exp/2021-07-03_E143_TwoPhotonDecay_ssanjari/analyzers/245/245MHz-2021.07.02.17.58.56.521.tiq'
    filename = '/lustre/ap/litv-exp/2021-07-03_E143_TwoPhotonDecay_ssanjari/ntcap/iq/IQ_2021-06-30_23-27-34_part3/0002757.iq.tdms', data_time = 10, skip_time = 2.3, binning = 4096
    '''
    
def get_energy_isomer_deltaf(charge, zz, nn, ref_iso, brho, gammat, delta_f, harmonic):
    '''
    Example: 72 Br, harmonic = 200, charge = 35, ref_iso = '72Br', ref_nuclei = f'{ref_iso}+{charge}', brho = 6.937117, gammat = 1.056, delta_f = 549
    Example: 72 Ge-> harmonic = 209, zz = 32, nn = 40, charge = 32, ref_iso = '72Ge', brho = 6.9303, gammat = 1.5, delta_f = 1.85e3
    '''
    imp = ImportData(ref_iso, charge, brho, gammat)
    imp.calculate_moqs([Particle(zz, nn, AMEData(), imp.ring)])
    
    imp._calculate_srrf()
    e_isomer = get_energy_isomer(delta_f, imp.moq[ref_nuclei], charge, imp.frequence_rel, gammat, harmonic)
    print(e_isomer)

def determine_deltaf_betweentwopartiucles():
    ring = Ring('ESR', 108.4) # 108.43 Ge
    p1 = Particle(36, 36, AMEData(), ring)
    p4 = Particle(36, 36, AMEData(), ring)
    p2 = Particle(35, 37, AMEData(), ring)
    p3 = Particle(35, 37, AMEData(), ring)
    p5 = Particle(34, 36, AMEData(), ring)
    p6 = Particle(21, 16, AMEData(), ring)
    p3.qq = 34
    p4.qq = 35
    p6.qq = 21
    moq1 = p1.get_ionic_moq_in_u()
    moq2 = p2.get_ionic_moq_in_u()
    moq3 = p3.get_ionic_moq_in_u()
    moq4 = p4.get_ionic_moq_in_u()
    moq5 = p5.get_ionic_moq_in_u()
    moq6 = p5.get_ionic_moq_in_u()
    deltamoq = abs(moq6 -moq4)
    moqratio = deltamoq / moq4
    # m2 = AMEData.to_mev(moq2 * p2.qq)
    # gamma = 
    deltaf = moqratio * 406.386 * 10**6 * (1 / 1.4**2)
    print(deltaf, moq1, moq2, moq3, moq4, moq5, moq6)

def calculate_moq_from_deltaf():
    deltaf = 12.5*10**3
    meassured_ref_frequency = 406.886 * 10**6
    gammat = 1.3956
    ref_moq = 2.054784866998157
    return deltaf / meassured_ref_frequency * ref_moq * gammat**2 + ref_moq

def candidates_by_moq_terminal(moq_of_unknown_particle, span = 0.001): #updated
    possible_particles = get_all_in_moq_window(moq_of_unknown_particle, span)
    return possible_particles

def candidates_by_moq_terminal_2(moq_of_unknown_particle, span = 0.001): #new
    possible_particles = get_all_in_moq_window(moq_of_unknown_particle, span)
    MyData = ImportData(refion, alphap)
    MyData._calculate_srrf(moqs = possible_particles, fref = fref)
    MyData._simulated_data(particles = True)
    return MyData.simulated_data_dict

def candidates_by_moq(moq, refion, fref, alphap, span = 0.001): #new
    possible_particles = get_all_in_moq_window(moq, span)
    print(possible_particles)
    #MyData = ImportData(refion, alphap)
    #MyData._calculate_srrf(moqs = possible_particles, fref = fref)
    #MyData._simulated_data(particles = True)
    #MyView = CreateGUI(refion, MyData.nuclei_names, 1, 0, True)
    #MyView._view()

#def candidates_by_harmonic(particles, fcen, fspan):
#    compute srf particles
    
#    multiply all of them by integers until they fall in fspan
#    print particle and frequency
    
def something2():
    ring = Ring('ESR', 108.4) # 108.43 Ge
    #Se70 = Particle(34, 36, AMEData(), ring)
    p = Particle(35, 37, AMEData(), ring)
    mass = AMEData.to_mev(p.get_ionic_moq_in_u() * p.qq)
    print(mass)

if __name__ == '__main__':
    #main()
    #determine_deltaf_betweentwopartiucles()
    #calculate_moq_from_deltaf()
    #something()
    #something2()
    #candidates_by_moq(2.2470182028797723, '72Ge+32', 407095000, 0.513868, 0.0005)
    #candidates_by_moq(2.054784866998157, '72Br+35', 407095000, 0.513868, 0.0015)
    candidates_by_moq(2.056321608595159, '70Se+34', 407095000, 0.513868, 0.0015)
