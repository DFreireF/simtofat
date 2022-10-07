from barion.ring import Ring
from barion.amedata import *
from barion.particle import *
from lisereader.reader import *
from pysimtof.importdata import *
from pysimtof.creategui import *
import numpy as np


def get_all_in_moq_window(moq_cen, moq_span, ame_data = None, ring = None):

    if not ame_data:
        ame = AMEData()
        ame_data = ame.ame_table
    if not ring:
        ring = Ring('ESR', 108.43)
    
    moq_in_window = dict()
    p = Particle(1, 1, ame, ring)
    particles = p.get_all_in_all()
    
    for particle in particles:
        moq = particle.get_ionic_moq_in_u()
        if (moq >= moq_cen - moq_span / 2 and moq <= moq_cen + moq_span / 2):
            particle_name = f'{particle.tbl_aa}{particle.tbl_name}+{particle.qq}'
            moq_in_window[particle_name] = moq
    return moq_in_window

def get_isomers_particle(particle):
    isomers = list()
    #read data from somewhere
    return isomers

def get_particles_for_barion_from_lise(lise_data, ame_data):
    return [Particle(lise[2], lise[3], ame_data, ring) for lise in lise_data for ame in ame_data if lise[0] == ame[6] and lise[1] == ame[5]]
    
def import_particles_from_lise(lisefile):
    
    lise_file = LISEreader(lisefile)
    lise_data = lise_file.get_info_all()
    return lise_data

def get_energy_isomer(delta_f, moq_gs, charge, fref, gammat):
    
    moq_isomer = moq_gs * ( -delta_f * gammat**2 / fref + 1)
    energy_respect_gs = (moq_isomer - moq_gs) * charge
    return AMEData.to_mev(energy_respect_gs)

def get_delta_frecuency_isomer(energy_isomer_respect_gs, moq_gs, charge, fref, gammat):
    
    #Introudce the energy in MeV
    return (fref * energy_isomer_respect_gs) / (AMEData.to_mev(charge * moq_gs) * gammat**2)

def get_frecuency_particle_ring(particle, f_ref, gammat, moq_ref):
    
    #f_ref is the real revolution frec (harmonic) of the ref. particle
    moq = particle.get_ionic_moq_in_u()
    return (1 - 1 / gammat / gammat * (moq-moq_ref) / moq_ref) * f_ref

def moq_of_particle(zz, nn, amedata, ring):
    return Particle(zz, nn, amedata, ring).get_ionic_moq_in_u()

def get_mass_particle(particle):
    return AMEData.to_mev(particle.get_ionic_moq_in_u() * particle.qq)

def correct_shift(xx, yy, zz, every = 7):
    # xx, yy, zz -> Spectrogram variables
    # f(t,p) ; t(p,f) ; p(t,f)
    freq_average_per_tframe = dict()
    for time_frame, _ in enumerate((zz[:,0])):
        freq_average = np.array([])
        i = 0
        while (i + every <= len(zz[0,:])):
            imin = i
            imax = i + every
            freq_average = np.append(freq_average, np.average(zz[time_frame, imin:imax]))
            freq_average_per_tframe[f'{time_frame}'] = freq_average
            i = i + 1
    ref = -1
    delta=0
    nzz = np.array([])
    for key in reversed(freq_average_per_tframe.keys()):
        max = freq_average_per_tframe[key][:].max()
        min = freq_average_per_tframe[key][:].min()
        if max >= 1.4 * min:
            index_max = np.argmax(freq_average_per_tframe[key][:])
            index_max_f = int((index_max + index_max + every - 1) / 2) # This has to be integer, if every is even (like 5,7,9)
            if ref == -1: 
                ref = index_max_f
            else:
                delta = ref - index_max_f
            nzz = np.append(nzz, np.roll(zz[int(key), :], delta))
        else: sys.exit('Empty?')
    nzz = np.reshape(nzz, np.shape(zz))
    return xx, yy, nzz
