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

def get_particles_for_barion_from_lise(lise_data, ame_data):
    return [Particle(lise[2], lise[3], AMEData(), ring) for lise in lise_data for ame in ame_data if lise[0] == ame[6] and lise[1] == ame[5]]
    
def import_particles_from_lise(lisefile):
    
    lise_file = LISEreader(lisefile)
    lise_data = lise_file.get_info_all()
    return lise_data

def get_energy_isomer(delta_f, moq_gs, charge, fref, gammat, harmonic):
    
    moq_isomer = moq_gs * ( -delta_f / harmonic * gammat**2 / fref + 1)
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
