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
            print(particle_name)
            moq_in_window[particle_name] = moq
            
    return moq_in_window

def simulation_moq_window(moq_cen, moq_span, ref_nuclei, ref_charge, brho, gammat, lise, harmonics, filename, data_time, skip_time, binning):
    moqs_in_window = get_all_in_moq_window(moq_cen, moq_span)
    mydata = ImportData(ref_nuclei, ref_charge, brho, gammat)
    mydata._set_secondary_args(lise, harmonics)
    mydata._set_tertiary_args(filename, data_time, skip_time, binning)
    mydata._exp_data()
    mydata._calculate_srrf(moqs_in_window)
    mydata._simulated_data(particles = True)
    mycanvas = CreateGUI(ref_nuclei, mydata.nuclei_names, 1, 0)
    mycanvas._view(mydata.exp_data, mydata.simulated_data_dict, filename)
    

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

def get_moq_of_particle(particle):
    return particle.get_ionic_moq_in_u()

def main():
    charge = 32
    ref_iso = '72Ge'
    brho = 6.9303
    gammat = 1.5
    delta_f = 1.85e3
    imp = ImportData(ref_iso, charge, brho, gammat)
    imp.calculate_moqs([Particle(32, 40, AMEData(), imp.ring)])
    imp._caluclate_srrf()
    e_isomer = get_energy_isomer(delta_f, imp.moq['72Ge+32'], charge, imp.frequence_rel, gammat, 209)
    print(e_isomer)

def moq_of_particle():
    moq = get_moq_of_particle(Particle(35, 37, AMEData(), Ring('ESR', 108.43)))
    print(moq)
    
def controller2():
    moq_cen = 2.054
    moq_span = 0.002
    ref_charge = 32
    ref_nuclei = '72Ge'
    brho = 6.937117
    gammat = 1.395
    lise = 'e.lpp'
    harmonics = [200]#208, 209, 210]
    filename = '/lustre/ap/litv-exp/2021-07-03_E143_TwoPhotonDecay_ssanjari/analyzers/410/410MHz-2021.07.02.17.58.30.204.tiq'
    data_time = 10
    skip_time = 2.5
    binning = 2048
    simulation_moq_window(moq_cen, moq_span, ref_nuclei, ref_charge, brho, gammat, lise, harmonics, filename, data_time, skip_time, binning)

if __name__ == '__main__':
    #moq_of_particle()
    #main()
    controller2()
    
    
