from barion.ring import Ring
from barion.amedata import *
from barion.particle import *
from lisereader.reader import *
from pysimtof.importdata import *
import numpy as np


def get_all_in_moq_window(moq_cen, moq_span, ame_data=None, ring=None):
    if ame_data:
        ame_data = ame_data
    else:
        ame = AMEData()
        ame_data = ame.ame_table
    if ring:
        ring = ring
    else:
        ring = Ring('ESR', 108.43)
    moq_window=np.array([])
    p=Particle(1, 1, ame_data, ring)
    particles=p.get_all_in_all()
    for particle in particles:
        moq=particle.get_ionic_moq_in_u()
        if (moq >= moq_cen - moq_span / 2 and moq <= moq_cen + moq_span / 2):
            moq_window = np.append(moq_window, particle)
    return moq_window

def import_particles_from_lise(lisefile):
    lise_file = LISEreader(lisefile)
    lise_data = lise_file.get_info_all()
    return lise_data

def get_energy_isomer(delta_f, moq_gs, charge, fref, gammat, harmonic):
    moq_isomer=moq_gs * ( -delta_f / harmonic * gammat**2 / fref + 1)
    energy_respect_gs = (moq_isomer - moq_gs) * charge
    return AMEData.to_mev(energy_respect_gs)

def get_delta_frecuency_isomer(energy_isomer_respect_gs, moq_gs, charge, fref, gammat):
    #Introudce the energy in MeV
    return (fref * energy_isomer_respect_gs) / (AMEData.to_mev(charge * moq_gs) * gammat**2)

def get_frecuency_particle_ring(particle, f_ref, gammat, moq_ref):
    #f_ref is the real revolution frec (harmonic) of the ref. particle
    moq = particle.get_ionic_moq_in_u()
    return (1 - 1 / gammat / gammat * (moq-moq_ref) / moq_ref) * f_ref

def main():
    charge = 32
    ref_iso = '72Ge'
    brho = 6.9303
    gammat = 1.395
    delta_f = 6e3
    imp = ImportData()
    imp.set_ref_ion(ref_iso, charge)
    imp._import('e.lpp')
    imp._calculate(brho, gammat, charge)
    e_isomer = get_energy_isomer(delta_f, imp.moq[imp.ref_ion], charge, imp.frequence_rel, gammat, 209)
    print(e_isomer)

if __name__ == '__main__':
    main()
    
    
