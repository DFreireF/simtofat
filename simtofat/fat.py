from barion.ring import Ring
from barion.amedata import *
from barion.particle import *
from lisereader.reader import *
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

def get_energy_isomer(delta_f, moq_mother, mother_charge, frel):
    moq_isomer=moq_mother*(delta_f*gammat**2/frel+1)
    mass_isomer = moq_isomer * mother_charge
    return mass_isomer - moq_mother * mother_charge
