import argparse, sys

parser = argparse.ArgumentParser(f"ldmx fire {sys.argv[0]}")
parser.add_argument("nevents", type=int)
arg = parser.parse_args()

nevents = arg.nevents

"""
Configure electron with:
 - 1 to 8 GeV energy
 - varying angle/momentum
 - upstream tagger position
"""

# https://github.com/LDMX-Software/SimCore/blob/f4202a00f58609b4b0e70485b5675bc21144b241/python/generators.py.in#L294
position = [ -299.2386690686212, 0.0, -6000.0 ]
momentum = [ 434.59663056485   , 0.0, 7988.698356992288]
energy_up = 8
energy_dn = 1

# defines starting direction
import math
momentum_mag = math.sqrt(sum(map(lambda x: x*x, momentum)))
unit_direction = list(map(lambda x: x/momentum_mag, momentum))
direction = "/gps/direction {} {} {}".format(*unit_direction)

ele_5deg20 = [
    "/gps/particle e-",
    "/gps/pos/type Plane",
    direction,
    # Linear energy
    "/gps/ene/type Lin",
    f"/gps/ene/min {energy_dn} GeV",
    f"/gps/ene/max {energy_up} GeV",
    "/gps/ene/gradient 0",
    "/gps/ene/intercept 1",
    # Square
    "/gps/pos/shape Square",
    "/gps/pos/centre -44. 0. -880. mm",
    "/gps/pos/halfx 0.01 mm",
    "/gps/pos/halfy 0.01 mm",
    # angles
    "/gps/ang/type iso",
    "/gps/ang/mintheta 3.1414 rad",  # about 20 deg
    "/gps/ang/maxtheta 3.1415 rad",  # about 5 deg
]
print(ele_5deg20)

from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process("electron")
p.libraries.append("libSimCore.so")
p.libraries.append("libHcal.so")
p.libraries.append("libEcal.so")

from LDMX.SimCore import generators
from LDMX.SimCore import simulator
from LDMX.Recon import pfReco
ecalPF = pfReco.pfEcalClusterProducer()
ecalPF.doSingleCluster = False                                                                                                                        
ecalPF.logEnergyWeight = True

particle_source = generators.gps("particle_source", ele_5deg20)
sim = simulator.simulator("mySim")
sim.setDetector("ldmx-det-v14-8gev", True)
sim.description = "Single electron particle source"
sim.beamSpotSmear = [20.0, 80.0, 0.0]  # mm
sim.generators.append(particle_source)

p.run = 1
p.outputFiles = [f"data/egps_upstreamtagger_{energy_up}_{energy_dn}_gev.root"]
p.maxEvents = nevents
p.logFrequency = 1


import LDMX.Ecal.EcalGeometry
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions
import LDMX.Ecal.digi as ecal_digi

p.sequence = [
    sim,
    ecal_digi.EcalDigiProducer(),
    ecal_digi.EcalRecProducer(),
    ecalPF,
]
