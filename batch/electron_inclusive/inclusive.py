import sys

nevents = int(sys.argv[1])
output_dir = str(sys.argv[2])
output_run = int(sys.argv[3])

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
    #direction,
    "/gps/direction 0 0 1",
    # Linear energy
    "/gps/ene/type Lin",
    f"/gps/ene/min {energy_dn} GeV",
    f"/gps/ene/max {energy_up} GeV",
    "/gps/ene/gradient 0",
    "/gps/ene/intercept 1",
    # Square
    "/gps/pos/shape Square",
    "/gps/pos/centre 0 0 0 mm",
    # "/gps/pos/centre -44. 0. -880. mm",
    "/gps/pos/halfx 0.01 mm",
    "/gps/pos/halfy 0.01 mm",
    # angles
    "/gps/ang/type iso",
    "/gps/ang/mintheta 2.8 rad", # about 20 deg
    "/gps/ang/maxtheta 3.05 rad", # about 5 deg
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

mySim = simulator.simulator("mySim")
# ldmx-det-v14-8gev? 
mySim.setDetector( 'ldmx-det-v14' )
from LDMX.SimCore import generators as gen
mySim.generators.append( gen.single_8gev_e_upstream_tagger() )
mySim.beamSpotSmear = [20.,80.,0.]
mySim.description = 'single electron'

p.run = output_run
p.outputFiles = [f"{output_dir}/{output_run}_egps_upstreamtagger_{energy_up}_{energy_dn}_gev.root"]
p.maxEvents = nevents
p.logFrequency = 1

from LDMX.TrigScint.trigScint import TrigScintDigiProducer
from LDMX.TrigScint.trigScint import TrigScintClusterProducer
from LDMX.TrigScint.trigScint import trigScintTrack
ts_digis = [
        TrigScintDigiProducer.pad1(),
        TrigScintDigiProducer.pad2(),
        TrigScintDigiProducer.pad3(),
        ]
for d in ts_digis :
    d.randomSeed = 1

import LDMX.Ecal.EcalGeometry
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions
import LDMX.Ecal.digi as ecal_digi

p.sequence = [
    mySim,
    ecal_digi.EcalDigiProducer(),
    ecal_digi.EcalRecProducer(),
    ecalPF,
    *ts_digis,
    TrigScintClusterProducer.pad1(),
    TrigScintClusterProducer.pad2(),
    TrigScintClusterProducer.pad3(),
]

p.keep = [
    "drop TrigScintScoringPlaneHits.*",
    "drop SimParticles.*",
    "drop Trigger.*",
    "drop TaggerSimHits.*",
    "drop RecoilSimHits.*",
    "drop HcalSimHits.*",
    "drop EcalSimHits.*",
    "drop TargetSimHits.*",
    "drop Magnet.*",
    "drop TriggerPad1SimHits.*",
    "drop TriggerPad2SimHits.*",
    "drop TriggerPad3SimHits.*",
    "drop EcalDigis.*",
    "drop EcalRecHits.*",
    "drop HcalScoringPlaneHits.*",
    "drop TrackerScoringPlaneHits.*",
    "drop trigScintDigisPad1.*",
    "drop trigScintDigisPad2.*",
]

