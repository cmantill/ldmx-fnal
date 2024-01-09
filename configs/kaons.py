from LDMX.Framework import ldmxcfg
from LDMX.SimCore import generators
from LDMX.SimCore import simulator
import sys

p=ldmxcfg.Process("egsim")
p.run = int(sys.argv[1])
sim = simulator.simulator( "kaon" )
sim.setDetector( "ldmx-det-v14" )
sim.description = "Kaon fired from upstream of the HCal"

kaon_gun = generators.gun('kaon_gun')
kaon_gun.particle = 'kaon0L'
kaon_gun.energy = 1.0 # GeV
kaon_gun.position = [0., 0., 840.]
kaon_gun.direction = [0., 0., 1.]
sim.generators.append(kaon_gun)

import LDMX.Ecal.EcalGeometry
import LDMX.Ecal.ecal_hardcoded_conditions
import LDMX.Hcal.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions

p.sequence=[  
        sim,
]

p.outputFiles=[sys.argv[2]]
p.maxEvents = 100
p.logFrequency = 1
