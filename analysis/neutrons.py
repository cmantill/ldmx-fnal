import numpy as np
import os,math,sys
import argparse
import pickle

import EventTree
from cppyy.gbl import ldmx
# import libDetDescr as DD
import ROOT as r

"""
Neutron analysis
"""

def main(arg):
    print(arg.input_files)
    trees_by_filename = dict.fromkeys(arg.input_files)
    for filename in trees_by_filename.keys():
        trees_by_filename[filename] = EventTree.EventTree(filename)

    var_names = [
            "event_nhcalrechit",  # number of rechits
            "event_nhcalrechit_1pe",  # number of rechits with 
            "event_nhcalrechit_5pe",  # number of rechits with
            "event_hcalrechit_sumenergy",  # summed energy of the event
            "event_hcalrechit_maxenergy",  # max energy of the event
            "event_hcalrechit_maxlayer",  # maximum layer reached in the event
            "event_hcalrechit_maxlayer_nhits",  # maximum layer with number of hits
            "event_hcalrechit_maxlayer_z",  # maximum layer reached in the z direction
            "event_hcalrechit_maxpe",  # maximum potential energy 
            "event_hcalrechit_maxpe_layer",  # layer where the maximum potential energy occured (?)
            "event_hcalrechit_sumpe",  # sum of the potential energies of the event
            "event_hcalrechit_isback", # determines if the hit occurs at the back of the Hcal
            "event_hcalrechit_nuniquelayers",  # number of unique layers interacted with in the event
            "hcalrechit_uniquelayer",  # unique layer of each hit
            "hcalrechit_uniquelayer_nhits",  # number of hits in each unique hit
            "hcalrechit_uniquelayer_nstrips",  # number of strips interacted with in each unique layer in each event
            "hcalrechit_uniquelayer_sumenergy",  # summed energy in each unique layer for each event

            "hcalrechit_x",  # x position of each rechit
            "hcalrechit_y",  # y position of each rechit
            "hcalrechit_z",  # z position of each rechit
            "hcalrechit_minz", # minimum z position of each rechit
            "event_hcalrechit_minPEover1_zval",  # z value of the lowest-PE recHit with over 1 PE of reco energy
            "event_hcalrechit_minPEover5_zval",  # z value of the lowest-PE recHit with over 5 PE of reco energy
            "hcalrechit_maxz", # maximum z position of each rechit
            "hcalrechit_energy",  # energy of each rechit
            "hcalrechit_layer",  # layer of each rechit
            "hcalrechit_section",  # section of each rechit
            "hcalrechit_strip",  # strip of each rechit
            "hcalrechit_pe",  # pe of each rechit
            "didithit", #did it or didn't it
            "misspercentage", # how many missed
            "event_nhcalrechit_1peallhits", # the array with hits over 1 pe and zeroes
            #"layer_phi", # theta of a layer using max x and y
            #"phivariable", #theta var to be used for plotting directly without processing
            #"thetaphi",
            "theta", #theta angle of scatter using atan over z
            "phi",  #phi is just the angle per layer aka only 2d
            #"symmtheta", #old symmetric theta, use clus_theta now
            "maxtheta", #max theta 
            "maxphi", #max phi
            "maxsymmtheta", #maximum symmetric theta, needs to be modified to use max clus theta instead of old symmtheta
            #"symthetaStd",
            #"maxthetaStd", 
            #"maxphiStd", 
            #"phiStd", 
            #"thetaStd", 
            #"symthetaStdAvg", 
            #"thetaStdAvg", 
            #"phiStdAvg",
            #"thetamean",
            "sum_e", #sum energy for weights
            "clus_x", #cluster x energy weight
            "clus_y", #cluster y energy weight
            "clus_z", #cluster z energy weight
            "clus_theta", #symmetric energy weighted theta! use this one
            "clusxuw", # unweighted avg xpos
            "clusyuw", # unweighted avg ypos
            "cluszuw", #unweighted avg zpos

            "event_neutron_e",  # energy of the neutron in each event
            "event_neutron_pz",  # momentum of the neutron in the z direction of each event
            "event_neutron_kine",  # neutron kinetic energy in each event
            "event_neutron_theta",  # angle of the event

            "event_necalrechit", # number of ecal hits
            "EcalScoringPlaneHits", #number of ecal scoring plane hits, using this in meantime while main code is debugged for 870 and 690 ecal hits
            "ecalrechit_x", # x position of ecal rechit
            "ecalrechit_y", # y position of ecal rechit
            "ecalrechit_z", # z position of ecal rechit
            "ecalrechit_energy", # energy of each ecal rechit
            "ecalrechitlayer", # layer of each ecal rechit
            "event_ecalminlayer", # minimum layer of each ecal event 
            "event_ecalmaxlayer", # maximum layer of each ecal event

            "event_ecalrechit_sumenergy", # sum energy per ecal event
            "event_ecalrechit_maxenergy", # max energy per ecal event
            "event_ecalrechit_avgenergy", # avg energy per ecal event
            "event_ecalrechit_maxz", # maximum z position per ecal event
    "event_ecalrechit_maxe_zposition", # z position of hit with the max energy
        "event_ecalrechit_minz", # minimum z position per ecal event with null binned to 0
        "event_ecalrechit_maxe_layer" # max energy layer per ecal event


    ]

    variables_by_filename = {}
    

    for filename,tree in trees_by_filename.items():
        count=0
        variables = dict.fromkeys(var_names)
        for key in variables.keys():
            variables[key] = []

        for ie,event in enumerate(tree):
            count +=1
            if count>1000: break
            if arg.max_events!=-1 and ie>=arg.max_events: continue

            # Hcal RecHits
            hits = dict.fromkeys([
                "energy","xpos","ypos","zpos",
                "layer","strip","section","pe", "xenergy", "yenergy", "zenergy", "theta", "radius", "phi", "symmtheta", "sum_e", "clus_x", "clus_y", "clus_z", "clus_theta", "clusxuw", "cluszuw", "clusyuw"
                ])
            for key in hits.keys(): hits[key] = []
            nhitsback=0
            ih = 0
            for ih,hit in enumerate(event.HcalRecHits):
                if hit.getSection()!=0: 
                    continue
                nhitsback+=1
                hits["energy"].append(hit.getEnergy())
                #print(hits["energy"])
                hits["xpos"].append(hit.getXPos())
                hits["ypos"].append(hit.getYPos())
                hits["zpos"].append(hit.getZPos())
                hits["pe"].append(hit.getPE())
                hits["xenergy"].append(hit.getXPos()*hit.getEnergy())
                hits["yenergy"].append(hit.getYPos()*hit.getEnergy())
                hits["zenergy"].append(hit.getZPos()*hit.getEnergy())

                hits["sum_e"] = np.sum(hit.getEnergy())
                hits["clus_x"] = np.sum(hit.getXPos() * hit.getEnergy()) / np.sum(hit.getEnergy())
                hits["clus_y"] = np.sum(hit.getYPos() * hit.getEnergy()) / np.sum(hit.getEnergy())
                hits["clus_z"] = np.sum(hit.getZPos() * hit.getEnergy()) / np.sum(hit.getEnergy())
                hits["clusxuw"] = (hit.getXPos()/len(hits))
                hits["clusyuw"] = (hit.getYPos()/len(hits))
                hits["cluszuw"] = (hit.getZPos()/len(hits))
                #hits["clus_theta"] = np.degrees(math.atan(np.sign((np.sum(hit.getYPos() * hit.getEnergy()) / np.sum(hit.getEnergy()))*(np.hypot(np.sum(hit.getXPos() * hit.getEnergy())/np.sum(hit.getEnergy()),np.sum(hit.getYPos() * hit.getEnergy()) / np.sum(hit.getEnergy())))/( np.sum(hit.getZPos() * hit.getEnergy()) / np.sum(hit.getEnergy())))))
                hits["clus_theta"] = np.sign(hits["clus_y"]) * np.hypot(hits["clus_x"],hits["clus_y"])/hits["clus_z"]
                
                #hit_id = hit.getID()
                #hit_hcalid = DD.HcalID(hit_id)
                #hits["section"].append(hit_hcalid.section())
                #hits["layer"].append(hit_hcalid.layer())
                #hits["strip"].append(hit_hcalid.strip())

                hits["section"].append(hit.getSection())
                hits["layer"].append(hit.getLayer())
                hits["strip"].append(hit.getStrip())
    
                hits["radius"].append(np.sqrt(np.square(hit.getXPos()*hit.getEnergy()) + np.square(hit.getYPos()*hit.getEnergy())))
                hits["theta"].append(math.degrees(math.atan((np.sqrt(np.square(hit.getXPos()*hit.getEnergy()) + np.square(hit.getYPos()*hit.getEnergy()))))/(np.max(hit.getZPos()*hit.getEnergy()))))
                #if hit.getYPos() < 0:
                #    hits["symmtheta"].append(-(math.degrees(math.atan((np.sqrt(np.square(hit.getXPos()*hit.getEnergy()) + np.square(hit.getYPos()*hit.getEnergy()))))/(np.max(hit.getZPos()*hit.getEnergy())))))
                #else:
                #    hits["symmtheta"].append(math.degrees(math.atan((np.sqrt(np.square(hit.getXPos()*hit.getEnergy()) + np.square(hit.getYPos()*hit.getEnergy()))))/(np.max(hit.getZPos()*hit.getEnergy()))))
                    #print(hits["symmtheta"])
                if hit.getYPos()!=0 and hit.getXPos()!=0:
                    hits["phi"].append(math.degrees(math.atan((hit.getYPos()*hit.getEnergy())/(hit.getXPos()*hit.getEnergy()))))
                #print(hits["theta"])


            nhits = ih

            variables["event_nhcalrechit"].append(nhits)
            variables["hcalrechit_x"].extend(hits["xpos"])
            variables["hcalrechit_y"].extend(hits["ypos"])
            variables["hcalrechit_z"].extend(hits["zpos"])


            variables["hcalrechit_energy"].extend(hits["energy"])
            variables["hcalrechit_layer"].extend(hits["layer"])
            variables["hcalrechit_strip"].extend(hits["strip"])
            variables["hcalrechit_section"].extend(hits["section"])
            variables["hcalrechit_pe"].extend(hits["pe"])




            for key,item in hits.items(): hits[key] = np.array(item)
            variables["event_nhcalrechit_1peallhits"].append((hits["pe"]>1).sum())
            #if nhitsback > 0 :
            if (hits["section"]==0).sum()>0:
                try:
                    variables["theta"].extend(hits["theta"])
                    variables["phi"].extend(hits["phi"])
                    #variables["symmtheta"].extend(hits["symmtheta"])
                    variables["maxtheta"].append(np.max(hits["theta"]))
                    #variables["maxtheta"].append(max(hits["theta"].min(), hits["theta"].max(), key=abs))
                    variables["maxphi"].append(np.max(hits["phi"]))
                    #variables["maxphi"].append(max(hits["phi"].min(), hits["phi"].max(), key=abs))
                    #variables["maxsymmtheta"].append(np.max(hits["symmtheta"]))

                    #variables["maxsymmtheta"].append(np.abs(hits["symmtheta"]).max())
                    variables["maxsymmtheta"].append(max(hits["symmtheta"].min(), hits["symmtheta"].max(), key=abs))
                    #variables["symthetaStd"].append(np.std(max(hits["symmtheta"].min(), hits["symmtheta"].max(), key=abs)))
                    #variables["maxthetaStd"].append(np.std(np.max(hits["theta"])))
                    #variables["maxphiStd"].append(np.std(np.max(hits["phi"])))
                    #variables["phiStd"].append(np.std(hits["phi"]))
                    #variables["thetaStd"].append(np.std(hits["theta"]))
                    #variables["symthetaStdAvg"].append(np.mean(np.std(max(hits["symmtheta"].min(), hits["symmtheta"].max(), key=abs))))
                    #variables["thetaStdAvg"].append(np.mean(np.std(hits["theta"])))
                    #variables["phiStdAvg"].append(np.mean(np.std(hits["phi"])))
                    #variables["thetamean"].append(np.mean(hits["theta"]))
                    variables["clusxuw"].append(np.mean(hits["xpos"]))
                    variables["clusyuw"].append(np.mean(hits["ypos"]))
                    variables["cluszuw"].append(np.mean(hits["zpos"]))
                    #variables["sum_e"] = np.sum([hit.getEnergy() for hit in hits])
                    #variables["clus_x"] = np.sum([hit.getXPos() * hit.getEnergy() for hit in hits]) / variables["sum_e"]
                    #variables["clus_y"] = np.sum([hit.getYPos() * hit.getEnergy() for hit in hits]) / variables["sum_e"]
                    #variables["clus_z"] = np.sum([hit.getZPos() * hit.getEnergy() for hit in hits]) / variables["sum_e"]
                    variables["clus_theta"].append(hits["clus_theta"])
                    variables["sum_e"].append(hits["sum_e"])
                    variables["clus_x"].append(hits["clus_x"])
                    variables["clus_y"].append(hits["clus_y"])
                    variables["clus_z"].append(hits["clus_z"])
                    #variables["stdangposx"].append(np.std(hits["clus_x"]))
                    #variables["stdangposy"].append(np.std(hits["clus_y"]))

                    #print(variables["clus_theta"])


                    #max(hits["symmtheta"].min(), hits["symmtheta"].max(), key=abs)
                    #print(variables["maxsymmtheta"])
                except ValueError:
                    pass

                variables["event_nhcalrechit_1pe"].append((hits["pe"]>1).sum())
                variables["event_nhcalrechit_5pe"].append((hits["pe"]>5).sum())



                variables["event_hcalrechit_sumenergy"].append(np.sum(hits["energy"]))
                variables["event_hcalrechit_maxenergy"].append(np.max(hits["energy"]))
                variables["hcalrechit_minz"].append(np.min(hits["zpos"]))

            else:
                variables["hcalrechit_minz"].append(0)

                variables["event_hcalrechit_sumenergy"].append(0)


            if (hits["pe"]>1).sum() > 0:
                pes = hits["pe"]
                highPE_indices = np.where(pes>1)
                minPE_Val = np.min(pes[highPE_indices])
                minPE_Index = np.argmin(pes[highPE_indices])
                minPE_Z = hits["zpos"][minPE_Index]
                #print(pes, highPE_indices, pes[highPE_indices], minPE_Val, minPE_Index, minPE_Z)
                variables["event_hcalrechit_minPEover1_zval"].append(minPE_Z)

            if (hits["pe"]>5).sum() > 0:
                pes = hits["pe"]
                highPE_indices2 = np.where(pes>5)
                minPE_Val2 = np.min(pes[highPE_indices2])
                minPE_Index2 = np.argmin(pes[highPE_indices2])
                minPE_Z2 = hits["zpos"][minPE_Index2]
                #print(pes, highPE_indices2, pes[highPE_indices2], minPE_Val2, minPE_Index2, minPE_Z2)
                variables["event_hcalrechit_minPEover5_zval"].append(minPE_Z2)


              #  variables["hcalrechit_maxz"].append(np.max(hits["zpos"]))
              #  maxlayer = np.max(hits["layer"])
              #  mask_maxlayer = hits["layer"]==maxlayer
              #  hits_maxlayer = np.where(mask_maxlayer)
              #  z_maxlayer = np.unique(hits["zpos"][mask_maxlayer])[0]


              #  variables["event_hcalrechit_maxlayer"].append(maxlayer)
              #  variables["event_hcalrechit_maxlayer_nhits"].append(len(list(hits_maxlayer)))
              #  variables["event_hcalrechit_maxlayer_z"].append(z_maxlayer)

              #  variables["event_hcalrechit_maxpe"].append(np.max(hits["pe"]))
              #  variables["event_hcalrechit_maxpe_layer"].append(np.unique(hits["layer"][hits["pe"]==np.max(hits["pe"])])[0])

              #  variables["event_hcalrechit_sumpe"].append(np.sum(hits["pe"]))



                nhitsback=0

                #radius = max(sqrt(((np.max(hits["ypos"]))**2 + (np.min(hits["xpos"]))**2), ((np.max(hits["ypos"]))**2 + (np.max(hits["xpos"]))**2),
                        #((np.min(hits["ypos"]))**2 + (np.min(hits["xpos"]))**2), ((np.min(hits["ypos"]))**2 + (np.max(hits["xpos"]))**2)))
                #variables["thetaPerHit"].append(math.degrees(math.atan(radius/np.max(hits["zpos"]))))


            isback = 0
            if (hits["section"]==0).sum()>0: 
                isback = 1
                variables["event_hcalrechit_isback"].append(isback)

            uniquelayer = np.unique(hits["layer"])
            sumenergy = []
            nhits = []
            strips = []
            nstrips = []
            for layer in uniquelayer:
                masklayer = (hits["layer"] == layer)
                sumenergy.append(np.sum(hits["energy"][masklayer]))
                nhits.append(masklayer.sum())                    
                strips.extend(list(hits["strip"][masklayer]))
                nstrips.append(len(list(hits["strip"][masklayer])))

             #   variables["event_hcalrechit_nuniquelayers"].append(len(list(uniquelayer)))
             #   variables["hcalrechit_uniquelayer"].extend(list(uniquelayer))
             #   variables["hcalrechit_uniquelayer_nhits"].extend(list(nhits))
             #   variables["hcalrechit_uniquelayer_nstrips"].extend(nstrips)
             #   variables["hcalrechit_uniquelayer_sumenergy"].extend(sumenergy)

            #else: 
             #   variables["hcalrechit_minz"].append(0)


            # "Truth" neutrons
            fill_neutron = False
            neutron_count = 0
            for id_particle in event.SimParticles:
                neutron_count+=1
                if fill_neutron:
                    continue
                if neutron_count>1:
                    neutron_count=0
                    continue
                part = id_particle.second
                if part.getPdgID()==2112 and not fill_neutron:
                    fill_neutron=True
                    energy = part.getEnergy()
                    momentum = part.getMomentum()
                    mass = part.getMass()

                    variables["event_neutron_e"].append(energy)
                    variables["event_neutron_pz"].append(momentum[2])
                    variables["event_neutron_kine"].append(energy - mass)
                    neutron = r.TLorentzVector()
                    neutron.SetPxPyPzE(momentum[0],momentum[1],momentum[2],energy)
                    variables["event_neutron_theta"].append(neutron.Theta())
                    neutron_count=0

            #ECal
            jh=0
            ehits = dict.fromkeys([
                "energy","xpos","ypos","zpos", "layer",
                ])
            for key in ehits.keys(): ehits[key] = []
            for jh,hit in enumerate(event.EcalRecHits):
                if hit.getEnergy()>0:
                    ehits["energy"].append(hit.getEnergy())
                ehits["xpos"].append(hit.getXPos())
                ehits["ypos"].append(hit.getYPos())
                ehits["zpos"].append(hit.getZPos())


                ehit_id = hit.getID()
                #ehit_ecalid = DD.EcalID(ehit_id)

                #ehit_ecalid = ehits.getLayerID(ehit.id)


                ehits["layer"].append((ehit_id >> 17) & 0xFF)

            variables["ecalrechitlayer"].extend(ehits["layer"])    

            nhits=jh


            #variables["event_necalrechit"].append(nhits)  #go back to this once we fix the rechits for 840 and 690
            #variables["ecalrechit_x"].extend(ehits["xpos"])
            #variables["ecalrechit_y"].extend(ehits["ypos"])
            #variables["ecalrechit_z"].extend(ehits["zpos"])


            for key,item in ehits.items(): ehits[key] = np.array(item)

            if nhits>0:
                #variables["ecalrechit_energy"].extend(ehits["energy"])
                #variables["event_ecalrechit_sumenergy"].append(np.sum(ehits["energy"]))
                #variables["event_ecalrechit_maxenergy"].append(np.max(ehits["energy"]))
                #variables["event_ecalrechit_avgenergy"].append(np.mean(ehits["energy"]))
                #variables["event_ecalrechit_maxz"].append(np.max(ehits["zpos"]))
                #variables["event_ecalrechit_maxe_zposition"].append(
                #        ehits["zpos"][ehits["energy"]==np.max(ehits["energy"])]
            #            )

                #variables["event_ecalrechit_maxe_layer"].append(
                #    ehits["layer"][ehits["energy"]==np.max(ehits["energy"])]
                #)

                variables["event_ecalminlayer"].append(np.min(ehits["layer"]))
                variables["event_ecalmaxlayer"].append(np.max(ehits["layer"]))
                variables["event_ecalrechit_minz"].append(np.min(ehits["zpos"]))

            else:
                variables["event_ecalrechit_minz"].append(0)


            #EcalScoringPlaneHits
            #kh=0
            #ephits = dict.fromkeys([
            #    "energy","xpos","ypos","zpos",
            #])
            #for key in ephits.keys(): ephits[key] = []
            #for kh,hit in enumerate(event.EcalScoringPlaneHits):
             #   ephits["energy"].append(hit.getEnergy())



            #nhits=kh


            #variables["EcalScoringPlaneHits"].append(nhits)  #go back to this once we fix the rechits for 840 and 69i0


        # print(variables)
        for key,arr in variables.items():
            variables[key] = np.array(arr)
            #misspercentarray = []
        #misspercentage = variables["event_nhcalrechit_1peallhits"][variables["event_nhcalrechit_1peallhits"]>0].size/variables["event_nhcalrechit_1peallhits"].size * 100

        #print(misspercentage)
        #misspercentarray.append(variables["event_nhcalrechit_1peallhits"][variables["event_nhcalrechit_1peallhits"]>0].size/variables["event_nhcalrechit_1peallhits"].size * 100)
        #print(misspercentarray)
        
        #print("completed file:" + filename)
        #print("here's the stdev x for " + filename)
        #print(np.mean(np.std(variables["clus_x"])))
        #print("here's the mean x for " + filename)
        #print(np.mean(variables["clus_x"]))
        #print("here's the stdev y for " + filename)
        #print(np.mean(np.std(variables["clus_y"])))
        #print("here's the mean y for " + filename)
        #print(np.mean(variables["clus_y"]))
        variables_by_filename[filename] = variables

    os.system(f"mkdir -p output/")
    with open(f"output/{arg.output}.pkl", "wb") as pickle_file:
        pickle.dump(variables_by_filename,pickle_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(f'ldmx python3 {sys.argv[0]}')
    parser.add_argument('input_files',nargs='+')
    parser.add_argument('--output',required=True,help='Output name of pickle file (without extension)')
    parser.add_argument('--max_events',default=-1,type=int)
    arg = parser.parse_args()

    main(arg)
