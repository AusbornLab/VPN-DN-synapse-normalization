from neuron import h, gui
from neuron.units import ms, mV
import time as clock
import os
import re
import numpy
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import random
import math
import seaborn as sns
import csv
import pandas as pd
from tkinter import Tk
import tkinter.filedialog as fd
from pathlib import Path

#2 = axon | 3 = dendrite | 0 = unlabeled | 1 = soma
h.load_file("stdrun.hoc")

def instantiate_swc(filename):
    ''' 
    Load swc file and instantiate it as cell
    Code source: https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3257
    '''

    # load helper library, included with Neuron
    h.load_file('import3d.hoc')

    # load data
    cell = h.Import3d_SWC_read()
    cell.input(filename)

    # # instantiate
    i3d = h.Import3d_GUI(cell,0)
    i3d.instantiate(None)
    
def change_Ra(ra=41.0265, electrodeSec=None, electrodeVal=None):
    for sec in h.allsec():
        sec.Ra = ra
    if electrodeSec is not None:
        electrodeSec.Ra = electrodeVal

def change_gLeak(gleak=0.000239761, erev=-58.75, electrodeSec=None, electrodeVal=None):
    for sec in h.allsec():
        sec.insert('pas')
        for seg in sec:
            seg.pas.g = gleak
            seg.pas.e = erev
    if electrodeSec is not None:
        for seg in electrodeSec:
            seg.pas.g = electrodeVal
            seg.pas.e = 0

def change_memCap(memcap=1.3765, electrodeSec=None, electrodeVal=None):
    for sec in h.allsec():
        sec.cm = memcap
    if electrodeSec is not None:
        electrodeSec.cm = electrodeVal

def nsegDiscretization(sectionListToDiscretize):
    #this function iterates over every section, calculates its spatial constant (lambda), and checks if the length of the segments within this section are less than 1/10  lambda
    #if true, nothing happens
    #if false, the section is broken up until into additional segments until the largest segment is no greater in length than 1/10 lambda

    #NOTE TO SELF, NOT IDEAL BC ONLY INCREASES, DOESNT DECREASE, WRITE CODE TO DECREMENT SEGMENT NUMBER IF APPLICABLE
    #NOTE need a check on if num seg > max allowaable num seg to avoid error
    #TODO FIX TO SAVE ON COMPUTATIONAL COMPLEXITY

    for sec in sectionListToDiscretize:
        #old code which may be useful for debugging
        #secMorphDict = sec.psection()['morphology']

        #by calling sec.psection, we can return the density mechanisms for the given section in a dictionary
        #this lets us access the gleak value for a given section, which we use to calculate the membrane resistance
        secDensityMechDict = sec.psection()['density_mechs']

        #using the membrane resistance, section diameter, and section axial resistance, we calculate the spatial constant (lambda) for a given section
        secLambda = math.sqrt( ( (1 / secDensityMechDict['pas']['g'][0]) * sec.diam) / (4*sec.Ra) )

        #debugging print statement
        #print("\nlength of section", sec, "is", (sec.L/sec.nseg)/secLambda, "lambda | sec.L:", sec.L, "| sec.nseg:", sec.nseg, "| lambda:", secLambda)

        #if the segment length of a section is greater than 1/10 lambda, then the section is broken into additional segments
        #the segment count per section must be odd to avoid issues with the way NEURON handles midpoints for internal calculations, so if the necessary number of segments to reach 
        #a max length of 1/10 lambda is even, an extra section is added
        if (sec.L/sec.nseg)/secLambda > 0.1:
            #debugging print statement
            #print("\nlength of section", sec, "is", (sec.L/sec.nseg)/secLambda, "lambda")

            numSeg_log2 = math.log2(sec.L/secLambda / 0.1)
            numSeg = math.ceil(2**numSeg_log2)
            if numSeg % 2 == 0:
                numSeg += 1
            sec.nseg = numSeg

            #debugging print statements
            #print("\nlength of section", sec, "is now", (sec.L/sec.nseg)/secLambda, "lambda")
            #print("fixed by using a total of", sec.nseg, "segments")
    return

def initializeModel(morph_file, neuron_name):

    cell = instantiate_swc(morph_file)

    allSections_nrn = h.SectionList()
    for sec in h.allsec():
        allSections_nrn.append(sec=sec)
    
    # Create a Python list from this SectionList
    # Select sections from the list by their index

    allSections_py = [sec for sec in allSections_nrn]    

    #extra sectionLists if necessary for visualization purposes
    colorR = h.SectionList()
    colorB = h.SectionList()
    colorK = h.SectionList()
    colorG = h.SectionList()
    colorV = h.SectionList()

    if neuron_name == "DNp01":
        sizIndex = 0
    elif neuron_name == "DNp02":
        sizIndex = 2
    elif neuron_name == "DNp03":
        sizIndex = 2
    elif neuron_name == "DNp04":
        sizIndex = 4
    elif neuron_name == "DNp06":
        sizIndex = 4

    axonList = h.SectionList()
    tetherList = h.SectionList()
    dendList = h.SectionList()

    for sec in allSections_py:
        if "soma" in sec.name():
            somaSection = sec
            colorG.append(somaSection)
        elif "axon" in sec.name():
            axonList.append(sec)
        elif "dend_11" in sec.name():
            tetherList.append(sec)
        elif "dend_6" in sec.name():
            axonEnd = sec
        elif "dend_12" in sec.name():
            sizSection = sec
        else:
            dendList.append(sec)
    i = 0
    # for sec in axonList:
    #     #if i <= 4:
    #     if i == sizIndex:
    #         sizSection = sec
    #     i += 1
    colorV.append(sizSection)

    allSections_py = createAxon(axonEnd, allSections_py, neuron_name)
    allSections_py, electrodeSec = createElectrode(somaSection, allSections_py, neuron_name)

    shape_window = h.PlotShape(h.SectionList(allSections_py))           # Create a shape plot
    shape_window.exec_menu('Show Diam')    # Show diameter of sections in shape plot
    shape_window.color_all(9)

    shape_window.color_list(axonList, 2)
    shape_window.color_list(colorG, 4)  
    shape_window.color_list(tetherList, 1)
    shape_window.color_list(dendList, 3)
    shape_window.color_list(colorV, 7)


    if neuron_name == "DNp01":
        erev = -66.6544#-72.5
        change_memCap(memcap=4.167)#3.5531)
        change_Ra(ra=17.6461)#30.4396)
        change_gLeak(gleak= 0.0011196588, erev=erev)#0.0012663288, erev=erev)
        # change_Ra(34)#33.2617)
        # change_gLeak(3.4e-4, erev=erev)#4.4e-9)    
        # change_memCap(1)#1.4261)
    elif neuron_name == "DNp02": 
        erev = -70.8#-70.812 TOOK FIRST 5000 DATAPOINTS FROM EACH HJ CURRENT STEP AND TOOK MEAN
        # change_Ra(91)
        # change_gLeak(0.0002415, erev=erev)  
        # change_memCap(1)
        change_memCap(memcap=1.595)
        change_Ra(ra=34.4495)
        change_gLeak(gleak=0.000147726, erev=erev)
    elif neuron_name == "DNp03":
        erev = -58.75#68
        # change_Ra(33)
        # change_gLeak(gleak=0.00034, erev=erev)
        # change_memCap(1)
        change_memCap(memcap=1.595)
        change_Ra(ra=34.4495)
        change_gLeak(gleak=0.000147726, erev=erev)
    elif neuron_name == "DNp04":
        erev=-72.5
        change_memCap(memcap=1.595)
        change_Ra(ra=34.4495)
        change_gLeak(gleak=0.000147726, erev=erev)
    elif neuron_name == "DNp06":
        erev = -60#-60.0406 SEE DNP02
        # change_Ra(91)
        # change_gLeak(0.0002415, erev=erev)  
        # change_memCap(1)
        change_memCap(memcap=1.595)
        change_Ra(ra=34.4495)
        change_gLeak(gleak=0.000147726, erev=erev)
    else:
        erev=-72.5
        change_Ra()
        change_gLeak()
        change_memCap()

    nsegDiscretization(allSections_py)
    

    return cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, shape_window, electrodeSec

def createAxon(axonEnd, pySectionList, neuron_name=None):
    equivCylAxon = h.Section()

    # print(axonEnd.n3d())
    #take last point and second to last and draw line through for dist = hiro/bhandawat or ams neutube
    lineStart_X = axonEnd.x3d(axonEnd.n3d()-2)
    lineStart_Y = axonEnd.y3d(axonEnd.n3d()-2)
    lineStart_Z = axonEnd.z3d(axonEnd.n3d()-2)
    lineEnd_X = axonEnd.x3d(axonEnd.n3d()-1)
    lineEnd_Y = axonEnd.y3d(axonEnd.n3d()-1)
    lineEnd_Z = axonEnd.z3d(axonEnd.n3d()-1)

    lineDir = numpy.array([lineEnd_X, lineEnd_Y, lineEnd_Z]) - numpy.array([lineStart_X, lineStart_Y, lineStart_Z])
    line_direction_norm = lineDir / numpy.linalg.norm(lineDir)

    # TODO: ADD OTHERS
    if neuron_name == "DNp01":
        # equivCylHeight = 184.656
        # equivCylDiam = 1.34*2
        equivCylHeight = 241.69#184.656
        equivCylDiam = 3.32*2#1.34*2*2
        # equivCylHeight = 360.346
        # equivCylDiam = 1.199*2
    elif neuron_name == "DNp02":
        equivCylHeight = 220.654
        equivCylDiam = 0.624*2
    elif neuron_name == "DNp03":
        # equivCylHeight = 180.14
        # equivCylDiam = 0.249*2
        equivCylHeight = 175.671
        equivCylDiam = 0.251*4
    elif neuron_name == "DNp04":
        #according to Namiki paper, DNp04 axon looks to be as long as DNp01 and as thick as DNp03, so we
        #use the DNp01 height variable and the DNp03 diam variable to create the approximate axon
        equivCylHeight = 360.346
        equivCylDiam = 0.249*2
    elif neuron_name == "DNp06":
        equivCylHeight = 222.766
        equivCylDiam = 0.392*2

    new_point = numpy.array([lineEnd_X, lineEnd_Y, lineEnd_Z]) + equivCylHeight * line_direction_norm
    ENDPT = [lineEnd_X, lineEnd_Y, lineEnd_Z]
    # print(ENDPT, new_point)

    x1 = axonEnd.x3d(axonEnd.n3d() - 1)
    y1 = axonEnd.y3d(axonEnd.n3d() - 1)
    z1 = axonEnd.z3d(axonEnd.n3d() - 1)
    d1 = axonEnd(1).diam
    # print(d1)


    equivCylAxon.pt3dadd(x1, y1, z1, equivCylDiam)
    equivCylAxon.pt3dadd(new_point[0], new_point[1], new_point[2], equivCylDiam)
    pySectionList.append(equivCylAxon)
    equivCylAxon.connect(axonEnd, 1)
    # print(equivCylAxon(0.5).area())

    return pySectionList

def createElectrode(somaSection, pySectionList, neuron_name=None):
    electrodeSec = h.Section()
    electrodeSec.L = 10
    electrodeSec.diam = 1
    electrodeSec.connect(somaSection, 0)

    pySectionList.append(electrodeSec)

    # equivCylAxon.connect(axonEnd, 1)
    # print(equivCylAxon(0.5).area())

    return pySectionList, electrodeSec

def write_csv(name,data):
    # print("writing to:", os.getcwd()+name)
    # print("we are in", os.getcwd())
    with open(name, 'a') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(data)

def returnQuints(dir_path):
    dirList = os.listdir(dir_path)
    numTriplets =  int(len([entry for entry in dirList if os.path.isfile(os.path.join(dir_path, entry))])/5)

    dirList.sort()
    quintList = []

    for triplet_set in range(numTriplets):
        print('reading data from simulation', triplet_set, '/', numTriplets)
        synSite_df = dirList[0+(triplet_set*5)]
        synVoltSIZ_df = dirList[1+(triplet_set*5)]
        synVoltSoma_df = dirList[2+(triplet_set*5)]
        synVoltSynpt_df = dirList[3+(triplet_set*5)]
        time_df = dirList[4+(triplet_set*5)]
        
        time_df = pd.read_csv(dir_path + '/' + time_df, header=None)
        synVoltSIZ_df = pd.read_csv(dir_path + '/' + synVoltSIZ_df, header=None)
        synVoltSoma_df = pd.read_csv(dir_path + '/' + synVoltSoma_df, header=None)
        synVoltSynpt_df = pd.read_csv(dir_path + '/' + synVoltSynpt_df, header=None)
        synSite_df = pd.read_csv(dir_path + '/' + synSite_df)
        
        quintList.append([time_df, synVoltSIZ_df, synVoltSoma_df, synVoltSynpt_df, synSite_df])
    
    return quintList

def loadSynMapDataframe(neuron_name=None):
    #select synapse data file using tkinter
    #TODO: make sure presynID is being read as a string
    Tk().withdraw()
    fd_title = "Select synapse map data file to use for synapse activation"
    syn_file = fd.askopenfilename(filetypes=[("csv file","*.csv")], initialdir=r"datafiles/morphologyData", title = fd_title)
    # synMap_df = pd.read_csv(syn_file)
    synMap_df = pd.read_csv(syn_file, dtype={'pre': str})

    if neuron_name == "DNp01":
        pass
    elif neuron_name == "DNp02":
        pass
    elif neuron_name == "DNp03":
        synMap_df = synMap_df.drop(synMap_df.columns[8:10], axis=1)
    elif neuron_name == "DNp04":
        synMap_df = synMap_df.drop(synMap_df.columns[8:10], axis=1)
    elif neuron_name == "DNp06":
        synMap_df = synMap_df.drop(synMap_df.columns[8:11], axis=1)
    # print(synMap_df.head())
    # print(synMap_df.dtypes)
    #quit(0)
    return synMap_df


#####
#Simulations
def activateSynapses(syn_df, MODE, interval = None, num_spike = None, somaSection=None, recordSection=None, sizSection=None, erev=-72.5, simTime=100, onsetTime=None, neuron_name=None, dirname=None, filename=None):
    #TODO: MAKE RECORDSECTION A LIST?
    #this function activates all synapses that exist within the dataframe (syn_df) that is passed to
    #whether the synapses are activated one at a time, all at once as a group, or sequentially (in either forward
    #or reverse order) is determind by the MODE parameter
    #which can be set to SINGLE, GROUP, SEQUENTIAL_F, or SEQUENTIAL_R
    #TODO: collapse SEQ_F and SEQ_R into one MODE

    #other parameters are as follows:
    #somaSection --- the section of the morphology that corresponds to the soma
    
    #recordSection --- the section of the morphology that is being recorded from. If no recordSection is passed to this function,
    #in the case of single synapse activation, the recordSection for each synapse will be the section to which that synapse is mapped to,
    #and in the case of group synapse activation, the recordSection will default to the soma. If no somaSection is provided, an error
    #will be raised

    #erev --- membrane resting/reversal potential (mV)

    #simTime --- simulation duration (ms) | onsetTime --- synapse activation start time (ms)

    #dirname --- TODO: FILL OUT

    #neuronName --- TODO: FILL OUT

    #filename --- TODO: FILL OUT

    if MODE == "SINGLE_EXP2":
        h.v_init = erev
        t_vec = h.Vector()
        t_vec.record(h._ref_t)
        seriesList = []
        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':
                pass
            else:
                synapses = {}
                syn_current_siz = {}
                syn_current_soma = {}
                syn_current_synsec = {}
                syn_current_synpt = {}

                seriesList.append(row)

                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break
                currSynRangeVar = row.loc['mappedSegRangeVar']

                stim = h.NetStim()
                stim.number = 1
                stim.start = 50
                stim.noise = 0
                exp2Gmax = 0.00027 
                exp2Rise = 0.2
                exp2Decay = 1.1           

                synapses['syn_{}'.format(0)] = h.Exp2Syn(currSynSec(currSynRangeVar))
                synapses['syn_{}'.format(0)].tau1 = exp2Rise
                synapses['syn_{}'.format(0)].tau2 = exp2Decay
                synapses['syn_{}'.format(0)].e = -10#-10


                ncstim = h.NetCon(stim, synapses['syn_{}'.format(0)])
                ncstim.delay = 0  
                ncstim.weight[0] = exp2Gmax
                ########################

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(0)] = h.Vector()
                if neuron_name == "DNp01":
                    # print("correct SIZ")
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.0551186)._ref_v)
                elif neuron_name == "DNp03":
                    # print(sizSection)
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.0551186)._ref_v)
                else:
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.5)._ref_v)

                syn_current_soma['v_vec_{}'.format(0)] = h.Vector()
                syn_current_soma['v_vec_{}'.format(0)].record(somaSection(0.5)._ref_v)

                syn_current_synpt['v_vec_{}'.format(0)] = h.Vector()
                syn_current_synpt['v_vec_{}'.format(0)].record(currSynSec(currSynRangeVar)._ref_v)
                
                #Run the simulation
                print("Round number {}/{}".format(index, syn_df.shape[0]))
                h.tstop = simTime
                h.run()


                synInfo_df = pd.DataFrame(seriesList)

                # Save Data 
                # Save simulated section and simulation data
                #TODO: SET UP FILE AUTOWRITE SO THAT WAY NEURON NAME IS PARAM IN FILE NAME

                if dirname is not None:
                    parent_dir = "datafiles/simulationData"
                    path = os.path.join(parent_dir, dirname)
                    try:
                        os.mkdir(path)
                    except OSError as error:
                        pass 

                    if filename is None:
                        outFilename = path + '/' + MODE
                        write_csv(outFilename + '_time.csv',t_vec)        
                        write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(0)])
                        write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(0)])
                        write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(0)])
                    else:
                        filename = path + '/' + filename
                        write_csv(filename + '_time.csv',t_vec)        
                        write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(0)])
                        write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(0)])
                        write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(0)])
                else:
                    if filename is None:
                        outFilename = 'datafiles/' + MODE
                        write_csv(outFilename + '_time.csv',t_vec)        
                        write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(0)])
                        write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(0)])
                        write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(0)])
                    else:
                        filename = 'datafiles/' + filename
                        write_csv(filename + '_time.csv',t_vec)        
                        write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(0)])
                        write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(0)])
                        write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(0)])


                synapses = None 
                syn_current = None
                tVec = None
                vVec = None
        if filename is None:
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
        else:
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

    elif MODE == "GROUP_EXP2":
        h.v_init = erev
        t_vec = h.Vector()
        t_vec.record(h._ref_t)
        seriesList = []
        syn_current_siz = {}
        syn_current_soma = {}
        syn_current_synpt = {}
        synpt_recording_list = []

        synapses = {}  # Using a dictionary to store synapses
        exp2Gmax = 0.00027 
        exp2Rise = 0.2
        exp2Decay = 1.1 

        for index, row in syn_df.iterrows():
            if not pd.isna(row['pre']):
                seriesList.append(row)

                # ONSET = onsetTime + numpy.random.normal(0, 0.5)  # Uncomment this line if you want to add random noise to ONSET
                currSynSecStr = row['mappedSection']

                # Find the corresponding section for the synapse
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break

                currSynRangeVar = row['mappedSegRangeVar']

                # Assign properties to the synapse
                syn = h.Exp2Syn(currSynSec(currSynRangeVar))
                syn.tau1 = exp2Rise
                syn.tau2 = exp2Decay
                syn.e = -10  # -10

                # Save the synapse in the dictionary
                synapses['syn_{}'.format(index)] = syn

                # Create recording vector to track voltage at the synapse location
                syn_current_synpt['v_vec_{}'.format(index)] = h.Vector()
                syn_current_synpt['v_vec_{}'.format(index)].record(currSynSec(currSynRangeVar)._ref_v)
                synpt_recording_list.append(syn_current_synpt['v_vec_{}'.format(index)])

        # Create NetStim and NetCon for activation
        nc = h.NetStim()
        nc.number = 1
        nc.start = onsetTime
        nc.noise = 0

        ncs = h.List()
        for i, syn in enumerate(synapses.values()):
            ncs.append(h.NetCon(nc, syn))
            ncs.object(i).weight[0] = exp2Gmax 
                # print(synInfo_df.head())

        syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
        if neuron_name == "DNp01":
            
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        elif neuron_name == "DNp03":
            print(sizSection.name())
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        else:
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.05)._ref_v)

        syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
        syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)
        h.tstop = simTime
        h.run()

        synInfo_df = pd.DataFrame(seriesList)

        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                outFilename = path + '/' + MODE
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = path + '/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                outFilename = 'datafiles/' + MODE
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None  

    elif MODE == "GROUP_EXP2_interval":
        h.v_init = erev
        t_vec = h.Vector()
        t_vec.record(h._ref_t)
        seriesList = []
        syn_current_siz = {}
        syn_current_soma = {}
        syn_current_synpt = {}
        synpt_recording_list = []

        synapses = {}  # Using a dictionary to store synapses
        exp2Gmax = 0.00027 
        exp2Rise = 0.2
        exp2Decay = 1.1 

        for index, row in syn_df.iterrows():
            if not pd.isna(row['pre']):
                seriesList.append(row)

                # ONSET = onsetTime + numpy.random.normal(0, 0.5)  # Uncomment this line if you want to add random noise to ONSET
                currSynSecStr = row['mappedSection']

                # Find the corresponding section for the synapse
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break

                currSynRangeVar = row['mappedSegRangeVar']

                # Assign properties to the synapse
                syn = h.Exp2Syn(currSynSec(currSynRangeVar))
                syn.tau1 = exp2Rise
                syn.tau2 = exp2Decay
                syn.e = -10  # -10

                # Save the synapse in the dictionary
                synapses['syn_{}'.format(index)] = syn

                # Create recording vector to track voltage at the synapse location
                syn_current_synpt['v_vec_{}'.format(index)] = h.Vector()
                syn_current_synpt['v_vec_{}'.format(index)].record(currSynSec(currSynRangeVar)._ref_v)
                synpt_recording_list.append(syn_current_synpt['v_vec_{}'.format(index)])

        # Create NetStim and NetCon for activation
        nc = h.NetStim()
        nc.number = num_spike
        nc.start = onsetTime
        nc.interval = interval
        nc.noise = 0

        ncs = h.List()
        for i, syn in enumerate(synapses.values()):
            ncs.append(h.NetCon(nc, syn))
            ncs.object(i).weight[0] = exp2Gmax 
                # print(synInfo_df.head())

        syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
        if neuron_name == "DNp01":
            
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        elif neuron_name == "DNp03":
            print(sizSection.name())
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        else:
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.05)._ref_v)

        syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
        syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)
        h.tstop = simTime
        h.run()

        synInfo_df = pd.DataFrame(seriesList)

        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                outFilename = path + '/' + MODE
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = path + '/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                outFilename = 'datafiles/' + MODE
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None  


    elif MODE == "GROUP":
        h.v_init = erev   # Resting membrane potential 

        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        synapses = {}
        syn_current_siz = {}
        syn_current_soma = {}
        syn_current_synpt = {}
        
        seriesList = []

        synpt_recording_list = []

        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
                pass
            else:
                seriesList.append(row)

                ONSET = onsetTime
                
                mu, sigma = 0, 0.5 # mean and standard deviation
                s = numpy.random.default_rng().normal(mu, sigma, 1)
                # ONSET += s#random.randint(-25, 25)
                # ONSET += 0#random.randint(-25, 25)
                # print(ONSET)

                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break
                currSynRangeVar = row.loc['mappedSegRangeVar']

                ###### FOR AMS FINAL SYNAPSE PROPERTY ASSIGNMENT (SINGLE ACT)

                synapses['syn_{}'.format(index)] = h.AlphaSynapse(currSynSec(currSynRangeVar))
                synapses['syn_{}'.format(index)].onset = ONSET       # Start time of synapse activation (ms)
                synapses['syn_{}'.format(index)].tau = 0.5        # Time to max conductance (ms)
                synapses['syn_{}'.format(index)].gmax = (5.5E-5)*4#40#1.5627E-5  # max conductance (umho)
                synapses['syn_{}'.format(index)].e = -10#-10

                #######################

                # Create recording vector to track current injection 
                # syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
                # syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

                # syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
                # syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

                # syn_current_synsec['v_vec_{}'.format(0)] = h.Vector()
                # syn_current_synsec['v_vec_{}'.format(0)].record(currSynSec(0.5)._ref_v)

                syn_current_synpt['v_vec_{}'.format(index)] = h.Vector()
                syn_current_synpt['v_vec_{}'.format(index)].record(currSynSec(currSynRangeVar)._ref_v)
                synpt_recording_list.append(syn_current_synpt['v_vec_{}'.format(index)])
                #print(synapses)
                # print(syn_current)
                #synAreaList.append(sec(0.5).diam)
                #syn_to_SIZ_dist = h.distance(sizSection(0.5), sec(synSegRangeVarList[i]))
                #synDistList.append(syn_to_SIZ_dist)
                # Run the simulation
                #TODO make sure printout is correct below
        
        #  syn_current_siz['v_vec_{}'.format(0)] = h.Vector()
        #         if neuron_name == "DNp01":
        #             syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.02)._ref_v)
        #         else:
        #             syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.5)._ref_v)


        syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
        if neuron_name == "DNp01":
            
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        elif neuron_name == "DNp03":
            print(sizSection.name())
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        else:
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.05)._ref_v)
        # syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

        syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
        syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)
        # print(synapses)
        # print("Round number {}/{}".format(index, syn_df.shape[0]))
        h.tstop = simTime
        h.run()

                # Save Data 
                # Save simulated section and simulation data
                #TODO: SET UP FILE AUTOWRITE SO THAT WAY NEURON NAME IS PARAM IN FILE NAME
        synInfo_df = pd.DataFrame(seriesList)
        # print(synInfo_df.head())

        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = path + '/' + MODE + '_' + datetime
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                    # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = path + '/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = 'datafiles/' + MODE + '_' + datetime
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                for recVec in synpt_recording_list:
                    write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None  

    elif MODE == "SEQUENTIAL_F":
        print('FORWARD')
        h.v_init = erev   # Resting membrane potential 

        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        seriesList = []

        synapses = {}
        syn_current_siz = {}
        syn_current_soma = {}

        ONSET = onsetTime

        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
                pass
            else:
                seriesList.append(row)

                ONSET += 1.5

                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break
                currSynRangeVar = row.loc['mappedSegRangeVar']

                synapses['syn_{}'.format(index)] = h.AlphaSynapse(currSynSec(currSynRangeVar))
                synapses['syn_{}'.format(index)].onset = ONSET       # Start time of synapse activation (ms)
                synapses['syn_{}'.format(index)].tau = 0.5        # Time to max conductance (ms)
                synapses['syn_{}'.format(index)].gmax = 1.5627E-5  # max conductance (umho)
                synapses['syn_{}'.format(index)].e = 10

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
                syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

                syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
                syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

                #print(synapses)
                # print(syn_current)
                #synAreaList.append(sec(0.5).diam)
                #syn_to_soma_dist = h.distance(somaSection(0.5), sec(currSynRangeVar))
                #print(ONSET, syn_to_soma_dist)
                #synDistList.append(syn_to_SIZ_dist)
                # Run the simulation
                #TODO make sure printout is correct below

        h.tstop = simTime
        h.run()

        synInfo_df = pd.DataFrame(seriesList)

        #Save data 
        #TODO: SET UP FILE AUTOWRITE SO THAT WAY NEURON NAME IS PARAM IN FILE NAME
        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = path + '/' + MODE + '_' + datetime
                # write_csv(outFilename + '_synapse_sites.csv', [sec])
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
                
            else:
                filename = path + '/' + filename
                # write_csv(filename + '_synapse_sites.csv', [sec])
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_siz['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = 'datafiles/' + MODE + '_' + datetime
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None  

    elif MODE == "SEQUENTIAL_F_EXP2":
        ONSET = onsetTime
        exp2Gmax = 0.00027  # (5.5E-5)*4 uS
        exp2Rise = 0.2
        exp2Decay = 1.1 
        print('FORWARD')
        h.v_init = erev   # Resting membrane potential 

        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        seriesList = []

        synapses = {}
        syn_current_siz = {}
        syn_current_soma = {}

        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':
                pass
            else:
                seriesList.append(row)

                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break
                currSynRangeVar = row.loc['mappedSegRangeVar']

                # Assign properties to the synapse
                syn = h.Exp2Syn(currSynSec(currSynRangeVar))
                syn.tau1 = exp2Rise
                syn.tau2 = exp2Decay
                syn.e = -10  # -10

                # Save the synapse in the dictionary
                synapses['syn_{}'.format(index)] = syn

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
                syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

                syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
                syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

        # # Create NetStim and NetCon for activation
        # nc = h.NetStim()
        # nc.number = 1
        # nc.noise = 0

        ncs = h.List()
        for i, syn in enumerate(synapses.values()):
            nc = h.NetStim()
            nc.number = 1
            nc.noise = 0
            nc.start = ONSET
            ncs.append(h.NetCon(nc, syn))
            ncs.object(i).weight[0] = exp2Gmax 
            ONSET += 1.5  # Increment onset time by 1.5 ms for each synapse in reverse order
        h.tstop = simTime
        h.run()


        synInfo_df = pd.DataFrame(seriesList)

        #Save data 
        #TODO: SET UP FILE AUTOWRITE SO THAT WAY NEURON NAME IS PARAM IN FILE NAME
        
        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = path + '/' + MODE + '_' + datetime
                # write_csv(outFilename + '_synapse_sites.csv', [sec])
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
                
            else:
                filename = path + '/' + filename
                # write_csv(filename + '_synapse_sites.csv', [sec])
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_siz['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = 'datafiles/' + MODE + '_' + datetime
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None

    elif MODE == "SEQUENTIAL_R_EXP2":
        ONSET = onsetTime
        exp2Gmax = 0.00027  # (5.5E-5)*4 uS
        exp2Rise = 0.2
        exp2Decay = 1.1 
        print('REVERSE')
        h.v_init = erev   # Resting membrane potential 

        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        seriesList = []

        synapses = {}
        syn_current_siz = {}
        syn_current_soma = {}

        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':
                pass
            else:
                seriesList.append(row)

                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break
                currSynRangeVar = row.loc['mappedSegRangeVar']

                # Assign properties to the synapse
                syn = h.Exp2Syn(currSynSec(currSynRangeVar))
                syn.tau1 = exp2Rise
                syn.tau2 = exp2Decay
                syn.e = -10  # -10

                # Save the synapse in the dictionary
                synapses['syn_{}'.format(index)] = syn

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
                syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

                syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
                syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

        # # Create NetStim and NetCon for activation
        # nc = h.NetStim()
        # nc.number = 1
        # nc.noise = 0

        ncs = h.List()
        for i, syn in enumerate(synapses.values()):
            nc = h.NetStim()
            nc.number = 1
            nc.noise = 0
            nc.start = ONSET
            ncs.append(h.NetCon(nc, syn))
            ncs.object(i).weight[0] = exp2Gmax 
            ONSET -= 1.5  # Increment onset time by 1.5 ms for each synapse in reverse order
        h.tstop = simTime
        h.run()

        synInfo_df = pd.DataFrame(seriesList)

        #Save data 
        #TODO: SET UP FILE AUTOWRITE SO THAT WAY NEURON NAME IS PARAM IN FILE NAME
        
        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = path + '/' + MODE + '_' + datetime
                # write_csv(outFilename + '_synapse_sites.csv', [sec])
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
                
            else:
                filename = path + '/' + filename
                # write_csv(filename + '_synapse_sites.csv', [sec])
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_siz['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = 'datafiles/' + MODE + '_' + datetime
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None

    elif MODE == "SEQUENTIAL_R":
        print('REVERSE')
        h.v_init = erev   # Resting membrane potential 

        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        seriesList = []

        synapses = {}
        syn_current_siz = {}
        syn_current_soma = {}

        ONSET = onsetTime

        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
                pass
            else:
                seriesList.append(row)

                ONSET -= 1.5

                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break
                currSynRangeVar = row.loc['mappedSegRangeVar']

                synapses['syn_{}'.format(index)] = h.AlphaSynapse(currSynSec(currSynRangeVar))
                #TODO MAKE ONSET TIME A VAR(?)
                synapses['syn_{}'.format(index)].onset = ONSET       # Start time of synapse activation (ms)
                synapses['syn_{}'.format(index)].tau = 0.5        # Time to max conductance (ms)
                synapses['syn_{}'.format(index)].gmax = 1.5627E-5  # max conductance (umho)
                synapses['syn_{}'.format(index)].e = 10

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
                syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

                syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
                syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

        h.tstop = simTime
        h.run()

        synInfo_df = pd.DataFrame(seriesList)

        #Save data 
        #TODO: SET UP FILE AUTOWRITE SO THAT WAY NEURON NAME IS PARAM IN FILE NAME
        
        if dirname is not None:
            parent_dir = "datafiles/simulationData"
            path = os.path.join(parent_dir, dirname)
            try:
                os.mkdir(path)
            except OSError as error:
                pass 

            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = path + '/' + MODE + '_' + datetime
                # write_csv(outFilename + '_synapse_sites.csv', [sec])
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
                
            else:
                filename = path + '/' + filename
                # write_csv(filename + '_synapse_sites.csv', [sec])
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_siz['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                datetime = clock.strftime("%Y%m%d-%H%M%S")
                outFilename = 'datafiles/' + MODE + '_' + datetime
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        

        synapses = None 
        syn_current = None
        tVec = None
        vVec = None  

    return tVec, vVec

def selectSynapses(subsetInfo, synMap_df):

    synSubset_df = pd.DataFrame()
    
    if subsetInfo[0] == "RAND":
        #TODO: if rand is nan, REDO GEN
        #random.seed(a=1)
        randN = subsetInfo[1]
        randSelectedSynIdxs = []
        randSelectedSyns = pd.DataFrame()

        if randN > synMap_df.shape[0]-1:
            # synSubset_df = synMap_df
            # randN = synMap_df.shape[0]-1
            synSubset_df = synMap_df
            print('TRYING TO SUBSET MORE SYNS THAN IN AVAILABLE POPULATION /// RETURNING WHOLE POPULATION AS SUBSET')
            # quit(1)
        else:
            uarray = numpy.random.choice(numpy.arange(0, synMap_df.shape[0]-1), replace=False, size=(1, randN))

            # print(type(uarray[0]), len(uarray[0]))
            # # # check that each item occurs only once
            # # print((numpy.bincount(uarray.ravel()) == 1).all())
            # # # True
            # quit(0)
            for element in uarray[0]:
                # print(element)
                randSelectedSyns = pd.concat([randSelectedSyns, synMap_df.iloc[[element]]])
            # quit(0)

            # for randSynCt in range(randN):
            #     # RANDVAR = random.randint(0, synMap_df.shape[0]-1)
            #     randSelectedSyns = pd.concat([randSelectedSyns, synMap_df.iloc[[RANDVAR]]])

            synSubset_df = randSelectedSyns

    elif subsetInfo[0] == "CLOSEST":
        closestN = subsetInfo[1]
        closestMethod = "PATH" #"EUCLID"
        STARTSYN = synMap_df.iloc[subsetInfo[2]]
        #print(STARTSYN)

        if str(STARTSYN.loc['pre']) == 'nan':#numpy.isnan(STARTSYN.loc['pre']):
                print("ERROR: RAND SYN IS NAN")
                raise IndexError
        else:
            for sec in h.allsec():
                if sec.name() == STARTSYN.loc["mappedSection"]:
                    STARTSYN_SEC = sec
                    break
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
                print("BOO")
                pass
            else:
                for sec in h.allsec():
                    if sec.name() == row.loc["mappedSection"]:
                        CURRSEC = sec
                        #print(row)
                        break
                #NOTE: IF PATH ON SAME SYN, IT IS 0, MAYBE DEFAULT TO EUCLID IF THIS IS THE CASE AS AN EXCEPTION?
                CURR_SYNDIST =  h.distance(STARTSYN_SEC(STARTSYN.loc['mappedSegRangeVar']), CURRSEC(row.loc['mappedSegRangeVar']))
                SYNDIST_DATA.append([row, CURR_SYNDIST])
                #CALC DIST FROM SYN TO OTHER SYNS USING SPECIFIED METHOD (START WITH PATH)
                #SORT, SELECT closestN closest syns and concat them to subsetdf
        SYNDIST_DATA.sort(key = lambda x: x[1])
        #print(SYNDIST_DATA[0:10])

        seriesList = []
        closestSelectedSyns = pd.DataFrame()
        for closestSynCt in range(closestN):
            #print(closestSynCt)
            closestSynX = SYNDIST_DATA[closestSynCt]
            #print(type(closestSynX[0]))
            seriesList.append(closestSynX[0])
            #print(SYNDIST_DATA[closestSynCt], SYNDIST_DATA[closestSynCt[0]])
            closestSelectedSyns = pd.concat([closestSelectedSyns, closestSynX[0]], axis=0)
            #print(closestSelectedSyns)
        synSubset_df = closestSelectedSyns
        #print(seriesList[0])
        cols = ['pre','pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z', 'type', 'mappedSection', 'mappedSegRangeVar']
        synSubset_df = pd.DataFrame(seriesList)#, axis=0)#columns=cols)
        #print(synSubset_df.head())
    elif subsetInfo[0] == "CLOSEST2SEC":
        closestN = subsetInfo[1]
        closestMethod = "PATH" #"EUCLID"

        STARTSEC = subsetInfo[2]
        STARTSEC_RV = subsetInfo[3]
        #print(STARTSYN)
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
                # print("BOO")
                pass
            else:
                for sec in h.allsec():
                    if sec.name() == row.loc["mappedSection"]:
                        CURRSEC = sec
                        #print(row)
                        break
                #NOTE: IF PATH ON SAME SYN, IT IS 0, MAYBE DEFAULT TO EUCLID IF THIS IS THE CASE AS AN EXCEPTION?
                CURR_SYNDIST =  h.distance(STARTSEC(STARTSEC_RV), CURRSEC(row.loc['mappedSegRangeVar']))
                SYNDIST_DATA.append([row, CURR_SYNDIST])
                #CALC DIST FROM SYN TO OTHER SYNS USING SPECIFIED METHOD (START WITH PATH)
                #SORT, SELECT closestN closest syns and concat them to subsetdf
        # test_row = SYNDIST_DATA[0]
        # print(test_row)
        # TEMP_VAOL = test_row[1]
        # quit(0)
        SYNDIST_DATA.sort(key = lambda x: x[1])
        # print(SYNDIST_DATA[-10:])
        # for item in SYNDIST_DATA:
        #     if item[1] == TEMP_VAOL:
        #         print(item)
        # quit(0)

        seriesList = []
        closestSelectedSyns = pd.DataFrame()
        for closestSynCt in range(closestN):
            #print(closestSynCt)
            closestSynX = SYNDIST_DATA[closestSynCt]
            #print(type(closestSynX[0]))
            seriesList.append(closestSynX[0])
            #print(SYNDIST_DATA[closestSynCt], SYNDIST_DATA[closestSynCt[0]])
            closestSelectedSyns = pd.concat([closestSelectedSyns, closestSynX[0]], axis=0)
            #print(closestSelectedSyns)
        synSubset_df = closestSelectedSyns
        #print(seriesList[0])
        cols = ['pre','pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z', 'type', 'mappedSection', 'mappedSegRangeVar']
        synSubset_df = pd.DataFrame(seriesList)#, axis=0)#columns=cols)
        #print(synSubset_df.head())
    elif subsetInfo[0] == "FARTHEST2SEC":
        farthestN = subsetInfo[1]
        STARTSEC = subsetInfo[2]
        STARTSEC_RV = subsetInfo[3]
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if str(row.loc['pre']) == 'nan':
                # print("BOO")
                pass
            else:
                for sec in h.allsec():
                    if sec.name() == row.loc["mappedSection"]:
                        CURRSEC = sec
                        break
                #NOTE: IF PATH ON SAME SYN, IT IS 0, MAYBE DEFAULT TO EUCLID IF THIS IS THE CASE AS AN EXCEPTION?
                
                CURR_SYNDIST =  h.distance(STARTSEC(STARTSEC_RV), CURRSEC(row.loc['mappedSegRangeVar']))
                # print(STARTSEC, CURRSEC, CURR_SYNDIST)

                SYNDIST_DATA.append([row, CURR_SYNDIST])

                #CALC DIST FROM SYN TO OTHER SYNS USING SPECIFIED METHOD (START WITH PATH)
                #SORT, SELECT closestN closest syns and concat them to subsetdf
        # quit(0)
        SYNDIST_DATA.sort(key=lambda x: x[1], reverse=True)
        
        for i in range(len(SYNDIST_DATA)):
            t = SYNDIST_DATA[i]
            # print(t[0].loc['mappedSection'], t[0].loc['mappedSegRangeVar'], t[1])

        seriesList = []
        farthestSelectedSyns = pd.DataFrame()
        for farthestSynCt in range(farthestN):
            #print(closestSynCt)
            farthestSynX = SYNDIST_DATA[farthestSynCt]
            #print(type(closestSynX[0]))
            seriesList.append(farthestSynX[0])
            #print(SYNDIST_DATA[closestSynCt], SYNDIST_DATA[closestSynCt[0]])
            farthestSelectedSyns = pd.concat([farthestSelectedSyns, farthestSynX[0]], axis=0)
            #print(closestSelectedSyns)
        synSubset_df = farthestSelectedSyns
        #print(seriesList[0])
        cols = ['pre','pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z', 'type', 'mappedSection', 'mappedSegRangeVar']
        synSubset_df = pd.DataFrame(seriesList)#, axis=0)#columns=cols)
        #print(synSubset_df.head())
    return synSubset_df

def activateClosestFarthest(neuron_name, erev, synMap_df, sizSection, somaSection):
    numSyn = 50
    if neuron_name == "DNp01":
        subsetInfo = ["CLOSEST2SEC", numSyn, sizSection, 0.0551186]
    else:
        subsetInfo = ["CLOSEST2SEC", numSyn, sizSection, 0.5]
    print(subsetInfo)
    SYNDF_C2S = selectSynapses(subsetInfo, synMap_df)
    print("Activating...")
    tVecRand, vVecRand = activateSynapses(SYNDF_C2S, MODE="GROUP", somaSection=somaSection, sizSection=sizSection, recordSection=sizSection,
        erev=erev, simTime=100, onsetTime=50, neuron_name=None, dirname='DNp01_final_sims_1113/close_vs_far/'+neuron_name+"_VPN_50close_vs_far_wRand", filename='closeGroup')
    
    if neuron_name == "DNp01":
        subsetInfo = ["FARTHEST2SEC", numSyn, sizSection, 0.0551186]
    else:
        subsetInfo = ["FARTHEST2SEC", numSyn, sizSection, 0.5]
    print(subsetInfo)
    SYNDF_F2S = selectSynapses(subsetInfo, synMap_df)
    print('Activating...')
    tVecRand, vVecRand = activateSynapses(SYNDF_F2S, MODE="GROUP", somaSection=somaSection, sizSection=sizSection, recordSection=sizSection,
        erev=erev, simTime=100, onsetTime=50, neuron_name=None, dirname='DNp01_final_sims_1113/close_vs_far/'+neuron_name+"_VPN_50close_vs_far_wRand", filename='farGroup')
    
    for i in range(50):
        print("Random trial", i+1, "/ 50")
        subsetInfo = ["RAND", numSyn]
        randGroup_filename = 'randGroup_' + str(i+1)
        SYNDF_R = selectSynapses(subsetInfo, synMap_df)
        tVecRand, vVecRand = activateSynapses(SYNDF_R, MODE="GROUP", somaSection=somaSection, sizSection=sizSection, recordSection=sizSection,
        erev=erev, simTime=100, onsetTime=50, neuron_name=None, dirname='DNp01_final_sims_1113/close_vs_far/'+neuron_name+"_VPN_50close_vs_far_wRand", filename=randGroup_filename)

def activateClosestFarthest_test(neuron_name, erev, synMap_df, sizSection, somaSection, numsyn= None):
    numSyn = numsyn
    numSyn_str = str(numsyn)
    if neuron_name == "DNp01":
        subsetInfo = ["CLOSEST2SEC", numSyn, sizSection, 0.0551186]
    elif neuron_name == "DNp03":
        subsetInfo = ["CLOSEST2SEC", numSyn, sizSection, 0.0551186]
    else:
        subsetInfo = ["CLOSEST2SEC", numSyn, sizSection, 0.5]
    print(subsetInfo)
    SYNDF_C2S = selectSynapses(subsetInfo, synMap_df)
    print("Activating...")
    tVecRand, vVecRand = activateSynapses(SYNDF_C2S, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, sizSection=sizSection, recordSection=sizSection,
        erev=erev, onsetTime=50, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/close_vs_far/'+neuron_name+"_VPN_"+numSyn_str+"_close_vs_far_wRand", filename='closeGroup')
    
    if neuron_name == "DNp01":
        subsetInfo = ["FARTHEST2SEC", numSyn, sizSection, 0.0551186]
    if neuron_name == "DNp03":
        subsetInfo = ["FARTHEST2SEC", numSyn, sizSection, 0.0551186]
    else:
        subsetInfo = ["FARTHEST2SEC", numSyn, sizSection, 0.5]
    print(subsetInfo)
    SYNDF_F2S = selectSynapses(subsetInfo, synMap_df)
    print('Activating...')
    tVecRand, vVecRand = activateSynapses(SYNDF_F2S, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, sizSection=sizSection, recordSection=sizSection,
        erev=erev, onsetTime= 50, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/close_vs_far/'+neuron_name+"_VPN_"+numSyn_str+"_close_vs_far_wRand", filename='farGroup')
    
    for i in range(50):
        print("Random trial", i+1, "/ 50")
        subsetInfo = ["RAND", numSyn]
        randGroup_filename = 'randGroup_' + str(i+1)
        SYNDF_R = selectSynapses(subsetInfo, synMap_df)
        tVecRand, vVecRand = activateSynapses(SYNDF_R, MODE="GROUP_EXP2", somaSection=somaSection, sizSection=sizSection, recordSection=sizSection,
        erev=erev, onsetTime= 50, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/close_vs_far/'+neuron_name+"_VPN_"+numSyn_str+"_close_vs_far_wRand", filename=randGroup_filename)


####
    
# def activateLinear250(neuron_name, erev, synMap_df, sizSection, somaSection):
#     for i in range(250):
#         print(i+1, "/ 250")

#         subsetInfo = ["RAND", i+1]
#         SYNDF = selectSynapses(subsetInfo, synMap_df)
#         ID_TEMP = SYNDF.axes[0].tolist()
#         tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname='DNp01_final_sims_1113/linear250/'+neuron_name+"_VPN_incRandSynCt250", filename=None)

def activateLinear250_test(neuron_name, erev, synMap_df, sizSection, somaSection):
    for i in range(250):
        print(i+1, "/ 250")

        subsetInfo = ["RAND", i+1]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        ID_TEMP = SYNDF.axes[0].tolist()
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, onsetTime=50, somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims_1113/linear250/'+neuron_name+"_VPN_only_250", filename=None)


# def activateAllSingle(neuron_name, erev, synMap_df, sizSection, somaSection):
#     tVec, vVec = activateSynapses(synMap_df, MODE="SINGLE", somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, onsetTime=50, neuron_name=neuron_name, dirname='DNp01_final_sims_1113/single/allSyn_synCon5pt5en5x4_'+neuron_name, filename=None)#allSynapsesSingleActivation_4siteRecord_wAxon", filename=None)


def activateAllSingleTEST(neuron_name, erev, synMap_df, sizSection, somaSection):
    #synMap_df = synMap_df.head(100), you can choose how many synapses to activate
    tVec, vVec = activateSynapses(synMap_df, MODE="SINGLE_EXP2", interval = None, num_spike= None, somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims_1113/singleTEST/allSyn_synCon5pt5en5x4_'+neuron_name, filename=None)#allSynapsesSingleActivation_4siteRecord_wAxon", filename=None)
    


# def activateRandomClose(neuron_name, erev, synMap_df, sizSection, somaSection):
#     for i in range(25):
#         print(i+1, "/ 25")

#         # subsetInfo = ["RAND", i+1]
#         subsetInfo = ["RAND", 50]
#         SYNDF = selectSynapses(subsetInfo, synMap_df)
#         ID_TEMP = SYNDF.axes[0].tolist()
#         tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+"_rand50record50", filename=None)
        
#         subsetInfo = ["CLOSEST", 50, ID_TEMP[0]]
#         SYNDFC = selectSynapses(subsetInfo, synMap_df)
#         tVecClose, vVecClose = activateSynapses(SYNDFC, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+"_close50record50", filename=None)

def activateRandomClose_TEST(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")
        synMap_df.reset_index(drop=True)
        subsetInfo = ["RAND", numSyn]
        #subsetInfo = ["RAND", numSyn]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        ID_TEMP = SYNDF.axes[0].tolist()
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_rand{numSyn}record{trials}", filename=None)
        
        subsetInfo = ["CLOSEST", numSyn, ID_TEMP[0]]
        SYNDFC = selectSynapses(subsetInfo, synMap_df)
        tVecClose, vVecClose = activateSynapses(SYNDFC, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_close{numSyn}record{trials}", filename=None)

def activateRandom_TEST(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")

        subsetInfo = ["RAND", numSyn]
        #subsetInfo = ["RAND", numSyn]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        ID_TEMP = SYNDF.axes[0].tolist()
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_rand{numSyn}record{trials}", filename=None)

def activate_close_syns(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")
        random_value = random.randint(1, 1100)
        subsetInfo = ["CLOSEST", numSyn, random_value]
        SYNDFC = selectSynapses(subsetInfo, synMap_df)
        tVecClose, vVecClose = activateSynapses(SYNDFC, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_close{numSyn}record{trials}", filename=None)

def activate_indiv_syns_rand_VPN_close(neuron_name, erev, sizSection, somaSection, syn_count =None, trials=None, simulation = None):
    if simulation == "Random":
        dir_path = f'datafiles/simulationData/{neuron_name}_rand{syn_count}record{trials}'
    elif simulation == "Close":
        dir_path = f'datafiles/simulationData/{neuron_name}_close{syn_count}record{trials}'
    elif simulation == "Partner":
        dir_path = 'datafiles/simulationData/'+neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+'_partner'

    quintList_rand = returnQuints(dir_path)
    #quintList_close = returnQuints(dir_path_close)

    for quint_ct, quint in enumerate(quintList_rand):
        synSite = quint[4]
        print(synSite)
        tVec, vVec = activateSynapses(synSite, MODE="SINGLE_EXP2", interval = None, num_spike= None, somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims_1113/{simulation}_single_syn_trials', filename=None)
        



def activate_random_syns(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")
        subsetInfo = ["RAND", numSyn]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_rand{numSyn}record{trials}", filename=None)
# def activateLinearVPN(neuron_name, erev, synMap_df, sizSection, somaSection):

#     synType_df = synMap_df[['type']]
#     synTypeList = []
#     for index, row in synType_df.iterrows():
#         #print(row)
#         synTypeList.append(row.loc['type'])
#     synTypeSet = set(synTypeList)
#     uniqueSynTypes = list(synTypeSet)
#     print(uniqueSynTypes)

#     for synType in uniqueSynTypes:
#         if isinstance(synType, float):
#             pass
#         else:
#             synSubset_df = filterSynapsesToCriteria([synType], synMap_df)

#             #CREATE DF OF ALL SYNS OF TYPE
#             for i in range(650):
#                 print(i+1, "/ 650 for", synType)
#                 # print(i, type(i))
#                 # print(i+1, type(i+1))
#                 subsetInfo = ["RAND", i+1]
#                 SYNDF = selectSynapses(subsetInfo, synSubset_df)
#                 print(SYNDF.head())
#                 tVec, vVec = activateSynapses(SYNDF, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+"_incRandSynCt650_VPNnoJitter_FIX", filename=synType+'_'+str(i+1)+'syn')

def activateLinearVPN_test(neuron_name, erev, synMap_df, sizSection, somaSection):
    synType_df = synMap_df[['type']]
    synTypeList = []
    for index, row in synType_df.iterrows():
        #print(row)
        synTypeList.append(row.loc['type'])
    synTypeSet = set(synTypeList)
    uniqueSynTypes = list(synTypeSet)
    print(uniqueSynTypes)
    for synType in uniqueSynTypes:
        if isinstance(synType, float):
            pass
        else:
            synSubset_df = filterSynapsesToCriteria([synType], synMap_df)
            synSubset_df_len = len(synSubset_df)
            #CREATE DF OF ALL SYNS OF TYPE
            for i in range(len(synSubset_df)):
                print(i+1, f"/ {synSubset_df_len} for", synType)
                # print(i, type(i))
                # print(i+1, type(i+1))
                subsetInfo = ["RAND", i+1]
                SYNDF = selectSynapses(subsetInfo, synSubset_df)
                print(SYNDF.head())
                tVec, vVec = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, onsetTime=50, somaSection=somaSection, sizSection=sizSection, recordSection=sizSection, erev=erev, simTime=75, neuron_name=None, dirname=neuron_name+"_incRandSynCt650_VPNnoJitter_FIX", filename=synType+'_'+str(i+1)+'syn')
   

def filterSynapsesToCriteria(subsetInfo, synMap_df):
    # print(subsetInfo, type(subsetInfo))
    filteredSyn_List = []
    for index, row in synMap_df.iterrows():
        if row.loc['type'] in subsetInfo:
            filteredSyn_List.append(row)
    filteredSyn_df = pd.DataFrame(filteredSyn_List)
    return filteredSyn_df

def plotMorphology(allSections_py, somaSection, sizSection, colorB):#axonList, sizSection, colorB):
    morph_shape_window = h.PlotShape(h.SectionList(allSections_py))           # Create a shape plot
    morph_shape_window.exec_menu('Show Diam')    # Show diameter of sections in shape plot
    morph_shape_window.color_all(9)

    #colorR for axon sections | colorB for test work | colorG for soma | colorK for presumed SIZ
    colorR = h.SectionList()
    # colorB = h.SectionList()
    colorG = h.SectionList()
    colorK = h.SectionList()

    colorG.append(somaSection)
    colorK.append(sizSection)
    
    # for sec in axonList:
    #     colorR.append(sec)

        # if len(sec.children()) < 1:
        #     colorB.append(sec)
        #    axonEnd = sec

    # morph_shape_window.color_list(colorR, 2)
    # morph_shape_window.color_list(colorG, 4)  
    # morph_shape_window.color_list(colorK, 1)
    morph_shape_window.color_list(colorB, 3)
    return morph_shape_window

def removeAxonalSynapses_REV(syn_df, dendList, axonList=None, tetherList=None, somaSection=None):
    # for index, row in syn_df.iterrows():
    #     print(row.loc['mappedSection'])
    synSecSeries = syn_df['mappedSection']
    tetherSynSites = syn_df[synSecSeries.str.startswith('dend_11')]
    somaSynSites = syn_df[synSecSeries.str.startswith('soma')]
    axonSynSites = syn_df[synSecSeries.str.startswith('axon')]
    dendSynSites = syn_df[synSecSeries.str.startswith('dend[')]
    return dendSynSites#, axonSynSites, tetherSynSites, somaSynSites


def activateLPLC2(neuron_name, LPLC2_synMap_df, erev, allSections_py, sizSection, somaSection, dirname=None, filename=None):
    #select neuron
    simTime = 80
    h.v_init = erev

    #define time vector
    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    seriesList = []
    
    synapseList = []

    spikeTimeCourse = [10.0550 , 10.0839, 10.2554, 10.5002, 10.6056, 10.6212,
    10.6302, 10.6419, 10.6509, 10.6594, 10.6776, 10.6837, 10.6893,
    10.6985, 10.7042, 10.7097, 10.7208, 10.7259, 10.7347, 10.7412,
    10.7466, 10.7529, 10.7576, 10.7634, 10.7698, 10.7768, 10.7819,
    10.7909, 10.8059, 10.8148, 10.8191, 10.8283, 10.8330]
    spikeTimeCourse = [x+40 for x in spikeTimeCourse]

    stim = h.NetStim()                                   #Create a netstim()
    vehicleStim = h.NetStim()


    stim.number = 1                                         #Number of Stimulations
    stim.start = 10                                       #Start time in ms

    stimDriverList = []

    for tVal in spikeTimeCourse:
        vehicleStim.number = 1                                         #Number of Stimulations
        vehicleStim.start = 0
        cStimDriver = h.NetCon(vehicleStim, stim)
        cStimDriver.delay = tVal
        cStimDriver.weight[0] = 1
        stimDriverList.append(cStimDriver)


    synct = 0
    IDX = 0
        

    h.v_init = erev   # Resting membrane potential 

    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    synapses = {}
    syn_current_siz = {}
    syn_current_soma = {}
    syn_current_synpt = {}
    
    seriesList = []

    synpt_recording_list = []


    
    # vsList.append(vs)
    # vs.play(stvec)

    for index, row in LPLC2_synMap_df.iterrows():
        if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
            pass
        else:
            seriesList.append(row)


            currSynSecStr = row.loc['mappedSection']
            for sec in h.allsec():
                if sec.name() == currSynSecStr:
                    currSynSec = sec
                    break
            currSynRangeVar = row.loc['mappedSegRangeVar']

            cSyn = h.Exp2Syn(currSynSec(currSynRangeVar))
            ncstim = h.NetCon(stim, cSyn)                           #Create a Netcon()
            ncstim.delay = 0                                     #Delay of stimulus in ms
            ncstim.weight[0] = 5.5E-5*4                                  #NetCon weight in nS.

            cSyn.tau1 = 0.5 #ms 0.2
            cSyn.tau2 = 0.5 #ms 1.1
            #cSyn.i = 0.004 #nA
            cSyn.e = -10 #mV
            synapseList.append([cSyn, row.loc["type"]])
    
    cSyn = synapseList[0]

    # vs = h.VecStim()
    # vs.play(stvec)
    # nc = h.NetCon(vs, cSyn)
    # #ncList.append(nc)
    # #nc.delay = synct*3
    # nc.weight[0] = 5.5e-5


    syn_current_synpt['v_vec_{}'.format(index)] = h.Vector()
    syn_current_synpt['v_vec_{}'.format(index)].record(currSynSec(currSynRangeVar)._ref_v)
    synpt_recording_list.append(syn_current_synpt['v_vec_{}'.format(index)])


    syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
    if neuron_name == "DNp01":
        
        syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
    elif neuron_name == "DNp03":
        # print(sizSection.name())
        syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
    else:
        syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.5)._ref_v)
    # syn_current_siz['v_vec_{}'.format(index)].record(recordSection(0.5)._ref_v)

    syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
    syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

    TEMPSOMA = h.Vector()
    TEMPSIZ = h.Vector()
    TEMPSOMA.record(somaSection(0.5)._ref_v)
    TEMPSIZ.record(sizSection(0.5)._ref_v)
    # print(synapses)
    # print("Round number {}/{}".format(index, syn_df.shape[0]))

    # h.tstop = simTime
    # h.run()
    # h.finitialize(-72.5 * mV)
    h.finitialize(erev * mV)
    h.continuerun(simTime * ms)

    synInfo_df = pd.DataFrame(seriesList)
    # print(synInfo_df.head())
    # print(list(syn_current_siz))

    plt.plot(t_vec, TEMPSOMA, 'g-')
    plt.plot(t_vec, TEMPSIZ, 'r-')
    plt.legend(['SOMA', 'SIZ'])
    plt.show()


    MODE = 'TESTLPLC2'

    # if dirname is not None:
    #     parent_dir = "datafiles/simulationData"
    #     path = os.path.join(parent_dir, dirname)
    #     try:
    #         os.mkdir(path)
    #     except OSError as error:
    #         pass 

    #     if filename is None:
    #         datetime = clock.strftime("%Y%m%d-%H%M%S")
    #         outFilename = path + '/' + MODE + '_' + datetime
    #         write_csv(outFilename + '_time.csv',t_vec)        
    #         write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
    #             # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
    #     else:
    #         filename = path + '/' + filename
    #         write_csv(filename + '_time.csv',t_vec)        
    #         write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
    #         # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
    # else:
    #     if filename is None:
    #         datetime = clock.strftime("%Y%m%d-%H%M%S")
    #         outFilename = 'datafiles/' + MODE + '_' + datetime
    #         write_csv(outFilename + '_time.csv',t_vec)        
    #         write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
    #         # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
    #     else:
    #         filename = 'datafiles/' + filename
    #         write_csv(filename + '_time.csv',t_vec)        
    #         write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
    #         # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

    # synapses = None 
    # syn_current = None
    # tVec = None
    # vVec = None  

#def activateFakeSpikeSyn(neuron_name, erev, synMap_df, sizSection, somaSection, electrodeSec):

    h.v_init = erev   # Resting membrane potential 

    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    fakeSyn = h.AlphaSynapse(sizSection(0.0551186))
    fakeSyn.onset = 10       # Start time of synapse activation (ms)
    fakeSyn.tau = 0.1        # Time to max conductance (ms)
    fakeSyn.gmax = (5.5E-5)*168000#40 #1.5627E-5  # max conductance (umho)  |  https://github.com/txzliu/connectome_to_function/blob/main/attrs_per_mEPSP.py
    fakeSyn.e = -10#-10

    fakeSynSIZ = h.Vector()
    fakeSynSIZ.record(sizSection(0.0551186)._ref_v)

    fakeSynSoma = h.Vector()
    fakeSynSoma.record(somaSection(0.5)._ref_v)

    fakeSynEsec = h.Vector()
    fakeSynEsec.record(electrodeSec(0.5)._ref_v)

    h.finitialize(erev * mV)
    h.continuerun(50 * ms)

    plt.plot(t_vec, fakeSynSIZ, 'r')
    plt.plot(t_vec, fakeSynSoma, 'g')
    plt.axhline(y=0, color='black', linestyle='-')
    plt.axhline(y=erev+10, color='black', linestyle='--')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    # plt.plot(t_vec, fakeSynEsec, 'b')
    plt.legend(['SIZ', 'Soma', '0 mV', '10mV from resting'])
    maxSIZ_depol = fakeSynSIZ.max()
    print(erev)
    print(fakeSynSoma.max() - erev)
    

    plt.title('Depolarization of DNp01 soma in response to SIZ "fake spike", SIZ depolarized to {} mV, SIZ depolarized to {} mV above resting'.format(round(maxSIZ_depol, 3), round(fakeSynSoma.max() - erev, 3)))
   

    plt.savefig('DNp01_finalSIZfakeSpike.svgz', format='svgz')
    plt.show()

def activateLPLC2_TEST(neuron_name, LPLC2_synMap_df, erev, allSections_py, sizSection, somaSection, dirname=None, filename=None, onsetTime = None, randomseed = None):
       # Set a seed for repeatability
    numpy.random.seed(randomseed)  # Change the seed value as needed

    # Randomly select a neuron based on the "pre" column
    selected_neuron = numpy.random.choice(LPLC2_synMap_df['pre'].dropna().unique())

    # Filter the DataFrame to get only the synapses of the selected neuron
    selected_neuron_df = LPLC2_synMap_df[LPLC2_synMap_df['pre'] == selected_neuron]
    print(selected_neuron)

    # Initialize simulation parameters and variables
    simTime = 250
    h.v_init = erev
    t_vec = h.Vector()
    t_vec.record(h._ref_t)
    seriesList = []

    # Stimulation Setup
    spikeTimeCourse = [10.0550, 10.0839, 10.2554, 10.5002, 10.6056, 10.6212,
                    10.6302, 10.6419, 10.6509, 10.6594, 10.6776, 10.6837, 10.6893,
                    10.6985, 10.7042, 10.7097, 10.7208, 10.7259, 10.7347, 10.7412,
                    10.7466, 10.7529, 10.7576, 10.7634, 10.7698, 10.7768, 10.7819,
                    10.7909, 10.8059, 10.8148, 10.8191, 10.8283, 10.8330]
    spikeTimeCourse = [x + 40 for x in spikeTimeCourse]

    stim = h.NetStim()
    vehicleStim = h.NetStim()

    stim.number = 1
    stim.start = spikeTimeCourse[0]

    stimDriverList = []

    for tVal in spikeTimeCourse:
        vehicleStim.number = 1
        vehicleStim.start = 0
        cStimDriver = h.NetCon(vehicleStim, stim)
        cStimDriver.delay = tVal
        cStimDriver.weight[0] = 1
        stimDriverList.append(cStimDriver)

    # Synapse setup for the selected neuron
    synapses = {}
    syn_current_siz = {}
    syn_current_soma = {}
    syn_current_synpt = {}
    synpt_recording_list = []

    # Synapse parameters
    exp2Gmax = 0.00027  # (5.5E-5)*4 uS
    exp2Rise = 0.2
    exp2Decay = 1.1

    for index, row in selected_neuron_df.iterrows():
        if not pd.isna(row['pre']):
            seriesList.append(row)

            currSynSecStr = row['mappedSection']

            for sec in h.allsec():
                if sec.name() == currSynSecStr:
                    currSynSec = sec
                    break

            currSynRangeVar = row['mappedSegRangeVar']

            # Assign properties to the synapse
            syn = h.Exp2Syn(currSynSec(currSynRangeVar))
            syn.tau1 = exp2Rise
            syn.tau2 = exp2Decay
            syn.e = -10  # -10

            # Save the synapse in the dictionary
            synapses['syn_{}'.format(index)] = syn

            # Create recording vector to track voltage at the synapse location
            syn_current_synpt['v_vec_{}'.format(index)] = h.Vector()
            syn_current_synpt['v_vec_{}'.format(index)].record(currSynSec(currSynRangeVar)._ref_v)
            synpt_recording_list.append(syn_current_synpt['v_vec_{}'.format(index)])

    # Create artificial spiker
    spiker = h.NetStim()
    spiker.number = len(spikeTimeCourse)
    spiker.start = spikeTimeCourse[0]

    netcon_list = []
    for synapse in synapses.values():
        nc = h.NetCon(spiker, synapse)
        nc.weight[0] = exp2Gmax  # Adjust weight as needed
        netcon_list.append(nc)

    TEMPSOMA = h.Vector()
    TEMPSIZ = h.Vector()
    TEMPSOMA.record(somaSection(0.5)._ref_v)
    TEMPSIZ.record(sizSection(0.5)._ref_v)
    # print(synapses)
    # print("Round number {}/{}".format(index, syn_df.shape[0]))

    # h.tstop = simTime
    # h.run()
    # h.finitialize(-72.5 * mV)
    h.finitialize(erev * mV)
    h.continuerun(simTime * ms)

    synInfo_df = pd.DataFrame(seriesList)
    # print(synInfo_df.head())
    # print(list(syn_current_siz))

    # # Plotting
    # plt.plot(t_vec, TEMPSOMA, 'g-')
    # plt.plot(t_vec, TEMPSIZ, 'r-')
    # plt.legend(['SOMA', 'SIZ'])
    # plt.show()
    # plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    print(syn_current_synpt.items())
    for index in syn_current_synpt:
        vector = syn_current_synpt[index]
        plt.plot(t_vec, vector, label=f'Synapse {index}')

    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage at Synapse')
    plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1), ncol=2)
    plt.show()

    # plt.xlabel('Time (ms)')
    # plt.ylabel('Voltage at Synapse')
    # plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1), ncol=2)
    # plt.show()


    MODE = 'TESTLPLC2'

    # if dirname is not None:
    #     parent_dir = "datafiles/simulationData"
    #     path = os.path.join(parent_dir, dirname)
    #     try:
    #         os.mkdir(path)
    #     except OSError as error:
    #         pass 

    #     if filename is None:
    #         datetime = clock.strftime("%Y%m%d-%H%M%S")
    #         outFilename = path + '/' + MODE + '_' + datetime
    #         write_csv(outFilename + '_time.csv',t_vec)        
    #         write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
    #             # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
    #     else:
    #         filename = path + '/' + filename
    #         write_csv(filename + '_time.csv',t_vec)        
    #         write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
    #         # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
    # else:
    #     if filename is None:
    #         datetime = clock.strftime("%Y%m%d-%H%M%S")
    #         outFilename = 'datafiles/' + MODE + '_' + datetime
    #         write_csv(outFilename + '_time.csv',t_vec)        
    #         write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
    #         # write_csv(outFilename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
    #     else:
    #         filename = 'datafiles/' + filename
    #         write_csv(filename + '_time.csv',t_vec)        
    #         write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
    #         write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
    #         for recVec in synpt_recording_list:
    #             write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
    #         # write_csv(filename + '_synaptic_currents_synpt.csv', syn_current_synpt['v_vec_{}'.format(index)])
    #         synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

    # synapses = None 
    # syn_current = None
    # tVec = None
    # vVec = None  

def activateFakeSpikeSyn(neuron_name, erev, synMap_df, sizSection, somaSection, electrodeSec):

    h.v_init = erev   # Resting membrane potential 

    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    fakeSyn = h.AlphaSynapse(sizSection(0.0551186))
    fakeSyn.onset = 10       # Start time of synapse activation (ms)
    fakeSyn.tau = 0.1        # Time to max conductance (ms)
    fakeSyn.gmax = (5.5E-5)*168000#40 #1.5627E-5  # max conductance (umho)  |  https://github.com/txzliu/connectome_to_function/blob/main/attrs_per_mEPSP.py
    fakeSyn.e = -10#-10

    fakeSynSIZ = h.Vector()
    fakeSynSIZ.record(sizSection(0.0551186)._ref_v)

    fakeSynSoma = h.Vector()
    fakeSynSoma.record(somaSection(0.5)._ref_v)

    fakeSynEsec = h.Vector()
    fakeSynEsec.record(electrodeSec(0.5)._ref_v)

    h.finitialize(erev * mV)
    h.continuerun(50 * ms)

    plt.plot(t_vec, fakeSynSIZ, 'r')
    plt.plot(t_vec, fakeSynSoma, 'g')
    plt.axhline(y=0, color='black', linestyle='-')
    plt.axhline(y=erev+10, color='black', linestyle='--')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    # plt.plot(t_vec, fakeSynEsec, 'b')
    plt.legend(['SIZ', 'Soma', '0 mV', '10mV from resting'])
    maxSIZ_depol = fakeSynSIZ.max()
    print(erev)
    print(fakeSynSoma.max() - erev)
    

    plt.title('Depolarization of DNp01 soma in response to SIZ "fake spike", SIZ depolarized to {} mV, SIZ depolarized to {} mV above resting'.format(round(maxSIZ_depol, 3), round(fakeSynSoma.max() - erev, 3)))
   

    plt.savefig('DNp01_finalSIZfakeSpike.svgz', format='svgz')
    plt.show()


def activateCurrentInjectionSIZ(neuron_name, erev, sizSection, somaSection, electrodeSec, continueRun=150, current=None, injDur=None, delay= None):

    h.v_init = erev   # Resting membrane potential 

    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    stimobj = h.IClamp(sizSection(0.5))
    stimobj.delay = delay
    stimobj.dur = injDur
    stimobj.amp = current
    stimobj.i = 0

    SIZsec = h.Vector()
    SIZsec.record(sizSection(0.0551186)._ref_v)
    # vfakeSynSIZ = fakeSynSIZ.to_python()
    # vfakeSynSIZ = numpy.array(vfakeSynSIZ)

    Somasec = h.Vector()
    Somasec.record(somaSection(0.5)._ref_v)

    Electrodesec = h.Vector()
    Electrodesec.record(electrodeSec(0.5)._ref_v)

    h.finitialize(erev * mV)
    h.continuerun(continueRun * ms)


    plt.plot(t_vec, SIZsec, 'r')
    plt.plot(t_vec, Somasec, 'g')
    plt.axhline(y=0, color='black', linestyle='-')
    plt.axhline(y=erev+10, color='black', linestyle='--')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    plt.plot(t_vec, Electrodesec, 'b')
    plt.legend(['SIZ', 'Soma', '0 mV', '10mV from resting', "ElectrodeSec"])
    maxSIZ_depol = SIZsec.max()
    print(erev)
    print(Somasec.max() - erev)
    maxSIZ_depol = SIZsec.max()
    print(f'Membrane potential at max SIZ depol: ' + str(maxSIZ_depol))
    print(f'Max SIZ depol: ' + str(maxSIZ_depol -erev))
    print(f'Membrane potential change at Soma: ' + str(Somasec.max() - erev))

    plt.title('Depolarization of {} soma in response to SIZ current injection, SIZ depolarized to {} mV, SOMA depolarized to {} mV above resting'.format(neuron_name, round(maxSIZ_depol, 3), round(Somasec.max() - erev, 3)))
    # plt.savefig('DNp01_finalSIZfakeSpike.svgz', format='svgz')
    plt.show()

def activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=150, current=None, injDur=None, delay= None):

    h.v_init = erev   # Resting membrane potential 

    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    stimobj = h.IClamp(somaSection(0.5))
    stimobj.delay = delay
    stimobj.dur = injDur
    stimobj.amp = current
    stimobj.i = 0

    SIZsec = h.Vector()
    SIZsec.record(sizSection(0.0551186)._ref_v)

    Somasec = h.Vector()
    Somasec.record(somaSection(0.5)._ref_v)

    Electrodesec = h.Vector()
    Electrodesec.record(electrodeSec(0.5)._ref_v)

    h.finitialize(erev * mV)
    h.continuerun(continueRun * ms)


    plt.plot(t_vec, SIZsec, 'r')
    plt.plot(t_vec, Somasec, 'g')
    plt.axhline(y=0, color='black', linestyle='-')
    plt.axhline(y=erev+10, color='black', linestyle='--')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    plt.plot(t_vec, Electrodesec, 'b')
    plt.legend(['SIZ', 'Soma', '0 mV', '10mV from resting', "ElectrodeSec"])
    maxSIZ_depol = SIZsec.max()
    print(f'Membrane potential at max SIZ depol: ' + str(maxSIZ_depol))
    print(f'Max SIZ depol: ' + str(maxSIZ_depol -erev))
    print(f'Membrane potential change at Soma: ' + str(Somasec.max() - erev))
    print(f'Ratio of Soma max amplitude/SIZ max amplitude: ' + str((Somasec.max() - erev)/(maxSIZ_depol -erev )) )
    

    plt.title('Depolarization of DNp01 SIZ in response to SOMA current injection, SIZ depolarized to {} mV, SIZ depolarized to {} mV above resting'.format(round(maxSIZ_depol, 3), round(SIZsec.max() - erev, 3)))
    plt.show()

def activate_VPNs_by_receptive_field(neuron_name, erev, synMap_df, sizSection, somaSection, numspike = None, interval = None, region = None, VPN_subtype = None):

    if region == "Dorsal":
        if VPN_subtype == "LC4":
            ids = ["720575940638496720", "720575940640612480", "720575940626916936", "720575940622939836", "720575940630484495"]
        elif VPN_subtype == "LPLC2":
            ids = ["720575940617289782", "720575940622680088", "720575940630831185", "720575940613470947", "720575940619398127"]
        elif VPN_subtype == "LC22":
            ids = ["720575940623299111", "720575940633150751", "720575940615071675", "720575940650869238", "720575940616040153"]
        elif VPN_subtype == "LPLC1":
            ids = ["720575940616881821", "720575940613304803", "720575940636983518", "720575940645750324", "720575940610902442"]
        elif VPN_subtype == "LPLC4":
            ids = ["720575940646414132", "720575940630303870", "720575940638338281", "720575940615081366", "720575940633176749"]
        elif VPN_subtype == "LC4 and LPLC2":
            ids = ["720575940617289782", "720575940622680088", "720575940630831185", "720575940613470947", "720575940619398127", "720575940638496720", "720575940640612480", "720575940626916936", "720575940622939836", "720575940630484495"]
        synMap_df = synMap_df[synMap_df['pre'].isin(ids)]
        print("Activating " + str(len(synMap_df)) + " number of synpases")
    elif region == "Ventral":
        if VPN_subtype == "LC4":
            ids = ["720575940615575007", "720575940613260959", "720575940628064260", "720575940649229433", "720575940617266722"]
        elif VPN_subtype == "LC22":
            ids = ["720575940623493373", "720575940618785995", "720575940632549346", "720575940620132999", "720575940628063387"]
        elif VPN_subtype == "LPLC2":
            ids = ["720575940622093546", "720575940622009284", "720575940623925774", "720575940630731847", "720575940619263937"]
        elif VPN_subtype == "LPLC1":
            ids = ["720575940628716294", "720575940625369884", "720575940646708788", "720575940627806992", "720575940632086573"]
        elif VPN_subtype == "LPLC4":
            ids = ["720575940609063236", "720575940616593714", "720575940635213431", "720575940623383432", "720575940634839964"]
        elif VPN_subtype == "LC4 and LPLC2":
            ids = ["720575940622093546", "720575940622009284", "720575940623925774", "720575940630731847", "720575940619263937","720575940615575007", "720575940613260959", "720575940628064260", "720575940649229433", "720575940617266722"]
        synMap_df = synMap_df[synMap_df['pre'].isin(ids)]
        print("Activating " + str(len(synMap_df)) + " number of synpases")
    elif region == "Anterior":
        if VPN_subtype == "LC4":
            ids = []
        elif VPN_subtype == "LPLC2":
            ids = []
        synMap_df = synMap_df[synMap_df['pre'].isin(ids)]
    elif region == "Posterior":
        if VPN_subtype == "LC4":
            ids = []
        elif VPN_subtype == "LPLC2":
            ids = []
        synMap_df = synMap_df[synMap_df['pre'].isin(ids)]

    if  numspike == None and interval == None:
        tVec, vVec = activateSynapses(synMap_df, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=None, sizSection=sizSection, erev=erev, simTime=100, onsetTime=50, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/Receptive_field_act/'+neuron_name+"_single_spike_"+region+"_"+VPN_subtype, filename=VPN_subtype)
    else:
       numspike_str = str(numspike)
       interval_str = str(interval)
       tVec, vVec = activateSynapses(synMap_df, MODE="GROUP_EXP2_interval", somaSection=somaSection, num_spike= numspike, interval= interval, recordSection=None, sizSection=sizSection, erev=erev, simTime=500, onsetTime=50, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/Receptive_field_act/'+neuron_name+"_"+numspike_str+"_spikes_"+region+interval_str+"_ms_interval", filename=VPN_subtype)

def activate_VPNs_by_receptive_field2(neuron_name, erev, synMap_df, sizSection, somaSection, numspike = None, interval = None, region = None, VPN_subtype = None):

    if region == "Dorsal":
        if VPN_subtype == "LC4":
            ids = ["720575940638496720", "720575940640612480", "720575940626916936", "720575940622939836", "720575940630484495"]
        elif VPN_subtype == "LPLC2":
            ids = ["720575940617289782", "720575940622680088", "720575940630831185", "720575940613470947", "720575940619398127"]
        elif VPN_subtype == "LC22":
            ids = ["720575940623299111", "720575940633150751", "720575940615071675", "720575940650869238", "720575940616040153"]
        elif VPN_subtype == "LPLC1":
            ids = ["720575940616881821", "720575940613304803", "720575940636983518", "720575940645750324", "720575940610902442"]
        elif VPN_subtype == "LPLC4":
            ids = ["720575940646414132", "720575940630303870", "720575940638338281", "720575940615081366", "720575940633176749"]
        elif VPN_subtype == "LC4 and LPLC2":
            ids = ["720575940617289782", "720575940622680088", "720575940630831185", "720575940613470947", "720575940619398127", "720575940638496720", "720575940640612480", "720575940626916936", "720575940622939836", "720575940630484495"]
        elif VPN_subtype == "LC4 and LPLC2_1":
             ids = ["720575940638496720", "720575940617289782"]
        elif VPN_subtype == "LC4_1":
             ids = ["720575940638496720",]
        elif VPN_subtype == "LPLC2_1":
             ids = ["720575940617289782"]
        synMap_df = synMap_df[synMap_df['pre'].isin(ids)]
        print("Activating " + str(len(synMap_df)) + " number of synpases")
    elif region == "Ventral":
        if VPN_subtype == "LC4":
            ids = ["720575940615575007", "720575940613260959", "720575940628064260", "720575940649229433", "720575940617266722"]
        elif VPN_subtype == "LC22":
            ids = ["720575940623493373", "720575940618785995", "720575940632549346", "720575940620132999", "720575940628063387"]
        elif VPN_subtype == "LPLC2":
            ids = ["720575940622093546", "720575940622009284", "720575940623925774", "720575940630731847", "720575940619263937"]
        elif VPN_subtype == "LPLC1":
            ids = ["720575940628716294", "720575940625369884", "720575940646708788", "720575940627806992", "720575940632086573"]
        elif VPN_subtype == "LPLC4":
            ids = ["720575940609063236", "720575940616593714", "720575940635213431", "720575940623383432", "720575940634839964"]
        elif VPN_subtype == "LC4 and LPLC2":
            ids = ["720575940622093546", "720575940622009284", "720575940623925774", "720575940630731847", "720575940619263937","720575940615575007", "720575940613260959", "720575940628064260", "720575940649229433", "720575940617266722"]
        synMap_df = synMap_df[synMap_df['pre'].isin(ids)]
        print("Activating " + str(len(synMap_df)) + " number of synpases")

    if  numspike == None and interval == None:
        tVec, vVec = activateSynapses(synMap_df, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=None, sizSection=sizSection, erev=erev, simTime=100, onsetTime=50, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/Receptive_field_act/'+neuron_name+"_single_spike_"+region+"_"+VPN_subtype, filename=VPN_subtype)


def VPN_paired_activations_by_neurons(neuron_name, synMap_df, sizSection, somaSection, erev, population = None):
    if population == None:
        synMap_VPN_df = synMap_df
    elif population =="LC4":
        synMap_VPN_df =  filterSynapsesToCriteria(['LC4'], synMap_df)
    elif population =="LC6":
        synMap_VPN_df =  filterSynapsesToCriteria(['LC6'], synMap_df)
    elif population =="LC22":
        synMap_VPN_df =  filterSynapsesToCriteria(['LC22'], synMap_df)
    elif population =="LPLC1":
        synMap_VPN_df =  filterSynapsesToCriteria(['LPLC1'], synMap_df)
    elif population =="LPLC2":
        synMap_VPN_df =  filterSynapsesToCriteria(['LPLC2'], synMap_df)
    elif population =="LPLC4":
        synMap_VPN_df =  filterSynapsesToCriteria(['LPLC4'], synMap_df)

# Simulation for activating all synapses that belong to a given VPN
    synID_list = []
    for index, row in synMap_VPN_df.iterrows():
        synID_list.append(row['pre'])

    synID_countDic = {}
    for id in synID_list:
        synID_countDic[id] = synID_countDic.get(id, 0) + 1

    filtered_ids = []
    for id, count in synID_countDic.items():
        if 1 <= count <= 50:
            filtered_ids.append(id)


    for num_selections in range(2, 20):
        for _ in range(5):  # Repeat 5 times
            selected_ids = random.sample(filtered_ids, min(len(filtered_ids), num_selections))
            filtered_df = synMap_df[synMap_df['pre'].isin(selected_ids)]
            syn_count = len(filtered_df)
            unique_ids = len(filtered_df['pre'].unique())
            print(filtered_df)
            print(f'neurons: {unique_ids} synapses: {syn_count}')
            if population == None:
                tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2", somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+f'_final_sims_1113/paired_VPN/'+neuron_name+"_partners", filename=f'syns_{syn_count}_neurons_{unique_ids}')
            else:
                VPN = population
                tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2", somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims_1113/paired_VPN/{neuron_name}_{VPN}_partners', filename=f'syns_{syn_count}_neurons_{unique_ids}')

def VPN_paired_activations(neuron_name, synMap_df, sizSection, somaSection, erev, shuffled = False):
    if shuffled == False:
        synMap_VPN_df = synMap_df
    # Simulation for activating all synapses that belong to a given VPN
        synID_list = []
        for index, row in synMap_VPN_df.iterrows():
            synID_list.append(row['pre'])

        synID_countDic = {}
        for id in synID_list:
            synID_countDic[id] = synID_countDic.get(id, 0) + 1

        filtered_ids = []
        for id, count in synID_countDic.items():
            if 1 <= count <= 40:
                filtered_ids.append(id)


        for id in filtered_ids:
            count = synID_countDic.get(id)
            print(id, count)
            filtered_df = synMap_df[synMap_df['pre'] == id]
            # print(filtered_df)
            tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2", somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_partner", filename='syns_'+str(count)+'_'+str(id))
    
    elif shuffled == True:
        synMap_VPN_df = synMap_df
    # Simulation for activating all synapses that belong to a given VPN
        synID_list = []
        for index, row in synMap_VPN_df.iterrows():
            synID_list.append(row['pre'])

        synID_countDic = {}
        for id in synID_list:
            synID_countDic[id] = synID_countDic.get(id, 0) + 1

        filtered_ids = []
        for id, count in synID_countDic.items():
            if 1 <= count <= 40:
                filtered_ids.append(id)


        for id in filtered_ids:
            count = synID_countDic.get(id)
            print(id, count)
            filtered_df = synMap_df[synMap_df['pre'] == id]
            # print(filtered_df)
            tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2", somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_partner_shuffled", filename='syns_'+str(count)+'_'+str(id))

def VPN_random_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, erev, numsyn_lim = None):
    numsyn_lim = numsyn_lim
    #Simulation for activating a random subset of synapses synapses that belong to any VPNS
    for j in range(1, numsyn_lim):
        subsetInfo = ["RAND", j]
        for i in range(10):

            SYNDF = selectSynapses(subsetInfo, synMap_VPN_df)
            tVec, vVec = activateSynapses(SYNDF, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand", filename='syns_'+str(j)+'_trial_'+str(i))
def VPN_random_subset_within_type_activations(neuron_name, synMap_df, sizSection, somaSection, erev):
#Simulation for selecting random synapses from a specific population. It is a repeat of the above simulations, but only selecting for LC4 or LPLC2 synapses. 
    #First subset the datframe to LC4 or LPLC2 or any VPN
    LPLC2_synMap_df = filterSynapsesToCriteria(['LPLC2'], synMap_df)
    LC4_synMap_df = filterSynapsesToCriteria(['LC4'], synMap_df)
    LC6_synMap_df = filterSynapsesToCriteria(['LC6'], synMap_df)
    LC22_synMap_df = filterSynapsesToCriteria(['LC22'], synMap_df)
    LPLC1_synMap_df = filterSynapsesToCriteria(['LPLC1'], synMap_df)
    LPLC4_synMap_df = filterSynapsesToCriteria(['LPLC4'], synMap_df)
    # for j in range(1,40):
    #     subsetInfo = ["RAND", j]
    #     for i in range(10):
    #         SYNDF_LC4 = selectSynapses(subsetInfo, LC4_synMap_df)
    #         tVec, vVec = activateSynapses(SYNDF_LC4, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand_LC4", filename='syns_'+str(j)+'_trial_'+str(i))
    # for j in range(1,40):
    #     subsetInfo = ["RAND", j]
    #     for i in range(10):
    #         SYNDF_LC6 = selectSynapses(subsetInfo, LC6_synMap_df)
    #         tVec, vVec = activateSynapses(SYNDF_LC6, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand_LC6", filename='syns_'+str(j)+'_trial_'+str(i))    
    for j in range(1, 40):
        subsetInfo = ["RAND", j]
        for i in range(10):
            SYNDF_LPLC2 = selectSynapses(subsetInfo, LPLC2_synMap_df)
            tVec, vVec = activateSynapses(SYNDF_LPLC2, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand_LPLC2", filename='syns_'+str(j)+'_trial_'+str(i))
    for j in range(1, 40):
        subsetInfo = ["RAND", j]
        for i in range(10):
            SYNDF_LC22 = selectSynapses(subsetInfo, LC22_synMap_df)
            tVec, vVec = activateSynapses(SYNDF_LC22, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand_LC22", filename='syns_'+str(j)+'_trial_'+str(i))
    for j in range(1, 40):
        subsetInfo = ["RAND", j]
        for i in range(10):
            SYNDF_LPLC1 = selectSynapses(subsetInfo, LPLC1_synMap_df)
            tVec, vVec = activateSynapses(SYNDF_LPLC1, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand_LPLC1", filename='syns_'+str(j)+'_trial_'+str(i))
    for j in range(1, 40):
        subsetInfo = ["RAND", j]
        for i in range(10):
            SYNDF_LPLC4 = selectSynapses(subsetInfo, LPLC4_synMap_df)
            tVec, vVec = activateSynapses(SYNDF_LPLC4, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_rand_LPLC4", filename='syns_'+str(j)+'_trial_'+str(i))

def shuffle_data(neuron_name, VPN, synSite_df, int_seed= None):
    synMap_df = synSite_df
    int_seed = int_seed
        
    if VPN == "LC4":
        mask = synMap_df['type'] == 'LC4'
    elif VPN == "LC6":
        mask = synMap_df['type'] == 'LC6'
    elif VPN == "LC22":
        mask = synMap_df['type'] == 'LC22'
    elif VPN == "LPLC1":
        mask = synMap_df['type'] == 'LPLC1'
    elif VPN == "LPLC2":
        mask = synMap_df['type'] == 'LPLC2'
    elif VPN == "LPLC4":
        mask = synMap_df['type'] == 'LPLC4'
    else:
        raise ValueError("Invalid VPN value")
        

    synMap_df = synMap_df[mask]

    # reset the index of the DataFrame to ensure that the index values match the row numbers in the DataFrame
    synMap_df = synMap_df.reset_index(drop=True)

    # drop rows with missing values in post_x, post_y, post_z columns
    synMap_df = synMap_df.dropna(subset=['post_x', 'post_y', 'post_z'])
    all_ids = synMap_df['pre'].tolist()

    groups = synMap_df.groupby(['post_x', 'post_y', 'post_z'])

# ensure that the unique ids are sorted to ensure reproducibility
    # set a random seed for reproducibility
    numpy.random.seed(int_seed)
    
    shuffled_ids = numpy.random.permutation(all_ids)

    synMap_df['neuron_id'] = shuffled_ids

    shuffled = synMap_df.drop(columns=['pre'])

    # Rename the 'neuron_id' column to 'pre'
    shuffled = shuffled.rename(columns={'neuron_id': 'pre'})

    return synMap_df, shuffled

def str2sec(sec_as_str):
    for sec in h.allsec():
        if sec.name() == sec_as_str:
            sec_as_sec = sec
            break
    return sec_as_sec

def calculate_syn_to_siz_dist(row, sizSection):
    # for index, row in synMap_VPN_df.iterrows():
        # print(row)
    synSec = str2sec(row['mappedSection'])
    syn_to_siz_dist = h.distance(synSec(row['mappedSegRangeVar']), sizSection(0.05))#sizSection(0.0551186))
    return syn_to_siz_dist

def calculate_syn_to_soma_dist(row, somaSection):
    # for index, row in synMap_VPN_df.iterrows():
        # print(row)
    synSec = str2sec(row['mappedSection'])
    syn_to_soma_dist = h.distance(synSec(row['mappedSegRangeVar']), somaSection(0.05))#sizSection(0.0551186))
    return syn_to_soma_dist


def activate_close_syns_by_avg_dist_to_SIZ(neuron_name, sizSection, erev, somaSection, synMap_df, syncount, trials):

    synMap_VPN_df = synMap_df
    if isinstance(synMap_VPN_df, tuple):
        synMap_VPN_df = synMap_VPN_df[0]

    synMap_VPN_df['synToSIZ_dist'] = synMap_VPN_df.apply(calculate_syn_to_siz_dist, axis=1, args=(sizSection,))

    mean_synToSIZ_dist = synMap_VPN_df['synToSIZ_dist'].mean()
    target_average_range = (mean_synToSIZ_dist - 5, mean_synToSIZ_dist + 5) 
    print(mean_synToSIZ_dist)
    print(target_average_range)
    print(synMap_VPN_df['synToSIZ_dist'].min())
    print(synMap_VPN_df['synToSIZ_dist'].max())
    selected_rows = synMap_VPN_df[
    (synMap_VPN_df['synToSIZ_dist'] >= target_average_range[0]) &
    (synMap_VPN_df['synToSIZ_dist'] <= target_average_range[1])
]

    for trial in range(1, trials):  # 5 trials
        # Loop through lengths
        for numSyn in range(1, syncount):  # 1 to 26
            # Remove any duplicate rows and then randomly select rows for the current trial and length
            subsetInfo = ["RAND", 50]
            selected_rows = selected_rows.reset_index(drop=True)
            SYNDF = selectSynapses(subsetInfo, selected_rows)
            SYNDF = SYNDF.reset_index(drop=True)
            x = random.randint(1,45)
            ID_TEMP = SYNDF.axes[0].tolist()
            subsetInfo = ["CLOSEST", numSyn, ID_TEMP[x]]
            SYNDF2 = selectSynapses(subsetInfo, selected_rows)
            # selected_subset = selected_rows.drop_duplicates().sample(n=length, replace=False)

            # Check if the number of selected rows is equal to the specified length
            if len(SYNDF2) == numSyn:
                tVecClose, vVecClose = activateSynapses(SYNDF2, MODE="GROUP_EXP2", interval=None, num_spike=None, somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+f"_rand_by_average_syn_dist_to_SIZ", filename=None)
                print(SYNDF2)
                print(SYNDF2['synToSIZ_dist'].mean())
                print(f"Syncount {numSyn}, Trial {trial} done")
            else:
                print(f"Warning: Length of selected rows ({len(SYNDF2)}) does not match the specified length ({numSyn}).")

def avg_dist_to_siz_per_VPN(neuron_name, sizSection, somaSection, synMap_df):
    synMap_VPN_df = synMap_df
    if isinstance(synMap_VPN_df, tuple):
        synMap_VPN_df = synMap_VPN_df[0]

    synMap_VPN_df['synToSIZ_dist'] = synMap_VPN_df.apply(calculate_syn_to_siz_dist, axis=1, args=(sizSection,))
    synMap_VPN_df['synToSoma_dist'] = synMap_VPN_df.apply(calculate_syn_to_soma_dist, axis=1, args=(somaSection,))
    synMap_VPN_df['type'] = pd.Categorical(synMap_VPN_df['type'])

    grouped_data = synMap_VPN_df.groupby('type')[['synToSIZ_dist', 'synToSoma_dist']].agg(['mean', 'count']).reset_index()

    # Rename the columns for clarity
    grouped_data = grouped_data.rename(columns={'synToSIZ_dist': 'average_syn_to_siz_dist', 
                                                'synToSoma_dist': 'average_syn_to_soma_dist',
                                                'count': 'synapse_count'})

    # Print or use the resulting DataFrame
    print(grouped_data)

def soma_to_siz_dis(neuron_name, sizSection, somaSection):
    soma_to_siz = h.distance(somaSection(0.5), sizSection(0.05))
    print(f"{neuron_name} soma distance to SIZ: {soma_to_siz}") 


def main():
    #TODO: remove synpt from GROUP for actsyn
    neuron_name = "DNp03"

    Tk().withdraw()
    fd_title = "Select morphology file to initialize"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"datafiles/morphologyData", title=fd_title)

    cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, shape_window, electrodeSec = initializeModel(morph_file, neuron_name)

    if neuron_name == "DNp01": 
        erev = -66.4298 - 0.2
        raVal = 212#350#266.1
        gleakVal = 1/2300#1/1600#1800
        cmVal = 0.7#0.77#0.8

    elif neuron_name == "DNp02":
        #DNp03 paramters, NEED TO BE DOUBLE CHECKED AND FIXED IF NECESSARY
        erev = -61.45
        raVal = 97.449523130925186
        gleakVal = 1/2400
        cmVal = 0.7

    elif neuron_name == "DNp03":
        raVal = 50#350#266.1
        gleakVal = 1/3150#1800
        cmVal = 0.8#0.77#0.8#350#266.1
        erev = -61.15

    elif neuron_name == "DNp04":
        raVal = 50#350#266.1
        gleakVal = 1/3150#1800
        cmVal = 0.8#0.77#0.8#350#266.1
        erev = -61.15

    elif neuron_name == "DNp06":
        raVal = 50#350#266.1
        gleakVal = 1/3150#1800
        cmVal = 0.8#0.77#0.8#350#266.1
        erev = -61.15

    elec_raVal = 235.6                  
    # elec_gleakVal = 0
    elec_cmVal = 6.4

    sealCon_8GOhm = 0.0003978
    # sealCon_2GOhm = 0.0016
    elec_gleakVal = sealCon_8GOhm


    change_Ra(ra=raVal, electrodeSec=electrodeSec, electrodeVal = elec_raVal)
    change_gLeak(gleak=gleakVal, erev=erev, electrodeSec=electrodeSec, electrodeVal = elec_gleakVal)
    change_memCap(memcap=cmVal, electrodeSec=electrodeSec, electrodeVal = elec_cmVal)
    
    ########################################



    # Loading of synapse data

    synMap_df = loadSynMapDataframe(neuron_name)

    synMap_df = removeAxonalSynapses_REV(synMap_df, dendList, axonList=axonList, tetherList=tetherList, somaSection=somaSection)

    synMap_VPN_df = filterSynapsesToCriteria(['LC4', 'LC6', 'LC22', 'LPLC1', 'LPLC2', 'LPLC4'], synMap_df)
    
    #soma_to_siz_dis(neuron_name, sizSection, somaSection)
    #VPN_paired_activations_by_neurons(neuron_name, synMap_VPN_df, sizSection, somaSection, erev, population = 'LPLC2')
    # avg_dist_to_siz_per_VPN(neuron_name, sizSection, somaSection, synMap_VPN_df)



    # activate_VPNs_by_receptive_field2(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numspike = None, interval = None, region = "Dorsal", VPN_subtype = "LC4_1")
    # activate_VPNs_by_receptive_field2(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numspike = None, interval = None, region = "Dorsal", VPN_subtype = "LC4 and LPLC2_1")
    # activate_VPNs_by_receptive_field2(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numspike = None, interval = None, region = "Dorsal", VPN_subtype = "LPLC2_1")
    # activate_VPNs_by_receptive_field(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numspike = 5, interval = 10, region = "Ventral", VPN_subtype = "LPLC2")
    # activate_VPNs_by_receptive_field(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numspike = 5, interval = 10, region = "Ventral", VPN_subtype = "LPLC4")
    
    
    #activate_VPNs_by_receptive_field(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numspike = None, interval = None, region = "Ventral", VPN_subtype = "LC4 and LPLC2")
    

    # plot_peak_depol_vs_freq(neuron_name, somaSection, erev, MODE='Depol_SIZ_and_SOMA', interval = None, numspikes = None, region = None, VPN_subtype = None)
    # plot_peak_depol_vs_freq(neuron_name, somaSection, erev, MODE='Frequency', interval = None, numspikes = None, region = None, VPN_subtype = None)
    # quit(0)
    ########################################

    ###Loading in a particular synapse dataframe for single synapse activations
    #activateAllSingle(neuron_name, erev, synMap_df, sizSection, somaSection)

    # Activating synapses of far/close/random trials, need to select the appropriate mapped synapse dataframe which is found in the simulation folder for each trial/neuron.

    #synMap_df_close = loadSynMapDataframe(neuron_name)
    #synMap_df_far = loadSynMapDataframe(neuron_name)
    # synMap_df_rand_trial_1 = loadSynMapDataframe(neuron_name)

    #tVec, vVec = activateSynapses(synMap_df_rand_trial_1, MODE="SINGLE_EXP2", somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname='DNp01_final_sims_1113/singleTEST/50_rand_trial_2'+neuron_name, filename=None)
    # quit(0)


    ####Loading the synapse dataframes from simulations previously done for partner, random and close.

    activate_indiv_syns_rand_VPN_close(neuron_name, erev, sizSection, somaSection, syn_count =15, trials=50, simulation = "Partner")
    #activate_indiv_syns_rand_VPN_close(neuron_name, erev, sizSection, somaSection, syn_count =15, trials=50, simulation = "Close")
    quit(0)


    ########################################
    # Current injection simulations 
    #DNp01
    #activateCurrentInjectionSIZ(neuron_name, erev, sizSection, somaSection, electrodeSec= electrodeSec, continueRun=50, current=40.5, injDur=.1, delay= 5)
    #DNp03
    #activateCurrentInjectionSIZ(neuron_name, erev, sizSection, somaSection, electrodeSec, continueRun=50, current=4.5, injDur=.1, delay= 5)

    # # Note with these simulations changing the amount of current, the peak amplitudes of the soma and SIZ change, but the ratio is the same (normalization/length constant phenomonon?)
    #activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.001, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.005, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.01, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.05, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.1, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.2, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=.5, injDur=5, delay= 5)
    #activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=5, delay= 5)


    #######################

    #Plotting of DNp01 current injection to somas with different injection durations
    #DNp01
    # soma_current_inject_durations = {
    # 'Duration (ms)': [0.1, 0.5, 1, 5, 10, 20],
    # 'Membrane potential at max SIZ depol': [-66.50284688742303, -66.10232979694965, -65.62344259565623, -63.609277710028046, -63.323741614418985, -60.01440210822973],
    # 'Max SIZ depol': [0.12695311257697028, 0.5274702030503562, 1.006357404343774, 3.0205222899719573, 3.306058385581018, 6.615397891770272],
    # 'Membrane potential change at Soma': [2.5104095443859364, 9.647389679798614, 15.537313048127643, 25.82902973413333, 26.27872688471438, 52.374995264797974],
    # 'Ratio of Soma max amplitude/SIZ max amplitude': [19.77430480772106, 18.289923533135774, 15.439160064867036, 8.551179979662765, 7.948657833547627, 7.917134558142578]
    # } 

#     #DNp03
#     data = {
#     'Duration (ms)': [0.1, 0.5, 1, 5, 10, 20],
#     'Membrane potential at max SIZ depol': [-59.93043930741353, -56.0512937430746, -51.38777990117752, -30.126957830152353, -24.791417761008105, -23.905902504349164],
#     'Max SIZ depol': [1.2195606925864695, 5.0987062569253965, 9.762220098822482, 31.023042169847646, 36.35858223899189, 37.244097495650834],
#     'Membrane potential change at Soma': [26.313105513456094, 93.88357022707838, 149.39840203537233, 249.7563463599723, 257.80021003823487, 258.90821601696814],
#     'Ratio of Soma max amplitude/SIZ max amplitude': [21.57588849280696, 18.413214155955654, 15.303732196469607, 8.050672303270085, 7.090491272285177, 6.951657669975814]
# }

#     # Create DataFrame
#     df = pd.DataFrame(data)

#     # df = pd.DataFrame(soma_current_inject_durations)

#     # Plotting
#     plt.plot(df['Duration (ms)'], df['Ratio of Soma max amplitude/SIZ max amplitude'], marker='o')
#     plt.xlabel('Injection Duration (ms)')
#     plt.ylabel('Ratio of Soma max amplitude / SIZ max amplitude')
#     plt.title('Ratio of Soma to SIZ vs current injection duration')
#     plt.show()


    #Plotting of DNp01 current injection to somas with different injection durations
    #DNp01
    # soma_current_inject_strenghts = {
    # 'Strength (pA)': [1, 5, 10, 50, 100, 200, 500, 1000, 5000, 10000, 20000],
    # 'Membrane potential at max SIZ depol': [-66.599539664746, -66.5875668146761, -66.57260053776832, -66.45287032250609,  -66.3032075534283, -66.00388201527271, -65.10590540080597, -63.609277710028046, -51.63625618380463, -36.66997927602547, -6.737425460465518],
    # 'Max SIZ depol': [0.03026033525399896, 0.04223318532390863,  0.05719946223167938, 0.17692967749391642, 0.32659244657170916, 0.6259179847272947, 1.523894599194037, 3.0205222899719573, 14.99354381619537, 29.95982072397453, 59.89237453953449],
    # 'Membrane potential change at Soma': [0.24231966683194628, 0.34476895639070904, 0.4728305683391767, 1.4973234639268185, 2.777939583411367, 5.339171822380479, 13.022868539287792, 25.82902973413333, 128.2783192928979, 256.33993124135293, 512.4631551382729],
    # 'Ratio of Soma max amplitude/SIZ max amplitude': [8.00783153253146, 8.163460883816686, 8.266346393678226, 8.462816895024876, 8.505829245507126, 8.530146045742232, 8.545780361827765, 8.551179979662765, 8.555570375186235, 8.556123669866418, 8.556400695050085]
    # } 

        #DNp03
#     soma_current_inject_strenghts = {
#     'Strength (pA)': [1, 5, 10, 50, 100, 200, 500, 1000],
#     'Membrane potential at max SIZ depol': [-60.843836676284276, -60.7208918638329, -60.5671535018546, -59.33724660602825, -57.799862986245316, -54.72509574667944, -45.500794027981755, -30.126957830152353],
#     'Max SIZ depol': [0.3061633237157224, 0.42910813616710186, 0.5828464981454005, 1.8127533939717466, 3.350137013754683, 6.424904253320555, 15.649205972018244, 31.023042169847646],
#     'Membrane potential change at Soma': [2.200415644676184, 3.1916305824751703, 4.430649254723896, 14.34279863271371, 26.732985355200988, 51.51335880017555, 125.85447913509947, 249.7563463599723],
#     'Ratio of Soma max amplitude/SIZ max amplitude': [7.187064792644156, 7.437823507574548, 7.60174294402057, 7.912162062644719, 7.9796692629116865, 8.017762875384816, 8.04222778843445, 8.050672303270085]
# }


    # df = pd.DataFrame(soma_current_inject_strenghts)

    # # Plotting
    # # plt.plot(df['Strength (pA)'], df['Ratio of Soma max amplitude/SIZ max amplitude'], marker='o')
    # # plt.ylim(7, 8.2)
    # # plt.xlim(0, 600)
    # # plt.xlabel('Current injection (pA)')
    # # plt.ylabel('Ratio of Soma max amplitude / SIZ max amplitude')
    # # plt.title('Ratio of Soma to SIZ vs current injection strength')
    # # plt.show()

    # plt.plot(df['Strength (pA)'], df['Max SIZ depol'], marker='o')
    # plt.ylim(0, 60)
    # plt.xlim(0, 20000)
    # plt.xlabel('Current injection (pA)')
    # plt.ylabel('SIZ max amplitude')
    # plt.title(f'{neuron_name} SIZ peak amplitude vs current injection strength')
    # plt.show()


    #######################

    # #Note that changing the injection duration will impact the ratio seen above. 
    
    #activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=.1, delay= 5)
    #quit(0)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=.5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=1, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=5, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=10, delay= 5)
    # activateCurrentInjectionSOMA(erev, sizSection, somaSection, electrodeSec, continueRun=50, current=1, injDur=20, delay= 5)

    # Synapse at SIZ simulation 
    #activateFakeSpikeSyn(neuron_name, erev, synMap_df, sizSection, somaSection, electrodeSec)

     ########################################

    #LPLC2_synMap_df = filterSynapsesToCriteria(['LPLC2'], synMap_df)
    #activateLPLC2(neuron_name, LPLC2_synMap_df, erev, allSections_py, sizSection, somaSection, dirname='73123_sims/LPLC2_test/HOLDING', filename=None)
    #activateLPLC2_TEST(neuron_name, LPLC2_synMap_df, erev, allSections_py, sizSection, somaSection, dirname=None, filename=None, onsetTime = 10, randomseed=7) # This function still needs to be fixed.
    #quit(0)
        
    # for id in filtered_ids:
    #     count = synID_countDic.get(id)
    #     print(id, count)
    #     filtered_df = synMap_df[synMap_df['pre'] == id]
    #     # print(filtered_df)
    #     # activateLPLC2(neuron_name, LPLC2_synMap_df, erev, sizSection, somaSection, dirname='73123_sims/LPLC2_test/HOLDING_TEST', filename='syns_'+str(count)+'_'+str(id))
    #     tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname='73123_sims/LPLC2_test/HOLDING', filename='syns_'+str(count)+'_'+str(id))
    #     quit(0)

    ##################
    
    #activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 2, trials=50)
    #activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 5, trials=50)
    #activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 10, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 15, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 20, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 40, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 70, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 80, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 90, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 100, trials=50)
    # activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 200, trials=50)
    # for x in range(3, 25):
    #     activate_close_syns(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn =x, trials=10)
    # for x in range(1, 25):
    #     activate_random_syns(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn =x, trials=10)
    # quit(0)
##################



    # synType_df = synMap_df[['type']]
    # synTypeList = []
    # for index, row in synType_df.iterrows():
    #     #print(row)
    #     synTypeList.append(row.loc['type'])
    # synTypeSet = set(synTypeList)
    # uniqueSynTypes = list(synTypeSet)
    # print(uniqueSynTypes)

    # for synType in uniqueSynTypes:
    #     if isinstance(synType, float):
    #         pass
    #     else:
    #         synSubset_df = filterSynapsesToCriteria([synType], synMap_df)
    #         print(synSubset_df.shape[0])
    # quit(0)

    # for i in range(50):
    #     subsetInfo = ["RAND", 10]
    #     SYNDF = selectSynapses(subsetInfo, synMap_df)
    #     tVec, vVec = activateSynapses(SYNDF, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname='73123_sims/random10/'+neuron_name+"_random10", filename=neuron_name+"_trial_"+str(i))
    # # quit(0)





    # for i in range(50):
    #     print("trial "+str(i+1))
    #     ID_TEMP = random.randint(0, synMap_df.shape[0])
    #     print(ID_TEMP)
    #     subsetInfo = ["CLOSEST", 10, ID_TEMP]
    #     SYNDFC = selectSynapses(subsetInfo, synMap_df)
    #     tVecClose, vVecClose = activateSynapses(SYNDFC, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname='73123_sims/random10/'+neuron_name+"_random10_CLUST", filename=neuron_name+"_trial_"+str(i))

    # quit(0)

    #Simulations for final figures
    # Random vs partner syn activations
    VPN_paired_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, erev)
    # VPN_random_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, erev, numsyn_lim = 40)
    #VPN_random_subset_within_type_activations(neuron_name, synMap_df, sizSection, somaSection, erev)

    # # Random vs close shunting figure 
    #activateRandomClose_TEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn=50, trials=40) # This currently works

    # #Activating all synapses once 
    # activateAllSingleTEST(neuron_name, erev, synMap_VPN_df, sizSection, somaSection) , #this one works

    # #Runs simulations of synapses in increasing number of synapses either up to 250, or for all VPN synapses within a VPN type
    # #activateLinear250_test(neuron_name, erev, synMap_VPN_df, sizSection, somaSection) #this one works 
    #activateLinearVPN_test(neuron_name, erev, synMap_VPN_df, sizSection, somaSection) # this one works

    #Runs simulations activating 50 random trials with any num of syns, activating the closest vs farthest 
    #activateClosestFarthest_test(neuron_name, erev, synMap_VPN_df,sizSection, somaSection, numsyn= 20) # this one works 
    
    # Shuffled simulation, the VPN set here decides which population to shuffeld synapses within.
    # VPN = "LC4"
    # synMap_df, shuffled_df = shuffle_data(neuron_name, VPN, synMap_VPN_df, int_seed= 21)
    # VPN_paired_activations(neuron_name, shuffled_df, sizSection, somaSection, erev, shuffled = True)
    # quit(0)


    #These should all be working with alphasyn
    #activateClosestFarthest(neuron_name, erev, synMap_VPN_df, sizSection, somaSection)
    # activateLinear250(neuron_name, erev, synMap_VPN_df, sizSection, somaSection)
    # activateAllSingle(neuron_name, erev, synMap_df, sizSection, somaSection)
    # activateRandomClose(neuron_name, erev, synMap_df, sizSection, somaSection)
    #activateLinearVPN(neuron_name, erev, synMap_df, sizSection, somaSection)

    # quit(0)



    # colorB_testList = returnAllNonDendriticTreeSections(neuron_name, allSections_py, somaSection)
    # window = plotMorphology(allSections_py, somaSection, sizSection, colorB_testList)
    # clock.sleep(1000)

    quit(0)


    colorB, SUBSET_DF = retrieveSynBranch(synMap_df, somaSection, "dend[598]")
    SUBSET_DF = SUBSET_DF.sample(n=10)
    #tVecF, vVecF = activateSynapses(SUBSET_DF, MODE="SEQUENTIAL_F", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname="7923_sims/sequential/"+neuron_name+"_seqAct_dend598", filename=None)
    #tVecR, vVecR = activateSynapses(SUBSET_DF, MODE="SEQUENTIAL_R", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname="113023_sims/sequential/"+neuron_name+"_seqAct_dend598", filename=None)
    tVecR, vVecR = activateSynapses(SUBSET_DF, MODE="SEQUENTIAL_R_EXP2", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname="113023_sims/sequential/"+neuron_name+"_seqAct_dend598_EXP2", filename=None)
    #tVecR, vVecR = activateSynapses(SUBSET_DF, MODE="SEQUENTIAL_F", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname="113023_sims/sequential/"+neuron_name+"_seqAct_dend598", filename=None)
    tVecR, vVecR = activateSynapses(SUBSET_DF, MODE="SEQUENTIAL_F_EXP2", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname="113023_sims/sequential/"+neuron_name+"_seqAct_dend598_EXP2", filename=None)
    quit(0)














    # SYNDF = synMap_df.loc[1:10]
    #TODO fix DNp06 directory name and reorg directory names in general tbh

    #TODO: PICK TRUE LEAVES FOR SEQ
    # print(synMap_df.head(100))
    #quit(0)
    colorB, SUBSET_DF = retrieveSynBranch(synMap_df, somaSection, "dend[598]")
    for sec in colorB:
        print(sec)
    quit(0)
    # print(SUBSET_DF)
    # print("BEFORE")
    # clock.sleep(10)


    ### REMOVES SYNAPSE ON AXON[0]

    SUBSET_DF2 = SUBSET_DF.iloc[:-11 , :]
    # print("AFTER")
    print(SUBSET_DF2)

    shape_window.color_list(colorB, 3)
    # print('zzzz')
    # clock.sleep(10000)
    # quit(0)
    # tVecF, vVecF = activateSynapses(SUBSET_DF2, MODE="SEQUENTIAL_F", somaSection=somaSection, recordSection=sizSection, erev=-68, simTime=200, onsetTime=75, neuron_name=None, filename="SEQ_FORWARD_1")
    # tVecR, vVecR = activateSynapses(SUBSET_DF2, MODE="SEQUENTIAL_R", somaSection=somaSection, recordSection=sizSection, erev=-68, simTime=200, onsetTime=75, neuron_name=None, filename="SEQ_REVERSE_1")
    
    tVecF, vVecF = activateSynapses(SUBSET_DF2, MODE="SEQUENTIAL_F", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname=neuron_name+"_seqAct_dend1505", filename=None)
    tVecR, vVecR = activateSynapses(SUBSET_DF2, MODE="SEQUENTIAL_R", somaSection=somaSection, recordSection=sizSection, erev=erev, simTime=600, onsetTime=300, neuron_name=None, dirname=neuron_name+"_seqAct_dend1505", filename=None)
    # quit(0)
    
    # plt.plot(tVecF, vVecF, 'g-')

    # quit(0)
        #plt.plot(tVecRand, vVecRand, 'r-')
        
        #KEY 2099
        #need to visualize synapses
        #activation time jitter (gaussian 0.5 ms?)

        # subsetInfo = ["CLOSEST", 10, ID_TEMP[0]]
        # SYNDFC = selectSynapses(subsetInfo, synMap_df)
        # print(SYNDFC.head())

        # tVecClose, vVecClose = activateSynapses_NEW(SYNDFC, MODE="GROUP", somaSection=somaSection, recordSection=sizSection, erev=-72.5, simTime=100, neuron_name=None, filename="jitter10")
        # plt.plot(tVecClose, vVecClose, 'g-')


        # plt.show()


    Rtime_df = pd.read_csv('rand50_t.csv', header=None)
    Rvolt_df = pd.read_csv('rand50_v.csv', header=None)
    Ctime_df = pd.read_csv('close50_t.csv', header=None)
    Cvolt_df = pd.read_csv('close50_v.csv', header=None)

    RtList = []
    RvList = []
    CtList = []
    CvList = []

    ct = 0
    for index, row in Rtime_df.iterrows():
        try:
            # ax1.plot(row, synCurr_df.iloc[ct])
            # ax2.plot(row, sizVolt_df.iloc[ct])
            RtList.append(row)
            RvList.append(Rvolt_df.iloc[ct])

            CtList.append(row)
            CvList.append(Cvolt_df.iloc[ct])
            ct +=1
            #print(ct)
        except IndexError:
            break
    # plt.show()
    print(Rvolt_df)
    print(Cvolt_df)
    # plt.plot(tVec, vSoma, 'm-')

    
    grid_RVC = gridspec.GridSpec(1, 2)
    fig = plt.figure()
    ax_R = fig.add_subplot(grid_RVC[0, 0]) # row 0, col 0
    ax_C = fig.add_subplot(grid_RVC[0, 1])
    for trial in range(len(RtList)):
        Rv = RvList[trial]#.to_python()
        Rv_np = numpy.array(Rv)
        RtVec = RtList[trial]#.to_python()
        RtVec_np = numpy.array(RtVec)
        ax_R.plot(RtVec_np, Rv_np, 'r-')#, linewidth=0.6)

        Cv = CvList[trial]#.to_python()
        Cv_np = numpy.array(Cv)
        CtVec = CtList[trial]#.to_python()
        CtVec_np = numpy.array(CtVec)
        ax_C.plot(CtVec_np, Cv_np, 'g-')#, linewidth=0.6)

        ax_R.set_ylabel('Voltage (mV)')
        ax_R.set_xlabel('Time (ms)')

        ax_C.set_ylabel('Voltage (mV)')
        ax_C.set_xlabel('Time (ms)')

        fig.tight_layout()
    plt.show()
    '''
    # tRand10, vRand10 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10, vRand10, 'm-')
    # syn_2, tRand10_2, vRand10_2 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_2, vRand10_2, 'g-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'b-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'b-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'b-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'b-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'r-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'r-')
    # syn_3, tRand10_3, vRand10_3 = selectSynapses(['RAND'], somaSection, neuron_name)
    # plt.plot(tRand10_3, vRand10_3, 'r-')
    # peak_vRand10_val = numpy.amax(numpy.array(vRand10.to_python()))
    # peak_vRand10_valIdx = numpy.argmax(numpy.array(vRand10.to_python()))
    # tRand10_np = numpy.array(tRand10.to_python())
    # peak_tRand10 = tRand10_np[peak_vRand10_valIdx]
    # print(peak_vRand10_val, peak_tRand10)
    

    # syn, tLPLC2, vLPLC2 = selectSynapses(['LPLC2'], somaSection, neuron_name)
    # plt.plot(tLPLC2, vLPLC2, 'm-')
    # peakLPLC2_val = numpy.amax(numpy.array(vLPLC2.to_python()))
    # peakLPLC2_valIdx = numpy.argmax(numpy.array(vLPLC2.to_python()))
    # tLPLC2_np = numpy.array(tLPLC2.to_python())
    # tLPLC2 = tLPLC2_np[peakLPLC2_valIdx]
    
    # syn, tLC4, vLC4 = selectSynapses(['LC4'], somaSection, neuron_name)
    # plt.plot(tLC4, vLC4, 'g-')
    # peakLC4_val = numpy.amax(numpy.array(vLC4.to_python()))
    # peakLC4_valIdx = numpy.argmax(numpy.array(vLC4.to_python()))
    # tLC4_np = numpy.array(tLC4.to_python())
    # tLC4 = tLC4_np[peakLC4_valIdx]

    # syn, tNA, vNA = selectSynapses(['na'], somaSection, neuron_name)
    # plt.plot(tNA, vNA, 'k-')
    # peakNA_val = numpy.amax(numpy.array(vNA.to_python()))
    # peakNA_valIdx = numpy.argmax(numpy.array(vNA.to_python()))
    # tNA_np = numpy.array(tNA.to_python())
    # tNA = tNA_np[peakNA_valIdx]

    # print(peakLC4_val, peakLPLC2_val, peakNA_val)
    # print(tLC4, tLPLC2, tNA, tLPLC2-tLC4)

    # plt.show()
    #quit(0)

    #TODO: WRITE SYNAREA, DIST, AND TYPE TO CSV AND READ AS DF SO YOU CAN PLOT WO HAVING TO SIM
    #synAreaList, synDistList, synTypeList = activateSynapses(sizSection, neuron_name)
    '''
    quit(0)
    #Misc dendrogram code
    from scipy.cluster.hierarchy import dendrogram, linkage

    # Initialize a dictionary to store distances
    synToSIZ_dist = {}
    for idx, synSite in synMap_df.iterrows():
        synSec_name = synSite['mappedSection']
        
        # Find the NEURON section object by name
        synSec = None
        for sec in h.allsec():
            if sec.name() == synSec_name:
                synSec = sec
                break
        
        if synSec is not None:
            # Use the found NEURON section object
            dist = h.distance(synSec(synSite['mappedSegRangeVar']), sizSection(0.1))
            synToSIZ_dist[idx] = dist
        else:
            print(f"Warning: Section '{synSec_name}' not found for synapse {idx + 1}. Skipping.")


    # Convert distances to a DataFrame for linkage
    distances_df = pd.DataFrame(list(synToSIZ_dist.values()), index=list(synToSIZ_dist.keys()), columns=['Distance'])
    # Perform hierarchical clustering on distances
    linked = linkage(distances_df, 'complete')

    # Plot the dendrogram horizontally
    dendrogram_plot = dendrogram(linked, orientation='left', labels=None, distance_sort='ascending', show_leaf_counts=True, )

    # Sort the DataFrame based on dendrogram leaves
    sorted_distances_df = distances_df.iloc[dendrogram_plot['leaves']]

    # Plot synapse points on the dendrogram with actual synapse distance on the x-axis
    for idx, dist in synToSIZ_dist.items():
        if idx in sorted_distances_df.index:
            leaf_pos = sorted_distances_df.index.get_loc(idx)
            plt.scatter(dist, leaf_pos, color='blue', marker='o')

    plt.title('Dendrogram with Synapse Distances to SIZ')
    plt.xlabel('Distance to SIZ')
    plt.yticks([])
    plt.show()
    quit(0)

main()
    
