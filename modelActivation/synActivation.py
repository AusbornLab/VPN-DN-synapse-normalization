from neuron import h, gui
from neuron.units import ms, mV
import os
import numpy
from matplotlib import pyplot as plt
import random
import math
import csv
import pandas as pd
from tkinter import Tk
import tkinter.filedialog as fd

#2 = axon | 3 = dendrite | 0 = unlabeled | 1 = soma
h.load_file("stdrun.hoc")

##########
#Initializing the model
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

    for sec in sectionListToDiscretize:
        secDensityMechDict = sec.psection()['density_mechs']

        #using the membrane resistance, section diameter, and section axial resistance, we calculate the spatial constant (lambda) for a given section
        secLambda = math.sqrt( ( (1 / secDensityMechDict['pas']['g'][0]) * sec.diam) / (4*sec.Ra) )

        if (sec.L/sec.nseg)/secLambda > 0.1:

            numSeg_log2 = math.log2(sec.L/secLambda / 0.1)
            numSeg = math.ceil(2**numSeg_log2)
            if numSeg % 2 == 0:
                numSeg += 1
            sec.nseg = numSeg
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
    lineStart_X = axonEnd.x3d(axonEnd.n3d()-2)
    lineStart_Y = axonEnd.y3d(axonEnd.n3d()-2)
    lineStart_Z = axonEnd.z3d(axonEnd.n3d()-2)
    lineEnd_X = axonEnd.x3d(axonEnd.n3d()-1)
    lineEnd_Y = axonEnd.y3d(axonEnd.n3d()-1)
    lineEnd_Z = axonEnd.z3d(axonEnd.n3d()-1)

    lineDir = numpy.array([lineEnd_X, lineEnd_Y, lineEnd_Z]) - numpy.array([lineStart_X, lineStart_Y, lineStart_Z])
    line_direction_norm = lineDir / numpy.linalg.norm(lineDir)

    if neuron_name == "DNp01":
        equivCylHeight = 241.69
        equivCylDiam = 3.32*2
    elif neuron_name == "DNp02":
        equivCylHeight = 220.654
        equivCylDiam = 0.624*2
    elif neuron_name == "DNp03":
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

    x1 = axonEnd.x3d(axonEnd.n3d() - 1)
    y1 = axonEnd.y3d(axonEnd.n3d() - 1)
    z1 = axonEnd.z3d(axonEnd.n3d() - 1)

    equivCylAxon.pt3dadd(x1, y1, z1, equivCylDiam)
    equivCylAxon.pt3dadd(new_point[0], new_point[1], new_point[2], equivCylDiam)
    pySectionList.append(equivCylAxon)
    equivCylAxon.connect(axonEnd, 1)

    return pySectionList

def createElectrode(somaSection, pySectionList, neuron_name=None):
    electrodeSec = h.Section()
    electrodeSec.L = 10
    electrodeSec.diam = 1
    electrodeSec.connect(somaSection, 0)

    pySectionList.append(electrodeSec)

    return pySectionList, electrodeSec

##########

# Helper functions
def write_csv(name,data):
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
    return synMap_df

def filterSynapsesToCriteria(subsetInfo, synMap_df):
    filteredSyn_List = []
    for index, row in synMap_df.iterrows():
        if row.loc['type'] in subsetInfo:
            filteredSyn_List.append(row)
    filteredSyn_df = pd.DataFrame(filteredSyn_List)
    return filteredSyn_df

def removeAxonalSynapses_REV(syn_df, dendList, axonList=None, tetherList=None, somaSection=None):
    synSecSeries = syn_df['mappedSection']
    tetherSynSites = syn_df[synSecSeries.str.startswith('dend_11')]
    somaSynSites = syn_df[synSecSeries.str.startswith('soma')]
    axonSynSites = syn_df[synSecSeries.str.startswith('axon')]
    dendSynSites = syn_df[synSecSeries.str.startswith('dend[')]
    return dendSynSites

def str2sec(sec_as_str):
    for sec in h.allsec():
        if sec.name() == sec_as_str:
            sec_as_sec = sec
            break
    return sec_as_sec

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

##########

#####
#Simulations
def activateSynapses(syn_df, MODE, interval = None, num_spike = None, somaSection=None, recordSection=None, sizSection=None, erev=-72.5, simTime=100, onsetTime=None, neuron_name=None, dirname=None, filename=None):
    #other parameters are as follows:
    #somaSection --- the section of the morphology that corresponds to the soma
    
    #recordSection --- the section of the morphology that is being recorded from. If no recordSection is passed to this function,
    #in the case of single synapse activation, the recordSection for each synapse will be the section to which that synapse is mapped to,
    #and in the case of group synapse activation, the recordSection will default to the soma. If no somaSection is provided, an error
    #will be raised

    #erev --- membrane resting/reversal potential (mV)

    #simTime --- simulation duration (ms) | onsetTime --- synapse activation start time (ms)
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
                synapses['syn_{}'.format(0)].e = -10


                ncstim = h.NetCon(stim, synapses['syn_{}'.format(0)])
                ncstim.delay = 0  
                ncstim.weight[0] = exp2Gmax
                ########################

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(0)] = h.Vector()
                if neuron_name == "DNp01":
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.0551186)._ref_v)
                elif neuron_name == "DNp03":
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
                syn.e = -10  

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
        tVec = None
        vVec = None  

    elif MODE == "GROUP_EXP2_interval": #Note this is not used within the paper, however is an extra simulation to use to activate neurons at a given time interval
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
                syn.e = -10  

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
        tVec = None
        vVec = None  

    return tVec, vVec

def selectSynapses(subsetInfo, synMap_df):

    synSubset_df = pd.DataFrame()
    
    if subsetInfo[0] == "RAND":
        randN = subsetInfo[1]
        randSelectedSynIdxs = []
        randSelectedSyns = pd.DataFrame()

        if randN > synMap_df.shape[0]-1:
           
            synSubset_df = synMap_df
            print('TRYING TO SUBSET MORE SYNS THAN IN AVAILABLE POPULATION /// RETURNING WHOLE POPULATION AS SUBSET')
        else:
            uarray = numpy.random.choice(numpy.arange(0, synMap_df.shape[0]-1), replace=False, size=(1, randN))

            for element in uarray[0]:
               
                randSelectedSyns = pd.concat([randSelectedSyns, synMap_df.iloc[[element]]])
            synSubset_df = randSelectedSyns

    elif subsetInfo[0] == "CLOSEST":
        closestN = subsetInfo[1]
        closestMethod = "PATH" 
        STARTSYN = synMap_df.iloc[subsetInfo[2]]
        #print(STARTSYN)

        if str(STARTSYN.loc['pre']) == 'nan':
                print("ERROR: RAND SYN IS NAN")
                raise IndexError
        else:
            for sec in h.allsec():
                if sec.name() == STARTSYN.loc["mappedSection"]:
                    STARTSYN_SEC = sec
                    break
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if str(row.loc['pre']) == 'nan':
                print("BOO")
                pass
            else:
                for sec in h.allsec():
                    if sec.name() == row.loc["mappedSection"]:
                        CURRSEC = sec
                        
                        break
                
                CURR_SYNDIST =  h.distance(STARTSYN_SEC(STARTSYN.loc['mappedSegRangeVar']), CURRSEC(row.loc['mappedSegRangeVar']))
                SYNDIST_DATA.append([row, CURR_SYNDIST])
                
        SYNDIST_DATA.sort(key = lambda x: x[1])
        

        seriesList = []
        closestSelectedSyns = pd.DataFrame()
        for closestSynCt in range(closestN):
            closestSynX = SYNDIST_DATA[closestSynCt]
            seriesList.append(closestSynX[0])
            closestSelectedSyns = pd.concat([closestSelectedSyns, closestSynX[0]], axis=0)
            
        synSubset_df = closestSelectedSyns
        
        cols = ['pre','pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z', 'type', 'mappedSection', 'mappedSegRangeVar']
        synSubset_df = pd.DataFrame(seriesList)
       
    elif subsetInfo[0] == "CLOSEST2SEC":
        closestN = subsetInfo[1]
        closestMethod = "PATH" 

        STARTSEC = subsetInfo[2]
        STARTSEC_RV = subsetInfo[3]
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if str(row.loc['pre']) == 'nan':
                
                pass
            else:
                for sec in h.allsec():
                    if sec.name() == row.loc["mappedSection"]:
                        CURRSEC = sec
                        break
                
                CURR_SYNDIST =  h.distance(STARTSEC(STARTSEC_RV), CURRSEC(row.loc['mappedSegRangeVar']))
                SYNDIST_DATA.append([row, CURR_SYNDIST])
        SYNDIST_DATA.sort(key = lambda x: x[1])
    
        seriesList = []
        closestSelectedSyns = pd.DataFrame()
        for closestSynCt in range(closestN):
            
            closestSynX = SYNDIST_DATA[closestSynCt]
           
            seriesList.append(closestSynX[0])
            
            closestSelectedSyns = pd.concat([closestSelectedSyns, closestSynX[0]], axis=0)
            
        synSubset_df = closestSelectedSyns
        
        cols = ['pre','pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z', 'type', 'mappedSection', 'mappedSegRangeVar']
        synSubset_df = pd.DataFrame(seriesList)
        
    elif subsetInfo[0] == "FARTHEST2SEC":
        farthestN = subsetInfo[1]
        STARTSEC = subsetInfo[2]
        STARTSEC_RV = subsetInfo[3]
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if str(row.loc['pre']) == 'nan':
                
                pass
            else:
                for sec in h.allsec():
                    if sec.name() == row.loc["mappedSection"]:
                        CURRSEC = sec
                        break
                
                
                CURR_SYNDIST =  h.distance(STARTSEC(STARTSEC_RV), CURRSEC(row.loc['mappedSegRangeVar']))
                

                SYNDIST_DATA.append([row, CURR_SYNDIST])

        SYNDIST_DATA.sort(key=lambda x: x[1], reverse=True)
        
        for i in range(len(SYNDIST_DATA)):
            t = SYNDIST_DATA[i]
            

        seriesList = []
        farthestSelectedSyns = pd.DataFrame()
        for farthestSynCt in range(farthestN):
            
            farthestSynX = SYNDIST_DATA[farthestSynCt]
            
            seriesList.append(farthestSynX[0])
           
            farthestSelectedSyns = pd.concat([farthestSelectedSyns, farthestSynX[0]], axis=0)
            
        synSubset_df = farthestSelectedSyns
        
        cols = ['pre','pre_x', 'pre_y', 'pre_z', 'post_x', 'post_y', 'post_z', 'type', 'mappedSection', 'mappedSegRangeVar']
        synSubset_df = pd.DataFrame(seriesList)
    return synSubset_df

def activateClosestFarthest(neuron_name, erev, synMap_df, sizSection, somaSection, numsyn= None):
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

def activateLinear250(neuron_name, erev, synMap_df, sizSection, somaSection):
    for i in range(250):
        print(i+1, "/ 250")

        subsetInfo = ["RAND", i+1]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        ID_TEMP = SYNDF.axes[0].tolist()
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, onsetTime=50, somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims_1113/linear250/'+neuron_name+"_VPN_only_250", filename=None)

def activateAllSingle(neuron_name, erev, synMap_df, sizSection, somaSection):
    #synMap_df = synMap_df.head(100), you can choose how many synapses to activate, if you do not partition then it will activate all synapses
    tVec, vVec = activateSynapses(synMap_df, MODE="SINGLE_EXP2", interval = None, num_spike= None, somaSection=somaSection, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims/singlesynactivation', filename=None)
    
def activateRandomClose(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")
        synMap_df.reset_index(drop=True)
        subsetInfo = ["RAND", numSyn]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        ID_TEMP = SYNDF.axes[0].tolist()
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_rand{numSyn}record{trials}", filename=None)
        
        subsetInfo = ["CLOSEST", numSyn, ID_TEMP[0]]
        SYNDFC = selectSynapses(subsetInfo, synMap_df)
        tVecClose, vVecClose = activateSynapses(SYNDFC, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_close{numSyn}record{trials}", filename=None)

def activateRandom(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")

        subsetInfo = ["RAND", numSyn]
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

def activate_random_syns(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None):
    numSyn= numSyn
    trials = trials
    for i in range(trials):
        print(i+1, f"/{trials}")
        subsetInfo = ["RAND", numSyn]
        SYNDF = selectSynapses(subsetInfo, synMap_df)
        tVecRand, vVecRand = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, somaSection=somaSection, onsetTime= 50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=None, dirname=neuron_name+f"_rand{numSyn}record{trials}", filename=None)

def activateLinearVPN(neuron_name, erev, synMap_df, sizSection, somaSection):
    synType_df = synMap_df[['type']]
    synTypeList = []
    for index, row in synType_df.iterrows():
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
            
            for i in range(len(synSubset_df)):
                print(i+1, f"/ {synSubset_df_len} for", synType)
                subsetInfo = ["RAND", i+1]
                SYNDF = selectSynapses(subsetInfo, synSubset_df)
                print(SYNDF.head())
                tVec, vVec = activateSynapses(SYNDF, MODE="GROUP_EXP2", interval = None, num_spike= None, onsetTime=50, somaSection=somaSection, sizSection=sizSection, recordSection=sizSection, erev=erev, simTime=75, neuron_name=None, dirname=neuron_name+"_incRandSynCt_VPN", filename=synType+'_'+str(i+1)+'syn')

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
            tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2", somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+"_partner", filename='syns_'+str(count)+'_'+str(id))
    
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
            tVec, vVec = activateSynapses(SYNDF, MODE="GROUP_EXP2", somaSection=somaSection, recordSection=sizSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+"_rand", filename='syns_'+str(j)+'_trial_'+str(i))


########
#Calculating average synapse distance to soma and SIZ 

def calculate_syn_to_siz_dist(row, sizSection):
    synSec = str2sec(row['mappedSection'])
    syn_to_siz_dist = h.distance(synSec(row['mappedSegRangeVar']), sizSection(0.05))
    return syn_to_siz_dist

def calculate_syn_to_soma_dist(row, somaSection):
    synSec = str2sec(row['mappedSection'])
    syn_to_soma_dist = h.distance(synSec(row['mappedSegRangeVar']), somaSection(0.05))
    return syn_to_soma_dist

def soma_to_siz_dis(neuron_name, sizSection, somaSection):
    soma_to_siz = h.distance(somaSection(0.5), sizSection(0.05))
    print(f"{neuron_name} soma distance to SIZ: {soma_to_siz}") 

########

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


def main():
    neuron_name = "DNp03"

    Tk().withdraw()
    fd_title = "Select morphology file to initialize"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"datafiles/morphologyData", title=fd_title)

    cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, shape_window, electrodeSec = initializeModel(morph_file, neuron_name)

    if neuron_name == "DNp01": 
        erev = -66.4298 - 0.2
        raVal = 212
        gleakVal = 1/2300
        cmVal = 0.7

    elif neuron_name == "DNp02":
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

    elif neuron_name == "DNp03":
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

    elif neuron_name == "DNp04":
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

    elif neuron_name == "DNp06":
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

    elec_raVal = 235.6                  
    elec_cmVal = 6.4

    sealCon_8GOhm = 0.0003978
    elec_gleakVal = sealCon_8GOhm


    change_Ra(ra=raVal, electrodeSec=electrodeSec, electrodeVal = elec_raVal)
    change_gLeak(gleak=gleakVal, erev=erev, electrodeSec=electrodeSec, electrodeVal = elec_gleakVal)
    change_memCap(memcap=cmVal, electrodeSec=electrodeSec, electrodeVal = elec_cmVal)
    
    ########################################

    # Loading of synapse data

    synMap_df = loadSynMapDataframe(neuron_name)

    synMap_df = removeAxonalSynapses_REV(synMap_df, dendList, axonList=axonList, tetherList=tetherList, somaSection=somaSection)

    synMap_VPN_df = filterSynapsesToCriteria(['LC4', 'LC6', 'LC22', 'LPLC1', 'LPLC2', 'LPLC4'], synMap_df)

    ########################################

    ###Loading in a particular synapse dataframe for single synapse activations
    ### Relation to Figure 9
    #activateAllSingle(neuron_name, erev, synMap_df, sizSection, somaSection)


    ########################################

    #Simulations for final figures
    ### Relation to Figure 10, Random vs partner syn activations

    #VPN_paired_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, erev)
    #VPN_random_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, erev, numsyn_lim = 40)


    # #Runs simulations of synapses in increasing number of synapses either up to 250, or for all VPN synapses within a VPN type
    #In relation to figure 11
    #activateLinearVPN(neuron_name, erev, synMap_VPN_df, sizSection, somaSection) 

    #In relation to figure 12
    #activateRandomClose(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn= 2, trials=50)
    for x in range(1, 25):
        activate_close_syns(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn =x, trials=10)
    for x in range(1, 25):
        activate_random_syns(neuron_name, erev, synMap_VPN_df, sizSection, somaSection, numSyn =x, trials=10)
    ##################
    

main()
    
