from neuron import h, gui
from neuron.units import ms, mV
import os
import numpy
import numpy as np
from matplotlib import pyplot as plt
import random
import math
import csv
import pandas as pd
from tkinter import Tk
from scipy.integrate import simpson
import tkinter.filedialog as fd
from itertools import combinations, product
from neuron import coreneuron
import traceback
import glob
import time
from random import randrange
from scipy.spatial.distance import cdist
import json
from pathlib import Path

#2 = axon | 3 = dendrite | 0 = unlabeled | 1 = soma
h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")
h.load_file('nrngui.hoc')
pc = h.ParallelContext()

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

def initializeModel(neuron_name, morph_file, coreneuron_flag = False, resting_ev=None, sizsection = None, capacitance = None, axial_resistivity = None, gleak_conductance= None, erev = None):
    """Initializes the model with the given parameters:
    resting_ev: membrane potential the cell sits at, which gets changed for simulations when initializing
    sizsection: The spike initiation zone of the model, highlighted in the .swc file with individual nodes labled "dend_12" 
    capacitance: Capacitance to give the cell 
    axial_resistivity = axial resitivity of the model
    gleak_conductance= leak conductance of the cell
    erev: reversal potential of the leak channel 
    """
    if coreneuron_flag == True:
        coreneuron.enable = True
        pc = h.ParallelContext()
        h.cvode.cache_efficient(1)
    else:
        pass
    cell = instantiate_swc(morph_file)

    allSections_nrn = h.SectionList()
    for sec in h.allsec():
        allSections_nrn.append(sec=sec)
    
    # Create a Python list from this SectionList
    # Select sections from the list by their index

    allSections_py = [sec for sec in allSections_nrn]

    #Define SIZ section
    if sizsection == None:
        sizIndex = 0
    else:
        sizIndex = sizsection
    axonList = h.SectionList()
    tetherList = h.SectionList()
    dendList = h.SectionList()
    preSIZ = h.SectionList()

    for sec in allSections_py:
        if "soma" in sec.name():
            somaSection = sec
        elif "axon" in sec.name():
            axonList.append(sec)
        elif "dend_6" in sec.name():
            axonEnd = sec
        elif "dend_11" in sec.name():
            tetherList.append(sec)
        elif "dend_12" in sec.name():
            sizSection = sec
        elif "dend_13" in sec.name():
            preSIZ.append(sec)
        else:
            dendList.append(sec)

    allSections_py = createAxon(axonEnd, allSections_py, neuron_name)
    allSections_py, electrodeSec = createElectrode(somaSection, allSections_py)
    if coreneuron_flag == True:
        h.finitialize(resting_ev * mV)
    # h.stdinit()

    change_memCap(memcap=capacitance)
    change_Ra(ra=axial_resistivity)
    change_gLeak(gleak= gleak_conductance, erev=erev)

    nsegDiscretization(allSections_py)
    
    #model parameters
    # num_sections = len(allSections_py)
    # total_length = sum([sec.L for sec in allSections_py])  # L is in micrometers by default in NEURON
    # avg_sec_length = total_length / num_sections if num_sections > 0 else 0

    # print(f"Number of Sections: {num_sections}")
    # print(f"Average Section Length (μm): {avg_sec_length:.2f}")
    # print(f"Total Length (μm): {total_length:.2f}")

    return cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, electrodeSec

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
    elif neuron_name == "DNp01_hemi":
        equivCylHeight = 441.69
        equivCylDiam = 3.63
    elif neuron_name == "DNp02":
        equivCylHeight = 220.654
        equivCylDiam = 0.624*2
    elif neuron_name == "DNp03":
        equivCylHeight = 175.671
        equivCylDiam = 0.251*4
    elif neuron_name == "DNp03_hemi":
        equivCylHeight = 375.671
        equivCylDiam = 1.3
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

def loadSynMapDataframe(neuron_name=None):
    Tk().withdraw()
    fd_title = "Select synapse map data file to use for synapse activation"
    syn_file = fd.askopenfilename(filetypes=[("csv file","*.csv")], initialdir=r"datafiles/morphologyData", title = fd_title)
    # synMap_df = pd.read_csv(syn_file)
    synMap_df = pd.read_csv(syn_file, dtype={'pre': str})

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

##########

#####
#Simulations
def activateSynapses(syn_df, MODE, somaSection=None, sizSection=None, erev=-72.5, simTime=100, onsetTime=None, neuron_name=None, dirname=None, filename=None):
    #other parameters are as follows:
    #somaSection --- the section of the morphology that corresponds to the soma
    #erev --- membrane resting/reversal potential (mV)
    #simTime --- simulation duration (ms) | onsetTime --- synapse activation start time (ms)
    #onsetTime --- synapse activation start time (ms)
    #neuron_name --- name of the neuron being simulated
    #dirname --- directory name to save simulation results
    #filename --- file name to save simulation results

    if MODE == "SINGLE_EXP2":
        pc = h.ParallelContext()
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
                if neuron_name == "DNp03_hemi":
                    exp2Gmax = 0.00008
                else:
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
                coreneuron.enable = True
                coreneuron.verbose = 0
                coreneuron.model_stats = True
                coreneuron.num_gpus = 1
                h.cvode.cache_efficient(1)
                h.finitialize(erev)
                pc.set_maxstep(10)


                ########################

                # Create recording vector to track current injection 
                syn_current_siz['v_vec_{}'.format(0)] = h.Vector()
                if neuron_name == "DNp01":
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.0551186)._ref_v)
                elif neuron_name == "DNp03":
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.0551186)._ref_v)
                else:
                    syn_current_siz['v_vec_{}'.format(0)].record(sizSection(0.0551186)._ref_v)

                syn_current_soma['v_vec_{}'.format(0)] = h.Vector()
                syn_current_soma['v_vec_{}'.format(0)].record(somaSection(0.5)._ref_v)

                syn_current_synpt['v_vec_{}'.format(0)] = h.Vector()
                syn_current_synpt['v_vec_{}'.format(0)].record(currSynSec(currSynRangeVar)._ref_v)
                
                #Run the simulation
                print("Round number {}/{}".format(index, syn_df.shape[0]))
                
                h.tstop = simTime
                pc.psolve(h.tstop)

                

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

    elif MODE == "GROUP_EXP2_no_dend_recording":
        pc = h.ParallelContext()
        h.dt = 0.1
        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        seriesList = []
        syn_current_siz = {}
        syn_current_soma = {}
        syn_current_synpt = {}
        synpt_recording_list = []
        synapses = {}  # Using a dictionary to store synapses
        if neuron_name == "DNp03_hemi":
            exp2Gmax = 0.00008
        else:
            exp2Gmax = 0.00027 
        exp2Rise = 0.2
        exp2Decay = 1.1 
        weights = {}

        for index, row in syn_df.iterrows():
            if not pd.isna(row['pre']):
                seriesList.append(row)

                currSynSecStr = row['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                        currSynSec = sec
                        break

                currSynRangeVar = row['mappedSegRangeVar']
                syn = h.Exp2Syn(currSynSec(currSynRangeVar))
                syn.tau1 = exp2Rise
                syn.tau2 = exp2Decay
                syn.e = -10  
                synapses[f'syn_{index}'] = syn

                # Record from synaptic site
                syn_current_synpt[f'v_vec_{index}'] = h.Vector()
                syn_current_synpt[f'v_vec_{index}'].record(currSynSec(currSynRangeVar)._ref_v)
                synpt_recording_list.append(syn_current_synpt[f'v_vec_{index}'])

        # Create shared NetStim
        nc = h.NetStim()
        nc.number = 1
        nc.start = onsetTime
        nc.noise = 0

        ncs = h.List()
        for i, syn in enumerate(synapses.values()):
            ncs.append(h.NetCon(nc, syn))
            ncs.object(i).weight[0] = exp2Gmax 


        # Record from SIZ
        syn_current_siz['v_vec_{}'.format(index)] = h.Vector()
        if neuron_name == "DNp01":
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        elif neuron_name == "DNp03":
            print(sizSection.name())
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.0551186)._ref_v)
        else:
            syn_current_siz['v_vec_{}'.format(index)].record(sizSection(0.05)._ref_v)

        # Record from soma
        syn_current_soma['v_vec_{}'.format(index)] = h.Vector()
        syn_current_soma['v_vec_{}'.format(index)].record(somaSection(0.5)._ref_v)

        coreneuron.enable = True
        coreneuron.verbose = 0
        coreneuron.model_stats = True
        coreneuron.num_gpus = 1
        h.cvode.cache_efficient(1)
        h.finitialize(erev)
        pc.set_maxstep(10)
        h.tstop = simTime
        pc.psolve(h.tstop)

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
                write_csv(outFilename + '_synaptic_currents_synpt.csv', [])
                # for recVec in synpt_recording_list:
                #     write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = path + '/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_synpt.csv', [])
                # for recVec in synpt_recording_list:
                #     write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)
        else:
            if filename is None:
                outFilename = 'datafiles/' + MODE
                write_csv(outFilename + '_time.csv',t_vec)        
                write_csv(outFilename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                write_csv(outFilename + '_synaptic_currents_synpt.csv', [])
                # for recVec in synpt_recording_list:
                #     write_csv(outFilename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(outFilename+'_synapse_sites.csv', index=False)
            else:
                filename = 'datafiles/' + filename
                write_csv(filename + '_time.csv',t_vec)        
                write_csv(filename + '_synaptic_currents_siz.csv', syn_current_siz['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_soma.csv', syn_current_soma['v_vec_{}'.format(index)])
                write_csv(filename + '_synaptic_currents_synpt.csv', [])
                # for recVec in synpt_recording_list:
                #     write_csv(filename + '_synaptic_currents_synpt.csv', recVec)
                synInfo_df.to_csv(filename + '_synapse_site.csv', index=False)

        synapses = None 
        tVec = None
        vVec = None

        pc.barrier()

    elif MODE == "PAIRWISE_EXP2_updated":
        """
        PAIRWISE_EXP2_updated:
        - Expects exactly two synapses in syn_df.
        - Stimulates both synapses sequentially, recording from the other synapse.
        - CoreNEURON compatible.
        - Saves results in the same format as V3.
        """

        pc = h.ParallelContext()

        if len(syn_df) != 2:
            raise ValueError(f"Expected exactly two synapses in syn_df, got {len(syn_df)}")

        syn1, syn2 = syn_df.iloc[0], syn_df.iloc[1]
        syn_pair_df = pd.DataFrame([syn1, syn2])

        # --- Save synapse info ---
        parent_dir = "datafiles/simulationData" if dirname else "datafiles"
        path = os.path.join(parent_dir, dirname) if dirname else parent_dir
        os.makedirs(path, exist_ok=True)
        outFilename = os.path.join(path, filename if filename else MODE)
        synapse_file = outFilename + '_synapse_sites.csv'
        if not os.path.exists(synapse_file):
            syn_pair_df.to_csv(synapse_file, index=False, header=True)
        else:
            syn_pair_df.to_csv(synapse_file, index=False, mode='a', header=False)

        # --- Synapse parameters ---
        exp2Gmax = 0.00027
        exp2Rise = 0.2
        exp2Decay = 1.1

        # --- Stimulate both synapses ---
        for stim_syn, record_syn in [(syn1, syn2), (syn2, syn1)]:

            # Determine NEURON sections and locations
            stim_sec = str2sec(stim_syn['mappedSection'])
            record_sec = str2sec(record_syn['mappedSection'])
            stim_loc = float(stim_syn['mappedSegRangeVar'])
            record_loc = float(record_syn['mappedSegRangeVar'])

            # --- Simulation vectors ---
            t_vec = h.Vector(); t_vec.record(h._ref_t)
            stim_voltage = h.Vector()
            recorded_voltage = h.Vector()

            # --- Setup synapse ---
            synapse = h.Exp2Syn(stim_sec(stim_loc))
            synapse.tau1 = exp2Rise
            synapse.tau2 = exp2Decay
            synapse.e = -10

            stim = h.NetStim()
            stim.number = 1
            stim.start = onsetTime
            stim.noise = 0

            ncstim = h.NetCon(stim, synapse)
            ncstim.delay = 0
            ncstim.weight[0] = exp2Gmax

            stim_voltage.record(stim_sec(stim_loc)._ref_v)
            recorded_voltage.record(record_sec(record_loc)._ref_v)

            # --- Run CoreNEURON simulation ---
            coreneuron.enable = True
            coreneuron.verbose = 0
            coreneuron.model_stats = True
            coreneuron.num_gpus = 1
            h.cvode.cache_efficient(1)
            h.finitialize(erev)
            pc.set_maxstep(10)
            h.tstop = simTime
            pc.psolve(h.tstop)

            # --- Extract results ---
            t = np.array(t_vec)
            v_stim = np.array(stim_voltage)
            v_record = np.array(recorded_voltage)

            baseline_idx = np.where(t < (onsetTime - 1))[0]
            if len(baseline_idx) == 0:
                raise ValueError("No time points before onsetTime - 1 ms for baseline")
            baseline_idx = baseline_idx[-1]

            baseline_stim = v_stim[baseline_idx]
            baseline_record = v_record[baseline_idx]
            resting_vm = np.mean([baseline_stim, baseline_record])

            v_stim = v_stim - baseline_stim + resting_vm
            v_record = v_record - baseline_record + resting_vm

            post_stim_mask = t >= onsetTime
            peak_stim = np.max(v_stim[post_stim_mask])
            peak_record = np.max(v_record[post_stim_mask])

            # --- Save summary ---
            results_file = outFilename + '_summary.csv'
            idx_stim = int(stim_syn['syn_index'])
            idx_record = int(record_syn['syn_index'])
            a, b = sorted([idx_stim, idx_record])
            pair_key = f"{a}-{b}"

            row = {
                "neuron_stim": str(stim_syn['pre']),
                "neuron_record": str(record_syn['pre']),
                "resting_vm": resting_vm,
                "peak_stim": peak_stim,
                "peak_record": peak_record,
                "stim_syn_index": idx_stim,
                "record_syn_index": idx_record,
                "pair_key": pair_key,
                "indices": f"({idx_stim}, {idx_record})"
            }
            df = pd.DataFrame([row])
            if not os.path.exists(results_file):
                df.to_csv(results_file, index=False, header=True)
            else:
                df.to_csv(results_file, index=False, mode='a', header=False)

            print(f"Saved summary results to {results_file}")

            # --- Cleanup ---
            del stim, synapse, recorded_voltage, stim_voltage, t_vec
            pc.barrier()

        vVec = None
        tVec = None

    elif MODE == "ALL_SYNAPSE_EXP2":
        """
        ALL_SYNAPSE_EXP2:
        - Iterates through all synapses in syn_df.
        - For each synapse, stimulates it once.
        - Records Vm from all synapses (including stimulated one) and the SIZ.
        - Saves all stim results into ONE CSV file.
        """

        pc = h.ParallelContext()

        # --- Use the folder path directly ---
        path = dirname if dirname else "datafiles"
        os.makedirs(path, exist_ok=True)

        outFilename = os.path.join(path, filename if filename else MODE)

        # --- Synapse parameters ---
        exp2Gmax = 0.00027
        exp2Rise = 0.2
        exp2Decay = 1.1

        all_results = []  # collect everything in one list

        # --- Loop over all synapses as stim site ---
        for stim_idx, stim_syn in syn_df.iterrows():
            print(f"[ALL_SYNAPSE_EXP2] Stimulating synapse {stim_idx} ({stim_syn['pre']})")

            stim_sec = str2sec(stim_syn['mappedSection'])
            stim_loc = float(stim_syn['mappedSegRangeVar'])

            # Create Exp2Syn at stimulation site
            synapse = h.Exp2Syn(stim_sec(stim_loc))
            synapse.tau1 = exp2Rise
            synapse.tau2 = exp2Decay
            synapse.e = -10

            # Setup stimulation
            stim = h.NetStim()
            stim.number = 1
            stim.start = onsetTime
            stim.noise = 0

            ncstim = h.NetCon(stim, synapse)
            ncstim.delay = 0
            ncstim.weight[0] = exp2Gmax

            # --- Recording setup ---
            t_vec = h.Vector(); t_vec.record(h._ref_t)
            v_traces = []
            for rec_idx, rec_syn in syn_df.iterrows():
                rec_sec = str2sec(rec_syn['mappedSection'])
                rec_loc = float(rec_syn['mappedSegRangeVar'])
                v = h.Vector()
                v.record(rec_sec(rec_loc)._ref_v)
                v_traces.append((rec_idx, rec_syn, v))

            # --- Add SIZ recording ---
            siz_v = h.Vector()
            if neuron_name == "DNp01":
                siz_v.record(sizSection(0.0551186)._ref_v)
            elif neuron_name == "DNp03":
                siz_v.record(sizSection(0.0551186)._ref_v)
            else:
                siz_v.record(sizSection(0.5)._ref_v)

            # --- Run simulation ---
            coreneuron.enable = True
            coreneuron.verbose = 0
            coreneuron.model_stats = True
            coreneuron.num_gpus = 1
            h.cvode.cache_efficient(1)
            h.finitialize(erev)
            pc.set_maxstep(10)
            h.tstop = simTime
            pc.psolve(h.tstop)
            # --- Extract results ---
            t = np.array(t_vec)
            post_stim_mask = t >= onsetTime
            baseline_mask = t < (onsetTime - 1)
            for rec_idx, rec_syn, v_vec in v_traces:
                v = np.array(v_vec)
                baseline_idx = np.argmin(np.abs(t - (onsetTime - 1)))  # closest time to onsetTime - 1 ms

                # Use the voltage at that exact time as baseline
                baseline_vm = v[baseline_idx]
                peak_vm = np.max(v[post_stim_mask]) if np.any(post_stim_mask) else v[-1]

                all_results.append({
                    "stim_syn_index": stim_idx,
                    "record_syn_index": rec_idx,
                    "stim_neuron": str(stim_syn['pre']),
                    "record_neuron": str(rec_syn['pre']),
                    "resting_vm": baseline_vm,
                    "peak_vm": peak_vm,
                    "peak_amp": peak_vm - baseline_vm,
                    "stim_section": stim_syn['mappedSection'],
                    "record_section": rec_syn['mappedSection'],
                    "stim_loc": stim_loc,
                    "record_loc": float(rec_syn['mappedSegRangeVar']),
                })

            # --- Append SIZ row ---
            siz_v_arr = np.array(siz_v)
            baseline_siz_idx = np.argmin(np.abs(t - (onsetTime - 1)))
            siz_baseline = siz_v_arr[baseline_siz_idx]
            siz_peak = np.max(siz_v_arr[post_stim_mask]) if np.any(post_stim_mask) else siz_v_arr[-1]

            all_results.append({
                "stim_syn_index": stim_idx,
                "record_syn_index": len(syn_df)+1,
                "stim_neuron": str(stim_syn['pre']),
                "record_neuron": "SIZ",
                "resting_vm": siz_baseline,
                "peak_vm": siz_peak,
                "peak_amp": siz_peak - siz_baseline,
                "stim_section": "na",
                "record_section": "siz",
                "stim_loc": "na",
                "record_loc": "na",
            })

            # Cleanup
            del synapse, stim, ncstim, t_vec, v_traces, siz_v
            pc.barrier()

        # --- Save ONE CSV with all stim results ---
        results_file = outFilename + "_stim_summary.csv"
        df = pd.DataFrame(all_results)
        df.to_csv(results_file, index=False, header=True)

        print(f"[ALL_SYNAPSE_EXP2] Saved all stim results to {results_file}")

        vVec = None
        tVec = None

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

def activateAllSingle(neuron_name, erev, synMap_df, sizSection, somaSection):
    #synMap_df = synMap_df.head(100), you can choose how many synapses to activate, if you do not partition then it will activate all synapses
    tVec, vVec = activateSynapses(synMap_df, MODE="SINGLE_EXP2", somaSection=somaSection, sizSection=sizSection, erev=erev, simTime=70, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims/singlesynactivation_fly_6', filename=None)
    
def activate_close_syns(neuron_name, erev, synMap_df, sizSection, somaSection, numSyn =None, trials=None, VPN = None, onsetTime=None):
    numSyn= numSyn
    trials = trials
    if VPN == None:
        VPN = "random"
    for i in range(trials):
        print(i+1, f"/{trials}")
        random_value = random.randint(0, len(synMap_df) - 1)
        subsetInfo = ["CLOSEST", numSyn, random_value]
        SYNDFC = selectSynapses(subsetInfo, synMap_df)
        tVecClose, vVecClose = activateSynapses(SYNDFC, MODE="GROUP_EXP2_no_dend_recording", somaSection=somaSection, sizSection=sizSection, erev=erev, simTime=100,  onsetTime=onsetTime, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims/CloseSims/{neuron_name}_close{numSyn}record{trials}', filename = f"{VPN}trial_{i+1:003d}" )

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
                tVec, vVec = activateSynapses(SYNDF, MODE="GROUP_EXP2_no_dend_recording", somaSection=somaSection, sizSection=sizSection, erev=erev, simTime=75, neuron_name=neuron_name, dirname=f'{neuron_name}_final_sims{neuron_name}_incRandSynCt_VPN_fly_6', filename=synType+'_'+str(i+1)+'syn')

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
            if 1 <= count <= 60:
                filtered_ids.append(id)


        for id in filtered_ids:
            count = synID_countDic.get(id)
            print(id, count)
            filtered_df = synMap_df[synMap_df['pre'] == id]
            # print(filtered_df)
            tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2_no_dend_recording", somaSection=somaSection, onsetTime=50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+"_partner", filename='syns_'+str(count)+'_'+str(id))
    
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
            if 1 <= count <= 60:
                filtered_ids.append(id)


        for id in filtered_ids:
            count = synID_countDic.get(id)
            print(id, count)
            filtered_df = synMap_df[synMap_df['pre'] == id]
            # print(filtered_df)
            tVecRand, vVecRand = activateSynapses(filtered_df, MODE="GROUP_EXP2_no_dend_recording", somaSection=somaSection, onsetTime=50, recordSection=sizSection, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims_1113/rand_vs_partner/'+neuron_name+"_partner_shuffled", filename='syns_'+str(count)+'_'+str(id))

def VPN_random_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, erev, numsyn_lim = None):
    numsyn_lim = numsyn_lim
    #Simulation for activating a random subset of synapses synapses that belong to any VPNS
    for j in range(1, numsyn_lim):
        subsetInfo = ["RAND", j]
        for i in range(10):

            SYNDF = selectSynapses(subsetInfo, synMap_VPN_df)
            tVec, vVec = activateSynapses(SYNDF, MODE="GROUP_EXP2_no_dend_recording", somaSection=somaSection, onsetTime= 50, sizSection=sizSection, erev=erev, simTime=100, neuron_name=neuron_name, dirname=neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+"_rand", filename='syns_'+str(j)+'_trial_'+str(i))


###############
#Simulations in relation to figure 10
#For running simulations with large synapse counts within shorter synapse spread ranges (example 10/20)
def select_syns_by_multiple_spreads(
    neuron_name,
    sizSection,
    erev,
    somaSection,
    synMap_df,
    syncount,
    trials,
    target_spreads=[10, 20, 30, 40, 50, 60, 70],
    spread_tolerance=2,
    VPN=None,
    max_runtime_hours=2,
    subset_size=100
):
    if VPN is None:
        VPN = "random"

    # Track trial progress for each target
    trial_counts = {t: 0 for t in target_spreads}
    selected_syn_groups = {t: [] for t in target_spreads}

    start_time = time.time()
    max_runtime = max_runtime_hours * 3600  # seconds

    while any(trial_counts[t] < trials for t in target_spreads):
        # Check runtime
        elapsed = time.time() - start_time
        if elapsed > max_runtime:
            print("\n Max runtime exceeded. Restarting with a new subset of synapses...\n")
            start_time = time.time()  # reset timer
            # just continue, nothing is reset — keeps trial_counts

        # Randomly select `syncount` synapses
        syn_subset = synMap_df.sample(n=subset_size)
        random_value = random.randint(0, len(syn_subset) - 1)
        subsetInfo = ["CLOSEST", syncount, random_value]
        # subsetInfo = ["RAND", syncount]
        syn_subset = selectSynapses(subsetInfo, syn_subset)
        synSecList = syn_subset["mappedSection"].tolist()
        synSegRangeVarList = syn_subset["mappedSegRangeVar"].tolist()

        # Compute all pairwise distances
        distances = []
        synapse_indices = range(len(synSecList))
        for i, j in combinations(synapse_indices, 2):
            synsec_1 = str2sec(synSecList[i])
            synsec_2 = str2sec(synSecList[j])
            dist = h.distance(
                synsec_1(synSegRangeVarList[i]),
                synsec_2(synSegRangeVarList[j])
            )
            distances.append(dist)

        avg_syn_spread = np.mean(distances)
        print(f"Avg synapse spread: {avg_syn_spread:.2f} elapsed time: {elapsed/60:.1f} min")

        # Check against all targets
        for target in target_spreads:
            if trial_counts[target] < trials:  # only consider if we still need trials
                if target <= avg_syn_spread <= (target + spread_tolerance):
                    # Save selection
                    trial_counts[target] += 1
                    selected_syn_groups[target].append(syn_subset)

                    print(f"\nRunning Simulation for target {target}, trial {trial_counts[target]}/{trials}")
                    tVecClose, vVecClose = activateSynapses(
                        syn_subset,
                        MODE="GROUP_EXP2_no_dend_recording",
                        somaSection=somaSection,
                        onsetTime=50,
                        sizSection=sizSection,
                        erev=erev,
                        simTime=100,
                        neuron_name=neuron_name,
                        dirname=f"{neuron_name}_final_sims/synapse_spread/{neuron_name}_{syncount}_syns_rand_by_syn_spread{target}",
                        filename=f"{VPN}_trial_{trial_counts[target]}"
                    )
                    print(f"Target {target} trial {trial_counts[target]} done")

                    # Print live counters for all targets
                    progress_str = " | ".join(
                        [f"{t}: {trial_counts[t]}/{trials}" for t in target_spreads]
                    )
                    print("Progress:", progress_str, "\n")

    print("All target spreads complete.")
    return selected_syn_groups

def run_with_resets(
    neuron_name,
    sizSection,
    erev,
    somaSection,
    synMap_df,
    syncount,
    trials,
    target_spreads,
    spread_tolerance,
    VPN=None,
    max_runtime_hours=2,
    subset_size=500
):
    trial_counts = {t: 0 for t in target_spreads}
    selected_syn_groups = {t: [] for t in target_spreads}

    while any(trial_counts[t] < trials for t in target_spreads):
        

        # run the V2 logic with this subset
        results = select_syns_by_multiple_spreads(
            neuron_name,
            sizSection,
            erev,
            somaSection,
            synMap_df,
            syncount,
            trials,
            target_spreads=target_spreads,
            spread_tolerance=spread_tolerance,
            VPN=VPN,
            max_runtime_hours=max_runtime_hours,
            subset_size=subset_size
        )

        # merge progress
        for t in target_spreads:
            trial_counts[t] = len(results[t]) + len(selected_syn_groups[t])
            selected_syn_groups[t].extend(results[t])

    return selected_syn_groups

####################
#Running electrotonic distance based synapse selection simulations
def run_all_synapse_exp2_in_folder(folder_path, sizSection, somaSection, erev=None, simTime=100, onsetTime=50):
    """
    Run ALL_SYNAPSE_EXP2 mode on all *_synapse_site.csv files in a given folder.
    
    Parameters
    ----------
    folder_path : str
        Path to the folder containing *_synapse_site.csv files.
    sizSection : hoc Section
        Section corresponding to the spike initiation zone (SIZ).
    somaSection : hoc Section
        Section corresponding to the soma.
    erev : float
        Reversal potential (default -72.5 mV).
    simTime : float
        Simulation duration (ms).
    onsetTime : float
        Synapse activation start time (ms).
    """

    # Find all synapse site CSVs in the folder
    syn_files = glob.glob(os.path.join(folder_path, "*_synapse_site.csv"))

    if not syn_files:
        print(f"No *_synapse_site.csv files found in {folder_path}")
        return

    for syn_file in syn_files:
        base_name = os.path.basename(syn_file).replace(".csv", "")
        print(f"\n[RUN_ALL] Processing {base_name}")

        # Load synapse dataframe
        syn_df = pd.read_csv(syn_file)

        # Run activation
        activateSynapses(
            syn_df,
            MODE="ALL_SYNAPSE_EXP2",
            somaSection=somaSection,
            sizSection=sizSection,
            erev=erev,
            simTime=simTime,
            onsetTime=onsetTime,
            dirname=folder_path,       # save in the same folder
            filename=base_name         # base name -> ensures _stim_summary is saved right
        )

########
#New simulation activating neurons base on receptive fields

def COM_neighborhood_activations(VPN, synMap_df, sizSection, somaSection, erev, neuron_name,
                                  n_trials=10, max_neighbors=9, random_seed=None):
    """
    For each k in 1..max_neighbors, randomly samples n_trials UNIQUE seed neurons,
    finds each seed's k nearest neighbors in DV/AP space, and activates all synapses
    for that group. Combinations are globally unique across all k levels — no two
    simulations activate the same set of neurons.

    The random_seed fully determines all sampling, making runs reproducible.
    Each simulation is tagged with its seed neuron + neighbor IDs for traceability.

    Output structure:
        {neuron_name}_final_sims/COM_neighborhood/
            {VPN}_k{k:02d}/
                trial{t:02d}_seed{seed_id}_k{k}_syns{N}.dat

    Parameters
    ----------
    VPN           : str   — e.g. "LC4", "LPLC2"
    synMap_df     : df    — synapse map with 'pre' column
    sizSection    : obj   — SIZ section for recording
    somaSection   : obj   — soma section
    erev          : float — reversal potential
    neuron_name   : str   — postsynaptic neuron name
    n_trials      : int   — number of unique combinations per k level (default 10)
    max_neighbors : int   — max k nearest neighbors (default 9 → groups of 2–10)
    random_seed   : int   — seed for full reproducibility (required for reuse)
    """

    # ── 1. Load receptive field COM data ──────────────────────────────────────
    vpn_file_map = {
        "LC4":   "datafiles/Receptive_field_files/LC4_receptive_field.xlsx",
        "LC4_hemi":   "datafiles/Receptive_field_files/LC4_hemibrain_centroids.xlsx",
        "LC6":   "datafiles/Receptive_field_files/LC6_receptive_field.xlsx",
        "LC22":  "datafiles/Receptive_field_files/LC22_receptive_field.xlsx",
        "LC22_hemi":  "datafiles/Receptive_field_files/LC22_hemibrain_centroids.xlsx",
        "LPLC1": "datafiles/Receptive_field_files/LPLC1_receptive_field.xlsx",
        "LPLC1_hemi": "datafiles/Receptive_field_files/LPLC1_hemibrain_centroids.xlsx",
        "LPLC2": "datafiles/Receptive_field_files/LPLC2_receptive_field.xlsx",
        "LPLC2_hemi": "datafiles/Receptive_field_files/LPLC2_hemibrain_centroids.xlsx",
        "LPLC4": "datafiles/Receptive_field_files/LPLC4_receptivefield_COM.xlsx",
        "LPLC4_hemi": "datafiles/Receptive_field_files/LPLC4_hemibrain_centroids.xlsx",
    }
    if VPN not in vpn_file_map:
        raise ValueError(f"Invalid VPN: '{VPN}'. Choose from {list(vpn_file_map.keys())}")

    VPN_COM = pd.read_excel(vpn_file_map[VPN], dtype={'updated_ids': str})

    # ── 2. Filter to neurons present in synMap_df ─────────────────────────────
    synMap_df = synMap_df.copy()
    synMap_df['pre'] = synMap_df['pre'].astype(str)
    VPN_COM['updated_ids'] = VPN_COM['updated_ids'].astype(str)

    available_ids = set(synMap_df['pre'].unique())
    VPN_COM_filt = VPN_COM[VPN_COM['updated_ids'].isin(available_ids)].reset_index(drop=True)

    n_neurons = len(VPN_COM_filt)
    if n_neurons < max_neighbors + 1:
        raise ValueError(
            f"Only {n_neurons} neurons with synapses available; "
            f"need at least {max_neighbors + 1}."
        )

    coords = VPN_COM_filt[['DV_raw_um', 'AP_raw_um']].values   # (N, 2)
    ids    = VPN_COM_filt['updated_ids'].values                  # (N,)

    # Pairwise Euclidean distances, shape (N, N)
    dist_matrix = cdist(coords, coords, metric='euclidean')

    # ── 3. Precompute each neuron's sorted neighbor list (excludes self) ───────
    # neighbor_rank[i] = list of neuron indices ordered nearest → farthest from i
    neighbor_rank = {}
    for i in range(n_neurons):
        dists = dist_matrix[i].copy()
        dists[i] = np.inf
        neighbor_rank[i] = list(np.argsort(dists))   # ascending distance, self excluded

    # ── 4. Build globally unique combination pool ──────────────────────────────
    # A combination is defined as a frozenset of neuron indices (seed + k neighbors).
    # We sample across ALL k levels so no two simulations share the same neuron set.
    rng = random.Random(random_seed)

    # For each k, we need n_trials unique seed indices whose k-neighbor groups
    # have not been used at ANY prior k level.
    used_combinations = set()   # frozenset of neuron index tuples, globally unique
    used_seeds        = set()   # seed indices already used at any k (optional strictness)

    # Simulation plan: list of dicts, sorted by k then trial
    sim_plan = []

    for k in range(1, max_neighbors + 1):
        group_size = k + 1

        # Candidate seed indices — shuffle for random ordering
        candidate_seeds = list(range(n_neurons))
        rng.shuffle(candidate_seeds)

        trials_found = 0
        for seed_i in candidate_seeds:
            if trials_found >= n_trials:
                break

            neighbor_indices = neighbor_rank[seed_i][:k]   # k nearest neighbors
            group_indices    = frozenset([seed_i] + list(neighbor_indices))

            # Skip if this exact combination was already used at a prior k
            if group_indices in used_combinations:
                continue

            used_combinations.add(group_indices)

            neighbor_ids = [ids[ni] for ni in neighbor_indices]
            group_ids    = [ids[seed_i]] + neighbor_ids
            distances    = [round(dist_matrix[seed_i, ni], 2) for ni in neighbor_indices]

            sim_plan.append({
                'k':             k,
                'group_size':    group_size,
                'trial':         trials_found + 1,
                'seed_idx':      seed_i,
                'seed_id':       ids[seed_i],
                'neighbor_ids':  neighbor_ids,
                'group_ids':     group_ids,
                'distances_um':  distances,
            })
            trials_found += 1

        if trials_found < n_trials:
            print(f"  WARNING: k={k} — only found {trials_found}/{n_trials} "
                  f"unique combinations (population may be small).")

    # ── 5. Log the full simulation plan ───────────────────────────────────────
    print(f"\nSimulation plan: {len(sim_plan)} total runs | "
          f"VPN={VPN} | seed={random_seed}\n")
    print(f"{'k':>4}  {'trial':>5}  {'group_size':>10}  "
          f"{'seed_id':>20}  {'n_syns_raw':>10}  {'distances_um'}")
    print("-" * 90)

    # ── 6. Execute simulations ─────────────────────────────────────────────────
    print(f"\nSimulation plan: {len(sim_plan)} total runs | VPN={VPN} | seed={random_seed}\n")
    print(f"{'k':>4}  {'trial':>5}  {'group_size':>10}  {'seed_id':>20}  {'n_syns':>10}")
    print("-" * 60)

    for sim in sim_plan:
        k         = sim['k']
        trial     = sim['trial']
        seed_id   = sim['seed_id']
        group_ids = sim['group_ids']

        group_synMap = synMap_df[synMap_df['pre'].isin(group_ids)].copy()
        n_syns = len(group_synMap)

        print(f"{k:>4}  {trial:>5}  {sim['group_size']:>10}  {seed_id:>20}  {n_syns:>10}")

        if group_synMap.empty:
            print(f"         ↳ WARNING: no synapses for group {group_ids}, skipping.")
            continue

        dirname  = (f"{neuron_name}_final_sims/COM_neighborhood/"
                    f"{VPN}_k{k:02d}")
        filename = (f"trial{trial:02d}_seed{seed_id}")

        tVec, vVec = activateSynapses(
            group_synMap,
            MODE="GROUP_EXP2_no_dend_recording",
            somaSection=somaSection,
            onsetTime=50,
            sizSection=sizSection,
            erev=erev,
            simTime=100,
            neuron_name=neuron_name,
            dirname=dirname,
            filename=filename
        )

    print(f"\n✓ Done. {len(sim_plan)} simulations complete.")
    print(f"  Reproduce exactly with: random_seed={random_seed}")


########

########
#New function for electrotonic distance between synapses, to replace physical distance/synapse spread
DEFAULT_RECEPTIVE_FIELD_DIR = Path("datafiles/Receptive_field_files")

def generate_neuron_pairs_csv(
    neuron_name: str,
    synMap_df: pd.DataFrame,
    VPN: str,
    json_path: str | Path | None = None,
    receptive_field_dir: str | Path = DEFAULT_RECEPTIVE_FIELD_DIR,
    output_root: str | Path = "datafiles/simulationData") -> pd.DataFrame:
    """
    Generate every unique unordered pair of VPN neurons present in both:

    1. ``synMap_df`` (using its ``updated_id`` column), and
    2. the selected VPN's receptive-field JSON (using ``updated_ids``).

    By default, the JSON is selected automatically from ``VPN`` as:

    datafiles/Receptive_field_files/{VPN}_receptive_field_with_polygons.json

    Supply ``json_path`` only when a particular VPN uses a different filename.

    The output is written exactly where ``select_syns_by_neuron_pairs`` expects:

    datafiles/simulationData/{neuron_name}_final_sims/
    {neuron_name}_pairwise_{VPN}/{neuron_name}_{VPN}_neuron_pairs.csv

    Returns
    -------
    pandas.DataFrame
        DataFrame with the columns ``neuron1`` and ``neuron2``.
    """
    if not neuron_name or not neuron_name.strip():
        raise ValueError("neuron_name cannot be empty.")
    if not VPN or not VPN.strip():
        raise ValueError("VPN cannot be empty.")

    neuron_name = neuron_name.strip()
    VPN = VPN.strip()
    receptive_field_dir = Path(receptive_field_dir)
    if json_path is None:
        json_path = receptive_field_dir / f"{VPN}_receptive_field_with_polygons.json"
    else:
        json_path = Path(json_path)
    output_root = Path(output_root)

    # if "updated_id" not in synMap_df.columns:
    #     raise ValueError("synMap_df must contain an 'updated_id' column.")
    if not json_path.exists():
        raise FileNotFoundError(f"Receptive-field JSON not found: {json_path}")

    with json_path.open("r", encoding="utf-8") as file:
        receptive_fields = json.load(file)

    if not isinstance(receptive_fields, list):
        raise ValueError("The receptive-field JSON must contain a list of records.")

    # Preserve the order in the JSON while removing missing and duplicate IDs.
    json_ids = []
    seen = set()
    for record in receptive_fields:
        updated_id = record.get("updated_ids")
        if updated_id is None:
            continue
        updated_id = str(updated_id).strip()
        if updated_id and updated_id not in seen:
            json_ids.append(updated_id)
            seen.add(updated_id)

    syn_map = synMap_df.copy()
    syn_map["updated_id"] = syn_map["pre"].astype(str).str.strip()

    # If a type column is available, restrict the synapse map to the requested
    # VPN before selecting neuron IDs.
    if "type" in syn_map.columns:
        syn_map = syn_map[syn_map["type"].astype(str) == VPN]

    synmap_ids = set(syn_map["updated_id"].dropna())
    matched_ids = [neuron_id for neuron_id in json_ids if neuron_id in synmap_ids]

    if len(matched_ids) < 2:
        raise ValueError(
            f"Only {len(matched_ids)} {VPN} neuron(s) matched between the JSON "
            "and synMap_df; at least 2 are required."
        )

    pairs_df = pd.DataFrame(
        combinations(matched_ids, 2),
        columns=["neuron1", "neuron2"],
    )

    output_dir = (
        output_root
        / f"{neuron_name}_final_sims"
        / f"{neuron_name}_pairwise_{VPN}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{neuron_name}_{VPN}_neuron_pairs.csv"
    pairs_df.to_csv(output_path, index=False)

    missing_from_synmap = [
        neuron_id for neuron_id in json_ids if neuron_id not in synmap_ids
    ]
    print(f"Matched {len(matched_ids)} {VPN} neurons.")
    print(f"Generated {len(pairs_df)} unique neuron pairs.")
    print(f"Saved neuron pairs to: {output_path.resolve()}")
    if missing_from_synmap:
        print(
            "Warning: "
            f"{len(missing_from_synmap)} JSON neuron(s) were absent from synMap_df "
            "and were excluded."
        )

    return pairs_df

def select_syns_by_neuron_pairs(
    neuron_name,
    sizSection,
    erev,
    somaSection,
    synMap_df,
    VPN=None,
    num_instances=10,
    instance_id=0,
):
    """
    Run pairwise synapse simulations for a neuron, resuming from previous runs
    and splitting workload across multiple instances.
    """
    dirname = f"{neuron_name}_final_sims/{neuron_name}_pairwise_{VPN}"
    # os.makedirs(os.path.join("datafiles", "simulationData", dirname), exist_ok=True)
    
    # Load neuron pairs
    pair_csv = os.path.join(
        "datafiles", "simulationData", dirname,
        f"{neuron_name}_{VPN}_neuron_pairs.csv"
    )
    if not os.path.exists(pair_csv):
        raise FileNotFoundError(f"Neuron pairs CSV not found: {pair_csv}")
    
    neuron_pairs_df = pd.read_csv(pair_csv, sep=None, engine="python", dtype=str)
    
    # Split workload by instance
    total_pairs = len(neuron_pairs_df)
    chunk_size = (total_pairs + num_instances - 1) // num_instances
    start_idx = instance_id * chunk_size
    end_idx = min(start_idx + chunk_size, total_pairs)
    subset_pairs_df = neuron_pairs_df.iloc[start_idx:end_idx]
    
    print(f"Instance {instance_id} will run {len(subset_pairs_df)} pairs out of {total_pairs}.")

    # Filter synMap for relevant neurons
    unique_neurons = pd.unique(subset_pairs_df[['neuron1','neuron2']].values.ravel())
    synMap_df = synMap_df.copy()
    synMap_df['pre'] = synMap_df['pre'].astype(str)
    synMap_df = synMap_df.reset_index(drop=True)
    synMap_df['syn_index'] = synMap_df.index.astype(int)
    filtered_synMap = synMap_df[synMap_df['pre'].isin(unique_neurons)]
    
    synapse_groups = {pre_id: filtered_synMap[filtered_synMap['pre']==pre_id]
                      for pre_id in filtered_synMap['pre'].unique()}
    
    # Loop through neuron pairs
    for neuron1, neuron2 in subset_pairs_df.to_numpy():
        neuron1 = str(neuron1)
        neuron2 = str(neuron2)
        if neuron1 not in synapse_groups or neuron2 not in synapse_groups:
            continue
        
        summary_file = os.path.join(
            "datafiles", "simulationData", dirname,
            f"neurons_{neuron1}_{neuron2}_synpair_summary.csv"
        )
        log_file = os.path.join(
            "datafiles", "simulationData", dirname,
            f"{neuron1}_{neuron2}_progress.log"
        )
        
        # Cartesian product of synapses
        synapses_neuron1 = synapse_groups[neuron1]
        synapses_neuron2 = synapse_groups[neuron2]
        all_pairs = [
            (int(s1['syn_index']), int(s2['syn_index']))
            for _, s1 in synapses_neuron1.iterrows()
            for _, s2 in synapses_neuron2.iterrows()
        ]
        
        # Load completed pairs
        completed_pairs = set()
        if os.path.exists(summary_file):
            summary_df = pd.read_csv(summary_file)
            for _, row in summary_df.iterrows():
                try:
                    idx1 = int(row['stim_syn_index'])
                    idx2 = int(row['record_syn_index'])
                    completed_pairs.add(tuple(sorted((idx1, idx2))))
                except Exception as e:
                    print(f"Warning reading summary row: {row}, error: {e}")
        
        # Only run remaining pairs
        remaining_pairs = [pair for pair in all_pairs if pair not in completed_pairs]
        print(remaining_pairs)
        print(f"[Instance {instance_id}] Neuron pair {neuron1}&{neuron2}: {len(remaining_pairs)} remaining pairs")
        
        for i, (idx1, idx2) in enumerate(remaining_pairs, start=1):
            try:
                syn1 = synapses_neuron1[synapses_neuron1['syn_index']==idx1].iloc[0]
                syn2 = synapses_neuron2[synapses_neuron2['syn_index']==idx2].iloc[0]
                syn_pair_df = pd.DataFrame([syn1, syn2])
                
                # Call your simulation function
                activateSynapses(
                    syn_pair_df,
                    MODE="PAIRWISE_EXP2_updated",
                    somaSection=somaSection,
                    sizSection=sizSection,
                    onsetTime=25,
                    erev=erev,
                    simTime=35,
                    neuron_name=neuron_name,
                    dirname=dirname,
                    filename=f"neurons_{neuron1}_{neuron2}_synpair"
                )
                
                # Log progress
                msg = f"[Instance {instance_id}] Simulated pair: ({idx1},{idx2}) Progress: {i}/{len(all_pairs)}"
                print(msg)
                with open(log_file, "a") as f:
                    f.write(msg + "\n")
                
            except Exception as e:
                err_msg = f"Error simulating pair ({idx1},{idx2}): {e}"
                print(err_msg)
                with open(log_file, "a") as f:
                    f.write(err_msg + "\n")
                    f.write(traceback.format_exc() + "\n")
                continue
        
        print(f"✅ [Instance {instance_id}] Completed all missing synapse pairs for neuron pair {neuron1}&{neuron2}\n")


########

def main():
    neuron_name = "DNp01"

    Tk().withdraw()
    fd_title = "Select morphology file to initialize"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"datafiles/morphologyData", title=fd_title)
    
    if neuron_name == "DNp01": 
        #old preprint passive properties based on average of 3 recordings, this was not corrected for LJP
        # erev = -66.2298
        # raVal = 212
        # gleakVal = 1/2300
        # cmVal = 0.7

        #new fit parameters based on avg of 7 recordings (see passivePropVisualizer.py for other parameters, for hemibrain models or manuscript)
        initial = -76.75
        erev = -77.5
        raVal = 250
        gleakVal = 1/2675
        cmVal = 0.8

        #Lower bound parameters from fitting fly 3
        # initial = -87.3
        # erev = -88.05
        # raVal = 250
        # gleakVal = 1/2820
        # cmVal = 0.50

        #Upper bound parameters from fitting fly 5
        # initial = -66.3
        # erev = -66.94
        # raVal = 300
        # gleakVal =  1/3275
        # cmVal = 0.79


    elif neuron_name == "DNp01_hemi": 
        initial = -76.75
        erev = -77.5
        raVal = 190
        gleakVal = 1/1900
        cmVal = 0.85
        

    elif neuron_name == "DNp02":
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

    elif neuron_name == "DNp03":
        #old preprint passive properties based on average of 3 recordings, this was not corrected for LJP
        initial = -61.15
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

        #new fit parameters based on avg of 8 recordings (see passivePropVisualizer.py for other parameters)
        # initial = -71.35
        # raVal = 55
        # gleakVal = 1/2725
        # cmVal = 1.1
        # erev = -73.71

        #Upper bound parameters from fitting fly 2
        # initial = -68.4
        # raVal = 50
        # gleakVal = 1/3150
        # cmVal = 0.8
        # erev = -70.70

        #Lower bound parameters from fitting fly 6
        # initial = -73.15
        # raVal = 30
        # gleakVal = 1/2700
        # cmVal = 2
        # erev = -74.8

    elif neuron_name == "DNp03_hemi":
        initial = -71.35
        raVal = 320
        gleakVal = 1/3755
        cmVal = 0.7
        erev = -73.73

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

    cell, allSections_py, allSections_nrn, somaSection, sizSection, erev, axonList, tetherList, dendList, electrodeSec = initializeModel(neuron_name, morph_file, coreneuron_flag = True, resting_ev= initial,  capacitance = cmVal, axial_resistivity = raVal, gleak_conductance= gleakVal, erev = erev)
    change_Ra(ra=raVal, electrodeSec=electrodeSec, electrodeVal = elec_raVal)
    change_gLeak(gleak=gleakVal, erev=erev, electrodeSec=electrodeSec, electrodeVal = elec_gleakVal)
    change_memCap(memcap=cmVal, electrodeSec=electrodeSec, electrodeVal = elec_cmVal)
    
    ########################################

    # Loading of synapse data

    synMap_df = loadSynMapDataframe(neuron_name)

    synMap_df = removeAxonalSynapses_REV(synMap_df, dendList, axonList=axonList, tetherList=tetherList, somaSection=somaSection)

    synMap_VPN_df = filterSynapsesToCriteria(['LC4', 'LC6', 'LC22', 'LPLC1', 'LPLC2', 'LPLC4'], synMap_df)

    #How to subset to a given population
    synMap_VPN_df_LPLC2 = filterSynapsesToCriteria(['LPLC2'], synMap_df)


    


    #######################################
  
    #Simulations for figure 4 supplement-2 
    #DNp01- LC4 and LPLC2
    #DNp03- LC4, LC22, LPLC1, LPLC4, LPLC2
    #Example using LPLC2 synapses
#     generate_neuron_pairs_csv(
#     neuron_name=neuron_name,
#     synMap_df=synMap_VPN_df_LPLC2,
#     VPN="LPLC2",
# )
#     select_syns_by_neuron_pairs(
#     neuron_name,
#     sizSection,
#     initial,
#     somaSection,
#     synMap_VPN_df_LPLC2,
#     VPN="LPLC2",
#     num_instances=1,
#     instance_id=0,
# )
    # quit()

    ########################################

    ###Loading in a particular synapse dataframe for single synapse activations
    ### Relation to Figure 8
    # activateAllSingle(neuron_name, initial, synMap_VPN_df, sizSection, somaSection)
    # quit()

    # # ########################################

    # # #Simulations for final figures
    # # ### Relation to Figure 9, Random vs partner syn activations

    # VPN_paired_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, initial, shuffled = False)

    # VPN_random_activations(neuron_name, synMap_VPN_df, sizSection, somaSection, initial, numsyn_lim = 60)



    # # #Runs simulations of synapses in increasing number of synapses either up to 250, or for all VPN synapses within a VPN type
    # #In relation to figure 9
    # activateLinearVPN(neuron_name, initial, synMap_VPN_df, sizSection, somaSection)

    #New simulations based on retinotopic activation of VPN neurons
    #simulations dependent on which modeling is being run, 

    # COM_neighborhood_activations(VPN= "LC4", synMap_df= synMap_df, sizSection= sizSection, somaSection= somaSection, erev= initial, neuron_name  = neuron_name, n_trials= 20, max_neighbors = 40, random_seed= 42)
    # COM_neighborhood_activations(VPN= "LPLC2", synMap_df= synMap_df, sizSection= sizSection, somaSection= somaSection, erev= initial, neuron_name  = neuron_name, n_trials= 20, max_neighbors = 40, random_seed= 42)
    # COM_neighborhood_activations(VPN= "LPLC1", synMap_df= synMap_df, sizSection= sizSection, somaSection= somaSection, erev= initial, neuron_name  = neuron_name, n_trials= 20, max_neighbors = 40, random_seed= 42)
    # COM_neighborhood_activations(VPN= "LPLC4", synMap_df= synMap_df, sizSection= sizSection, somaSection= somaSection, erev= initial, neuron_name  = neuron_name, n_trials= 20, max_neighbors = 40, random_seed= 42)
    # COM_neighborhood_activations(VPN= "LC22", synMap_df= synMap_df, sizSection= sizSection, somaSection= somaSection, erev= initial, neuron_name  = neuron_name, n_trials= 20, max_neighbors = 40, random_seed= 42)


    #In relation to figure 10
    syn_list_DNp01 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 22, 25, 29, 31] 
    syn_list_DNp03 = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 29, 32] 
    syn_list_DNp01_hemi = [5, 11, 12, 15, 16, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 36, 39, 47] 
    syn_list_DNp03_hemi = [6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 27, 30, 38] 


    #Runs close simulations and then runs electrotonic simulations for those activated synapses
    # trial_count = 150
    # for x in syn_list_DNp01:
    #     activate_close_syns(neuron_name, initial, synMap_VPN_df, sizSection, somaSection, numSyn =x, trials=150, onsetTime=30)

    #     folder_path = f"datafiles/simulationData/{neuron_name}_final_sims/CloseSims/{neuron_name}_close{x}record{trial_count}" #for all close syns
    #     run_all_synapse_exp2_in_folder(
    #         folder_path,
    #         sizSection,
    #         somaSection,
    #         erev=initial,
    #         simTime=55,
    #         onsetTime=30
    #     )

    #for Partner neurons/randomly selected synapses
    # folder_path = f"datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner" #change partner to rand for random syns
    # run_all_synapse_exp2_in_folder(
    #     folder_path,
    #     sizSection,
    #     somaSection,
    #     erev=initial,
    #     simTime=55,
    #     onsetTime=30
    # )


    ##################
    #Runs simulations based on the synapse spread targets and the tolerance from that target spread, and how many trials to try and run for each target spread
    #Note to run the below, the list cannot have 1 as a synapse count
    for x in syn_list_DNp01:
        run_with_resets(
        neuron_name,
        sizSection,
        erev,
        somaSection,
        synMap_VPN_df,
        syncount = x,
        trials = 50,
        target_spreads = [10, 20, 30, 40, 50, 60, 70],
        spread_tolerance = 10,
        VPN="LPLC2", #if none, changes label to "Random", otherwise VPN label is for individual populations
        max_runtime_hours=0.008,
        subset_size=100
        )
        
    #Electrotonic synapse activations:
    target_spread = [10, 20, 30, 40, 50, 60, 70]
    for x in target_spread:
        for y in syn_list_DNp01:
            folder_path = f"datafiles/simulationData/{neuron_name}_final_sims/synapse_spread/{neuron_name}_{y}_syns_rand_by_syn_spread{x}" #for all randomly selected syns
            try:
                run_all_synapse_exp2_in_folder(
                    folder_path,
                    sizSection,
                    somaSection,
                    erev=initial,
                    simTime=55,
                    onsetTime=30
                )
            except Exception as e:
                print(f"[WARNING] Skipping neuron {y} with spread {x}: {e}")
                continue

    

main()
    
