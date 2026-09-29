from neuron import h, gui
from neuron.units import ms, mV
import time as clock
import os
from itertools import combinations
import numpy
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import math
import numpy as np
import pandas as pd
from tkinter import Tk
import tkinter.filedialog as fd
import random
from scipy import stats
from matplotlib.lines import Line2D
from collections import defaultdict
import re
import seaborn as sns
from scipy.stats import levene, shapiro, linregress, ttest_1samp
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

#Functions are seperated into particular blocks of code and grouped to help troubleshoot.

#2 = axon | 3 = dendrite | 0 = unlabeled | 1 = soma
h.load_file("stdrun.hoc")


#########
#Initialization of models
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
    
def change_Ra(ra=28.08):
    for sec in h.allsec():
        sec.Ra = ra

def change_gLeak(gleak=5e-5, erev=-72.5):
    for sec in h.allsec():
        sec.insert('pas')
        for seg in sec:
            seg.pas.g = gleak
            seg.pas.e = erev

def change_memCap(memcap=1):
    for sec in h.allsec():
        sec.cm = memcap

def nsegDiscretization(sectionListToDiscretize):
    #this function iterates over every section, calculates its spatial constant (lambda), and checks if the length of the segments within this section are less than 1/10  lambda
    #if true, nothing happens

    for sec in sectionListToDiscretize:
        #by calling sec.psection, we can return the density mechanisms for the given section in a dictionary
        #this lets us access the gleak value for a given section, which we use to calculate the membrane resistance
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

def initializeModel(neuron_name):
    Tk().withdraw()
    fd_title = "Select morphology file to initialize"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"/morphologyData", title=fd_title)

    cell = instantiate_swc(morph_file)

    allSections_nrn = h.SectionList()
    for sec in h.allsec():
        allSections_nrn.append(sec=sec)
    
    # Create a Python list from this SectionList
    # Select sections from the list by their index

    allSections_py = [sec for sec in allSections_nrn]    



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


    # allSections_py = createAxon(axonEnd, allSections_py, neuron_name)

    if neuron_name == "DNp01":
        change_Ra(34)#33.2617)
        change_gLeak(3.4e-4, erev=-66.4298 - 0.2)#4.4e-9)    
        change_memCap(1)#1.4261)
        erev = -66.4298 - 0.2
    elif neuron_name == "DNp01_hemi":
        erev = -66.5
        change_Ra(60)
        gleakVal = 1/2000
        change_gLeak(1/2000, erev=-66.5) 
        change_memCap(1)
    elif neuron_name == "DNp02":# or neuron_name == "DNp06":
        change_Ra(91)
        change_gLeak(0.0002415, erev=-70.8)  
        change_memCap(1)
        erev=-70.8
    elif neuron_name == "DNp03":
        change_Ra(50)
        change_gLeak(gleak=1/3150, erev=-61.15)
        change_memCap(1)
        erev=-61.15
    elif neuron_name == "DNp06":
        change_Ra(91)
        change_gLeak(0.0002415, erev=-60)  
        change_memCap(1)
        erev = -60#-60.0406 SEE DNP02
    else:
        change_Ra()
        change_gLeak()
        change_memCap()
        erev=-72.5

    nsegDiscretization(allSections_py)

    return cell, allSections_py, allSections_nrn, somaSection, erev, axonList, tetherList, dendList, sizSection#, erev, axonList

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
        equivCylHeight = 241.69
        equivCylDiam = 3.32*2
    elif neuron_name == "DNp01_hemi":
        equivCylHeight = 241.69
        equivCylDiam = 3.32*2
    elif neuron_name == "DNp02":
        # equivCylHeight = 220.654/5
        # equivCylDiam = 0.624*2
        equivCylHeight = 171.521
        equivCylDiam = 0.658*2
    elif neuron_name == "DNp03":
        # equivCylHeight = 180.14
        # equivCylDiam = 0.249*2
        equivCylHeight = 175.671
        equivCylDiam = 0.251*2
    elif neuron_name == "DNp03_hemi":
        equivCylHeight = 375.671
        equivCylDiam = 1.3
    elif neuron_name == "DNp04":
        #according to Namiki paper, DNp04 axon looks to be as long as DNp01 and as thick as DNp03, so we
        #use the DNp01 height variable and the DNp03 diam variable to create the approximate axon
        equivCylHeight = 360.346/10
        equivCylDiam = 0.249*2
    elif neuron_name == "DNp06":
        equivCylHeight = 222.766/5
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

#########

#########
#Helper functions for visualizations/and analysis

def plotTraces(t_df, v_df, colorCode=None, FIGURE=None, gridPos=None, axList=None):


    tVecList = []
    vSynList = []
    # synSiteList = []

    ct = 0

    if isinstance(t_df, pd.Series):
        tVecList.append(t_df)
        vSynList.append(v_df)
    else:
        for index, row in t_df.iterrows():
            try:
                tVecList.append(row)
                vSynList.append(v_df.iloc[ct])
                # synSiteList.append(synSite_df.iloc[ct])
                ct +=1
            except IndexError:
                break
    if FIGURE is None:
        FIGURE = plt.figure()
        if gridPos is not None:
            axList = []
            gs_fig = gridspec.GridSpec(gridPos[0], 1)
            for x in range(gridPos[0]):
                axList.append(FIGURE.add_subplot(gs_fig[x, :]))
                # ax2_fig = fig.add_subplot(gs_fig3[1, :])
    else:
        FIGURE = plt.figure(FIGURE)
    for trial in range(len(tVecList)):
        vSyn = vSynList[trial]#.to_python()
        vSyn_np = numpy.array(vSyn)
        tVec = tVecList[trial]#.to_python()
        tVec_np = numpy.array(tVec)
        if colorCode is None:
            if gridPos is not None:
                axList[gridPos[1]].plot(tVec_np, vSyn_np)
            else:
                plt.plot(tVec_np, vSyn_np)#, 'k-')
        else:
            if colorCode[0] == '#':
                axList[gridPos[1]].plot(tVec_np, vSyn_np, '-', color=colorCode)
            else:
                if gridPos is not None:
                    axList[gridPos[1]].plot(tVec_np, vSyn_np, colorCode)
                else:
                    plt.plot(tVec_np, vSyn_np, colorCode)

    if axList is None:
        return FIGURE
    else:
        return FIGURE, axList
    
def plotSynapseLocations(syn_df, colorKey=2, syn_shape_window=None, sizeArg=0.5, sizSection=None):

    if syn_shape_window is None:
        syn_shape_window = h.Shape()
    else:
        syn_shape_window = syn_shape_window

    syn_shape_window.color_all(9)

    if isinstance(syn_df, pd.Series):
        currSynSecStr = syn_df.loc['mappedSection']
        for sec in h.allsec():
            if sec.name() == currSynSecStr:
                currSynSec = sec
                break
        synTemp = h.AlphaSynapse(currSynSec(syn_df.loc['mappedSegRangeVar']))
        syn_shape_window.point_mark(synTemp, colorKey, "O", sizeArg)
    else:
        for index, row in syn_df.iterrows():
            if str(row.loc['pre']) == 'nan':#numpy.isnan(row.loc['pre']):
                pass
            else:
                currSynSecStr = row.loc['mappedSection']
                for sec in h.allsec():
                    if sec.name() == currSynSecStr:
                            currSynSec = sec
                            break
                synTemp = h.AlphaSynapse(currSynSec(row.loc['mappedSegRangeVar']))
                syn_shape_window.point_mark(synTemp, colorKey, "O", sizeArg)

    if sizSection is not None:
        synTemp = h.AlphaSynapse(sizSection(0.05))
        syn_shape_window.point_mark(synTemp, 7, "O", 12)

    syn_shape_window.exec_menu("Shape Plot")  # Ensure window appears
    return syn_shape_window

def str2sec(sec_as_str):
    for sec in h.allsec():
        if sec.name() == sec_as_str:
            sec_as_sec = sec
            break
    return sec_as_sec

def filterSynapsesToCriteria(subsetInfo, synMap_df):
    filteredSyn_List = []
    indexList = []
    for index, row in synMap_df.iterrows():
        if row.loc['type'] in subsetInfo:
            # print(row)
            filteredSyn_List.append(row)
            indexList.append(index)
    filteredSyn_df = pd.DataFrame(filteredSyn_List)
    # quit(0)
    return filteredSyn_df, indexList

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

def removeAxonalSynapses(syn_df, dendList, axonList=None, tetherList=None, somaSection=None, sizSection=None):
    synSecSeries = syn_df['mappedSection']
    tetherSynSites = syn_df[synSecSeries.str.startswith('dend_11')]
    somaSynSites = syn_df[synSecSeries.str.startswith('soma')]
    axonSynSites = syn_df[synSecSeries.str.startswith('axon')]
    dendSynSites = syn_df[synSecSeries.str.startswith('dend[')]

    if sizSection is not None:
        dendwindow = plotSynapseLocations(dendSynSites, colorKey=3, sizSection=sizSection)
    else:
        dendwindow = plotSynapseLocations(dendSynSites, colorKey=3)
    return dendSynSites

##########

#Loading of data functions based on simulation outputs
def iterQuintsV3(dir_path):
    all_files = os.listdir(dir_path)

    # Full list of allowed presynaptic types
    match = re.compile(
        r'(LCs_4|LPLCs_4|LC4|LC4_2|LC4_3|LC4_4|LPLC1|LPLC1_2|LPLC1_3|LPLC1_4|LPLC2|LPLC2_2|LPLC2_3|LPLC2_4|LPLC4|LPLC4_2|LPLC4_3|LPLC4_4|LC22|LC22_2|LC22_3|LC22_4|LC6|random|random2)_trial_(\d+)_(.*)\.csv'
    )

    trial_dict = {}

    for fname in all_files:
        m = match.match(fname)
        if not m:
            continue
        presyn, trial_num, kind = m.groups()
        key = (presyn, int(trial_num))
        trial_dict.setdefault(key, {})[kind] = fname

    for key in sorted(trial_dict.keys()):
        presyn, trial_num = key
        files = trial_dict[key]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping trial {key} — missing required files.")
            continue

        try:
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)

            yield (time_df, vSIZ, vSoma, synSite_df)

        except Exception as e:
            print(f"❌ Error loading trial {key}: {e}")
            continue

def iterQuintsV4(dir_path):
    all_files = os.listdir(dir_path)

    pattern = re.compile(
        r"syns_(\d+)(?:_trial_(\d+)|_(\d+))_(synapse_site|synaptic_currents_siz|synaptic_currents_soma|time)\.csv"
    )

    group_dict = {}

    for fname in all_files:
        m = pattern.match(fname)
        if not m:
            continue

        group_id, trial_str, alt_id, kind = m.groups()
        # Use trial number if available, else use alt_id
        unique_id = trial_str if trial_str is not None else alt_id
        key = (group_id, int(unique_id))
        group_dict.setdefault(key, {})[kind] = fname

    for key in sorted(group_dict.keys()):
        group_id, unique_id = key
        files = group_dict[key]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping syns_{group_id}_{unique_id} — missing required files.")
            continue

        try:
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)

            yield (time_df, vSIZ, vSoma, synSite_df)

        except Exception as e:
            print(f"❌ Error reading files for syns_{group_id}_{unique_id}: {e}")
            continue

def iterQuintsV5(dir_path):
    """
    Handles filenames like: LC4(lVLPT8)_R_25syn_synapse_site.csv
    """
    all_files = os.listdir(dir_path)

    pattern = re.compile(
        r'(.+?)_(\d+)syn_(synapse_site|synaptic_currents_siz|synaptic_currents_soma|time)\.csv'
    )

    group_dict = {}
    for fname in all_files:
        m = pattern.match(fname)
        if not m:
            continue
        presyn_name, group_num, kind = m.groups()
        key = (presyn_name, int(group_num))
        group_dict.setdefault(key, {})[kind] = fname

    for key in sorted(group_dict.keys(), key=lambda x: (x[0], x[1])):
        presyn_name, group_num = key
        files = group_dict[key]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping {presyn_name}_{group_num}syn — missing required files.")
            continue

        try:
            print(f"Loading: {presyn_name}_{group_num}syn", flush=True)
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)

            yield (time_df, vSIZ, vSoma, synSite_df)

        except Exception as e:
            print(f"❌ Error loading {presyn_name}_{group_num}syn: {e}")
            continue

def iterQuintsV6(dir_path):
    """
    Handles filenames like: syns_1_720575940608088259_synapse_site.csv
    """
    all_files = os.listdir(dir_path)

    pattern = re.compile(
        r'(syns_\d+_\d+)_(synapse_site|synaptic_currents_siz|synaptic_currents_soma|time)\.csv'
    )

    group_dict = {}
    for fname in all_files:
        m = pattern.match(fname)
        if not m:
            continue
        group_id, kind = m.groups()
        key = group_id
        group_dict.setdefault(key, {})[kind] = fname

    for key in sorted(group_dict.keys()):
        files = group_dict[key]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping {key} — missing required files.")
            continue

        try:
            print(f"Loading: {key}", flush=True)
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)

            yield (time_df, vSIZ, vSoma, synSite_df)

        except Exception as e:
            print(f"❌ Error loading {key}: {e}")
            continue

def iterQuintsV7(dir_path):
    """
    Handles filenames like:
    randomtrial_144_synapse_site.csv
    randomtrial_144_synaptic_currents_siz.csv
    randomtrial_144_synaptic_currents_soma.csv
    randomtrial_144_time.csv
    """
    all_files = os.listdir(dir_path)

    # Regex for "name trialnum kind"
    pattern = re.compile(
        r'([A-Za-z0-9]+)trial_(\d+)_(synapse_site|synaptic_currents_siz|synaptic_currents_soma|time)\.csv'
    )

    group_dict = {}
    for fname in all_files:
        m = pattern.match(fname)
        if not m:
            continue
        presyn_name, trial_num, kind = m.groups()
        key = (presyn_name, int(trial_num))
        group_dict.setdefault(key, {})[kind] = fname

    for key in sorted(group_dict.keys(), key=lambda x: (x[0], x[1])):
        presyn_name, trial_num = key
        files = group_dict[key]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping {presyn_name}trial_{trial_num} — missing required files: {set(required) - set(files.keys())}")
            continue

        try:
            print(f"Loading: {presyn_name}trial_{trial_num}", flush=True)
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)

            yield (time_df, vSIZ, vSoma, synSite_df)

        except Exception as e:
            print(f"❌ Error loading {presyn_name}trial_{trial_num}: {e}")
            continue

def iterQuintsV8(dir_path):
    all_files = os.listdir(dir_path)

    # Regex for presynaptic files
    match = re.compile(
        r'(LCs_4|LPLCs_4|LC4|LC4_2|LC4_3|LC4_4|LPLC1|LPLC1_2|LPLC1_3|LPLC1_4|LPLC2|LPLC2_2|LPLC2_3|LPLC2_4|LPLC4|LPLC4_2|LPLC4_3|LPLC4_4|LC22|LC22_2|LC22_3|LC22_4|LC6|random|random2|random3||random4|random_1|random_2|random_3|)_trial_(\d+)_(.*)\.csv'
    )

    trial_dict = {}

    for fname in all_files:
        m = match.match(fname)
        if not m:
            continue
        presyn, trial_num, kind = m.groups()
        key = (presyn, int(trial_num))
        trial_dict.setdefault(key, {})[kind] = fname

    for key in sorted(trial_dict.keys()):
        presyn, trial_num = key
        files = trial_dict[key]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time', 'synapse_site_stim_summary']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping trial {key} — missing required files.")
            continue

        try:
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)
            stim_summary_df = pd.read_csv(os.path.join(dir_path, files['synapse_site_stim_summary']))

            yield (time_df, vSIZ, vSoma, synSite_df, stim_summary_df)

        except Exception as e:
            print(f"❌ Error loading trial {key}: {e}")
            continue

def iterQuintsV9(dir_path):
    all_files = os.listdir(dir_path)

    match = re.compile(r'randomtrial_(\d+)_(.*)\.csv')

    trial_dict = {}

    for fname in all_files:
        m = match.match(fname)
        if not m:
            continue
        trial_num, kind = m.groups()
        key = int(trial_num)
        trial_dict.setdefault(key, {})[kind] = fname

    for trial_num in sorted(trial_dict.keys()):
        files = trial_dict[trial_num]

        required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time', 'synapse_site_stim_summary']
        if not all(r in files for r in required):
            print(f"⚠️ Skipping trial {trial_num} — missing required files.")
            continue

        try:
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)
            stim_summary_df = pd.read_csv(os.path.join(dir_path, files['synapse_site_stim_summary']))

            yield (time_df, vSIZ, vSoma, synSite_df, stim_summary_df)

        except Exception as e:
            print(f"❌ Error loading trial {trial_num}: {e}")
            continue

def iterQuintsV10(dir_path):
    """
    Partner iterator like V9 style, keeping multiple files separate:
    - Skips empty '_synaptic_currents_synpt.csv' files
    - Groups each quint of files: (time_df, vSIZ_np, vsoma_np, synSite_df, stim_summary_df)
    """

    all_files = os.listdir(dir_path)

    # Skip intentionally empty synaptic currents files
    all_files = [f for f in all_files if not f.endswith("_synaptic_currents_synpt.csv")]

    # Regex to parse partner files: syns_<trial>_<id>_<kind>.csv
    pattern = re.compile(r"syns_(\d+)_(\d+)_(.*)\.csv")

    # Group files by trial number
    trial_dict = {}
    for fname in all_files:
        m = pattern.match(fname)
        if not m:
            continue
        trial_num, unique_id, kind = m.groups()
        trial_num = int(trial_num)
        trial_dict.setdefault(trial_num, {}).setdefault(kind, []).append(fname)

    for trial_num in sorted(trial_dict.keys()):
        files = trial_dict[trial_num]

        required_types = ['time', 'synaptic_currents_siz', 'synaptic_currents_soma', 
                          'synapse_site', 'synapse_site_stim_summary']
        if not all(r in files for r in required_types):
            print(f"⚠️ Skipping partner trial {trial_num} — missing required file types.")
            continue

        # There may be multiple files per type; iterate over them by index
        max_count = max(len(files[r]) for r in required_types)
        for i in range(max_count):
            try:
                # Pick the i-th file of each type, or the first if fewer exist
                def pick_file(kind):
                    return files[kind][i] if i < len(files[kind]) else files[kind][0]

                time_df = pd.read_csv(os.path.join(dir_path, pick_file('time')), header=None).astype(float)
                vSIZ_np = pd.read_csv(os.path.join(dir_path, pick_file('synaptic_currents_siz')), header=None).to_numpy()
                vsoma_np = pd.read_csv(os.path.join(dir_path, pick_file('synaptic_currents_soma')), header=None).to_numpy()
                synSite_df = pd.read_csv(os.path.join(dir_path, pick_file('synapse_site')))
                stim_summary_df = pd.read_csv(os.path.join(dir_path, pick_file('synapse_site_stim_summary')))

                yield (time_df, vSIZ_np, vsoma_np, synSite_df, stim_summary_df)

            except Exception as e:
                print(f"⚠️ Skipping partner trial {trial_num} group {i} due to error: {e}")
                continue

def iterQuintsV11(dir_path):
    """
    Iterates over simulation trials where filenames follow:
    syns_<count>_trial_<n>_<kind>.csv
    Required kinds: synapse_site, synaptic_currents_siz,
                    synaptic_currents_soma, time, synapse_site_stim_summary
    """

    all_files = os.listdir(dir_path)

    # Regex for files like: syns_1_trial_0_synaptic_currents_siz.csv
    match = re.compile(r'syns_(\d+)_trial_(\d+)_(.*)\.csv')

    trial_dict = {}

    for fname in all_files:
        m = match.match(fname)
        if not m:
            continue
        syn_count, trial_num, kind = m.groups()
        key = (int(syn_count), int(trial_num))
        trial_dict.setdefault(key, {})[kind] = fname

    for key in sorted(trial_dict.keys()):
        syn_count, trial_num = key
        files = trial_dict[key]

        required = [
            'synapse_site',
            'synaptic_currents_siz',
            'synaptic_currents_soma',
            'time',
            'synapse_site_stim_summary'
        ]
        if not all(r in files for r in required):
            print(f"⚠️ Skipping trial {key} — missing required files.")
            continue

        try:
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']), header=None).values
            vSoma = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df = pd.read_csv(os.path.join(dir_path, files['time']), header=None)
            stim_summary_df = pd.read_csv(os.path.join(dir_path, files['synapse_site_stim_summary']))

            yield (time_df, vSIZ, vSoma, synSite_df, stim_summary_df)

        except Exception as e:
            print(f"❌ Error loading trial {key}: {e}")
            continue

def iterQuintsNeighborhood(dir_path):
    """
    Handles filenames from COM_neighborhood_activations output:
        trial{02d}_seed{id}_{kind}.csv

    Yields: (time_df, vSIZ_np, vSoma_np, synSite_df, k, trial_num, seed_id)
    The k value and seed_id are parsed from dir_path and filename respectively.
    """
    all_files = os.listdir(dir_path)

    # Parse k from the dir_path itself: .../{VPN}_k{k:02d}
    k_match = re.search(r'_k(\d+)$', dir_path.rstrip('/'))
    k_val = int(k_match.group(1)) if k_match else -1

    pattern = re.compile(
        r'trial(\d+)_seed(\d+)_(synapse_site|synaptic_currents_siz|synaptic_currents_soma|time)\.csv'
    )

    group_dict = {}
    for fname in all_files:
        m = pattern.match(fname)
        if not m:
            continue
        trial_num, seed_id, kind = m.groups()
        key = (int(trial_num), seed_id)
        group_dict.setdefault(key, {})[kind] = fname

    required = ['synapse_site', 'synaptic_currents_siz', 'synaptic_currents_soma', 'time']

    for key in sorted(group_dict.keys()):
        trial_num, seed_id = key
        files = group_dict[key]

        if not all(r in files for r in required):
            print(f"⚠️ Skipping trial{trial_num:02d}_seed{seed_id} — "
                  f"missing: {set(required) - set(files.keys())}")
            continue

        try:
            print(f"Loading: k={k_val} trial={trial_num:02d} seed={seed_id}", flush=True)
            synSite_df = pd.read_csv(os.path.join(dir_path, files['synapse_site']))
            vSIZ      = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_siz']),  header=None).values
            vSoma     = pd.read_csv(os.path.join(dir_path, files['synaptic_currents_soma']), header=None).values
            time_df   = pd.read_csv(os.path.join(dir_path, files['time']), header=None)

            yield (time_df, vSIZ, vSoma, synSite_df, k_val, trial_num, seed_id)

        except Exception as e:
            print(f"❌ Error loading trial{trial_num:02d}_seed{seed_id}: {e}")
            continue

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

#########
#Plotting of single VPN activated synapses 

def plotallsynsexp2syn(neuron_name):
    if neuron_name == "DNp01":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synapse_sites.csv', dtype={'pre': str})
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_soma.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_time.csv', header=None)

        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_synpt_df_exp2, 'k', FIGURE=None, gridPos=[3, 0], axList=None)
        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_siz_df_exp2, 'b', FIGURE=TF, gridPos=[3,1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_soma_df_exp2, 'g', FIGURE=TF, gridPos=[3,2], axList=TF_AXLIST)

    elif neuron_name == "DNp03":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synapse_sites.csv', dtype={'pre': str})
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_soma.csv', header=None)
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_time.csv', header=None)

        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_synpt_df_exp2, 'k', FIGURE=None, gridPos=[3, 0], axList=None)
        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_siz_df_exp2, 'b', FIGURE=TF, gridPos=[3,1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_soma_df_exp2, 'g', FIGURE=TF, gridPos=[3,2], axList=TF_AXLIST)
        
        TF_AXLIST[0].set_ylim([-61, -58.75])
        TF_AXLIST[1].set_ylim([-61, -60.5])
        TF_AXLIST[2].set_ylim([-59.3, -59])

    
    TF_AXLIST[1].set_xlim([47.5, 60])
    TF_AXLIST[0].set_xlim([47.5, 60])
    TF_AXLIST[2].set_xlim([47.5, 60])

    TF_AXLIST[0].set_title('Recording at Dendrite')
    TF_AXLIST[1].set_title('Recording at SIZ')
    TF_AXLIST[2].set_title('Recording at SOMA')


    TF_AXLIST[0].set_xlabel('Time (ms)')
    TF_AXLIST[1].set_xlabel('Time (ms)')
    TF_AXLIST[0].set_ylabel('Voltage (mV)')
    TF_AXLIST[1].set_ylabel('Voltage (mV)')
    TF_AXLIST[2].set_ylabel('Voltage (mV)')

    TF.tight_layout()
    for ax in TF_AXLIST:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.show()

def plotallsynsexp2syn_by_type(neuron_name):
    if neuron_name == "DNp01":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synapse_sites.csv', dtype={'pre': str})
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_soma.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_time.csv', header=None)
        
        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_synpt_df_exp2, 'g', FIGURE=None, gridPos=[3, 0], axList=None)
        TF, TF_AXLIST = plotTraces(time_df_exp2[0:519], synCurr_siz_df_exp2[0:519], 'b', FIGURE=TF, gridPos=[3,1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[520:], synCurr_siz_df_exp2[520:], 'y', FIGURE=TF, gridPos=[3,1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[0:519], synCurr_soma_df_exp2[0:519], 'b', FIGURE=TF, gridPos=[3, 2], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[520:], synCurr_soma_df_exp2[520:], 'y', FIGURE=TF, gridPos=[3, 2], axList=TF_AXLIST)
        
        # TF_AXLIST[1].set_ylim([erev-0.001, SIZMAX+0.02])
        TF_AXLIST[2].set_ylim([-66.42, -66.36])

    elif neuron_name == "DNp03":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synapse_sites.csv', dtype={'pre': str})
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_soma.csv', header=None)
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_time.csv', header=None)

        TF, TF_AXLIST = plotTraces(time_df_exp2, synCurr_synpt_df_exp2, 'k', FIGURE=None, gridPos=[3, 0], axList=None)

        TF, TF_AXLIST = plotTraces(time_df_exp2[0:510], synCurr_siz_df_exp2[0:510], 'b', FIGURE=TF, gridPos=[3,1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[511:730], synCurr_siz_df_exp2[511:730], 'g', FIGURE=TF, gridPos=[3,1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[731:1825], synCurr_siz_df_exp2[731:1825], 'p', FIGURE=TF, gridPos=[3, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[1826:1842], synCurr_siz_df_exp2[1826:1842], 'y', FIGURE=TF, gridPos=[3, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[1843:], synCurr_siz_df_exp2[1843:], 'r', FIGURE=TF, gridPos=[3, 1], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df_exp2[0:510], synCurr_soma_df_exp2[0:510], 'b', FIGURE=TF, gridPos=[3,2], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[511:730], synCurr_soma_df_exp2[511:730], 'g', FIGURE=TF, gridPos=[3,2], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[731:1825], synCurr_soma_df_exp2[731:1825], 'p', FIGURE=TF, gridPos=[3, 2], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[1826:1842], synCurr_soma_df_exp2[1826:1842], 'y', FIGURE=TF, gridPos=[3, 2], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df_exp2[1843:], synCurr_soma_df_exp2[1843:], 'r', FIGURE=TF, gridPos=[3, 2], axList=TF_AXLIST)
        
        TF_AXLIST[0].set_ylim([-61, -58.75])
        TF_AXLIST[1].set_ylim([-61, -60.5])
        TF_AXLIST[2].set_ylim([-59.3, -59])
        

    TF_AXLIST[0].set_xlim([47.5, 60])
    TF_AXLIST[1].set_xlim([47.5, 60])
    TF_AXLIST[2].set_xlim([47.5, 60])

    TF_AXLIST[0].set_title('Recording at Dendrite')
    TF_AXLIST[2].set_title('Recording at SIZ')
    TF_AXLIST[1].set_title('Recording at Soma')

    TF_AXLIST[0].set_xlabel('Time (ms)')
    TF_AXLIST[1].set_xlabel('Time (ms)')
    TF_AXLIST[2].set_xlabel('Time (ms)')
    TF_AXLIST[0].set_ylabel('Voltage (mV)')
    TF_AXLIST[1].set_ylabel('Voltage (mV)')
    TF_AXLIST[2].set_ylabel('Voltage (mV)')
    TF.tight_layout()
    for ax in TF_AXLIST:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()

##########

# Recreating figures from paper

#############
#Retinotopy set up and analysis Figure 3 Supplemental 1, Figure 4D-E, Figure 4 Supplemental 1

def calculate_distance(neuron_name, VPN, synMap_df):
    synSecList = []
    synTypeList = []
    synpreList = []
    synSegRangeVarList = []

    if VPN not in ["LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4", "LC4_hemi", "LC22_hemi","LPLC1_hemi", "LPLC2_hemi", "LPLC4_hemi"]:
        raise ValueError("Invalid VPN value")

    if VPN in ["LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4"]:
        mask = synMap_df['type'] == VPN
        synMap_df = synMap_df[mask]
    elif VPN == "LC4_hemi":
        mask = synMap_df['type'] == 'LC4'
        synMap_df = synMap_df[mask]
    elif VPN == "LPLC2_hemi":
        mask = synMap_df['type'] == 'LPLC2'
        synMap_df = synMap_df[mask]
    elif VPN == "LC22_hemi":
        mask = synMap_df['type'] == 'LC22'
        synMap_df = synMap_df[mask]
    elif VPN == "LPLC1_hemi":
        mask = synMap_df['type'] == 'LPLC1'
        synMap_df = synMap_df[mask]
    elif VPN == "LPLC4_hemi":
        mask = synMap_df['type'] == 'LPLC4'
        synMap_df = synMap_df[mask]

    for index, row in synMap_df.iterrows():
        TEMP = row.loc["mappedSection"]
        for sec in h.allsec():
            if sec.name() == TEMP:
                synSecList.append(sec)
                synTypeList.append(row.loc["type"])
                synpreList.append(row.loc["pre"])
                synSegRangeVarList.append(row.loc["mappedSegRangeVar"])

  # Calculate distances for all unique synapse pairs
    synapse_indices = range(len(synSecList))
    syn_pairs = combinations(synapse_indices, 2)  # Get all unique pairs of synapses

    synapse_id_1 = []
    synapse_id_2 = []
    distances = []

    for i, j in syn_pairs:
        dist = h.distance(synSecList[i](synSegRangeVarList[i]), synSecList[j](synSegRangeVarList[j]))
        distances.append(dist)
        synapse_id_1.append(i)
        synapse_id_2.append(j)

    # Create a DataFrame to hold the information
    dist_df = pd.DataFrame({
        'Synapse_ID_1': synapse_id_1,
        'Synapse_ID_2': synapse_id_2,
        'Section_1': [synSecList[i] for i in synapse_id_1],
        'Section_2': [synSecList[i] for i in synapse_id_2],
        'Pre_ID_1': [synpreList[i] for i in synapse_id_1],
        'Pre_ID_2': [synpreList[i] for i in synapse_id_2],
        'Distance': distances,
    })
    return dist_df

def calculate_average_distance_between_synapses_per_neuron_pair(neuron_name, VPN, synMap_df):
    if VPN not in ["LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4", "LC4_hemi", "LC22_hemi","LPLC1_hemi", "LPLC2_hemi", "LPLC4_hemi"]:
        raise ValueError("Invalid VPN value")

    if VPN in ["LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4"]:
        mask = synMap_df['type'] == VPN
        synMap_df = synMap_df[mask]
    elif VPN == "LC4_hemi":
        mask = synMap_df['type'] == 'LC4'
        synMap_df = synMap_df[mask]
    elif VPN == "LPLC2_hemi":
        mask = synMap_df['type'] == 'LPLC2'
        synMap_df = synMap_df[mask]
    elif VPN == "LC22_hemi":
        mask = synMap_df['type'] == 'LC22'
        synMap_df = synMap_df[mask]
    elif VPN == "LPLC1_hemi":
        mask = synMap_df['type'] == 'LPLC1'
        synMap_df = synMap_df[mask]
    elif VPN == "LPLC4_hemi":
        mask = synMap_df['type'] == 'LPLC4'
        synMap_df = synMap_df[mask]


    synSecList = []
    synpreList = []
    synSegRangeVarList = []

    for index, row in synMap_df.iterrows():
        TEMP = row.loc["mappedSection"]
        for sec in h.allsec():
            if sec.name() == TEMP:
                synSecList.append(sec)
                synpreList.append(row.loc["pre"])
                synSegRangeVarList.append(row.loc["mappedSegRangeVar"])

    # Group synapses by 'pre' id
    synapse_groups = {}
    for i, pre_id in enumerate(synpreList):
        if pre_id in synapse_groups:
            synapse_groups[pre_id].append(i)
        else:
            synapse_groups[pre_id] = [i]

    # Generate unique pairs of 'pre' ids
    pre_id_pairs = list(combinations(synapse_groups.keys(), 2))

    # Generate unique pairs of synapse indices for each 'pre' id pair
    synapse_pairs = []
    distances = []
    for pre_id_pair in pre_id_pairs:
        for i in synapse_groups[pre_id_pair[0]]:
            for j in synapse_groups[pre_id_pair[1]]:
                synapse_pairs.append((i, j))
                dist = h.distance(synSecList[i](synSegRangeVarList[i]), synSecList[j](synSegRangeVarList[j]))
                distances.append(dist)

    # Compute the average distance for each unique neuron pair
    avg_distance_per_neuron_pair = {}
    for pre_id_pair in pre_id_pairs:
        dist_list = [distances[i] for i, (syn_i, syn_j) in enumerate(synapse_pairs) if pre_id_pair[0] in synpreList[syn_i] and pre_id_pair[1] in synpreList[syn_j]]
        avg_distance = sum(dist_list) / len(dist_list)
        avg_distance_per_neuron_pair[pre_id_pair] = avg_distance

    # Convert the results to a DataFrame
    neuron_syn_pair_dists = pd.DataFrame(list(avg_distance_per_neuron_pair.keys()), columns=['Pre_ID_1', 'Pre_ID_2'])
    neuron_syn_pair_dists['Average_Distance'] = list(avg_distance_per_neuron_pair.values())
    return neuron_syn_pair_dists

def nearest_neighbor(neuron_name, VPN, sizSection, synMap_df):
    synSecList = []
    synTypeList = []
    synpreList = []
    synSegRangeVarList = []

    
    if VPN == "LC4" or VPN == "LC4_hemi":
        mask = synMap_df['type'] == 'LC4'
    elif VPN == "LC6" or VPN == "LC6_hemi":
        mask = synMap_df['type'] == 'LC6'
    elif VPN == "LC22" or VPN == "LC22_hemi":
        mask = synMap_df['type'] == 'LC22'
    elif VPN == "LPLC1" or VPN == "LPLC1_hemi":
        mask = synMap_df['type'] == 'LPLC1'
    elif VPN == "LPLC2" or VPN == "LPLC2_hemi":
        mask = synMap_df['type'] == 'LPLC2'
    elif VPN == "LPLC4" or VPN == "LPLC4_hemi":
        mask = synMap_df['type'] == 'LPLC4'
    else:
        raise ValueError("Invalid VPN value")
        

    synMap_df = synMap_df[mask]
    for index, row in synMap_df.iterrows():
        synMap_df = synMap_df.reset_index(drop=True) 
        TEMP = row.loc["mappedSection"]
        for sec in h.allsec():
            if sec.name() == TEMP:
                synSecList.append(sec)
                synTypeList.append(row.loc["type"])
                synpreList.append(row.loc["pre"])
                synSegRangeVarList.append(row.loc["mappedSegRangeVar"])
    synDistList = []
    synDistSIZ = []
    synDistSIZ2 = []
    synpreIdList = [] # Add a list to store pre_id with minimum distance
    synpreIdList2 = [] # Add a list to store the pre_id of the nearest neighbor
    secpreIdlist = [] # Add a list to store the section 
    secpreIdlist2 = [] # Add a list to store the section of nearest neighbor

    t0 = clock.time()  # Start measuring time
    for i, sec in enumerate(synSecList):
        # Calculate distance to other synapses in synSecList
        distances = []
        dist_SIZ_1 = []
        dist_SIZ_2 = []
        neuron_id = []
        neuron_id_2 = []
        sec_neuron = []
        sec_neuron_2 = []
        for j, sec2 in enumerate(synSecList):
            if j != i:  # Skip the current synapse itself
                dist = h.distance(sec(synSegRangeVarList[i]), sec2(synSegRangeVarList[j]))
                syn_to_SIZ_dist = h.distance(sizSection(0.05), sec(synSegRangeVarList[i]))
                syn_to_SIZ_dist2 = h.distance(sizSection(0.05), sec(synSegRangeVarList[j]))
                dist_SIZ_1.append(syn_to_SIZ_dist)
                dist_SIZ_2.append(syn_to_SIZ_dist2)
                distances.append(dist)
                neuron_id.append(synpreList[i])  # Store pre_id for i-th synapse
                neuron_id_2.append(synpreList[j])  # Store pre_id for j-th synapse
                sec_neuron.append(synSecList[i])  # Store pre_id for i-th synapse
                sec_neuron_2.append(synSecList[j])  # Store pre_id for j-th synapse
                
        min_dist = np.min(distances) if distances else np.inf  # Get the minimum distance, or set to infinity if no distances are calculated
        synDistList.append(min_dist)
        min_pre_id = neuron_id[np.argmin(distances)] if distances else None  # Get the pre_id of the current synapse with minimum distance
        min_pre_id_2 = neuron_id_2[np.argmin(distances)] if distances else None  # Get the pre_id of nearest synapse
        synpreIdList.append(min_pre_id)  # Store pre_id for i-th synapse with minimum distance
        synpreIdList2.append(min_pre_id_2)  # Store pre_id for j-th synapse with minimum distance
        min_pre_sec = sec_neuron[np.argmin(distances)] if distances else None  # Get the pre_id of the current synapse with minimum distance
        min_pre_sec_2 = sec_neuron_2[np.argmin(distances)] if distances else None  # Get the pre_id of nearest synapse

        min_dist_SIZ = dist_SIZ_1[np.argmin(distances)] if distances else None  # Get the pre_id of the current synapse with minimum distance
        min_dist_SIZ_2 = dist_SIZ_2[np.argmin(distances)] if distances else None  # Get the pre_id of nearest synapse

        synDistSIZ.append(min_dist_SIZ)
        synDistSIZ2.append(min_dist_SIZ_2)
        secpreIdlist.append(min_pre_sec)
        secpreIdlist2.append(min_pre_sec_2)
        
    print("Total time taken:", clock.time() - t0, "seconds wall time")  # Print total time taken for the function to finish
    
    nearest_neighbor_df = pd.DataFrame({'pre_id' : synpreIdList,
                                'nn_pre_id' : synpreIdList2,
                                'distance_um' : synDistList,
                                'pre_sec' : secpreIdlist,
                                'nn_sec' : secpreIdlist,
                                'pre_dist_SIZ' : synDistSIZ,
                                'nn_dist_SIZ' : synDistSIZ2},     
                                columns=['pre_id','nn_pre_id', 'distance_um','pre_sec','nn_sec', 'pre_dist_SIZ','nn_dist_SIZ'])
    nearest_neighbor_df['same_pre'] = nearest_neighbor_df['pre_id'] == nearest_neighbor_df['nn_pre_id']
    nearest_neighbor_df['color'] = np.where(nearest_neighbor_df['pre_id'] == nearest_neighbor_df['nn_pre_id'], 'Blue', 'Red')
    nearest_neighbor_df['SIZ_dist_diff'] = abs(nearest_neighbor_df['pre_dist_SIZ'] - nearest_neighbor_df['nn_dist_SIZ'])
    
    avg_distance_per_neuron = nearest_neighbor_df.groupby('pre_id')['distance_um'].mean()
    avg_distance_ids = np.unique(nearest_neighbor_df['pre_id'])
    avg_distance_per_neuron_df = pd.DataFrame({'pre_id' : avg_distance_ids,
                                'avg_dist_nn_um' : avg_distance_per_neuron},
                                columns=['pre_id', 'avg_dist_nn_um'])
    
    return nearest_neighbor_df, avg_distance_per_neuron_df # Return dataframe with all information of distance and neighbor pair ids

def NN_COM_merge(VPN, nearest_neighbor_df, dist_df, neuron_pair_dist):
    nn_df = nearest_neighbor_df
    all_dist_df = dist_df
    neuron_syn_pair_dists = neuron_pair_dist

    if VPN == "LC4":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC4_hemi":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_hemibrain_centroids.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC22_hemi":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_hemibrain_centroids.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC6":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC6_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC6_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC22":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC1":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC2":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC2_hemi":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_hemibrain_centroids.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC1_hemi":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_hemibrain_centroids.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC4_hemi":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_hemibrain_centroids.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC4":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_receptive_field.xlsx", dtype={'updated_ids': str})
    else:
        raise ValueError("Invalid VPN value")
    # load COM data frame with 'updated_ids' as strs
    nn_and_COM_col_coded = pd.merge(VPN_COM, nn_df, left_on="updated_ids", right_on="pre_id")
    nn_and_COM_col_coded = pd.merge(VPN_COM, nn_and_COM_col_coded, left_on="updated_ids", right_on="nn_pre_id")
    nn_and_COM_col_coded['col'] = np.where(nn_and_COM_col_coded['same_pre'], 'blue', 'red')

    distances = []

    #loop through the rows of the dataframe and calculate the distance between the X and Y vectors
    for index, row in  nn_and_COM_col_coded.iterrows():
        x_vec = np.array([row['DV_raw_um_x'], row['AP_raw_um_x']])
        y_vec = np.array([row['DV_raw_um_y'], row['AP_raw_um_y']])
        distance = np.linalg.norm(x_vec - y_vec)
        distances.append(distance)

    # add the distances to a new column in the dataframe
    nn_and_COM_col_coded['Distance_COM'] = distances

    #
    all_syns_dist_df = pd.merge(VPN_COM, all_dist_df, left_on="updated_ids", right_on="Pre_ID_1")
    all_syns_dist_df = pd.merge(VPN_COM, all_syns_dist_df, left_on="updated_ids", right_on="Pre_ID_2")
    

    distances_2 = []
    #loop through the rows of the dataframe and calculate the distance between the X and Y vectors
    for index, row in  all_syns_dist_df.iterrows():
        x_vec = np.array([row['DV_raw_um_x'], row['AP_raw_um_x']])
        y_vec = np.array([row['DV_raw_um_y'], row['AP_raw_um_y']])
        distance_2 = np.linalg.norm(x_vec - y_vec)
        distances_2.append(distance_2)

    all_syns_dist_df['Distance_COM'] = distances_2
    
    neuron_syn_pair_dists = pd.merge(VPN_COM, neuron_syn_pair_dists, left_on="updated_ids", right_on="Pre_ID_1")
    neuron_syn_pair_dists = pd.merge(VPN_COM, neuron_syn_pair_dists, left_on="updated_ids", right_on="Pre_ID_2")

    distances_3 = []
    #loop through the rows of the dataframe and calculate the distance between the X and Y vectors
    for index, row in  neuron_syn_pair_dists.iterrows():
        x_vec = np.array([row['DV_raw_um_x'], row['AP_raw_um_x']])
        y_vec = np.array([row['DV_raw_um_y'], row['AP_raw_um_y']])
        distance_3 = np.linalg.norm(x_vec - y_vec)
        distances_3.append(distance_3)

    neuron_syn_pair_dists['Distance_COM'] = distances_3

    return nn_and_COM_col_coded, all_syns_dist_df, neuron_syn_pair_dists

def NN_stats(original_data=None, shuffled_data = None):
    # Your single number
    single_number = original_data

    # Your group of numbers
    group_of_numbers = shuffled_data
    shapiro_stat, shapiro_p_value = shapiro(shuffled_data)

    # Print the results
    print(f'Shapiro-Wilk statistic: {shapiro_stat}')
    print(f'P-value: {shapiro_p_value}')

    # Calculate the mean of the group
    group_mean = np.mean(group_of_numbers)

    # Perform a one-sample t-test
    t_statistic, p_value = ttest_1samp(group_of_numbers, single_number, alternative='less')

    # Print the results
    print(f'Single number: {single_number}')
    print(f'Mean of the group: {group_mean}')
    print(f'T-statistic: {t_statistic}')
    print(f'P-value: {p_value}')

    # Check for statistical significance (common alpha level is 0.05)
    if p_value < 0.05:
        print('The difference is statistically significant.')
    else:
        print('The difference is not statistically significant.')

def bar_whisker_plot_shuffled_vs_original_data(dataset= None):
    #Data  shown here was retrieved from the NN analysis above which generates csv's with the following values. 
    if dataset == "FAFBv630":
        data1 = [44.56, 36.88, 51.02, 49.8, 43.75, 57.27, 40, 58.52, 54.13, 46.93, 65.66, 44.72, 53.85, 45.02]
        data2 = [2.154, 2.455, 3.03, 3.03, 1.998, 6.385, 2.008, 1.914, 2.174, 1.213, 3.03, 1.352, 7.915, 2.005]
    elif dataset == "FAFBv783":
        data1 = [17.38, 16.81, 18.10, 19.13, 9.09, 38.71, 14.39, 30.78, 24.72, 15.49, 29.31, 18.76, 33.33, 16.70]
        data2 = [2.294, 1.421, 2.43, 2.323, 1.818, 6.647, 1.912, 1.809, 2.279, 1.233, 3.047, 1.267, 5.011, 2.023]
    elif dataset == "FAFBv783_updated":
        data1 = [14.92, 9.92, 13.63, 18.40, 28.57, 29.47, 12.37, 23.03, 18.45, 10.68, 34.82, 17.62, 20.45, 11.32]
        data2 = [1.892, 1.244, 1.988, 1.911, 0.793, 5.725, 2.019, 1.919, 1.967, 1.076, 2.321, 0.912, 3.636, 1.688]
    elif dataset == "hemibrain":
        data1 = [13.97, 9.47]
        data2 = [1.543,  1.372]

    statistic1, p_value1 = shapiro(data1)
    statistic2, p_value2 = shapiro(data2)

    # Print the results
    print("Shapiro-Wilk Test for data1:")
    print(f"Statistic: {statistic1}, p-value: {p_value1}")
    if p_value1 > 0.05:
        print("The data1 follows a normal distribution.")
    else:
        print("The data1 does not follow a normal distribution.")

    print("\nShapiro-Wilk Test for data2:")
    print(f"Statistic: {statistic2}, p-value: {p_value2}")
    if p_value2 > 0.05:
        print("The data2 follows a normal distribution.")
    else:
        print("The data2 does not follow a normal distribution.")

    statistic, p_value_mw = stats.mannwhitneyu(data1, data2)

    # Create box and whisker plot with connecting lines
    fig, ax = plt.subplots()
    bp = ax.boxplot([data1, data2], vert=True, labels=['Original Data', 'Randomized'], medianprops={'linewidth': 2})
    ax.scatter(np.ones_like(data1), data1, color='blue', label='Original Data')
    ax.scatter(2 * np.ones_like(data2), data2, color='black', label='Randomized')

    # Draw connecting lines
    for i in range(len(data1)):
        ax.plot([1, 2], [data1[i], data2[i]], color='gray', linestyle='--', alpha=0.25)

    # Calculate and print standard deviation
    std_col1 = np.std(data1)
    std_col2 = np.std(data2)
    print(f"Standard Deviation - Original Data: {std_col1:.2f}, Randomized Data: {std_col2:.2f}")

    # Calculate and print statistics
    stats_col1 = stats.describe(data1)
    stats_col2 = stats.describe(data2)
    print("\nStatistics - Original Data:")
    print(stats_col1)
    print("\nStatistics - Randomized Data:")
    print(stats_col2)

    # Perform Mann-Whitney U test and print results
    print(f"\nMann-Whitney U test between Column 1 and Column 2:")
    print(f"U-statistic: {statistic:.4f}")
    print(f"P-value: {p_value_mw}")

    # Display the plot
    ax.set_ylabel('Percentage of NN')
    plt.legend()
    plt.show()

def plot_retinotopy(neuron_syn_pair_dists, VPN, neuron_name):
    fig, ax1 = plt.subplots(figsize=(12, 12))

    # Scatter plot
    ax1.scatter(neuron_syn_pair_dists['Average_Distance'], neuron_syn_pair_dists['Distance_COM'], s=10, c='black')
    # Add the regression line
    m, b = np.polyfit(neuron_syn_pair_dists['Average_Distance'], neuron_syn_pair_dists['Distance_COM'], 1)
    plt.plot(neuron_syn_pair_dists['Average_Distance'], m*neuron_syn_pair_dists['Average_Distance'] + b)

    # Set labels and title
    ax1.set_xlabel(f'Average Synapse Distance for all {VPN} pairs (um)')
    ax1.set_ylabel(f'Distance between COMs for all {VPN} pairs (um)')
    ax1.set_title(f'{VPN} Dendrite Retinotopy for all {VPN} pairs to {neuron_name}')

    x = neuron_syn_pair_dists['Average_Distance']
    y = neuron_syn_pair_dists['Distance_COM']

    # Fit a linear regression line
    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    # Calculate the predicted values
    predicted_values = slope * x + intercept
    plt.plot(x, predicted_values, color='red', label='Linear Regression Line')

    # Print R^2 value
    r_squared = r_value ** 2
    print("R^2 value:", r_squared)
    print("p value:", p_value)
    # Show the plot
    plt.show()

def COM_receptive_field_split_and_visualization(VPN):
    if VPN == "LC4":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC6":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC6_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC6_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC22":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC1":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC2":
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_receptive_field.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC4":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
        # VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_receptive_field.xlsx", dtype={'updated_ids': str})
    else:
        raise ValueError("Invalid VPN value")
    ventral = VPN_COM[VPN_COM['DV_norm'] > 0.5]['DV_norm'].count()
    posterior = VPN_COM[VPN_COM['AP_norm'] < 0.5]['AP_norm'].count()

    dorsal = VPN_COM[VPN_COM['DV_norm'] < 0.5]['DV_norm'].count()
    anterior = VPN_COM[VPN_COM['AP_norm'] > 0.5]['AP_norm'].count()

    print(f"'Ventral': {ventral}, 'Dorsal': {dorsal}")
    print(f"'Anterior': {anterior}, 'Posterior': {posterior}")
    data = {
    'VPN': ['LC4_DV', 'LC4_AP', 'LC6_DV', 'LC6_AP', 'LC22_DV', 'LC22_AP',
            'LPLC1_DV', 'LPLC1_AP', 'LPLC2_DV', 'LPLC2_AP', 'LPLC4_DV', 'LPLC4_AP'],
    'Dorsal': [34, np.nan, 37, np.nan, 21, np.nan, 33, np.nan, 60, np.nan, 31, np.nan],
    'Ventral': [20, np.nan, 28, np.nan, 22, np.nan, 35, np.nan, 48, np.nan, 25, np.nan],
    'Anterior': [np.nan, 28, np.nan, 38, np.nan, 25, np.nan, 27, np.nan, 57, np.nan, 34],
    'Posterior': [np.nan, 26, np.nan, 27, np.nan, 18, np.nan, 41, np.nan, 51, np.nan, 22]
    }
    df = pd.DataFrame(data)

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # Bar width
    bar_width = 0.5

    # Bar positions
    positions = np.arange(len(df['VPN']))

    # Plotting dorsal/ventral subplot
    ax1.bar(positions - bar_width/2, df['Dorsal'], width=bar_width, label='Dorsal', color = 'red')
    ax1.bar(positions + bar_width/2, df['Ventral'], width=bar_width, label='Ventral', color='blue')

    ax1.set_xticks(positions)
    ax1.set_xticklabels(df['VPN'], rotation=45, ha='right')
    ax1.set_xlabel('VPN')
    ax1.set_ylabel('Neuron Count')
    ax1.set_title('Dorsal and Ventral Relationships')
    ax1.legend()

    # Plotting anterior/posterior subplot
    ax2.bar(positions - bar_width/2, df['Anterior'], width=bar_width, label='Anterior', color = 'purple')
    ax2.bar(positions + bar_width/2, df['Posterior'], width=bar_width, label='Posterior', color = 'cyan')

    ax2.set_xticks(positions)
    ax2.set_xticklabels(df['VPN'], rotation=45, ha='right')
    ax2.set_xlabel('VPN')
    ax2.set_title('Anterior and Posterior Relationships')
    ax2.legend()

    # Adjust layout for better visualization
    plt.tight_layout()

    # Show the plot
    plt.show()

def nearest_neighbor_analysis(neuron_name, synMap_df):
    same_pre_nn = synMap_df[synMap_df['pre_id'] == synMap_df['nn_pre_id']]
    diff_pre_nn = synMap_df[synMap_df['pre_id'] != synMap_df['nn_pre_id']]

    # calculate the percentage of synapses with the same pre id and different pre id
    same_pre_percent = len(same_pre_nn) / len(synMap_df) * 100
    diff_pre_percent = len(diff_pre_nn) / len(synMap_df) * 100

    print(f"Percentage of synapses with the same pre id as their nearest neighbor: {same_pre_percent:.2f}%")
    print(f"Percentage of synapses with a different pre id from their nearest neighbor: {diff_pre_percent:.2f}%")
    #############

def plot_vpn_dn_comparison_bargraph_with_stats(original_data, combined_data, labels=None, show=True):
    """
    Plots a grouped bar graph comparing original values to shuffled mean values,
    adds error bars (std dev), runs statistical tests, and annotates p-values and significance.
    """
    if len(original_data) != len(combined_data):
        raise ValueError("original_data and combined_data must have the same length.")

    if labels is None:
        labels = [f"VPN-DN {i+1}" for i in range(len(original_data))]

    means = [np.mean(data) for data in combined_data]
    std_devs = [np.std(data) for data in combined_data]
    p_values = []
    significance_marks = []

    # Statistical tests
    for orig, shuffled in zip(original_data, combined_data):
        _, shapiro_p = shapiro(shuffled)
        _, p = ttest_1samp(shuffled, orig, alternative='less')
        p_values.append(p)

        if p < 0.0001:
            significance_marks.append("****")
        elif p < 0.001:
            significance_marks.append("***")
        elif p < 0.01:
            significance_marks.append("**")
        elif p < 0.05:
            significance_marks.append("*")
        else:
            significance_marks.append("")

    # Plotting
    x = np.arange(len(original_data))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(10, len(labels)*0.6), 6))
    ax.bar(x - width/2, original_data, width, label='Original', color='blue')
    ax.bar(x + width/2, means, width, yerr=std_devs, capsize=5, label='Shuffled Mean ± SD', color='black')

    ax.set_ylabel('NN belonging to the same neuron (%)')
    ax.set_title('Original vs Shuffled Mean per VPN-DN')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', linestyle='--', alpha=0.5)

    # Annotate p-values and asterisks
    for i, (p, mark, orig, mean, std) in enumerate(zip(p_values, significance_marks, original_data, means, std_devs)):
        max_height = max(orig, mean + std)
        ax.text(i, max_height + 1, f"p={p:.3g}\n{mark}", ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.close()

def VPN_per_syn_count(dir_path, neuron_name):
    # Load the CSV
    file_path = f'{dir_path}/{neuron_name}_Partner_syn_spread_vs_dist_to_siz_data.csv'
    df = pd.read_csv(file_path)
    
    # Count number of neurons per synapse count
    synapse_counts = df['Synct'].value_counts().sort_index()

    # Compute average synapse count
    avg_synct = df['Synct'].mean()

    # Plot
    plt.figure(figsize=(8, 5))
    bars = plt.bar(synapse_counts.index, synapse_counts.values, color='skyblue', edgecolor='black')

    # Add count labels above bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, height + 0.5, str(int(height)), 
                 ha='center', va='bottom', fontsize=10)

    # Add red line for average synapse count
    plt.axvline(avg_synct, color='red', linestyle='--', linewidth=2, label=f'Avg Synct = {avg_synct:.2f}')

    # Labels and title
    plt.xlabel('Number of Synapses (Synct)', fontsize=12)
    plt.ylabel('Number of Neurons', fontsize=12)
    plt.title(f'Neuron Count per Synapse Number ({neuron_name})', fontsize=14)
    plt.xticks(synapse_counts.index)
    plt.legend()

    plt.show()


#Figure 4D,E
def syn_spread_attenuation(neuron_name, VPN=None):
    """
    Compute electrotonic distance between synapse pairs using pre-computed
    peaks from summary CSVs, and morphological synapse spread from synapse site CSVs.
    Works with flat folder structure (CSV pairs directly inside dir_path).
    Prints progress as neuron pairs are processed.
    Fits linear regression lines on plots and shows equation + R^2 with full precision.
    """

    # --- Locate simulation directory ---
    dir_path = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_pairwise_{VPN}"
    if not os.path.exists(dir_path):
        raise FileNotFoundError(f"Directory not found: {dir_path}")

    # Collect CSVs
    site_files = [f for f in os.listdir(dir_path) if f.endswith("_synpair_synapse_sites.csv")]
    summary_files = [f for f in os.listdir(dir_path) if f.endswith("_synpair_summary.csv")]

    # Extract common keys (neuron ID pairs) before "_synpair_"
    site_keys = {f.split("_synpair_")[0] for f in site_files}
    summary_keys = {f.split("_synpair_")[0] for f in summary_files}
    common_keys = sorted(site_keys & summary_keys)
    total_pairs = len(common_keys)

    print(f"Found {total_pairs} neuron pairs in {dir_path}")

    electrotonic_distances_per_VPN_pair = []
    dist_to_syn_per_VPN_pair = []
    neuron_pairs = []

    # --- Loop with progress tracking ---
    for idx, key in enumerate(common_keys, start=1):
        print(f"[{idx}/{total_pairs}] Processing neuron pair {key}...")

        site_file = os.path.join(dir_path, f"{key}_synpair_synapse_sites.csv")
        summary_file = os.path.join(dir_path, f"{key}_synpair_summary.csv")

        site_df = pd.read_csv(site_file)
        summary_df = pd.read_csv(summary_file)

        # Identify neuron pair
        unique_neurons = site_df['pre'].unique()
        if len(unique_neurons) != 2:
            print(f"  Skipping {key}: does not contain exactly 2 unique neurons")
            continue
        neuron_pairs.append((unique_neurons[0], unique_neurons[1]))

        # --- Electrotonic distance ---
        pair_e_distances = []
        for _, row in summary_df.iterrows():
            peak_stim = abs(row['peak_stim'] - row['resting_vm'])
            peak_rec = abs(row['peak_record'] - row['resting_vm'])
            if peak_stim > 0 and peak_rec > 0:
                log_ratio = np.log(peak_stim / peak_rec)
                pair_e_distances.append(log_ratio)

        avg_e_dist = np.nanmean(pair_e_distances) if pair_e_distances else np.nan
        electrotonic_distances_per_VPN_pair.append(avg_e_dist)

        # --- Synapse spread ---
        dist_to_syn = []
        for i in range(0, len(site_df), 2):
            if i + 1 >= len(site_df):
                continue
            row_i, row_j = site_df.iloc[i], site_df.iloc[i+1]

            synsec_1 = str2sec(row_i["mappedSection"])
            synsec_2 = str2sec(row_j["mappedSection"])
            segVar1, segVar2 = row_i["mappedSegRangeVar"], row_j["mappedSegRangeVar"]

            try:
                d = h.distance(synsec_1(segVar1), synsec_2(segVar2))
                dist_to_syn.append(d)
            except Exception:
                continue

        avg_dist_per_pair = np.mean(dist_to_syn) if dist_to_syn else np.nan
        dist_to_syn_per_VPN_pair.append(avg_dist_per_pair)

    # --- COM receptive field ---
    vpn_file_map = {
        "LC4": "LC4_receptive_field.xlsx",
        "LC6": "LC6_receptive_field.xlsx",
        "LC22": "LC22_receptive_field.xlsx",
        "LPLC1": "LPLC1_receptive_field.xlsx",
        "LPLC2": "LPLC2_receptive_field.xlsx",
        "LPLC4": "LPLC4_receptivefield_COM.xlsx"
    }
    if VPN not in vpn_file_map:
        raise ValueError("Invalid VPN value")

    VPN_COM = pd.read_excel(
        f"datafiles/Receptive_field_files/{vpn_file_map[VPN]}",
        dtype={'updated_ids': str}
    )

    distances = []
    for pre_id, post_id in neuron_pairs:
        pre_row = VPN_COM[VPN_COM['updated_ids'] == str(pre_id)]
        post_row = VPN_COM[VPN_COM['updated_ids'] == str(post_id)]
        if not pre_row.empty and not post_row.empty:
            x_vec = np.array([pre_row['DV_raw_um'].values[0], pre_row['AP_raw_um'].values[0]])
            y_vec = np.array([post_row['DV_raw_um'].values[0], post_row['AP_raw_um'].values[0]])
            distances.append(np.linalg.norm(x_vec - y_vec))
        else:
            distances.append(np.nan)

    # --- Helper function for regression plotting ---
    def plot_with_regression(ax, x, y, xlabel, ylabel, title):
        x = np.array(x)
        y = np.array(y)
        mask = ~np.isnan(x) & ~np.isnan(y)
        x, y = x[mask], y[mask]

        ax.plot(x, y, 'o', color='black')

        if len(np.unique(x)) > 1 and len(np.unique(y)) > 1:
            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            ax.plot(x, slope * x + intercept, 'r-')
            eqn_text = f"y={slope}x+{intercept}\n$R^2$={r_value**2}"
        else:
            eqn_text = "Regression not possible"

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.text(0.05, 0.95, eqn_text,
                transform=ax.transAxes, verticalalignment='top', color='red')

    # --- Plotting ---
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    plot_with_regression(
        axs[0],
        dist_to_syn_per_VPN_pair,
        electrotonic_distances_per_VPN_pair,
        "Average synapse spread (µm)",
        "Electrotonic distance",
        f"Electrotonic distance vs. Synapse spread {neuron_name}-{VPN}"
    )

    plot_with_regression(
        axs[1],
        electrotonic_distances_per_VPN_pair,
        distances,
        "Electrotonic distance",
        "Euclidean distance between COM points",
        f"Electrotonic distance vs. Visual space {neuron_name}-{VPN}"
    )

    fig.tight_layout()

    # Save the combined two-panel figure
    fig.savefig(
        os.path.join(
            dir_path,
            f"{neuron_name}_{VPN}_electrotonic_distance_analysis.svg"
        ),
        format="svg",
        bbox_inches="tight"
    )


    # --- Third figure ---
    fig3, ax3 = plt.subplots(figsize=(8, 6))

    plot_with_regression(
        ax3,
        dist_to_syn_per_VPN_pair,
        distances,
        "Average synapse spread (µm)",
        "Euclidean distance between COM points",
        f"Synapse spread vs. Visual space {neuron_name}-{VPN}"
    )

    fig3.tight_layout()

    # Save the third figure
    fig3.savefig(
        os.path.join(
            dir_path,
            f"{neuron_name}_{VPN}_synapse_spread_vs_visual_space.svg"
        ),
        format="svg",
        bbox_inches="tight"
    )

    plt.show()

#This looks as synapses of a given type and calculates the surface area they encompass on DN Dendrites
def calculate_branch_surface_area_and_export_swc(synSite_df, dendList, target_type, swc_filename=None, gap_allowance=0):
    """
    Calculates the surface area of dendritic branches containing synapses from a given type,
    with the option to include intervening sections up to a specified gap length,
    and optionally exports those branches to an SWC file.

    Args:
        synSite_df (pd.DataFrame): Synapse information with 'mappedSection' and 'type'.
        dendList (list): List of dendritic section objects.
        target_type (str): Presynaptic type to filter for (e.g. "LPLC2").
        swc_filename (str, optional): Output filename for the SWC export. If None, no file is written.
        gap_allowance (int, optional): Maximum number of intervening sections allowed to keep.
                                       Default = 0 (only sections with synapses).
                                       Example: gap_allowance=1 will include sections
                                       that fall in gaps of size 1 between synapses.

    Returns:
        float: Total surface area (in µm²).
        list: List of section objects included.
    """
    dendList = list(dendList)

    # === Step 1: Filter for synapses from the given type ===
    df_type = synSite_df[synSite_df['type'] == target_type]
    print(len(df_type), "total synapses in synSite_df")
    if df_type.empty:
        print(f"No synapses found for type '{target_type}'")
        return 0.0, []

    # === Step 2: Extract section indices from 'mappedSection' ===
    mapped_sec_ids = df_type['mappedSection'].dropna().apply(
        lambda s: int(re.search(r'\[(\d+)\]', s).group(1))
    )
    mapped_sec_ids = sorted(mapped_sec_ids.unique())

    # === Step 3: Expand spans with gap_allowance ===
    full_section_indices = set(mapped_sec_ids)
    for i in range(len(mapped_sec_ids) - 1):
        start, end = mapped_sec_ids[i], mapped_sec_ids[i + 1]
        gap = end - start - 1
        if gap > 0 and gap <= gap_allowance:
            # include all sections in the gap
            full_section_indices.update(range(start, end + 1))

    full_section_indices = sorted(full_section_indices)

    # === Step 4: Map indices to section objects ===
    selected_sections = []
    for i in full_section_indices:
        if i < len(dendList):
            selected_sections.append(dendList[i])
        else:
            print(f"⚠️ Section index {i} out of bounds")

    # === Step 5: Calculate surface area ===
    total_area = 0.0
    for sec in selected_sections:
        sec.push()
        area = 0
        for seg in sec:
            area += h.area(seg.x)
        total_area += area
        h.pop_section()
    print(f"total area: {total_area}")
    syn_count = int(len(df_type))
    density = syn_count/total_area
    print(f"density: {density}")

    # === Step 6: Optional SWC export ===
    if swc_filename:  # only write if a filename is provided
        swc_lines = []
        swc_id = 1
        for sec in selected_sections:
            h.define_shape()  # ensure 3D points are defined
            n3d = int(h.n3d(sec=sec))
            if n3d < 2:
                continue  # skip sections with no geometry

            for i in range(n3d):
                x = h.x3d(i, sec=sec)
                y = h.y3d(i, sec=sec)
                z = h.z3d(i, sec=sec)
                diam = h.diam3d(i, sec=sec)
                parent_id = swc_id - 1 if i > 0 else -1
                swc_lines.append(
                    f"{swc_id} 3 {x:.3f} {y:.3f} {z:.3f} {diam/2:.3f} {parent_id}"
                )
                swc_id += 1

        with open(swc_filename, 'w') as f:
            f.write("# SWC file of selected dendritic branches for type: " + target_type + "\n")
            f.write("\n".join(swc_lines))

        print(f"✅ Exported {len(swc_lines)} points to SWC: {swc_filename}")

    return total_area, selected_sections

def plot_nearest_neighbor_section_gaps(synSite_df, target_type=None, show_plot=True):
    """
    For each synapse, find its 1st, 2nd, and 3rd nearest neighbors (by section index)
    and count how many section indices are between them. Plot histogram distributions.

    Args:
        synSite_df (pd.DataFrame): Must have 'mappedSection' (e.g. "Section[123]") and 'type'.
        target_type (str, optional): Filter for synapses of this presynaptic type. 
                                     If None, use all synapses.
        show_plot (bool): Whether to display the histogram plot.

    Returns:
        dict of pd.Series: Gap counts for 1st, 2nd, and 3rd nearest neighbors.
    """
    df = synSite_df.copy()

    # --- Step 1: optional filtering ---
    if target_type is not None:
        df = df[df["type"] == target_type]
    if df.empty:
        print(f"No synapses found for type {target_type}")
        return {}

    # --- Step 2: extract integer section IDs ---
    sec_ids = df['mappedSection'].dropna().apply(
        lambda s: int(re.search(r'\[(\d+)\]', s).group(1))
    ).sort_values().to_numpy()

    if len(sec_ids) < 2:
        print("Not enough synapses for neighbor calculation")
        return {}

    # --- Step 3: compute neighbor gaps ---
    gaps_dict = {1: [], 2: [], 3: []}

    for sec in sec_ids:
        # Compute distances to all *other* synapses
        distances = np.abs(sec_ids - sec)
        distances = distances[distances != 0]  # remove self
        sorted_gaps = np.sort(distances) - 1   # sections between = distance - 1

        # Store first 3 neighbors (clamp at 0 for adjacency)
        for k in [1, 2, 3]:
            if len(sorted_gaps) >= k:
                gaps_dict[k].append(max(sorted_gaps[k-1], 0))

    gap_series_dict = {k: pd.Series(v) for k, v in gaps_dict.items()}

    # --- Step 4: histogram plots ---
    if show_plot:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
        for idx, k in enumerate([1, 2, 3]):
            series = gap_series_dict[k]
            if series.empty:
                continue
            bins = np.arange(series.max()+2) - 0.5
            axes[idx].hist(series, bins=bins, density=True,
                           edgecolor="black", alpha=0.7)
            axes[idx].set_xticks(range(int(series.max())+1))
            axes[idx].set_xlabel(f"Sections between {k}ᵗʰ neighbor")
            axes[idx].set_title(f"{k}ᵗʰ Nearest Neighbor")
            axes[idx].yaxis.set_major_formatter(
                plt.FuncFormatter(lambda y, _: f"{y*100:.0f}%")
            )

        axes[0].set_ylabel("Fraction of synapses (%)")
        fig.suptitle(f"Nearest neighbor section gaps ({target_type if target_type else 'All'})")
        plt.tight_layout()
        plt.show()

    return gap_series_dict

def plot_density_of_VPN_synsapses(csv_filename=None):
    """
    Plot synapse density per VPN neuron across all gaps in subplots,
    and save CSV with all data.
    """

    # === Hemibrain datasets ===
    hemibrain = {
        'DNp01': {
            'LC4': {0:  (2276, 4069.511812919178),
                    1:  (2276, 4282.966550454949),
                    2:  (2276, 4359.165045348288),
                    3:  (2276, 4565.127390189379),
                    5:  (2276, 4650.061220908717),
                    10: (2276, 4767.9772272870705),
                    15: (2276, 4861.093247981264)},
            'LPLC2': {0:  (1471, 1033.0515233120038),
                      1:  (1471, 1067.516037281448),
                      2:  (1471, 1204.044803113463),
                      3:  (1471, 1204.044803113463),
                      5:  (1471, 1469.2485644232079),
                      10: (1471, 1567.6714973928931),
                      15: (1471, 1567.6714973928931)}
        },
        'DNp03': {
            'LC4': {0:  (1000, 624.1071655186298),
                    1:  (1000, 719.9430074710078),
                    2:  (1000, 786.4483028460771),
                    3:  (1000, 829.2614310965629),
                    5:  (1000, 859.0808032094642),
                    10: (1000, 939.7238589980507),
                    15: (1000, 1003.3511017075788)},
            'LC22': {0:  (12, 7.683389844974828),
                     1:  (12, 8.018816530641567),
                     2:  (12, 8.018816530641567),
                     3:  (12, 8.57523828427733),
                     5:  (12, 8.57523828427733),
                     10: (12, 8.57523828427733),
                     15: (12, 8.57523828427733)},
            'LPLC1': {0:  (1375, 571.300931117936),
                      1:  (1375, 726.9734987816724),
                      2:  (1375, 794.9432694347323),
                      3:  (1375, 846.7839497729068),
                      5:  (1375, 965.6978667352806),
                      10: (1375, 1066.2787968080386),
                      15: (1375, 1146.1459822514814)},
            'LPLC4': {0:  (1006, 456.9673493109073),
                      1:  (1006, 513.9160920433771),
                      2:  (1006, 555.5588012388206),
                      3:  (1006, 600.5765066138594),
                      5:  (1006, 668.8605531745524),
                      10: (1006, 845.7731013148197),
                      15: (1006, 873.7227517049608)},
            'LPLC2': {0:  (0, 0),
                      1:  (0, 0),
                      2:  (0, 0),
                      3:  (0, 0),
                      5:  (0, 0),
                      10: (0, 0),
                      15: (0, 0)}
        }
    }

    # === FAFB datasets ===
    fafb = {
        'DNp01': {
            'LC4': {0:  (1401, 2191.9033),
                    1:  (1401, 2317.822),
                    2:  (1401, 2402.430),
                    3:  (1401, 2435.508),
                    5:  (1401, 2494.512),
                    10: (1401, 2592.412),
                    15: (1401, 2669.618)},
            'LPLC2': {0:  (615, 279.915),
                      1:  (615, 297.549),
                      2:  (615, 298.690),
                      3:  (615, 298.690),
                      5:  (615, 320.496),
                      10: (615, 326.149),
                      15: (615, 370.633)}
        },
        'DNp03': {
            'LC4': {0:  (712, 768.712),
                    1:  (712, 846.996),
                    2:  (712, 922.976),
                    3:  (712, 980.662),
                    5:  (712, 1033.397),
                    10: (712, 1093.706),
                    15: (712, 1199.664)},
            'LC22': {0:  (190, 219.652),
                     1:  (190, 233.318),
                     2:  (190, 261.262),
                     3:  (190, 275.238),
                     5:  (190, 302.793),
                     10: (190, 381.800),
                     15: (190, 398.352)},
            'LPLC1': {0:  (1269, 1216.482),
                      1:  (1269, 1377.142),
                      2:  (1269, 1499.105),
                      3:  (1269, 1617.578),
                      5:  (1269, 1784.122),
                      10: (1269, 1982.700),
                      15: (1269, 2247.848)},
            'LPLC4': {0:  (1103, 881.192),
                      1:  (1103, 984.014),
                      2:  (1103, 1045.567),
                      3:  (1103, 1101.629),
                      5:  (1103, 1139.075),
                      10: (1103, 1250.313),
                      15: (1103, 1312.374)},
            'LPLC2': {0:  (14, 31.142),
                      1:  (14, 31.846),
                      2:  (14, 31.846),
                      3:  (14, 31.846),
                      5:  (14, 38.404),
                      10: (14, 38.404),
                      15: (14, 38.404)}
        }
    }

    colors = {'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
              'LC22': '#FFF000', 'LPLC4': '#00FF00'}

    def get_marker(source):
        return 'o' if source == 'Hemibrain' else 's'

    # Determine all unique gaps across datasets
    all_gaps = sorted(set(
        list(fafb['DNp01']['LC4'].keys()) + list(fafb['DNp03']['LC4'].keys())
    ))

    n_gaps = len(all_gaps)
    fig, axes = plt.subplots(1, n_gaps, figsize=(4*n_gaps, 6), sharey=True)
    if n_gaps == 1:
        axes = [axes]

    rows = []

    for idx, g in enumerate(all_gaps):
        ax = axes[idx]
        for dataset_name, dataset in [('FAFB', fafb), ('Hemibrain', hemibrain)]:
            for neuron_name, neuron_data in dataset.items():
                for region, gap_dict in neuron_data.items():
                    if g not in gap_dict:
                        continue
                    synapses, area = gap_dict[g]
                    density = synapses / area if area != 0 else 0
                    color = colors.get(region, 'gray')
                    marker = get_marker(dataset_name)
                    facecolor = 'none' if neuron_name == 'DNp01' else color
                    ax.scatter(region, density, marker=marker, s=100,
                               facecolors=facecolor, edgecolors=color,
                               label=f'{dataset_name} {neuron_name} {region}' if idx==0 else None)

                    # Save all data
                    rows.append({
                        'dataset': dataset_name,
                        'neuron': neuron_name,
                        'region': region,
                        'gap': g,
                        'synapse_count': synapses,
                        'surface_area_um2': area,
                        'density_syn_per_um2': density
                    })

        ax.set_title(f'Gap {g}')
        ax.set_xlabel('Region')
        ax.grid(True, axis='y')
        if idx == 0:
            ax.set_ylabel('Synapse Density (synapses/µm²)')

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

    if csv_filename:
        df = pd.DataFrame(rows)
        df.to_csv(f'datafiles/simulationData/{csv_filename}', index=False)
        print(f"✅ Synapse density for all gaps saved to {csv_filename}")


#Figure 8
def plotAllSynsAndhighlightsynsFig8(neuron_name):
    if neuron_name == "DNp01":
        synSite_df = pd.read_csv("datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_time.csv', header=None)

    # if neuron_name == "DNp01":
    #     synSite_df = pd.read_csv("datafiles/simulationData/DNp01_final_sims/singlesynactivation_fly_5/SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
    #     synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation_fly_5/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
    #     synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation_fly_5/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
    #     time_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation_fly_5/SINGLE_EXP2_time.csv', header=None)

    elif neuron_name == "DNp01_hemi":
        synSite_df = pd.read_csv("datafiles/simulationData/DNp01_hemi_final_sims/singlesynactivation/SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp01_hemi_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp01_hemi_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp01_hemi_final_sims/singlesynactivation/SINGLE_EXP2_time.csv', header=None)

    # elif neuron_name == "DNp03":
    #     synSite_df = pd.read_csv("datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
    #     synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
    #     synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
    #     time_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_time.csv', header=None)

    if neuron_name == "DNp03":
        synSite_df = pd.read_csv("datafiles/simulationData/DNp03_final_sims/singlesynactivation_fly_6/SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation_fly_6/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation_fly_6/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation_fly_6/SINGLE_EXP2_time.csv', header=None)
        
    elif neuron_name == "DNp03_hemi":
        synSite_df = pd.read_csv("datafiles/simulationData/DNp03_hemi_final_sims/singlesynactivation/SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp03_hemi_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp03_hemi_final_sims/singlesynactivation/SINGLE_EXP2_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp03_hemi_final_sims/singlesynactivation/SINGLE_EXP2_time.csv', header=None)
        
    synSite_df, indexList = filterSynapsesToCriteria(['LC4','LC6', 'LC22', 'LPLC1','LPLC2', 'LPLC4'], synSite_df)
    i = 0
    for index, row in synSite_df.iterrows():
        print(row, indexList[i])
        i += 1
    synCurr_siz_df = synCurr_siz_df.iloc[indexList]
    synCurr_synpt_df = synCurr_synpt_df.iloc[indexList]
    time_df = time_df.iloc[indexList]
    
    SYNPTMAX = max(synCurr_synpt_df.max(axis=1))
    maxSort = synCurr_synpt_df.max(axis=1).sort_values(ascending=False)

    SIZMAX = synCurr_siz_df.iloc[:, 2002:].to_numpy().max()
    TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=None, gridPos=[2,1], axList=None)

    TESTWINDOW = plotSynapseLocations(synSite_df, colorKey=1, sizeArg=6)

    # clock.sleep(1200)
    if neuron_name == "DNp01":
        RSYN_IDX_1 = 200
        RSYN_IDX_2 = 39
        RSYN_IDX_3 = 481
        RSYN_IDX_4 = 1625
        RSYN_IDX_5 = 2014
        RSYN_IDX_6 = 746
    elif neuron_name == "DNp01_hemi":
        RSYN_IDX_1 = 501
        RSYN_IDX_2 = 39
        RSYN_IDX_4 = 1853
        RSYN_IDX_5 = 2270
    elif neuron_name == "DNp02":
        RSYN_IDX_1 = 3855
        RSYN_IDX_2 = 3834
        RSYN_IDX_3 = 2852
    elif neuron_name == "DNp03":
        RSYN_IDX_1 = 5
        RSYN_IDX_2 = 150
        RSYN_IDX_4 = 530
        RSYN_IDX_5 = 535
        RSYN_IDX_7 = 920
        RSYN_IDX_8 = 935
        RSYN_IDX_10 = 1828
        RSYN_IDX_11 = 1832
        RSYN_IDX_13 = 1948
        RSYN_IDX_14 = 2085
    elif neuron_name == "DNp03_hemi":
        #LC4
        RSYN_IDX_1 = 830
        RSYN_IDX_2 = 1335
       #LC22
        RSYN_IDX_4 = 1203
        RSYN_IDX_5 = 1212
       #LPLC1
        RSYN_IDX_7 = 2609
        RSYN_IDX_8 = 1362
        #LPLC4
        RSYN_IDX_10 = 2808
        RSYN_IDX_11 = 1276
       
        # RSYN_IDX_13 = 1948
        # RSYN_IDX_14 = 2085
    elif neuron_name == "DNp04":
        RSYN_IDX_1 = 292
        RSYN_IDX_2 = 153
    elif neuron_name == "DNp06":
        RSYN_IDX_1 = 2695
        RSYN_IDX_2 = 15276
    else:
        RSYN_IDX_1 = random.randint(0, synSite_df.shape[0])
        RSYN_IDX_2 = random.randint(0, synSite_df.shape[0])
        RSYN_IDX_3 = random.randint(0, synSite_df.shape[0])

    if neuron_name == "DNp01":
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_1], colorKey=3, sizeArg=4)#syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_2], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=4)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_4], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_5], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)

        TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_siz_df.iloc[RSYN_IDX_4], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_siz_df.iloc[RSYN_IDX_5], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_siz_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_siz_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_synpt_df, '#899499', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_synpt_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_synpt_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_synpt_df.iloc[RSYN_IDX_4], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_synpt_df.iloc[RSYN_IDX_5], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        
        TF_AXLIST[1].set_xlim([48, 60])
        TF_AXLIST[1].set_ylim([-77.44, SIZMAX+0.02])
        TF_AXLIST[0].set_xlim([48, 60])
        TF_AXLIST[0].set_ylim([-77.6, SYNPTMAX+0.15])
        TF_AXLIST[0].set_title('Recording at Synapse Site')
        TF_AXLIST[1].set_title('Recording at SIZ')

        TF_AXLIST[0].set_xlabel('Time (ms)')
        TF_AXLIST[1].set_xlabel('Time (ms)')
        TF_AXLIST[0].set_ylabel('Voltage (mV)')
        TF_AXLIST[1].set_ylabel('Voltage (mV)')

    elif neuron_name == "DNp01_hemi":
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_1], colorKey=3, sizeArg=4)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_2], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=4)

        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_4], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_5], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)

        TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_siz_df.iloc[RSYN_IDX_4], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_siz_df.iloc[RSYN_IDX_5], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_siz_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_siz_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_synpt_df, '#899499', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_synpt_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_synpt_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_synpt_df.iloc[RSYN_IDX_4], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_synpt_df.iloc[RSYN_IDX_5], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF_AXLIST[1].set_xlim([48, 60])
        TF_AXLIST[1].set_ylim([-77.48, SIZMAX+0.02])
        TF_AXLIST[0].set_xlim([48, 60])
        TF_AXLIST[0].set_ylim([-77.48, SYNPTMAX+0.06])
        TF_AXLIST[0].set_title('Recording at Synapse Site')
        TF_AXLIST[1].set_title('Recording at SIZ')

        TF_AXLIST[0].set_xlabel('Time (ms)')
        TF_AXLIST[1].set_xlabel('Time (ms)')
        TF_AXLIST[0].set_ylabel('Voltage (mV)')
        TF_AXLIST[1].set_ylabel('Voltage (mV)')


    elif neuron_name =="DNp03":
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_1], colorKey=3, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_2], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_4], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_5], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_7], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_8], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_10], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_11], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_13], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_14], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_siz_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_siz_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_siz_df.iloc[RSYN_IDX_4], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_siz_df.iloc[RSYN_IDX_5], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_7], synCurr_siz_df.iloc[RSYN_IDX_7], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_8], synCurr_siz_df.iloc[RSYN_IDX_8], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_siz_df.iloc[RSYN_IDX_10], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_siz_df.iloc[RSYN_IDX_11], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_13], synCurr_siz_df.iloc[RSYN_IDX_13], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_14], synCurr_siz_df.iloc[RSYN_IDX_14], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_synpt_df, '#899499', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_synpt_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_synpt_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_synpt_df.iloc[RSYN_IDX_4], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_synpt_df.iloc[RSYN_IDX_5], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_7], synCurr_synpt_df.iloc[RSYN_IDX_7], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_8], synCurr_synpt_df.iloc[RSYN_IDX_8], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_synpt_df.iloc[RSYN_IDX_10], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_synpt_df.iloc[RSYN_IDX_11], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_13], synCurr_synpt_df.iloc[RSYN_IDX_13], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_14], synCurr_synpt_df.iloc[RSYN_IDX_14], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF_AXLIST[1].set_xlim([48, 60])
        TF_AXLIST[1].set_ylim([-74.52, SIZMAX+0.04])
        TF_AXLIST[0].set_xlim([48, 60])
        TF_AXLIST[0].set_ylim([-74.6, SYNPTMAX+0.15])

        TF_AXLIST[0].set_title('Recording at Synapse Site')
        TF_AXLIST[1].set_title('Recording at SIZ')

        TF_AXLIST[0].set_xlabel('Time (ms)')
        TF_AXLIST[1].set_xlabel('Time (ms)')
        TF_AXLIST[0].set_ylabel('Voltage (mV)')
        TF_AXLIST[1].set_ylabel('Voltage (mV)')

    elif neuron_name =="DNp03_hemi":
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_1], colorKey=3, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_2], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_4], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_5], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_7], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_8], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_10], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_11], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_siz_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_siz_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_siz_df.iloc[RSYN_IDX_4], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_siz_df.iloc[RSYN_IDX_5], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_7], synCurr_siz_df.iloc[RSYN_IDX_7], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_8], synCurr_siz_df.iloc[RSYN_IDX_8], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_siz_df.iloc[RSYN_IDX_10], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_siz_df.iloc[RSYN_IDX_11], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_synpt_df, '#899499', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_synpt_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_synpt_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_synpt_df.iloc[RSYN_IDX_4], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_synpt_df.iloc[RSYN_IDX_5], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_7], synCurr_synpt_df.iloc[RSYN_IDX_7], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_8], synCurr_synpt_df.iloc[RSYN_IDX_8], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_synpt_df.iloc[RSYN_IDX_10], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_synpt_df.iloc[RSYN_IDX_11], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF_AXLIST[1].set_xlim([48, 60])
        TF_AXLIST[1].set_ylim([-73.5, SIZMAX+0.04])
        TF_AXLIST[0].set_xlim([48, 60])
        TF_AXLIST[0].set_ylim([-73.6, SYNPTMAX+0.04])
        TF_AXLIST[0].set_title('Recording at Synapse Site')
        TF_AXLIST[1].set_title('Recording at SIZ')

        TF_AXLIST[0].set_xlabel('Time (ms)')
        TF_AXLIST[1].set_xlabel('Time (ms)')
        TF_AXLIST[0].set_ylabel('Voltage (mV)')
        TF_AXLIST[1].set_ylabel('Voltage (mV)')
        
    np.random.seed(20)  # fixed seed for reproducibility

    for syn_type, group in synSite_df.groupby("type"):
        idx_min = group.index.min()
        idx_max = group.index.max()
        total = len(group)
        
        # sample up to 5 random indices from this group
        rand_idx = np.random.choice(group.index, size=min(15, total), replace=False)
        
        print(f"Type {syn_type}: index range {idx_min} – {idx_max}, total {total} synapses")
        print(f"  Random 5 indices: {sorted(rand_idx.tolist())}")
        plt.show()

def plot_peak_depol_at_SIZ_vs_dist_to_SIZ(neuron_name, sizSection):
    # === Color mapping for neuron types ===
    color_dict = {
        "LC4": "#0000FF",    # blue
        "LPLC4": "#00FF00",  # green
        "LC22": "#FFF000",   # yellow
        "LPLC2": "#FFA500",  # orange
        "LPLC1": "#FF0D00"   # red
    }

    # === Load data based on neuron ===
    if neuron_name == "DNp01":
        base = "datafiles/simulationData/DNp01_final_sims/singlesynactivation/"
        # base = "datafiles/simulationData/DNp01_final_sims/singlesynactivation_fly_5/"
    elif neuron_name == "DNp01_hemi":
        base = "datafiles/simulationData/DNp01_hemi_final_sims/singlesynactivation/"
    elif neuron_name == "DNp03":
        base = "datafiles/simulationData/DNp03_final_sims/singlesynactivation/"
        # base = "datafiles/simulationData/DNp03_final_sims/singlesynactivation_fly_6/"
    elif neuron_name == "DNp03_hemi":
        base = "datafiles/simulationData/DNp03_hemi_final_sims/singlesynactivation/"
    else:
        raise ValueError(f"Unknown neuron: {neuron_name}")

    synSite_df = pd.read_csv(base + "SINGLE_EXP2_synapse_sites.csv", dtype={'pre': str})
    synCurr_siz_df = pd.read_csv(base + "SINGLE_EXP2_synaptic_currents_siz.csv", header=None)
    synCurr_synpt_df = pd.read_csv(base + "SINGLE_EXP2_synaptic_currents_synpt.csv", header=None)
    time_df = pd.read_csv(base + "SINGLE_EXP2_time.csv", header=None)

    t_vec = time_df.values.flatten()

    dist_to_SIZ = []
    peak_depol_at_SIZ = []
    peak_depol_at_syn = []
    electrotonic_distances = []
    colors = []

    # === First figure: distance to SIZ ===
    gs_fig = gridspec.GridSpec(3, 1)
    TEST_IND = plt.figure(figsize=(7, 10))
    ax1_fig = TEST_IND.add_subplot(gs_fig[0, :])
    ax2_fig = TEST_IND.add_subplot(gs_fig[1, :], sharey=ax1_fig)
    ax3_fig = TEST_IND.add_subplot(gs_fig[2, :])

    # --- Loop over synapses ---
    for index, row in synSite_df.iterrows():
        neuron_type = row['type']  # column indicating neuron type
        color = color_dict.get(neuron_type, '#000000')  # fallback black
        colors.append(color)

        idx = synSite_df[(synSite_df['mappedSection'] == row['mappedSection']) &
                         (synSite_df['mappedSegRangeVar'] == row['mappedSegRangeVar'])].index[0]

        rest_SIZ = np.array(synCurr_siz_df.iloc[idx][1996])
        rest_syn = np.array(synCurr_synpt_df.iloc[idx][1996])

        synSec = str2sec(row['mappedSection'])
        distToSIZ = h.distance(sizSection(0.05), synSec(row['mappedSegRangeVar']))
        dist_to_SIZ.append(distToSIZ)

        trace_syn = np.array(synCurr_synpt_df.iloc[idx])
        trace_siz = np.array(synCurr_siz_df.iloc[idx])

        peak_depol_syn = np.max(trace_syn[2002:]) - rest_syn
        peak_depol_SIZ = np.max(trace_siz[2002:]) - rest_SIZ
        peak_depol_at_syn.append(peak_depol_syn)
        peak_depol_at_SIZ.append(peak_depol_SIZ)

        # Electrotonic distance
        if peak_depol_syn > 0:
            log_ratio = np.log(peak_depol_syn / peak_depol_SIZ)
        else:
            log_ratio = np.nan
        electrotonic_distances.append(log_ratio)

        # Plot solid points with assigned color
        ax1_fig.plot(distToSIZ, peak_depol_syn, 'o', color=color, markeredgecolor='none')
        ax2_fig.plot(distToSIZ, peak_depol_SIZ, 'o', color=color, markeredgecolor='none')
        ax3_fig.plot(distToSIZ, peak_depol_SIZ, 'o', color=color, markeredgecolor='none')

    # --- Linear regressions for subplots ---
    def add_regression(ax, x, y):
        if len(x) > 1:
            coeffs = np.polyfit(x, y, 1)
            fit_line = np.polyval(coeffs, x)
            y_pred = fit_line
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            R2 = 1 - ss_res / ss_tot
            ax.plot(np.sort(x), np.polyval(coeffs, np.sort(x)), 'r-', lw=2,
                    label=f"y={coeffs[0]:.3f}x+{coeffs[1]:.2f}, R²={R2:.2f}")
            ax.legend()

    x_arr = np.array(dist_to_SIZ)
    add_regression(ax2_fig, x_arr, np.array(peak_depol_at_SIZ))
    add_regression(ax3_fig, x_arr, np.array(peak_depol_at_SIZ))

    ax1_fig.set_xlabel('Distance from Synapse to SIZ (µm)')
    ax2_fig.set_xlabel('Distance from Synapse to SIZ (µm)')
    ax1_fig.set_ylabel('Depolarization (mV)')
    ax2_fig.set_ylabel('Depolarization (mV)')
    ax1_fig.set_title(f'Individual Synapse Depolarizations at Synapse Location ({neuron_name})')
    ax2_fig.set_title(f'Individual Synapse Depolarizations at SIZ ({neuron_name})')
    ax3_fig.set_xlabel('Distance from Synapse to SIZ (µm)')
    ax3_fig.set_ylabel('Depolarization (mV)')
    ax3_fig.set_title(f'Individual Synapse Depolarizations at SIZ ({neuron_name})')

    TEST_IND.tight_layout()

    # === Second figure: electrotonic distance ===
    gs2 = gridspec.GridSpec(3, 1)
    fig2 = plt.figure(figsize=(7, 10))
    ax4 = fig2.add_subplot(gs2[0, :])
    ax5 = fig2.add_subplot(gs2[1, :], sharey=ax4)
    ax6 = fig2.add_subplot(gs2[2, :])

    ax4.scatter(electrotonic_distances, peak_depol_at_syn, c=colors, edgecolors='none', s=50)
    ax4.set_ylabel("Depolarization at Synapse (mV)")
    ax4.set_title(f"Depolarization vs Electrotonic Distance ({neuron_name})")

    ax5.scatter(electrotonic_distances, peak_depol_at_SIZ, c=colors, edgecolors='none', s=50)
    ax5.set_ylabel("Depolarization at SIZ (mV)")

    ax6.scatter(electrotonic_distances, peak_depol_at_SIZ, c=colors, edgecolors='none', s=50)
    ax6.set_xlabel("Electrotonic Distance")
    ax6.set_ylabel("Depolarization at SIZ (mV)")
    ax6.set_ylim(bottom=0)

    fig2.tight_layout()

    # === Third figure: electrotonic vs physical distance ===
    fig3 = plt.figure(figsize=(18, 5))
    ax_dend = fig3.add_subplot(1, 3, 1)
    ax_phys = fig3.add_subplot(1, 3, 2)
    ax_elec = fig3.add_subplot(1, 3, 3)

    ax_dend.scatter(x_arr, peak_depol_at_syn, c=colors, edgecolors='none', s=50)
    add_regression(ax_dend, x_arr, np.array(peak_depol_at_syn))
    ax_dend.set_xlabel("Distance from Synapse to SIZ (µm)")
    ax_dend.set_ylabel("Depolarization at Dendrite (mV)")
    ax_dend.set_title(f"Dendritic Depol vs Physical Distance ({neuron_name})")

    ax_phys.scatter(x_arr, peak_depol_at_SIZ, c=colors, edgecolors='none', s=50)
    add_regression(ax_phys, x_arr, np.array(peak_depol_at_SIZ))
    ax_phys.set_xlabel("Distance from Synapse to SIZ (µm)")
    ax_phys.set_ylabel("Depolarization at SIZ (mV)")
    ax_phys.set_title(f"SIZ Depol vs Physical Distance ({neuron_name})")

    ax_elec.scatter(x_arr, electrotonic_distances, c=colors, edgecolors='none', s=50)
    add_regression(ax_elec, x_arr, np.array(electrotonic_distances))
    ax_elec.set_xlabel("Distance from Synapse to SIZ (µm)")
    ax_elec.set_ylabel("Electrotonic Distance (log ratio)")
    ax_elec.set_title(f"Electrotonic vs Physical Distance ({neuron_name})")

    fig3.tight_layout()

    fig4, (ax7, ax8) = plt.subplots(1, 2, figsize=(16, 6))

    # attach "pre" IDs, depol, type, and electrotonic distance
    depol_data = pd.DataFrame({
        "pre": synSite_df["pre"].astype(str),   # make sure dtype is str
        "peak_depol_SIZ": peak_depol_at_SIZ,
        "electrotonic_distance": electrotonic_distances,
        "type": synSite_df["type"]
    })

    # enforce an explicit order for categories (preserves file order)
    order = depol_data["pre"].unique().tolist()

    # --- Subplot 1: boxplot of depol grouped by pre ---
    sns.boxplot(
        data=depol_data,
        x="pre",
        y="peak_depol_SIZ",
        order=order,
        ax=ax7,
        color="lightgray",
        showfliers=False
    )

    sns.stripplot(
        data=depol_data,
        x="pre",
        y="peak_depol_SIZ",
        order=order,
        hue="type",
        palette=color_dict,
        jitter=True,
        dodge=False,
        alpha=0.8,
        size=5,
        ax=ax7
    )

    # compute group averages (aligned to the 'order')
    avg_depol = depol_data.groupby("pre")["peak_depol_SIZ"].mean().reindex(order).values
    xpos = np.arange(len(order))
    ax7.scatter(
        xpos,
        avg_depol,
        color="red",
        edgecolors="black",
        linewidths=0.6,
        zorder=11,
        s=80,
        label="Neuron average"
    )

    ax7.set_xlabel("Presynaptic Neuron (pre ID)")
    ax7.set_ylabel("Depolarization at SIZ (mV)")
    ax7.set_title(f"Synapse Groups by Presynaptic Neuron ({neuron_name})")
    ax7.tick_params(axis='x', rotation=90)

    # remove legend here (we'll add one on the right)
    ax7.legend([], [], frameon=False)

    # --- Subplot 2: boxplot of electrotonic distance grouped by pre ---
    sns.boxplot(
        data=depol_data,
        x="pre",
        y="electrotonic_distance",
        order=order,
        ax=ax8,
        color="lightgray",
        showfliers=False
    )

    sns.stripplot(
        data=depol_data,
        x="pre",
        y="electrotonic_distance",
        order=order,
        hue="type",
        palette=color_dict,
        jitter=True,
        dodge=False,
        alpha=0.8,
        size=5,
        ax=ax8
    )

    # compute group averages (skip NaNs automatically)
    avg_elec = depol_data.groupby("pre")["electrotonic_distance"].mean().reindex(order).values
    # only plot non-NaN averages
    valid_mask = ~np.isnan(avg_elec)
    ax8.scatter(
        xpos[valid_mask],
        avg_elec[valid_mask],
        color="red",
        edgecolors="black",
        linewidths=0.6,
        zorder=11,
        s=80,
        label="Neuron average"
    )

    ax8.set_xlabel("Presynaptic Neuron (pre ID)")
    ax8.set_ylabel("Electrotonic Distance (log ratio)")
    ax8.set_title(f"Electrotonic Distance by Presynaptic Neuron ({neuron_name})")
    ax8.tick_params(axis='x', rotation=90)

    # build legend on right subplot: type entries + red average marker
    handles, labels = ax8.get_legend_handles_labels()
    avg_handle = Line2D([0], [0], marker='o', color='red', markeredgecolor='black',
                        markersize=8, linestyle='', label='Neuron average')
    # add avg_handle to front so it shows in legend
    ax8.legend(handles + [avg_handle], labels + ['Neuron average'], title="Type",
            bbox_to_anchor=(1.05, 1), loc='upper left')

    fig4.tight_layout()

    plt.show()


#Figure9 
def plot_partner_VPN_activations(neuron_name, MODE= 'all/'):
    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_partner'
        quintList = returnQuints(dir_path)
    elif neuron_name == "DNp01_hemi":
        dir_path = 'datafiles/simulationData/DNp01_hemi_final_sims/rand_vs_partner/'+neuron_name+'_partner'
        quintList = iterQuintsV4(dir_path)
    elif neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_partner'
        quintList = returnQuints(dir_path)
    elif neuron_name == "DNp03_hemi":
        dir_path = 'datafiles/simulationData/DNp03_hemi_final_sims/rand_vs_partner/'+neuron_name+'_partner'
        quintList = iterQuintsV4(dir_path)
    fig = plt.figure()

    LC4_0 = True
    LC4_R = True
    LC6_0 = True
    LC22_0 = True
    LPLC1_0 = True
    LPLC1_R = True
    LPLC4_R = True
    LPLC2_R = True
    LPLC2_0 = True
    LPLC4_0 = True
    Rand_0 = True
    if neuron_name in ["DNp01", "DNp03"]:
        for quint_ct in range(len(quintList)):
            quint = quintList[quint_ct]
            vSIZ_np = numpy.array(quint[1])
            peakV_siz = numpy.max(vSIZ_np)
            Vrest_siz = vSIZ_np[0, 1996]
            synct = quint[4].shape
            synct = synct[0]
            
            JIT = numpy.random.default_rng().normal(0, 0.1, 1)
            
            if MODE == "VPN/":
                if quint[4]['type'].iloc[0] == "LC4":
                    color = '#0000FF'
                    if LC4_0:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color, label="LC4")#, markersize=2)
                        LC4_0 = False
                    else: 
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color)
                elif quint[4]['type'].iloc[0] == "LC6":
                    color = '#FF00E9'
                    if LC6_0:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color, label="LC6")#, markersize=2)
                        LC6_0 = False
                    else: 
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color)
                elif quint[4]['type'].iloc[0] == "LC22":
                    color = '#FFF000'
                    if LC22_0:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color, label="LC22")#, markersize=2)
                        LC22_0 = False
                    else: 
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color)
                elif quint[4]['type'].iloc[0] == "LPLC2":
                    color = '#FFA500'
                    if LPLC2_0:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color, label="LPLC2")
                        LPLC2_0 = False
                    else:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color)
                elif quint[4]['type'].iloc[0] == "LPLC1":
                    color = '#FF0D00'
                    if LPLC1_0:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color, label="LPLC1")
                        LPLC1_0 = False
                    else:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color)
                elif quint[4]['type'].iloc[0] == "LPLC4":
                    color = '#00FF00'
                    if LPLC4_0:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color, label="LPLC4")
                        LPLC4_0 = False
                    else:
                        plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'o', color=color)
            else:
                plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'ks')


    elif neuron_name in ["DNp01_hemi", "DNp03_hemi"]:
        for k, quint in enumerate(quintList):
            print(f"Processing group {k+1}")
            time_df, vSIZ_np, vsoma_np, synSite = quint

            peakV_siz = np.max(vSIZ_np)
            Vrest_siz = vSIZ_np[0, 1996]

            synct = synSite.shape[0]
            JIT = numpy.random.default_rng().normal(0, 0.1, 1)

            if MODE == "VPN/":
                syn_type = synSite['type'].iloc[0] if not synSite.empty else "unknown"

                color_map = {
                    "LC4": '#0000FF',
                    "LC6": '#FF00E9',
                    "LC22": '#FFF000',
                    "LPLC2": '#FFA500',
                    "LPLC1": '#FF0D00',
                    "LPLC4": '#00FF00'
                }

                label_flags = {
                    "LC4": LC4_0,
                    "LC4": LC4_R,
                    "LC6": LC6_0,
                    "LC22": LC22_0,
                    "LPLC2": LPLC2_0,
                    "LPLC1": LPLC1_0,
                    "LPLC4": LPLC4_0,
                    "LPLC4": LPLC4_R,
                    "LPLC2": LPLC2_R,
                    "LPLC1": LPLC1_R
                }

                color = color_map.get(syn_type, 'black')

                if label_flags.get(syn_type, False):
                    plt.plot(synct + JIT + 0.25, peakV_siz - Vrest_siz, 'o', color=color, label=syn_type)
                    exec(f"{syn_type}_0 = False")  # update flag dynamically
                else:
                    plt.plot(synct + JIT + 0.25, peakV_siz - Vrest_siz, 'o', color=color)
            else:
                plt.plot(synct + JIT + 0.25, peakV_siz - Vrest_siz, 'ks')
        
    if neuron_name == "DNp01":
        plt.xlim([0, 26])
    elif neuron_name == "DNp03":
        plt.xlim([0, 40])
    elif neuron_name == "DNp03_hemi":
        plt.xlim([0, 40])
    plt.ylabel("Peak Depolarization at SIZ (mV)")
    plt.xlabel("Number of Synapses")
    plt.title("Recording at SIZ")
    plt.suptitle("Voltage Response of {} to Activations of Partnered VPN Synapses".format(neuron_name))
    plt.tight_layout()
    plt.legend()
    plt.show()

def plot_rand_vs_partner_SIZ(neuron_name, MODE='all/'):
    if neuron_name in ["DNp01", "DNp01_hemi"]:
        suffixes = ["_partner", "_rand"]
    elif neuron_name in ["DNp03", "DNp03_hemi"]:
        suffixes = ["_partner", "_rand"]
    else:
        raise ValueError("Invalid neuron_name")

    color_map = {
        "LC4": '#0000FF',
        "LC6": '#FF00E9',
        "LC22": '#FFF000',
        "LPLC2": '#FFA500',
        "LPLC1": '#FF0D00',
        "LPLC4": '#00FF00'
    }

    fig, ax = plt.subplots(figsize=(8,6))
    plotted_labels = set()  # track which labels have already been added to the legend

    # First, find max synapse count for partner
    max_syn_partner = 0
    partner_dir = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner'
    # partner_dir = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner_fly_6'
    for quint in iterQuintsV4(partner_dir):
        time_df, vSIZ_np, vsoma_np, synSite = quint
        synct = synSite.shape[0]
        if synct > max_syn_partner:
            max_syn_partner = synct

    for suffix in suffixes:
        if "partner" in suffix:
            dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner'
            # dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner_fly_6'
        else:
            dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_rand'
            # dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_rand_fly_6'

        quintList = iterQuintsV4(dir_path)

        for k, quint in enumerate(quintList):
            print(f"Processing group {k+1}")
            time_df, vSIZ_np, vsoma_np, synSite = quint

            idx_30 = np.argmin(np.abs(time_df - 30))
            peakV_siz = np.max(vSIZ_np[0, idx_30:])
            Vrest_siz = vSIZ_np[0, idx_30]
            synct = synSite.shape[0]
            JIT = np.random.default_rng().normal(0, 0.15, 1)

            if suffix == "_partner" and MODE == "VPN/":
                syn_type = synSite['type'].iloc[0] if not synSite.empty else "unknown"
                color = color_map.get(syn_type, 'black')
                label = syn_type if syn_type not in plotted_labels else None
                if label is not None:
                    plotted_labels.add(syn_type)
                ax.plot(synct + JIT, peakV_siz - Vrest_siz, 'o', color=color, label=label)

            elif suffix == "_rand":
                # Only plot random points with syn_count <= max_syn_partner
                if synct <= max_syn_partner:
                    label = 'Random' if 'Random' not in plotted_labels else None
                    if label is not None:
                        plotted_labels.add('Random')
                    ax.plot(synct + JIT, peakV_siz - Vrest_siz, 'ko', label=label)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel("Peak Depolarization at SIZ (mV)")
    ax.set_xlabel("Number of Synapses")
    ax.set_title(f"Recording at SIZ of {neuron_name}")
    ax.legend()

    # Set x-limits
    xlims = {
    "DNp01": max_syn_partner + 0.5,
    "DNp03": max_syn_partner + 0.5,
    "DNp03_hemi": max_syn_partner + 0.5,
    "DNp01_hemi": max_syn_partner + 0.5
}
    ax.set_xlim([0, xlims.get(neuron_name, max_syn_partner)])

    plt.tight_layout()
    plt.show()

def plot_rand_vs_partner_SOMA(neuron_name, MODE='all/'):
    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp01_hemi":
        dir_path = 'datafiles/simulationData/DNp01_hemi_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp03_hemi":
        dir_path = 'datafiles/simulationData/DNp03_hemi_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    quintList = returnQuints(dir_path)

    fig = plt.figure()
    
    LC4_0 = True
    LC6_0 = True
    LC22_0 = True
    LPLC1_0 = True
    LPLC2_0 = True
    LPLC4_0 = True
    

    for quint_ct in range(len(quintList)):
        quint = quintList[quint_ct]
        vSoma_np = numpy.array(quint[2])
        peakV_soma = numpy.max(vSoma_np)
        Vrest_soma = vSoma_np[0, 1996]
        synct = quint[4].shape
        synct = synct[0]
        JIT = numpy.random.default_rng().normal(0, 0.1, 1)
        
        if MODE == "VPN/":
            if quint[4]['type'].iloc[0] == "LC4":
                color = '#0000FF'
                if LC4_0:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color, label="LC4")#, markersize=2)
                    LC4_0 = False
                else: 
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color)
            elif quint[4]['type'].iloc[0] == "LC6":
                color = '#FF00E9'
                if LC6_0:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color, label="LC6")#, markersize=2)
                    LC6_0 = False
                else: 
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color)
            elif quint[4]['type'].iloc[0] == "LC22":
                color = '#FFF000'
                if LC22_0:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color, label="LC22")#, markersize=2)
                    LC22_0 = False
                else: 
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color)
            elif quint[4]['type'].iloc[0] == "LPLC2":
                color = '#FFA500'
                if LPLC2_0:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color, label="LPLC2")
                    LPLC2_0 = False
                else:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color)
            elif quint[4]['type'].iloc[0] == "LPLC1":
                color = '#FF0D00'
                if LPLC1_0:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color, label="LPLC1")
                    LPLC1_0 = False
                else:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color)
            elif quint[4]['type'].iloc[0] == "LPLC4":
                color = '#00FF00'
                if LPLC4_0:
                    plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'o', color=color, label="LPLC4")
                    LPLC4_0 = False
                else:
                    plt.plot(synct+JIT+0.25,peakV_soma-Vrest_soma, 'o', color=color)
        else:
            plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'ks')

    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_rand'
    if neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_rand'
    quintList = returnQuints(dir_path)

    Rand_0 = True

    for quint_ct in range(len(quintList)):
        quint = quintList[quint_ct]
        vSoma_np = numpy.array(quint[2])
        peakV_soma = numpy.max(vSoma_np)
        Vrest_soma = vSoma_np[0, 1996]
        synct = quint[4].shape
        synct = synct[0]
        if synct >= 40:
            break

        JIT = numpy.random.default_rng().normal(0, 0.1, 1)
        if Rand_0:
            plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'ks', label='Random')#, markersize=2)
            Rand_0 = False
        else:
            plt.plot(synct+JIT+0.25, peakV_soma-Vrest_soma, 'ks')

    plt.gcf().get_axes()[0].spines['top'].set_visible(False)
    plt.gcf().get_axes()[0].spines['right'].set_visible(False)
    plt.xlim([0, 40])
    plt.ylabel("Peak Depolarization at SOMA (mV)")
    plt.xlabel("Number of Synapses")
    plt.title("Recording at SOMA")
    # plt.suptitle("Voltage Response of {} to Activations of Partnered VPN Synapses".format(neuron_name))
    plt.suptitle("Voltage Response of {} to Activations of Partnered vs Random Synapses".format(neuron_name))
    plt.tight_layout()
    plt.legend()
    plt.show()


def Plot_DN_VPN_syns(neuron_name, location=None, csv_filename=None):
    data = []  # List to collect rows for DataFrame
    if neuron_name in ["DNp01", "DNp03", "DNp01_hemi", "DNp03_hemi"]:
        dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_incRandSynCt_VPN'
        # dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_incRandSynCt_VPN_fly_5'
        quintList = list(iterQuintsV5(dir_path))
        print(f"Total quints loaded: {len(quintList)}", flush=True)
        total_quints = len(quintList)

        for k, quint in enumerate(quintList):
            print(f"Processing group {k+1} / {total_quints}", flush=True)
            time_df, vSIZ_np, vSoma_np, synSite = quint

            # find index closest to 30 ms
            idx_30 = np.argmin(np.abs(time_df - 30))

            # peak depolarizations
            peakV_siz = np.max(vSIZ_np[0, idx_30:])
            Vrest_siz = vSIZ_np[0, idx_30]
            peakV_soma = np.max(vSoma_np[0, idx_30:])
            Vrest_soma = vSoma_np[0, idx_30]

            synct = synSite.shape[0]
            syn_type = synSite['type'].iloc[0] if not synSite.empty and 'type' in synSite.columns else "na"

            # Store for group plotting later
            data.append({
                'syn_count': synct,
                'deflection_soma': peakV_soma - Vrest_soma,
                'deflection_siz': peakV_siz - Vrest_siz,
                'syn_type': syn_type
            })

    else:
        raise ValueError(f"Unknown neuron name: {neuron_name}")

    # === Group-level Plot ===
    df = pd.DataFrame(data)

    if csv_filename:
        df.to_csv(f'datafiles/simulationData/{neuron_name}_final_sims/{csv_filename}', index=False)
        print(f"✅ Data saved to CSV: {csv_filename}")

    color_map = {
        'na': '#000000',
        'LC4': '#0000FF',
        'LPLC2': '#FFA500',
        'LPLC1': '#FF0D00',
        'LC22': '#FFF000',
        'LPLC4': '#00FF00',
        'LC6': '#FF00E9'
    }

    plt.figure(figsize=(8,6))
    for syn_type, group in df.groupby('syn_type'):
        color = color_map.get(syn_type, '#888888')  # fallback color
        if location == "Soma":
            plt.scatter(group['syn_count'], group['deflection_soma'], label=syn_type, color=color, s=50)
        elif location == "SIZ":
            plt.scatter(group['syn_count'], group['deflection_siz'], label=syn_type, color=color, s=50)
        elif location == "Both":
            plt.scatter(group['syn_count'], group['deflection_soma'], label=f"{syn_type} Soma", color=color, marker='o', s=50)
            plt.scatter(group['syn_count'], group['deflection_siz'], label=f"{syn_type} SIZ", color=color, marker='s', s=50)

    plt.xlabel('Number of synapses activated')
    plt.ylabel('Voltage deflection from resting membrane potential (mV)')
    plt.title(f'Depolarization of {neuron_name} as a Function of Synapse Count')
    plt.legend()
    plt.tight_layout()
    plt.show()

def shunting_rate_vs_density_with_exp_curve(neuron_list, dataset_type="FAFB"):
    """
    Fits exponential depolarization vs synapse count, extracts tau, and plots:
        1) DNp01: Depol vs syn_count
        2) DNp03: Depol vs syn_count
        3) Shunting tau vs density (linear regression, LC22/LPLC2 omitted only here)
    Correctly uses density for the chosen dataset (FAFB or Hemibrain).
    """


    color_map = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    density_df = pd.read_csv("datafiles/simulationData/VPN-DNs_synapse_density.csv")
    results = []

    def exp_rise(x, A, tau, C):
        return A * (1 - np.exp(-x / tau)) + C

    fig, axs = plt.subplots(1, 3, figsize=(20, 6))
    ax_dnp01, ax_dnp03, ax_tau = axs

    for neuron_name in neuron_list:
        depol_csv = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
        df = pd.read_csv(depol_csv)

        # Keep LC22/LPLC2 in depol subplots but drop them for DNp03 density/tau calculation later
        drop_for_tau = ["LC22", "LPLC2"] if "DNp03" in neuron_name else []

        for vpn, vpn_df in df.groupby("syn_type"):
            x = vpn_df["syn_count"].values
            y = vpn_df["deflection_siz"].values
            if len(x) < 3:
                continue
            x0, y0 = x[0], y[0]

            def exp_rise_constrained(x, A, tau):
                C = y0 - A * (1 - np.exp(-x0 / tau))
                return A * (1 - np.exp(-x / tau)) + C

            p0 = [max(y) - y0, np.mean(x)]
            try:
                popt, _ = curve_fit(exp_rise_constrained, x, y, p0=p0, maxfev=5000)
                A, tau = popt
                C = y0 - A * (1 - np.exp(-x0 / tau))
            except RuntimeError:
                tau = np.nan
                A, C = max(y), y[0]

            # Lookup density for correct dataset
            if dataset_type.lower() == "fafb":
                sub_density = density_df[
                    (density_df["neuron"] == neuron_name) &
                    (density_df["region"] == vpn) &
                    (density_df["dataset"].str.lower() == "fafb")
                ]
            else:  # Hemibrain
                hemi_name = neuron_name.replace("_hemi","")
                sub_density = density_df[
                    (density_df["neuron"] == hemi_name) &
                    (density_df["region"] == vpn) &
                    (density_df["dataset"].str.lower() == "hemibrain")
                ]

            density = sub_density["density_syn_per_um2"].values[0] if not sub_density.empty else np.nan
            print(f"Neuron: {neuron_name}, VPN: {vpn}, Dataset: {dataset_type}, Density used: {density:.3f} syn/µm²")

            results.append({
                "neuron": neuron_name, "vpn": vpn, "tau": tau, "density": density, "A": A
            })

            # --- Depol subplots ---
            x_fit = np.linspace(min(x), max(x), 200)
            y_fit = exp_rise(x_fit, A, tau, C)
            color = color_map.get(vpn, "gray")
            marker_shape = "o" if "DNp01" in neuron_name else "s"
            ax = ax_dnp01 if "DNp01" in neuron_name else ax_dnp03
            ax.scatter(x, y, color=color, marker=marker_shape, s=80, edgecolors="k",
                       label=f"{vpn} ({neuron_name})")
            ax.plot(x_fit, y_fit, color=color, linestyle="--")
            y_tau = exp_rise(tau, A, tau, C)
            ax.plot(tau, y_tau, marker="X", color=color, markersize=10)
            ax.text(tau, y_tau, f"τ={tau:.1f}", fontsize=9, ha='left', va='bottom')

    # --- Tau vs density subplot (omit LC22/LPLC2 only for DNp03) ---
    res_df = pd.DataFrame(results)
    tau_df = res_df[~((res_df['vpn'].isin(drop_for_tau)) & res_df['neuron'].str.contains("DNp03"))]

    for _, row in tau_df.iterrows():
        color = color_map.get(row["vpn"], "gray")
        shape = "o" if "DNp01" in row["neuron"] else "s"
        ax_tau.scatter(row["density"], row["tau"], s=100, marker=shape, color=color, edgecolors="k",
                       label=f"{row['vpn']} ({row['neuron']})")

    # Linear regression and R²
    valid = tau_df.dropna(subset=["density", "tau"])
    if not valid.empty:
        X = valid["density"].values.reshape(-1,1)
        y = valid["tau"].values
        lr = LinearRegression().fit(X, y)
        x_line = np.linspace(min(X), max(X), 100)
        y_line = lr.predict(x_line.reshape(-1,1))
        ax_tau.plot(x_line, y_line, color='k', linestyle='--', label='Linear fit')
        r2 = lr.score(X, y)
        ax_tau.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax_tau.transAxes,
                    fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

    # --- Formatting ---
    for ax, title in zip([ax_dnp01, ax_dnp03, ax_tau],
                         ['DNp01 VPNs', 'DNp03 VPNs', 'Shunting τ vs Density']):
        ax.set_title(title)
        ax.grid(alpha=0.5, linestyle='--')

    ax_dnp01.set_xlabel("Synapse count"); ax_dnp01.set_ylabel("Depol at SIZ (mV)")
    ax_dnp03.set_xlabel("Synapse count"); ax_dnp03.set_ylabel("Depol at SIZ (mV)")
    ax_tau.set_xlabel("Density (synapses / µm²)"); ax_tau.set_ylabel("Shunting τ (ms)")

    # Combined legend for all
    handles, labels = [], []
    for ax in [ax_dnp01, ax_dnp03, ax_tau]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h); labels.extend(l)
    by_label = dict(zip(labels, handles))
    ax_tau.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05,1), loc='upper left')

    plt.tight_layout()
    plt.show()
    return res_df

def shunting_rate_vs_density(neuron_list, gap=0, dataset_type="FAFB"):
    """
    For each neuron in neuron_list:
    - Computes Depol vs synapse count
    - Computes shunting tau for all VPN-DNs across all gaps
    - Generates:
        1) Figure: Depol vs synapse count (DNp01, DNp03) + Shunting τ vs density at specified gap
        2) Figure: Shunting τ vs density for each unique gap in the CSV
    """
    color_map = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    # Load density data
    density_df = pd.read_csv("datafiles/simulationData/VPN-DNs_synapse_density.csv")
    density_df = density_df[density_df['dataset'].str.lower() == dataset_type.lower()]

    results = []

    # --- Compute τ for all neurons, VPNs, and gaps ---
    for neuron_name in neuron_list:
        if dataset_type == "FAFB":
            depol_csv = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
            df = pd.read_csv(depol_csv)
        elif dataset_type == "hemibrain":
            depol_csv = f"datafiles/simulationData/{neuron_name}_hemi_final_sims/{neuron_name}_hemi_VPN_population_depol_data.csv"
            df = pd.read_csv(depol_csv)

        for vpn, vpn_df in df.groupby("syn_type"):
            # Get all gaps for this neuron and VPN
            sub_density_all_gaps = density_df[
                (density_df["neuron"] == neuron_name) &
                (density_df["region"] == vpn)
            ]

            for _, sub_density_row in sub_density_all_gaps.iterrows():
                gap_val = sub_density_row["gap"]
                density = sub_density_row["density_syn_per_um2"]

                x = vpn_df["syn_count"].values
                y = vpn_df["deflection_siz"].values
                x0, y0 = x[0], y[0]

                def exp_rise_constrained(x, A, tau):
                    C = y0 - A * (1 - np.exp(-x0 / tau))
                    return A * (1 - np.exp(-x / tau)) + C

                p0 = [max(y) - y0, np.mean(x)]
                try:
                    popt, _ = curve_fit(exp_rise_constrained, x, y, p0=p0, maxfev=5000)
                    A, tau = popt
                except RuntimeError:
                    tau = np.nan
                    A = max(y)
                C = y0 - A * (1 - np.exp(-x0 / (tau if not np.isnan(tau) else 1)))

                results.append({
                    "neuron": neuron_name,
                    "vpn": vpn,
                    "tau": tau,
                    "density": density,
                    "A": A,
                    "gap": gap_val
                })

    res_df = pd.DataFrame(results)

    # --- First Figure: Depol vs Synapse Count + Shunting τ at specified gap ---
    fig, axs = plt.subplots(1, 3, figsize=(20, 6))
    ax_dnp01, ax_dnp03, ax_all = axs

    # Plot Depol vs Synapse Count
    for neuron_name in neuron_list:
        if dataset_type == "FAFB":
            depol_csv = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
            df = pd.read_csv(depol_csv)
        elif dataset_type == "hemibrain":
            depol_csv = f"datafiles/simulationData/{neuron_name}_hemi_final_sims/{neuron_name}_hemi_VPN_population_depol_data.csv"
            df = pd.read_csv(depol_csv)
        for vpn, vpn_df in df.groupby("syn_type"):
            x = vpn_df["syn_count"].values
            y = vpn_df["deflection_siz"].values
            color = color_map.get(vpn, "gray")
            marker_shape = "o" if neuron_name == "DNp01" else "s"
            ax = ax_dnp01 if neuron_name == "DNp01" else ax_dnp03
            ax.scatter(x, y, color=color, marker=marker_shape, s=80, edgecolors="k",
                       label=f"{vpn} ({neuron_name})")

    # Third subplot: Shunting τ vs density at the provided gap
    df_plot = res_df[res_df['gap'] == gap]
    df_plot = df_plot[~((df_plot['neuron'] == "DNp03") & (df_plot['vpn'].isin(["LPLC2"])))]
    df_plot = df_plot[~((df_plot['neuron'] == "DNp03") & (df_plot['vpn'].isin(["LC22"])))]

    for _, row in df_plot.iterrows():
        color = color_map.get(row["vpn"], "gray")
        shape = "o" if row["neuron"] == "DNp01" else "s"
        ax_all.scatter(row["density"], row["tau"], s=100, marker=shape, color=color,
                       edgecolors="k", label=f"{row['vpn']} ({row['neuron']})")

    valid = df_plot.dropna(subset=["density", "tau"])
    if not valid.empty:
        X = valid["density"].values.reshape(-1, 1)
        y = valid["tau"].values
        lr = LinearRegression().fit(X, y)
        x_line = np.linspace(min(X), max(X), 100)
        y_line = lr.predict(x_line.reshape(-1, 1))
        ax_all.plot(x_line, y_line, color='k', linestyle='--', label="Linear fit")
        r2 = lr.score(X, y)
        ax_all.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax_all.transAxes,
                    fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

    for ax, title in zip([ax_dnp01, ax_dnp03, ax_all], ['DNp01 VPNs', 'DNp03 VPNs', f'All VPN-DNs (gap={gap})']):
        ax.set_xlabel("Synapse count / Density (τ subplot)")
        ax.set_ylabel("Depol at SIZ (mV) / Shunting τ (ms)")
        ax.grid(alpha=0.5, linestyle='--')

    handles, labels = [], []
    for ax in [ax_dnp01, ax_dnp03, ax_all]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    ax_all.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

    # --- Second Figure: Shunting τ vs Density for all gaps ---
    unique_gaps = sorted(res_df['gap'].unique())
    n_gaps = len(unique_gaps)
    n_cols = 2
    n_rows = math.ceil(n_gaps / n_cols)

    fig2, axs2 = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))

    # Flatten axs2 to a 1D array for easy indexing
    axs2 = axs2.flatten()
    if len(unique_gaps) == 1:
        axs2 = [axs2]

    for i, gap_val in enumerate(unique_gaps):
        ax = axs2[i]
        df_plot = res_df[res_df['gap'] == gap_val]
        df_plot = df_plot[~((df_plot['neuron'] == "DNp03") & (df_plot['vpn'].isin(["LPLC2","LC22"])))]

        for _, row in df_plot.iterrows():
            color = color_map.get(row["vpn"], "gray")
            shape = "o" if row["neuron"] == "DNp01" else "s"
            ax.scatter(row["density"], row["tau"], s=100, marker=shape, color=color,
                       edgecolors="k", label=f"{row['vpn']} ({row['neuron']})")

        valid = df_plot.dropna(subset=["density","tau"])
        if not valid.empty:
            X = valid["density"].values.reshape(-1, 1)
            y = valid["tau"].values
            lr = LinearRegression().fit(X, y)
            x_line = np.linspace(min(X), max(X), 100)
            y_line = lr.predict(x_line.reshape(-1,1))
            ax.plot(x_line, y_line, color='k', linestyle='--', label="Linear fit")
            r2 = lr.score(X, y)
            ax.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax.transAxes,
                    fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

        ax.set_xlabel("Density (synapses / µm²)")
        ax.set_ylabel("Shunting τ (ms)")
        ax.set_title(f"Gap = {gap_val}")
        ax.grid(alpha=0.5, linestyle='--')

    handles, labels = [], []
    for ax in axs2:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    axs2[-1].legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    plt.show()

    return res_df

# New plotting for COM neighbor activations
def Plot_DN_VPN_neighborhood(neuron_name, VPN, location="SIZ", csv_filename=None,
                              max_neighbors=9):
    """
    Reads all k-level subdirectories from COM_neighborhood_activations output,
    plots peak SIZ (or soma) depolarization vs synapse count, colored by k.

    Parameters
    ----------
    neuron_name   : str  — e.g. "DNp01"
    VPN           : str  — e.g. "LC4", "LPLC2"
    location      : str  — "SIZ", "Soma", or "Both"
    csv_filename  : str  — optional; saves collected data to CSV
    max_neighbors : int  — highest k to look for (default 9)
    """

    base_path = f'datafiles/simulationData/{neuron_name}_final_sims/COM_neighborhood'
    data = []

    for k in range(1, max_neighbors + 1):
        dir_path = os.path.join(base_path, f'{VPN}_k{k:02d}')

        if not os.path.isdir(dir_path):
            print(f"⚠️ Directory not found, skipping: {dir_path}")
            continue

        quint_iter = list(iterQuintsNeighborhood(dir_path))
        print(f"k={k}: loaded {len(quint_iter)} trials", flush=True)

        for quint in quint_iter:
            time_df, vSIZ_np, vSoma_np, synSite_df, k_val, trial_num, seed_id = quint

            # Baseline index at onset (50 ms) — use 30 ms window before onset as rest
            idx_onset = np.argmin(np.abs(time_df.values.flatten() - 50))

            peakV_siz  = np.max(vSIZ_np[0,  idx_onset:])
            Vrest_siz  = vSIZ_np[0, idx_onset]
            peakV_soma = np.max(vSoma_np[0, idx_onset:])
            Vrest_soma = vSoma_np[0, idx_onset]

            syn_count = synSite_df.shape[0]

            data.append({
                'k':               k_val,
                'trial':           trial_num,
                'seed_id':         seed_id,
                'syn_count':       syn_count,
                'deflection_siz':  peakV_siz  - Vrest_siz,
                'deflection_soma': peakV_soma - Vrest_soma,
                'group_size':      k_val + 1,
            })

    if not data:
        print("No data loaded. Check your directory paths.")
        return

    df = pd.DataFrame(data)

    # ── Optional CSV export ────────────────────────────────────────────────────
    if csv_filename:
        save_path = os.path.join(
            f'datafiles/simulationData/{neuron_name}_final_sims', csv_filename
        )
        df.to_csv(save_path, index=False)
        print(f"✅ Data saved to: {save_path}")

    # ── Color map: one color per k level ──────────────────────────────────────
    cmap = plt.cm.viridis
    k_vals = sorted(df['k'].unique())
    norm   = plt.Normalize(vmin=min(k_vals), vmax=max(k_vals))
    color_map = {k: cmap(norm(k)) for k in k_vals}

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 6))

    for k_val, group in df.groupby('k'):
        color = color_map[k_val]
        label = f'k={k_val} ({k_val+1} cells)'

        if location == "SIZ":
            ax.scatter(group['syn_count'], group['deflection_siz'],
                       label=label, color=color, s=60, alpha=0.85)
        elif location == "Soma":
            ax.scatter(group['syn_count'], group['deflection_soma'],
                       label=label, color=color, s=60, alpha=0.85)
        elif location == "Both":
            ax.scatter(group['syn_count'], group['deflection_siz'],
                       label=f'{label} SIZ',  color=color, marker='s', s=60, alpha=0.85)
            ax.scatter(group['syn_count'], group['deflection_soma'],
                       label=f'{label} Soma', color=color, marker='o', s=60, alpha=0.85,
                       facecolors='none', linewidths=1.5)

    ax.set_xlabel('Number of synapses activated', fontsize=12)
    ax.set_ylabel('Voltage deflection from resting potential (mV)', fontsize=12)
    ax.set_title(f'{neuron_name} — {VPN} neighborhood activation\n'
                 f'Peak depolarization at {location} vs synapse count', fontsize=12)
    ax.legend(title='Neighbor count', bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=9, title_fontsize=10)
    plt.tight_layout()
    plt.show()

    return df

def Plot_DN_VPN_neighborhood_grid(neuron_name, VPN_list, location="SIZ", 
                                   csv_filename=None, max_neighbors=9):
    """
    Reads COM_neighborhood_activations output for each VPN in VPN_list and
    plots them as side-by-side subplots sharing the same axes scale.

    Parameters
    ----------
    neuron_name   : str       — e.g. "DNp01"
    VPN_list      : list[str] — e.g. ["LC4", "LPLC2"] or ["LC4", "LPLC1", "LPLC4"]
    location      : str       — "SIZ", "Soma", or "Both"
    csv_filename  : str       — optional; saves combined data to CSV
    max_neighbors : int       — highest k to look for (default 9)
    """

    base_path = f'datafiles/simulationData/{neuron_name}_final_sims/COM_neighborhood'

    # ── 1. Load all data upfront ───────────────────────────────────────────────
    all_data = {}   # VPN -> DataFrame

    for VPN in VPN_list:
        data = []

        for k in range(1, max_neighbors + 1):
            dir_path = os.path.join(base_path, f'{VPN}_k{k:02d}')

            if not os.path.isdir(dir_path):
                print(f"⚠️  [{VPN}] k={k}: directory not found, skipping.")
                continue

            quint_iter = list(iterQuintsNeighborhood(dir_path))
            print(f"[{VPN}] k={k}: loaded {len(quint_iter)} trials", flush=True)

            for quint in quint_iter:
                time_df, vSIZ_np, vSoma_np, synSite_df, k_val, trial_num, seed_id = quint

                idx_onset  = np.argmin(np.abs(time_df.values.flatten() - 50))
                peakV_siz  = np.max(vSIZ_np[0,  idx_onset:])
                Vrest_siz  = vSIZ_np[0, idx_onset]
                peakV_soma = np.max(vSoma_np[0, idx_onset:])
                Vrest_soma = vSoma_np[0, idx_onset]

                data.append({
                    'VPN':             VPN,
                    'k':               k_val,
                    'trial':           trial_num,
                    'seed_id':         seed_id,
                    'syn_count':       synSite_df.shape[0],
                    'deflection_siz':  peakV_siz  - Vrest_siz,
                    'deflection_soma': peakV_soma - Vrest_soma,
                    'group_size':      k_val + 1,
                })

        if data:
            all_data[VPN] = pd.DataFrame(data)
        else:
            print(f"⚠️  [{VPN}] No data loaded — skipping from plot.")

    if not all_data:
        print("No data loaded for any VPN. Check your directory paths.")
        return None

    # ── 2. Optional CSV export ─────────────────────────────────────────────────
    combined_df = pd.concat(all_data.values(), ignore_index=True)

    if csv_filename:
        save_path = os.path.join(
            f'datafiles/simulationData/{neuron_name}_final_sims', csv_filename
        )
        combined_df.to_csv(save_path, index=False)
        print(f"✅ Combined data saved to: {save_path}")

    # ── 3. Shared color map across all subplots (viridis by k) ────────────────
    all_k   = sorted(combined_df['k'].unique())
    cmap    = plt.cm.viridis
    norm    = plt.Normalize(vmin=min(all_k), vmax=max(all_k))
    color_map = {k: cmap(norm(k)) for k in all_k}

    # Shared axis limits for fair comparison across VPNs
    if location == "SIZ":
        y_col = 'deflection_siz'
    elif location == "Soma":
        y_col = 'deflection_soma'
    else:  # Both — use whichever is larger for ylim
        y_col = None

    if y_col:
        y_max = combined_df[y_col].max() * 1.1
        y_min = min(0, combined_df[y_col].min() * 1.1)
    else:
        y_max = max(combined_df['deflection_siz'].max(),
                    combined_df['deflection_soma'].max()) * 1.1
        y_min = min(0, combined_df[['deflection_siz', 'deflection_soma']].min().min() * 1.1)

    x_max = combined_df['syn_count'].max() * 1.05
    x_min = 0

    # ── 4. Build subplot grid ──────────────────────────────────────────────────
    n_vpns   = len(all_data)
    fig_w    = max(5 * n_vpns, 8)
    fig, axes = plt.subplots(1, n_vpns, figsize=(fig_w, 5.5),
                             sharey=True, sharex=True)

    # Ensure axes is always iterable even for a single VPN
    if n_vpns == 1:
        axes = [axes]

    for ax, VPN in zip(axes, [v for v in VPN_list if v in all_data]):
        df = all_data[VPN]

        for k_val, group in df.groupby('k'):
            color = color_map[k_val]
            label = f'k={k_val} ({k_val+1} cells)'

            if location == "SIZ":
                ax.scatter(group['syn_count'], group['deflection_siz'],
                           label=label, color=color, s=55, alpha=0.85)
            elif location == "Soma":
                ax.scatter(group['syn_count'], group['deflection_soma'],
                           label=label, color=color, s=55, alpha=0.85)
            elif location == "Both":
                ax.scatter(group['syn_count'], group['deflection_siz'],
                           label=f'{label} SIZ',  color=color,
                           marker='s', s=55, alpha=0.85)
                ax.scatter(group['syn_count'], group['deflection_soma'],
                           label=f'{label} Soma', color=color,
                           marker='o', s=55, alpha=0.85,
                           facecolors='none', linewidths=1.5)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_title(VPN, fontsize=13, fontweight='bold')
        ax.set_xlabel('Synapses activated', fontsize=10)
        ax.tick_params(axis='both', labelsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Shared y-axis label on leftmost subplot only
    axes[0].set_ylabel('Voltage deflection from rest (mV)', fontsize=10)

    # ── 5. Shared colorbar legend (k levels) ──────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical',
                        fraction=0.02, pad=0.02, aspect=30)
    cbar.set_label('k (# neighbors)', fontsize=10)
    cbar.set_ticks(all_k)
    cbar.set_ticklabels([f'k={k}' for k in all_k], fontsize=8)

    fig.suptitle(
        f'{neuron_name} — COM neighborhood activation  |  {location}',
        fontsize=13, y=1.02
    )

    plt.tight_layout()
    plt.show()

    return combined_df

def shunting_rate_vs_density_with_neighborhood(neuron_list, VPN_config, gap=0,
                                                dataset_type="FAFB", max_neighbors=9):
    """
    Identical to shunting_rate_vs_density, but overlays neighborhood simulation
    data (from COM_neighborhood_activations) onto the DNp01 and DNp03 depol subplots.

    Neighborhood points are colored by k (viridis), while the existing
    VPN population data retains its original VPN color coding.
    No exponential fitting is applied to the neighborhood data.

    Parameters
    ----------
    neuron_list  : list[str]  — e.g. ["DNp01", "DNp03"]
    VPN_config   : dict       — maps each neuron to its neighborhood VPN list
                                e.g. {"DNp01": ["LC4", "LPLC2"],
                                      "DNp03": ["LC4", "LPLC1", "LPLC4"]}
    gap          : float      — gap value for Figure 1 right subplot
    dataset_type : str        — "FAFB" or "hemibrain"
    max_neighbors: int        — highest k to load (default 9)
    """

    color_map = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    # k colormap — viridis, shared across both depol subplots
    cmap   = plt.cm.viridis
    norm_k = plt.Normalize(vmin=1, vmax=max_neighbors)
    k_color = {k: cmap(norm_k(k)) for k in range(1, max_neighbors + 1)}

    # ── Load density data (identical to original) ─────────────────────────────
    density_df = pd.read_csv("datafiles/simulationData/VPN-DNs_synapse_density.csv")
    density_df = density_df[density_df['dataset'].str.lower() == dataset_type.lower()]

    results = []

    # ── Compute τ for all neurons, VPNs, and gaps (identical to original) ─────
    for neuron_name in neuron_list:
        if dataset_type == "FAFB":
            depol_csv = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
        elif dataset_type == "hemibrain":
            depol_csv = f"datafiles/simulationData/{neuron_name}_hemi_final_sims/{neuron_name}_hemi_VPN_population_depol_data.csv"

        df = pd.read_csv(depol_csv)

        for vpn, vpn_df in df.groupby("syn_type"):
            sub_density_all_gaps = density_df[
                (density_df["neuron"] == neuron_name) &
                (density_df["region"] == vpn)
            ]
            for _, sub_density_row in sub_density_all_gaps.iterrows():
                gap_val = sub_density_row["gap"]
                density = sub_density_row["density_syn_per_um2"]

                x = vpn_df["syn_count"].values
                y = vpn_df["deflection_siz"].values
                x0, y0 = x[0], y[0]

                def exp_rise_constrained(x, A, tau):
                    C = y0 - A * (1 - np.exp(-x0 / tau))
                    return A * (1 - np.exp(-x / tau)) + C

                p0 = [max(y) - y0, np.mean(x)]
                try:
                    popt, _ = curve_fit(exp_rise_constrained, x, y, p0=p0, maxfev=5000)
                    A, tau = popt
                except RuntimeError:
                    tau = np.nan
                    A = max(y)
                C = y0 - A * (1 - np.exp(-x0 / (tau if not np.isnan(tau) else 1)))

                results.append({
                    "neuron": neuron_name, "vpn": vpn, "tau": tau,
                    "density": density, "A": A, "gap": gap_val
                })

    res_df = pd.DataFrame(results)

    # ── Load neighborhood data for overlay ────────────────────────────────────
    neighborhood_data = {}   # neuron_name -> DataFrame

    for neuron_name in neuron_list:
        base_path = (f"datafiles/simulationData/{neuron_name}_final_sims/COM_neighborhood"
                     if dataset_type == "FAFB" else
                     f"datafiles/simulationData/{neuron_name}_hemi_final_sims/COM_neighborhood")

        rows = []
        for VPN in VPN_config.get(neuron_name, []):
            for k in range(1, max_neighbors + 1):
                dir_path = os.path.join(base_path, f'{VPN}_k{k:02d}')
                if not os.path.isdir(dir_path):
                    continue
                for quint in iterQuintsNeighborhood(dir_path):
                    time_df, vSIZ_np, vSoma_np, synSite_df, k_val, trial_num, seed_id = quint
                    idx_onset = np.argmin(np.abs(time_df.values.flatten() - 50))
                    rows.append({
                        'k':              k_val,
                        'syn_count':      synSite_df.shape[0],
                        'deflection_siz': np.max(vSIZ_np[0, idx_onset:]) - vSIZ_np[0, idx_onset],
                    })
        if rows:
            neighborhood_data[neuron_name] = pd.DataFrame(rows)
        else:
            print(f"⚠️  No neighborhood data found for {neuron_name}.")

    # ── FIGURE 1 (identical layout to original + neighborhood overlay) ─────────
    fig, axs = plt.subplots(1, 3, figsize=(20, 6))
    ax_dnp01, ax_dnp03, ax_all = axs
    neuron_ax = {neuron_list[0]: ax_dnp01, neuron_list[1]: ax_dnp03}

    # Original VPN population scatter (unchanged from original)
    for neuron_name in neuron_list:
        if dataset_type == "FAFB":
            depol_csv = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
        elif dataset_type == "hemibrain":
            depol_csv = f"datafiles/simulationData/{neuron_name}_hemi_final_sims/{neuron_name}_hemi_VPN_population_depol_data.csv"
        df = pd.read_csv(depol_csv)

        ax = neuron_ax[neuron_name]
        for vpn, vpn_df in df.groupby("syn_type"):
            color  = color_map.get(vpn, "gray")
            marker = "o" if neuron_name == "DNp01" else "s"
            ax.scatter(vpn_df["syn_count"].values, vpn_df["deflection_siz"].values,
                       color=color, marker=marker, s=80, edgecolors="k",
                       label=f"{vpn} ({neuron_name})", zorder=3)

    # Neighborhood overlay on top (colored by k, semi-transparent, no edge)
    for neuron_name, nd_df in neighborhood_data.items():
        ax = neuron_ax[neuron_name]
        for k_val, k_group in nd_df.groupby('k'):
            ax.scatter(k_group['syn_count'], k_group['deflection_siz'],
                       color=k_color[k_val], s=40, alpha=0.6,
                       marker='^', zorder=2,
                       label=f'neighborhood k={k_val}' if neuron_name == neuron_list[0] else '_nolegend_')

        # Colorbar per depol subplot
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm_k)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.03, aspect=25)
        cbar.set_label('k (neighbors)', fontsize=8)
        cbar.set_ticks(range(1, max_neighbors + 1))
        cbar.set_ticklabels([str(k) for k in range(1, max_neighbors + 1)], fontsize=7)

    # Third subplot: τ vs density (identical to original)
    df_plot = res_df[res_df['gap'] == gap]
    df_plot = df_plot[~((df_plot['neuron'] == "DNp03") & (df_plot['vpn'].isin(["LPLC2", "LC22"])))]

    for _, row in df_plot.iterrows():
        color = color_map.get(row["vpn"], "gray")
        shape = "o" if row["neuron"] == "DNp01" else "s"
        ax_all.scatter(row["density"], row["tau"], s=100, marker=shape,
                       color=color, edgecolors="k",
                       label=f"{row['vpn']} ({row['neuron']})")

    valid = df_plot.dropna(subset=["density", "tau"])
    if not valid.empty:
        X      = valid["density"].values.reshape(-1, 1)
        y_tau  = valid["tau"].values
        lr     = LinearRegression().fit(X, y_tau)
        x_line = np.linspace(X.min(), X.max(), 100)
        y_line = lr.predict(x_line.reshape(-1, 1))
        ax_all.plot(x_line, y_line, color='k', linestyle='--', label="Linear fit")
        r2 = lr.score(X, y_tau)
        ax_all.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax_all.transAxes,
                    fontsize=12, verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.7))

    for ax, title in zip([ax_dnp01, ax_dnp03, ax_all],
                          ['DNp01 VPNs', 'DNp03 VPNs', f'All VPN-DNs (gap={gap})']):
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Synapse count / Density (τ subplot)")
        ax.set_ylabel("Depol at SIZ (mV) / Shunting τ (ms)")
        ax.grid(alpha=0.5, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    handles, labels = [], []
    for ax in [ax_dnp01, ax_dnp03, ax_all]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    ax_all.legend(by_label.values(), by_label.keys(),
                  bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.show()

    # ── FIGURE 2: τ vs density for all gaps (identical to original) ───────────
    unique_gaps = sorted(res_df['gap'].unique())
    n_gaps = len(unique_gaps)
    n_cols = 2
    n_rows = math.ceil(n_gaps / n_cols)

    fig2, axs2 = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows))
    axs2 = axs2.flatten()
    if len(unique_gaps) == 1:
        axs2 = [axs2]

    for i, gap_val in enumerate(unique_gaps):
        ax = axs2[i]
        df_plot = res_df[res_df['gap'] == gap_val]
        df_plot = df_plot[~((df_plot['neuron'] == "DNp03") & (df_plot['vpn'].isin(["LPLC2", "LC22"])))]

        for _, row in df_plot.iterrows():
            color = color_map.get(row["vpn"], "gray")
            shape = "o" if row["neuron"] == "DNp01" else "s"
            ax.scatter(row["density"], row["tau"], s=100, marker=shape,
                       color=color, edgecolors="k",
                       label=f"{row['vpn']} ({row['neuron']})")

        valid = df_plot.dropna(subset=["density", "tau"])
        if not valid.empty:
            X      = valid["density"].values.reshape(-1, 1)
            y_tau  = valid["tau"].values
            lr     = LinearRegression().fit(X, y_tau)
            x_line = np.linspace(X.min(), X.max(), 100)
            y_line = lr.predict(x_line.reshape(-1, 1))
            ax.plot(x_line, y_line, color='k', linestyle='--', label="Linear fit")
            r2 = lr.score(X, y_tau)
            ax.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax.transAxes,
                    fontsize=12, verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.7))

        ax.set_xlabel("Density (synapses / µm²)")
        ax.set_ylabel("Shunting τ (ms)")
        ax.set_title(f"Gap = {gap_val}")
        ax.grid(alpha=0.5, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for j in range(n_gaps, len(axs2)):
        axs2[j].set_visible(False)

    handles, labels = [], []
    for ax in axs2[:n_gaps]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    axs2[-1].legend(by_label.values(), by_label.keys(),
                    bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.show()

    return res_df

def shunting_rate_vs_density_with_neighborhood_fit(neuron_list, VPN_config, gap=0,
                                                    dataset_type="FAFB", max_neighbors=9):
    color_map = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    cmap    = plt.cm.viridis
    norm_k  = plt.Normalize(vmin=1, vmax=max_neighbors)
    k_color = {k: cmap(norm_k(k)) for k in range(1, max_neighbors + 1)}

    density_df = pd.read_csv("datafiles/simulationData/VPN-DNs_synapse_density.csv")
    density_df = density_df[density_df['dataset'].str.lower() == dataset_type.lower()]

    # ── Shared exponential fit helper ─────────────────────────────────────────
    def fit_tau(x, y):
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]

        if len(x) < 3:
            return np.nan, np.nan

        x0, y0 = x[0], y[0]

        def exp_rise_constrained(x, A, tau):
            C = y0 - A * (1 - np.exp(-x0 / tau))
            return A * (1 - np.exp(-x / tau)) + C

        p0 = [max(y) - y0, np.mean(x)]
        try:
            popt, _ = curve_fit(exp_rise_constrained, x, y, p0=p0, maxfev=5000)
            A, tau = popt
        except RuntimeError:
            tau = np.nan
            A   = max(y)
        return A, tau

    # ── 1. Fit τ from original population data ────────────────────────────────
    results = []

    for neuron_name in neuron_list:
        depol_csv = (
            f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
            if dataset_type == "FAFB" else
            f"datafiles/simulationData/{neuron_name}_hemi_final_sims/{neuron_name}_hemi_VPN_population_depol_data.csv"
        )
        df = pd.read_csv(depol_csv)

        for vpn, vpn_df in df.groupby("syn_type"):
            sub_density_all_gaps = density_df[
                (density_df["neuron"] == neuron_name) &
                (density_df["region"] == vpn)
            ]
            for _, sub_density_row in sub_density_all_gaps.iterrows():
                x = vpn_df["syn_count"].values.astype(float)
                y = vpn_df["deflection_siz"].values.astype(float)
                A, tau = fit_tau(x, y)
                results.append({
                    "neuron":  neuron_name,
                    "vpn":     vpn,
                    "tau":     tau,
                    "density": sub_density_row["density_syn_per_um2"],
                    "A":       A,
                    "gap":     sub_density_row["gap"],
                    "source":  "population",
                })

    # ── 2. Load neighborhood data & fit τ per (neuron, VPN) ───────────────────
    neighborhood_data = {}
    nd_results        = []

    for neuron_name in neuron_list:
        base_path = (
            f"datafiles/simulationData/{neuron_name}_final_sims/COM_neighborhood"
            if dataset_type == "FAFB" else
            f"datafiles/simulationData/{neuron_name}_hemi_final_sims/COM_neighborhood"
        )

        rows = []
        for VPN in VPN_config.get(neuron_name, []):
            vpn_rows = []
            for k in range(1, max_neighbors + 1):
                dir_path = os.path.join(base_path, f'{VPN}_k{k:02d}')
                if not os.path.isdir(dir_path):
                    continue
                for quint in iterQuintsNeighborhood(dir_path):
                    time_df, vSIZ_np, vSoma_np, synSite_df, k_val, trial_num, seed_id = quint
                    idx_onset  = np.argmin(np.abs(time_df.values.flatten() - 50))
                    deflection = np.max(vSIZ_np[0, idx_onset:]) - vSIZ_np[0, idx_onset]

                    if not np.isfinite(deflection):
                        print(f"⚠️  Skipping trial {trial_num} seed {seed_id} "
                              f"(k={k_val}, {VPN}) — non-finite deflection.")
                        continue

                    entry = {
                        'vpn':            VPN,
                        'k':              k_val,
                        'syn_count':      synSite_df.shape[0],
                        'deflection_siz': deflection,
                    }
                    rows.append(entry)
                    vpn_rows.append(entry)

            # Fit τ to this VPN's full neighborhood cloud (all k pooled)
            if vpn_rows:
                vpn_nd_df = pd.DataFrame(vpn_rows).sort_values('syn_count')
                x_nd = vpn_nd_df['syn_count'].values.astype(float)
                y_nd = vpn_nd_df['deflection_siz'].values.astype(float)
                A_nd, tau_nd = fit_tau(x_nd, y_nd)

                if np.isnan(tau_nd):
                    print(f"⚠️  τ fit failed for {neuron_name} {VPN} neighborhood — "
                          f"will be excluded from density plot.")

                sub_density_all_gaps = density_df[
                    (density_df["neuron"] == neuron_name) &
                    (density_df["region"] == VPN)
                ]
                for _, sub_density_row in sub_density_all_gaps.iterrows():
                    nd_results.append({
                        "neuron":  neuron_name,
                        "vpn":     VPN,
                        "tau":     tau_nd,
                        "density": sub_density_row["density_syn_per_um2"],
                        "A":       A_nd,
                        "gap":     sub_density_row["gap"],
                        "source":  "neighborhood",
                    })

        if rows:
            neighborhood_data[neuron_name] = pd.DataFrame(rows)
        else:
            print(f"⚠️  No neighborhood data found for {neuron_name}.")

    res_df    = pd.DataFrame(results)
    nd_res_df = pd.DataFrame(nd_results) if nd_results else pd.DataFrame(
                    columns=["neuron","vpn","tau","density","A","gap","source"])
    combined_res_df = pd.concat([res_df, nd_res_df], ignore_index=True)

    # ── FIGURE 1 ───────────────────────────────────────────────────────────────
    fig, axs = plt.subplots(1, 3, figsize=(20, 6))
    ax_dnp01, ax_dnp03, ax_all = axs
    neuron_ax = {neuron_list[0]: ax_dnp01, neuron_list[1]: ax_dnp03}

    # Original population scatter
    for neuron_name in neuron_list:
        depol_csv = (
            f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_VPN_population_depol_data.csv"
            if dataset_type == "FAFB" else
            f"datafiles/simulationData/{neuron_name}_hemi_final_sims/{neuron_name}_hemi_VPN_population_depol_data.csv"
        )
        df = pd.read_csv(depol_csv)
        ax = neuron_ax[neuron_name]
        for vpn, vpn_df in df.groupby("syn_type"):
            color  = color_map.get(vpn, "gray")
            marker = "o" if neuron_name == "DNp01" else "s"
            ax.scatter(vpn_df["syn_count"].values, vpn_df["deflection_siz"].values,
                       color=color, marker=marker, s=80, edgecolors="k",
                       label=f"{vpn} ({neuron_name})", zorder=3)

    # Neighborhood overlay (colored by k)
    for neuron_name, nd_df in neighborhood_data.items():
        ax = neuron_ax[neuron_name]
        for k_val, k_group in nd_df.groupby('k'):
            ax.scatter(k_group['syn_count'], k_group['deflection_siz'],
                       color=k_color[k_val], s=40, alpha=0.6, marker='^', zorder=2,
                       label=f'neighborhood k={k_val}' if neuron_name == neuron_list[0] else '_nolegend_')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm_k)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.03, aspect=25)
        cbar.set_label('k (neighbors)', fontsize=8)
        cbar.set_ticks(range(1, max_neighbors + 1))
        cbar.set_ticklabels([str(k) for k in range(1, max_neighbors + 1)], fontsize=7)

    # Right subplot: population τ + neighborhood τ
    df_gap = combined_res_df[combined_res_df['gap'] == gap].copy()
    df_gap = df_gap[~((df_gap['neuron'] == "DNp03") & (df_gap['vpn'].isin(["LPLC2", "LC22"])))]

    for _, row in df_gap.iterrows():
        color = color_map.get(row["vpn"], "gray")
        if row["source"] == "population":
            marker = "o" if row["neuron"] == "DNp01" else "s"
            ax_all.scatter(row["density"], row["tau"], s=100, marker=marker,
                           color=color, edgecolors="k",
                           label=f"{row['vpn']} ({row['neuron']})")
        else:
            ax_all.scatter(row["density"], row["tau"], s=100, marker="^",
                           color=color, edgecolors="k", alpha=0.7,
                           label=f"{row['vpn']} ({row['neuron']}) neighborhood")

    valid = df_gap.dropna(subset=["density", "tau"])
    if not valid.empty:
        X      = valid["density"].values.reshape(-1, 1)
        y_tau  = valid["tau"].values
        lr     = LinearRegression().fit(X, y_tau)
        x_line = np.linspace(X.min(), X.max(), 100)
        y_line = lr.predict(x_line.reshape(-1, 1))
        ax_all.plot(x_line, y_line, color='k', linestyle='--', label="Linear fit (combined)")
        r2 = lr.score(X, y_tau)
        ax_all.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax_all.transAxes,
                    fontsize=12, verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.7))

    for ax, title in zip([ax_dnp01, ax_dnp03, ax_all],
                          ['DNp01 VPNs', 'DNp03 VPNs', f'All VPN-DNs (gap={gap})']):
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Synapse count / Density (τ subplot)")
        ax.set_ylabel("Depol at SIZ (mV) / Shunting τ (ms)")
        ax.grid(alpha=0.5, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    handles, labels = [], []
    for ax in [ax_dnp01, ax_dnp03, ax_all]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    ax_all.legend(by_label.values(), by_label.keys(),
                  bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.show()

    # ── FIGURE 2: τ vs density for every gap ──────────────────────────────────
    unique_gaps = sorted(combined_res_df['gap'].unique())
    n_gaps  = len(unique_gaps)
    n_cols  = 2
    n_rows  = math.ceil(n_gaps / n_cols)

    fig2, axs2 = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows))
    axs2 = axs2.flatten()
    if len(unique_gaps) == 1:
        axs2 = [axs2]

    for i, gap_val in enumerate(unique_gaps):
        ax   = axs2[i]
        df_g = combined_res_df[combined_res_df['gap'] == gap_val].copy()
        df_g = df_g[~((df_g['neuron'] == "DNp03") & (df_g['vpn'].isin(["LPLC2", "LC22"])))]

        for _, row in df_g.iterrows():
            color = color_map.get(row["vpn"], "gray")
            if row["source"] == "population":
                marker = "o" if row["neuron"] == "DNp01" else "s"
                ax.scatter(row["density"], row["tau"], s=100, marker=marker,
                           color=color, edgecolors="k",
                           label=f"{row['vpn']} ({row['neuron']})")
            else:
                ax.scatter(row["density"], row["tau"], s=100, marker="^",
                           color=color, edgecolors="k", alpha=0.7,
                           label=f"{row['vpn']} ({row['neuron']}) neighborhood")

        valid_g = df_g.dropna(subset=["density", "tau"])
        if not valid_g.empty:
            X_g    = valid_g["density"].values.reshape(-1, 1)
            y_g    = valid_g["tau"].values
            lr_g   = LinearRegression().fit(X_g, y_g)
            x_line = np.linspace(X_g.min(), X_g.max(), 100)
            y_line = lr_g.predict(x_line.reshape(-1, 1))
            ax.plot(x_line, y_line, color='k', linestyle='--', label="Linear fit (combined)")
            r2 = lr_g.score(X_g, y_g)
            ax.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax.transAxes,
                    fontsize=12, verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.7))

        ax.set_xlabel("Density (synapses / µm²)")
        ax.set_ylabel("Shunting τ (ms)")
        ax.set_title(f"Gap = {gap_val}")
        ax.grid(alpha=0.5, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for j in range(n_gaps, len(axs2)):
        axs2[j].set_visible(False)

    handles, labels = [], []
    for ax in axs2[:n_gaps]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    axs2[-1].legend(by_label.values(), by_label.keys(),
                    bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.show()

    return combined_res_df

def shunting_rate_vs_density_neighborhood_only(
    neuron_list,
    VPN_config,
    gap=0,
    dataset_type="FAFB",
    max_neighbors=9,
):
    color_map = {
        'LC4': '#0000FF',
        'LPLC2': '#FFA500',
        'LPLC1': '#FF0D00',
        'LC22': '#FFF000',
        'LPLC4': '#00FF00',
        'LC6': '#FF00E9'
    }

    density_df = pd.read_csv("datafiles/simulationData/VPN-DNs_synapse_density.csv")
    density_df = density_df[density_df['dataset'].str.lower() == dataset_type.lower()].copy()

    def canonical_vpn_name(vpn_name):
        return vpn_name.replace("_hemi", "")

    def fit_tau(x, y):
        """Fit a constrained exponential rise. Returns A, tau, x0, y0 so the
        fitted curve can be reconstructed later for plotting."""
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]

        if len(x) < 3:
            return np.nan, np.nan, np.nan, np.nan

        order = np.argsort(x)
        x = x[order]
        y = y[order]

        x0, y0 = x[0], y[0]

        def exp_rise_constrained(x, A, tau):
            C = y0 - A * (1 - np.exp(-x0 / tau))
            return A * (1 - np.exp(-x / tau)) + C

        p0 = [max(y) - y0, np.mean(x)]
        try:
            popt, _ = curve_fit(exp_rise_constrained, x, y, p0=p0, maxfev=5000)
            A, tau = popt
        except RuntimeError:
            tau = np.nan
            A = np.nan

        return A, tau, x0, y0

    def print_linfit_stats(label, x, y):
        """Run linregress and print R2, slope, equation, and p-value to terminal."""
        res = linregress(x, y)
        r2 = res.rvalue ** 2
        print(
            f"\n[Linear fit: {label}]\n"
            f"  Equation : tau = {res.slope:.5g} * density + {res.intercept:.5g}\n"
            f"  Slope    : {res.slope:.5g}\n"
            f"  Intercept: {res.intercept:.5g}\n"
            f"  R^2      : {r2:.5f}\n"
            f"  p-value  : {res.pvalue:.5g}\n"
        )
        return res, r2

    neighborhood_data = {}
    nd_results = []
    fit_curves = {}  # (neuron_name, vpn_label) -> (A, tau, x0, y0, xmin, xmax)

    for neuron_name in neuron_list:
        base_path = (
            f"datafiles/simulationData/{neuron_name}_final_sims/COM_neighborhood"
            if dataset_type.lower() == "fafb"
            else f"datafiles/simulationData/{neuron_name}_hemi_final_sims/COM_neighborhood"
        )

        rows = []

        for VPN in VPN_config.get(neuron_name, []):
            vpn_rows = []
            vpn_label = canonical_vpn_name(VPN)

            for k in range(1, max_neighbors + 1):
                dir_path = os.path.join(base_path, f"{VPN}_k{k:02d}")
                if not os.path.isdir(dir_path):
                    continue

                for quint in iterQuintsNeighborhood(dir_path):
                    time_df, vSIZ_np, vSoma_np, synSite_df, k_val, trial_num, seed_id = quint
                    idx_onset = np.argmin(np.abs(time_df.values.flatten() - 50))
                    deflection = np.max(vSIZ_np[0, idx_onset:]) - vSIZ_np[0, idx_onset]

                    if not np.isfinite(deflection):
                        print(
                            f"⚠️ Skipping trial {trial_num} seed {seed_id} "
                            f"(k={k_val}, {VPN}) — non-finite deflection."
                        )
                        continue

                    entry = {
                        "neuron": neuron_name,
                        "vpn": vpn_label,
                        "vpn_dir": VPN,
                        "k": k_val,
                        "syn_count": synSite_df.shape[0],
                        "deflection_siz": deflection,
                    }
                    rows.append(entry)
                    vpn_rows.append(entry)

            if vpn_rows:
                vpn_nd_df = pd.DataFrame(vpn_rows).sort_values("syn_count")
                x_nd = vpn_nd_df["syn_count"].values.astype(float)
                y_nd = vpn_nd_df["deflection_siz"].values.astype(float)
                A_nd, tau_nd, x0_nd, y0_nd = fit_tau(x_nd, y_nd)

                if np.isnan(tau_nd):
                    print(
                        f"⚠️ τ fit failed for {neuron_name} {VPN} neighborhood — "
                        f"will be excluded from density plot."
                    )
                else:
                    fit_curves[(neuron_name, vpn_label)] = (
                        A_nd, tau_nd, x0_nd, y0_nd, x_nd.min(), x_nd.max()
                    )

                sub_density = density_df[
                    (density_df["neuron"] == neuron_name) &
                    (density_df["region"] == vpn_label)
                ]

                if sub_density.empty:
                    print(
                        f"⚠️ No density match for neuron={neuron_name}, "
                        f"VPN={VPN}, canonical={vpn_label}, dataset={dataset_type}"
                    )

                for _, sub_density_row in sub_density.iterrows():
                    nd_results.append({
                        "neuron": neuron_name,
                        "vpn": vpn_label,
                        "tau": tau_nd,
                        "density": sub_density_row["density_syn_per_um2"],
                        "A": A_nd,
                        "gap": sub_density_row["gap"],
                        "source": "neighborhood",
                    })

        if rows:
            neighborhood_data[neuron_name] = pd.DataFrame(rows)
        else:
            print(f"⚠️ No neighborhood data found for {neuron_name}.")

    nd_res_df = pd.DataFrame(nd_results)
    if nd_res_df.empty:
        print("⚠️ No neighborhood τ results were generated.")
        return nd_res_df

    # ── FIGURE 1 ─────────────────────────────────────────────────────────────
    fig, axs = plt.subplots(1, 3, figsize=(20, 6))
    ax_dnp01, ax_dnp03, ax_all = axs
    neuron_ax = {neuron_list[0]: ax_dnp01, neuron_list[1]: ax_dnp03}

    for neuron_name, nd_df in neighborhood_data.items():
        ax = neuron_ax[neuron_name]
        for vpn, vpn_group in nd_df.groupby("vpn"):
            color = color_map.get(vpn, "gray")
            marker = "o" if neuron_name == "DNp01" else "s"
            ax.scatter(
                vpn_group["syn_count"],
                vpn_group["deflection_siz"],
                color=color,
                marker=marker,
                s=45,
                alpha=0.75,
                edgecolors="k",
                linewidths=0.4,
                label=f"{vpn} ({neuron_name})",
                zorder=3,
            )

            # plot the fitted exponential-rise curve for this population
            fit = fit_curves.get((neuron_name, vpn))
            if fit is not None:
                A_f, tau_f, x0_f, y0_f, xmin_f, xmax_f = fit
                x_curve = np.linspace(xmin_f, xmax_f, 200)
                C_f = y0_f - A_f * (1 - np.exp(-x0_f / tau_f))
                y_curve = A_f * (1 - np.exp(-x_curve / tau_f)) + C_f
                ax.plot(
                    x_curve, y_curve,
                    color=color, linewidth=2, alpha=0.85, zorder=2,
                )

    df_gap = nd_res_df[nd_res_df["gap"] == gap].copy()
    df_gap = df_gap[~((df_gap["neuron"] == "DNp03") & (df_gap["vpn"].isin(["LPLC2", "LC22"])))]

    for _, row in df_gap.iterrows():
        color = color_map.get(row["vpn"], "gray")
        marker = "o" if row["neuron"] == "DNp01" else "s"
        ax_all.scatter(
            row["density"],
            row["tau"],
            s=100,
            marker=marker,
            color=color,
            edgecolors="k",
            alpha=0.85,
            label=f"{row['vpn']} ({row['neuron']})",
        )

    valid = df_gap.dropna(subset=["density", "tau"])
    if not valid.empty and len(valid) >= 2:
        x_v = valid["density"].values
        y_v = valid["tau"].values
        res, r2 = print_linfit_stats(f"Neighborhood τ vs density (gap={gap})", x_v, y_v)
        x_line = np.linspace(x_v.min(), x_v.max(), 100)
        y_line = res.slope * x_line + res.intercept
        ax_all.plot(
            x_line, y_line,
            color="k", linestyle="--", linewidth=2,
            label="Linear fit (neighborhood)"
        )
        ax_all.text(
            0.05, 0.95,
            f"$R^2$ = {r2:.3f}\nslope = {res.slope:.3g}\np = {res.pvalue:.3g}",
            transform=ax_all.transAxes,
            fontsize=11,
            verticalalignment="top",
            bbox=dict(facecolor="white", alpha=0.7),
        )

    for ax, title in zip(
        [ax_dnp01, ax_dnp03, ax_all],
        ["DNp01 neighborhood VPNs", "DNp03 neighborhood VPNs", f"Neighborhood τ vs density (gap={gap})"]
    ):
        ax.set_title(title, fontsize=12)
        ax.grid(alpha=0.5, linestyle="--")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax_dnp01.set_xlabel("Neighborhood synapse count")
    ax_dnp01.set_ylabel("Depol at SIZ (mV)")
    ax_dnp03.set_xlabel("Neighborhood synapse count")
    ax_dnp03.set_ylabel("Depol at SIZ (mV)")
    ax_all.set_xlabel("Density (synapses / µm²)")
    ax_all.set_ylabel("Shunting τ (ms)")

    handles, labels = [], []
    for ax in [ax_dnp01, ax_dnp03, ax_all]:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    by_label = dict(zip(labels, handles))
    ax_all.legend(by_label.values(), by_label.keys(),
                  bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

    plt.tight_layout()
    plt.show()

    # ── FIGURE 2 ─────────────────────────────────────────────────────────────
    unique_gaps = sorted(nd_res_df["gap"].dropna().unique())
    n_gaps = len(unique_gaps)
    n_cols = 2
    n_rows = math.ceil(n_gaps / n_cols) if n_gaps > 0 else 1

    fig2, axs2 = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows))
    axs2 = np.atleast_1d(axs2).flatten()

    for i, gap_val in enumerate(unique_gaps):
        ax = axs2[i]
        df_g = nd_res_df[nd_res_df["gap"] == gap_val].copy()
        df_g = df_g[~((df_g["neuron"] == "DNp03") & (df_g["vpn"].isin(["LPLC2", "LC22"])))]

        for _, row in df_g.iterrows():
            color = color_map.get(row["vpn"], "gray")
            marker = "o" if row["neuron"] == "DNp01" else "s"
            ax.scatter(
                row["density"], row["tau"],
                s=100, marker=marker,
                color=color, edgecolors="k", alpha=0.85,
                label=f"{row['vpn']} ({row['neuron']})"
            )

        valid_g = df_g.dropna(subset=["density", "tau"])
        if not valid_g.empty and len(valid_g) >= 2:
            x_g = valid_g["density"].values
            y_g = valid_g["tau"].values
            res_g, r2_g = print_linfit_stats(f"Neighborhood τ vs density (gap={gap_val})", x_g, y_g)
            x_line = np.linspace(x_g.min(), x_g.max(), 100)
            y_line = res_g.slope * x_line + res_g.intercept
            ax.plot(
                x_line, y_line,
                color="k", linestyle="--", linewidth=2,
                label="Linear fit (neighborhood)"
            )
            ax.text(
                0.05, 0.95,
                f"$R^2$ = {r2_g:.3f}\nslope = {res_g.slope:.3g}\np = {res_g.pvalue:.3g}",
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment="top",
                bbox=dict(facecolor="white", alpha=0.7),
            )

        ax.set_xlabel("Density (synapses / µm²)")
        ax.set_ylabel("Shunting τ (ms)")
        ax.set_title(f"Gap = {gap_val}")
        ax.grid(alpha=0.5, linestyle="--")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for j in range(n_gaps, len(axs2)):
        axs2[j].set_visible(False)

    if n_gaps > 0:
        handles, labels = [], []
        for ax in axs2[:n_gaps]:
            h, l = ax.get_legend_handles_labels()
            handles.extend(h)
            labels.extend(l)
        by_label = dict(zip(labels, handles))
        axs2[min(n_gaps - 1, len(axs2) - 1)].legend(
            by_label.values(), by_label.keys(),
            bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8
        )

    plt.tight_layout()
    plt.show()

    return nd_res_df


#Figure 10
#Must run this first before running the below, and run partner and random seperately.
def plot_syn_spread_vs_dist_to_siz_partner_vs_rand(neuron_name, sizSection=None, plots = None):
    if plots == "Rand":
        dir_path = 'datafiles/simulationData/'+neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+'_rand'
        quintList = iterQuintsV6(dir_path)

    elif plots =="Partner":
        dir_path = 'datafiles/simulationData/'+neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+'_partner'
        quintList = iterQuintsV6(dir_path)
        print(quintList)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))
    synct_list = []
    peakV_siz_list= []
    peakV_soma_list =[]
    avg_dist_list = []
    avg_syn_spread_list = []
    synTypeList = []
    synpreList = []
    k = 0

    for k, quint in enumerate(quintList, start=1):
            print(f"Processing group {k}", flush=True)
            time_df, vSIZ_np, vSoma_np, synSite = quint
            idx_30 = np.argmin(np.abs(time_df - 30))
            # resting potential at 30 ms
            Vrest_siz = vSIZ_np[0, idx_30]
            peakV_siz = np.max(vSIZ_np[0, idx_30:])
            k+=1
            print(k)
            synSegRangeVarList = []
            synSecList = []
            peakV_soma = np.max(vSoma_np[0, idx_30:])
            Vrest_soma = vSoma_np[0, idx_30]

            synct = synSite.shape[0]
            synct_list.append(synct)

            distances_to_SIZ = []
            for _, row in synSite.iterrows():
                synSec = str2sec(row['mappedSection'])
                distToSIZ = h.distance(sizSection(0.05), synSec(row['mappedSegRangeVar']))
                distances_to_SIZ.append(distToSIZ)
                 # Append synapse information to lists
                synSecList.append(row.loc["mappedSection"])
                synTypeList.append(row.loc["type"])
                synpreList.append(row.loc["pre"])
                synSegRangeVarList.append(row.loc["mappedSegRangeVar"])

            avg_dist = numpy.mean(distances_to_SIZ)
            avg_dist_list.append(avg_dist)
            peakV_siz_list.append(peakV_siz- Vrest_siz)
            peakV_soma_list.append(peakV_soma- Vrest_soma)

            # Calculate distances for all unique synapse pairs
            synapse_indices = range(len(synSecList))
            syn_pairs = combinations(synapse_indices, 2)  # Get all unique pairs of synapses
            distances = []

            for i, j in syn_pairs:
                # Retrieve the sections and segments for each synapse
                synsec_1 = str2sec(synSecList[i])
                synsec_2 = str2sec(synSecList[j])

                # Calculate the distance between synapses
                dist = h.distance(synsec_1(synSegRangeVarList[i]), synsec_2(synSegRangeVarList[j]))
                distances.append(dist)
            avg_syn_spread = numpy.mean(distances)
            avg_syn_spread_list.append(avg_syn_spread)
    
            
    df = pd.DataFrame({
        'Synct': synct_list,
        'Peak_Depol_SIZ': peakV_siz_list,
        'Peak_Depol_SOMA': peakV_soma_list,
        'avg_dist_to_SIZ': avg_dist_list,
        'Avg_syn_spread': avg_syn_spread_list
    })

    csv_filename = f'{neuron_name}_{plots}_syn_spread_vs_dist_to_siz_data.csv'
    dir_path = f'datafiles/simulationData/{neuron_name}_final_sims'
    full_path = os.path.join(dir_path, csv_filename)
    df.to_csv(full_path, index=False)

    ax1.plot(df['Avg_syn_spread'], df['Peak_Depol_SIZ'], 'x', color='black')
    ax1.set_xlabel("Average synapse spread")
    ax1.set_ylabel("Peak depolarization at SIZ")

    ax1.legend()

    ax2.plot(df['avg_dist_to_SIZ'], df['Avg_syn_spread'], 'x', color='red')
    ax2.set_xlabel("Average distance of synapes to SIZ")
    ax2.set_ylabel("Average synapse spread")
    ax2.legend()

    ax3.plot(df['avg_dist_to_SIZ'], df['Peak_Depol_SIZ'], 'x', color='blue')
    ax3.set_xlabel("Average distance of synapses to SIZ")
    ax3.set_ylabel("Peak depolarization at SIZ")
    ax3.legend()

    plt.show()

#Must run this first before running the below, and run partner and random seperately.
def plot_syn_spread_vs_dist_to_siz(neuron_name, erev, MODE='all/', sizSection=None, plots = None,synnum = None,trials = None):
    synnum = synnum
    trials = trials
    if plots == "Rand":
        dir_path_rand = 'datafiles/simulationData/'+neuron_name+f'_rand{synnum}record{trials}'
        quintList = iterQuintsV6(dir_path_rand)
    elif plots =="Close":
        dir_path_close = 'datafiles/simulationData/'+neuron_name+f'_close{synnum}record{trials}'
        quintList = iterQuintsV7(dir_path_close)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))
    synct_list = []
    peakV_siz_list= []
    avg_dist_list = []
    avg_syn_spread_list = []
    
    synTypeList = []
    synpreList = []
    k = 0

    for k, quint in enumerate(quintList, start=1):
            print(f"Processing group {k}", flush=True)
            time_df, vSIZ_np, vSoma_np, synSite = quint
            idx_30 = np.argmin(np.abs(time_df - 30))
            # resting potential at 30 ms
            Vrest_siz = vSIZ_np[0, idx_30]
            peakV_siz = np.max(vSIZ_np[0, idx_30:])

            synSegRangeVarList = []
            synSecList = []
            k+=1
            print(k)
            synct = synSite.shape[0]
            synct_list.append(synct)

            distances_to_SIZ = []
            for _, row in synSite.iterrows():
                synSec = str2sec(row['mappedSection'])
                distToSIZ = h.distance(sizSection(0.055), synSec(row['mappedSegRangeVar']))
                distances_to_SIZ.append(distToSIZ)
                 # Append synapse information to lists
                synSecList.append(row.loc["mappedSection"])
                synTypeList.append(row.loc["type"])
                synpreList.append(row.loc["pre"])
                synSegRangeVarList.append(row.loc["mappedSegRangeVar"])

            avg_dist = numpy.mean(distances_to_SIZ)
            avg_dist_list.append(avg_dist)
            peakV_siz_list.append(peakV_siz- Vrest_siz)

            # Calculate distances for all unique synapse pairs
            synapse_indices = range(len(synSecList))
            syn_pairs = combinations(synapse_indices, 2)  # Get all unique pairs of synapses
            distances = []

            for i, j in syn_pairs:
                # Retrieve the sections and segments for each synapse
                synsec_1 = str2sec(synSecList[i])
                synsec_2 = str2sec(synSecList[j])

                # Calculate the distance between synapses
                dist = h.distance(synsec_1(synSegRangeVarList[i]), synsec_2(synSegRangeVarList[j]))
                distances.append(dist)
            avg_syn_spread = numpy.mean(distances)
            avg_syn_spread_list.append(avg_syn_spread)
    
            
    df = pd.DataFrame({
        'Synct': synct_list,
        'Peak_Depol_SIZ': peakV_siz_list,
        'avg_dist_to_SIZ': avg_dist_list,
        'Avg_syn_spread': avg_syn_spread_list
    })

    csv_filename = f'{neuron_name}_{plots}_syn_spread_vs_dist_to_siz_data_{trials}_trials_{synnum}_synapses.csv'
    dir_path = 'datafiles/simulationData'
    full_path = os.path.join(dir_path, csv_filename)
    df.to_csv(full_path, index=False)

    ax1.plot(df['Avg_syn_spread'], df['Peak_Depol_SIZ'], 'o', color='black')
    ax1.set_xlabel("Average synapse spread")
    ax1.set_ylabel("Peak depolarization at SIZ")

    ax1.legend()

    ax2.plot(df['avg_dist_to_SIZ'], df['Avg_syn_spread'], 'o', color='red')
    ax2.set_xlabel("Average distance of synapes to SIZ")
    ax2.set_ylabel("Average synapse spread")
    ax2.legend()

    ax3.plot(df['avg_dist_to_SIZ'], df['Peak_Depol_SIZ'], 'o', color='blue')
    ax3.set_xlabel("Average distance of synapses to SIZ")
    ax3.set_ylabel("Peak depolarization at SIZ")
    ax3.legend()

    # plt.savefig(f'{plots}_syn_spread_vs_dist_to_siz_vs_SIZ_depol_{trials}_trials_{synnum}_synapses.png')
    plt.show()

def plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_physical(neuron_name, sizSection=None, Dist=None, syncount=None):
    """
    Plots synapse spread vs distance to SIZ for a given neuron.
    Skips the folder if it does not exist.
    """
    Dist = str(Dist)
    syncount = str(syncount)
    dir_path = f'datafiles/simulationData/{neuron_name}_{syncount}_syns_rand_by_syn_spread{Dist}'

    # Skip if folder doesn't exist
    if not os.path.exists(dir_path):
        print(f"⚠️ Folder not found: {dir_path}. Skipping this pair.")
        return

    quintList = iterQuintsV3(dir_path)  # V3 iterator

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))
    
    synct_list = []
    peakV_siz_list = []
    peakV_soma_list = []
    avg_dist_list = []
    avg_syn_spread_list = []
    synTypeLabels = []
    spread_by_type = defaultdict(list)

    color_map = {
        'LC4': '#0000FF',
        'LPLC2': '#FFA500',
        'LPLC1': '#FF0D00',
        'LC22': '#FFF000',
        'LPLC4': '#00FF00',
        'LC6': '#FF00E9',
        "Mixed": "black"
    }

    for k, quint in enumerate(quintList):
        print(f"Processing group {k+1}", flush=True)
        time_df, vSIZ_np, vsoma_np, synSite = quint

        idx_30 = np.argmin(np.abs(time_df - 30))
        Vrest_siz = vSIZ_np[0, idx_30]
        peakV_siz = np.max(vSIZ_np[0, idx_30:])
        peakV_soma = np.max(vsoma_np[0, idx_30:])
        Vrest_soma = vsoma_np[0, idx_30]
        synSecList, synSegRangeVarList = [], []
        synct = synSite.shape[0]
        synct_list.append(synct)

        distances_to_SIZ = []
        syn_types = set()

        for _, row in synSite.iterrows():
            synSec = str2sec(row['mappedSection'])
            distToSIZ = h.distance(sizSection(0.05), synSec(row['mappedSegRangeVar']))
            distances_to_SIZ.append(distToSIZ)

            synSecList.append(row['mappedSection'])
            synSegRangeVarList.append(row['mappedSegRangeVar'])
            syn_types.add(row['type'])

        avg_dist = np.mean(distances_to_SIZ)
        avg_dist_list.append(avg_dist)
        peakV_siz_list.append(peakV_siz - Vrest_siz)
        peakV_soma_list.append(peakV_soma - Vrest_soma)

        syn_type = syn_types.pop() if len(syn_types) == 1 else "Mixed"
        synTypeLabels.append(syn_type)

        # Compute synapse spread
        distances = []
        for i, j in combinations(range(len(synSecList)), 2):
            sec1, sec2 = str2sec(synSecList[i]), str2sec(synSecList[j])
            dist = h.distance(sec1(synSegRangeVarList[i]), sec2(synSegRangeVarList[j]))
            distances.append(dist)

        avg_syn_spread = np.mean(distances)
        avg_syn_spread_list.append(avg_syn_spread)
        spread_by_type[syn_type].append(avg_syn_spread)

        color = color_map.get(syn_type, "black")
        ax1.plot(avg_syn_spread, peakV_siz - Vrest_siz, 'x', color=color)
        ax2.plot(avg_dist, avg_syn_spread, 'x', color=color)
        ax3.plot(avg_dist, peakV_siz - Vrest_siz, 'x', color=color)

    # DataFrame
    df = pd.DataFrame({
        'Synct': synct_list,
        'Peak_Depol_SIZ': peakV_siz_list,
        'Peak_Depol_SOMA': peakV_soma_list,
        'avg_dist_to_SIZ': avg_dist_list,
        'Avg_syn_spread': avg_syn_spread_list,
        'Type': synTypeLabels
    })

    # Plot formatting
    ax1.set_xlabel("Average synapse spread"); ax1.set_ylabel("Peak depolarization at SIZ")
    ax1.set_title("Peak SIZ depol vs Synapse Spread")

    ax2.set_xlabel("Average distance to SIZ"); ax2.set_ylabel("Average synapse spread")
    ax2.set_title("Spread vs Distance to SIZ")

    ax3.set_xlabel("Average distance to SIZ"); ax3.set_ylabel("Peak depolarization at SIZ")
    ax3.set_title("Depol vs Distance to SIZ")

    plt.tight_layout()
    # plt.show()

    # Print average spread per type
    print("\nAverage Synapse Spread per Synapse Type:")
    for syn_type, spreads in spread_by_type.items():
        print(f"{syn_type}: {np.mean(spreads):.2f}")

    # === Save DataFrame to CSV ===
    save_name = f"{neuron_name}_{syncount}_syn_spread{Dist}.csv"
    df.to_csv(save_name, index=False)
    print(f"✅ Data saved to {save_name}")

    return df

def plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=None, Dist=None, syncount=None, close_analysis = False, partner_analysis = False, trials = None):
    """
    Extended FINAL version:
    - Original synapse spread vs SIZ plots (3 subplots)
    - Additional electrotonic distance analysis with new 3 subplots
    """
    
    if close_analysis == True:
        syncount = int(syncount)
        dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/CloseSims/{neuron_name}_close{syncount}record{trials}'
        quintList = iterQuintsV9(dir_path) 
    elif partner_analysis == True:
        dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner'
        quintList = iterQuintsV10(dir_path) 
    else:
        Dist = str(Dist)
        syncount = int(syncount)
        dir_path = f'datafiles/simulationData/{neuron_name}_final_sims/synapse_spread/{neuron_name}_{syncount}_syns_rand_by_syn_spread{Dist}'
        quintList = iterQuintsV8(dir_path) 

    if not os.path.exists(dir_path):
        print(f"⚠️ Folder not found: {dir_path}. Skipping this pair.")
        return

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))

    synct_list, peakV_siz_list, peakV_soma_list = [], [], []
    avg_dist_list, avg_syn_spread_list = [], []
    avg_elec_syns_list, avg_elec_SIZ_list = [], []
    synTypeLabels = []
    spread_by_type = defaultdict(list)

    color_map = {
        'LC4': '#0000FF',
        'LPLC2': '#FFA500',
        'LPLC1': '#FF0D00',
        'LC22': '#FFF000',
        'LPLC4': '#00FF00',
        'LC6': '#FF00E9',
        "Mixed": "black"
    }

    for k, quint in enumerate(quintList):
        print(f"Processing group {k+1}", flush=True)
        time_df, vSIZ_np, vsoma_np, synSite, stim_summary_df = quint

        # --- Existing spread & depol calculations ---
        idx_30 = np.argmin(np.abs(time_df - 30))
        Vrest_siz = vSIZ_np[0, idx_30]
        peakV_siz = np.max(vSIZ_np[0, idx_30:])
        peakV_soma = np.max(vsoma_np[0, idx_30:])
        Vrest_soma = vsoma_np[0, idx_30]
        synct = synSite.shape[0]
        synct_list.append(synct)

        distances_to_SIZ, syn_types = [], set()
        synSecList, synSegRangeVarList = [], []

        for _, row in synSite.iterrows():
            synSec = str2sec(row['mappedSection'])
            distToSIZ = h.distance(sizSection(0.05), synSec(row['mappedSegRangeVar']))
            distances_to_SIZ.append(distToSIZ)
            synSecList.append(row['mappedSection'])
            synSegRangeVarList.append(row['mappedSegRangeVar'])
            syn_types.add(row['type'])

        avg_dist = np.mean(distances_to_SIZ)
        avg_dist_list.append(avg_dist)
        peakV_siz_list.append(peakV_siz - Vrest_siz)
        peakV_soma_list.append(peakV_soma - Vrest_soma)

        syn_type = syn_types.pop() if len(syn_types) == 1 else "Mixed"
        synTypeLabels.append(syn_type)

        distances = []
        for i, j in combinations(range(len(synSecList)), 2):
            sec1, sec2 = str2sec(synSecList[i]), str2sec(synSecList[j])
            dist = h.distance(sec1(synSegRangeVarList[i]), sec2(synSegRangeVarList[j]))
            distances.append(dist)

        avg_syn_spread = np.mean(distances)
        avg_syn_spread_list.append(avg_syn_spread)
        spread_by_type[syn_type].append(avg_syn_spread)

        color = color_map.get(syn_type, "black")
        ax1.plot(avg_syn_spread, peakV_siz - Vrest_siz, 'x', color=color)
        ax2.plot(avg_dist, avg_syn_spread, 'x', color=color)
        ax3.plot(avg_dist, peakV_siz - Vrest_siz, 'x', color=color)

        # --- New electrotonic distance analysis ---
        ratios_syns, ratios_SIZ = [], []

        # find the numeric index for SIZ
        siz_index = synct +1  # or len(syn_df), whatever you used when assigning SIZ

        # iterate over unique stim_syn_index
        for stim_idx in stim_summary_df["stim_syn_index"].unique():
            group = stim_summary_df[stim_summary_df["stim_syn_index"] == stim_idx]

            # row where stim == record (same synapse)
            same_row = group[group["stim_syn_index"] == group["record_syn_index"]]
            if same_row.empty:
                continue
            stim_peak = same_row["peak_amp"].values[0]
            if stim_peak == 0:
                continue  # skip to avoid log(0)

            # log(stim / other synapses)
            syn_only = group[group["record_syn_index"] != siz_index]  # exclude SIZ
            for _, row in syn_only.iterrows():
                if row["peak_amp"] != 0:
                    ratios_syns.append(np.log(stim_peak / row["peak_amp"]))

            # log(stim / SIZ
            siz_row = group[group["record_syn_index"] == siz_index]
            if not siz_row.empty and siz_row["peak_amp"].values[0] != 0:
                ratios_SIZ.append(np.log(stim_peak / siz_row["peak_amp"].values[0]))

        # average the ln ratios
        avg_elec_syns = np.mean(ratios_syns) if ratios_syns else np.nan
        avg_elec_SIZ = np.mean(ratios_SIZ) if ratios_SIZ else np.nan
        avg_elec_syns_list.append(avg_elec_syns)
        avg_elec_SIZ_list.append(avg_elec_SIZ)
    # --- DataFrame ---
    df = pd.DataFrame({
        'Synct': synct_list,
        'Peak_Depol_SIZ': peakV_siz_list,
        'Peak_Depol_SOMA': peakV_soma_list,
        'avg_dist_to_SIZ': avg_dist_list,
        'Avg_syn_spread': avg_syn_spread_list,
        'avg_electrotonic_dist_syns': avg_elec_syns_list,
        'avg_electrotonic_dist_SIZ': avg_elec_SIZ_list,
        'Type': synTypeLabels
    })

    # --- First figure (spread-based) ---
    ax1.set_xlabel("Average synapse spread"); ax1.set_ylabel("Peak depolarization at SIZ")
    ax1.set_title("Peak SIZ depol vs Synapse Spread")

    ax2.set_xlabel("Average distance to SIZ"); ax2.set_ylabel("Average synapse spread")
    ax2.set_title("Spread vs Distance to SIZ")

    ax3.set_xlabel("Average distance to SIZ"); ax3.set_ylabel("Peak depolarization at SIZ")
    ax3.set_title("Depol vs Distance to SIZ")

    plt.tight_layout()

    # --- Second figure (electrotonic-based) ---
    fig2, (bx1, bx2, bx3) = plt.subplots(3, 1, figsize=(8, 15))

    for syn_type, color in color_map.items():
        mask = df["Type"] == syn_type
        bx1.plot(df.loc[mask, "avg_electrotonic_dist_syns"], df.loc[mask, "Peak_Depol_SIZ"], 'x', color=color)
        bx2.plot(df.loc[mask, "avg_electrotonic_dist_SIZ"], df.loc[mask, "avg_electrotonic_dist_syns"], 'x', color=color)
        bx3.plot(df.loc[mask, "avg_electrotonic_dist_SIZ"], df.loc[mask, "Peak_Depol_SIZ"], 'x', color=color)

    bx1.set_xlabel("Avg Electrotonic Distance to Synapses"); bx1.set_ylabel("Peak depol at SIZ")
    bx1.set_title("Peak SIZ depol vs Electrotonic Spread")

    bx2.set_xlabel("Avg Electrotonic Distance to SIZ"); bx2.set_ylabel("Electrotonic distance to synapses")
    bx2.set_title("Electrotonic Spread vs Dist to SIZ")

    bx3.set_xlabel("Avg Electrotonic Distance to SIZ"); bx3.set_ylabel("Peak depol at SIZ")
    bx3.set_title("Depol vs Electrotonic Dist to SIZ")

    plt.tight_layout()
    plt.show()

    # --- Save ---
    if close_analysis == True:
        save_name = f"datafiles/simulationData/{neuron_name}_final_sims/CloseSims/{neuron_name}_close{syncount}_syn_spread_final.csv"  
    elif partner_analysis == True:
        save_name = f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_partner_syn_spread_final.csv" 
    else:
        save_name = f"datafiles/simulationData/{neuron_name}_final_sims/synapse_spread/{neuron_name}_{syncount}_syn_spread{Dist}_final.csv"
    df.to_csv(save_name, index=False)
    print(f"✅ Data saved to {save_name}")

    return df

def plot_partner_close_all_electro_bins(neuron_name, n_e_bins=7, export_csv=True, csv_path=None):

    type_colors = {
        'LC4': '#0000FF',
        'LPLC2': '#FFA500',
        'LPLC1': '#FF0D00',
        'LC22': '#FFF000',
        'LPLC4': '#00FF00',
        'LC6': '#FF00E9'
    }

    def try_read_csv(candidates):
        for p in candidates:
            if os.path.exists(p):
                try:
                    return pd.read_csv(p)
                except Exception as e:
                    print(f"⚠️ Could not read {p}: {e}")
        return pd.DataFrame()

    # --- Define syn counts ---
    if neuron_name == "DNp01":
        all_syncounts = [1,2,3,4,5,6,7,8,9,11,13,22,25,29,31]
    elif neuron_name == "DNp03":
        all_syncounts = [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,29,32]
    elif neuron_name == "DNp01_hemi":
        all_syncounts = [5,11,12,15,16,18,19,20,21,22,23,25,26,27,28,29,30,31,32,35,36,39,47]
    elif neuron_name == "DNp03_hemi":
        all_syncounts = [6,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,27,30,38]
    else:
        raise ValueError("Unknown neuron_name")

    # --- Load Random + Close data ---
    def load_for_synct(sc):
        parts = []
        for Dist in range(10, 71, 10):
            fname = f"{neuron_name}_{sc}_syn_spread{Dist}_final.csv"
            df = try_read_csv([
                os.path.join("datafiles", "simulationData", f"{neuron_name}_final_sims", "synapse_spread", fname),
                fname
            ])
            if not df.empty:
                df = df.copy()
                df["Dist"] = Dist
                df["Group"] = "All"
                df["Synct"] = sc
                parts.append(df)

        close_fname = f"{neuron_name}_close{sc}_syn_spread_final.csv"
        close_df = try_read_csv([
            os.path.join("datafiles", "simulationData", f"{neuron_name}_final_sims", "CloseSims", close_fname),
            close_fname
        ])
        if not close_df.empty:
            close_df = close_df.copy()
            close_df["Group"] = "Close"
            close_df["Synct"] = sc
            parts.append(close_df)

        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    dfs = [load_for_synct(sc) for sc in all_syncounts]
    dfs = [df for df in dfs if not df.empty]
    if not dfs:
        print("❌ No data found.")
        return None

    # --- Load VPNs ---
    partner_fname = f"{neuron_name}_partner_syn_spread_final.csv"
    partner_df = try_read_csv([
        os.path.join("datafiles", "simulationData", f"{neuron_name}_final_sims", partner_fname),
        partner_fname
    ])
    if not partner_df.empty:
        partner_df = partner_df.copy()
        partner_df["Group"] = "VPNs"

    all_data = pd.concat(dfs + ([partner_df] if not partner_df.empty else []), ignore_index=True)
    if "avg_electrotonic_dist_syns" not in all_data.columns or "Peak_Depol_SIZ" not in all_data.columns:
        print("❌ Missing required columns.")
        return all_data

    all_data["Peak_Depol_norm"] = all_data["Peak_Depol_SIZ"] / all_data["Synct"]

    # --- Electrotonic bins ---
    global_e_vals = all_data["avg_electrotonic_dist_syns"].dropna()
    e_bins = np.linspace(0, float(global_e_vals.max()), n_e_bins + 1)
    e_labels = [f"{e_bins[i]:.2f}-{e_bins[i+1]:.2f}" for i in range(len(e_bins)-1)]
    all_data.loc[:, "Electro_bin"] = pd.cut(
        all_data["avg_electrotonic_dist_syns"], bins=e_bins,
        labels=e_labels, include_lowest=True
    )

    all_data.loc[:, "Combined_Group"] = np.where(
        all_data["Group"].isin(["All","Close"]),
        "Random+Close",
        all_data["Type"]
    )

    
        # --- Figure setup ---
    fig = plt.figure(figsize=(14,5))  # wider overall
    gs = fig.add_gridspec(1,2, width_ratios=[3,1], wspace=0.35)  # left subplot wider

    # --- Left: boxplots + scatter ---
    ax1 = fig.add_subplot(gs[0,0])
    vpn_types = sorted(all_data[all_data["Group"]=="VPNs"]["Type"].dropna().unique())
    box_groups = ["Random+Close"] + vpn_types
    colors = {"Random+Close":"black"}
    colors.update(type_colors)

    total_bin_width = 1.2  # allocate slightly more width per bin
    n_groups = len(box_groups)
    box_width = total_bin_width / n_groups * 0.8  # shrink boxes to increase spacing between them
    bin_spacing = 2.5  # increase horizontal spacing between bins
    bin_positions = np.arange(len(e_labels)) * bin_spacing  # spread bins further apart

    for i, bin_label in enumerate(e_labels):
        for j, grp in enumerate(box_groups):
            if grp=="Random+Close":
                sub = all_data[(all_data["Electro_bin"]==bin_label) & (all_data["Group"].isin(["All","Close"]))] 
            else:
                sub = all_data[(all_data["Electro_bin"]==bin_label) & 
                            (all_data["Group"]=="VPNs") & (all_data["Type"]==grp)]
            if sub.empty: continue

            # position each box inside bin with extra spacing
            bin_center = bin_positions[i]
            pos = bin_center - total_bin_width/2 + j*(total_bin_width/n_groups) + (total_bin_width/n_groups)/2
            edgecolor = 'black' if grp=="Random+Close" else colors[grp]

            # boxplot
            ax1.boxplot(sub["Peak_Depol_norm"], positions=[pos], widths=box_width,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor='none', edgecolor=edgecolor),
                        medianprops=dict(color=edgecolor),
                        whiskerprops=dict(color=edgecolor),
                        capprops=dict(color=edgecolor))

            # scatter points
            q1, q3 = sub["Peak_Depol_norm"].quantile([0.25,0.75])
            iqr = q3 - q1
            non_outliers = sub[(sub["Peak_Depol_norm"] >= q1-1.5*iqr) & (sub["Peak_Depol_norm"] <= q3+1.5*iqr)]
            jitter = (np.random.rand(len(non_outliers)) - 0.5) * box_width * 1.2  # slightly wider jitter

            if grp=="Random+Close":
                sub_random = non_outliers[non_outliers["Group"]=="All"]
                sub_close = non_outliers[non_outliers["Group"]=="Close"]
                ax1.scatter(np.full(len(sub_random), pos)+jitter[:len(sub_random)],
                            sub_random["Peak_Depol_norm"], color="black", s=25, edgecolor=None, alpha=1.0)
                ax1.scatter(np.full(len(sub_close), pos)+jitter[:len(sub_close)],
                            sub_close["Peak_Depol_norm"], color="gray", s=25, edgecolor=None, alpha=1.0)
            else:
                ax1.scatter(np.full(len(non_outliers), pos)+jitter,
                            non_outliers["Peak_Depol_norm"], color=colors[grp], s=25, edgecolor=None, alpha=1.0)

    ax1.set_xticks(bin_positions)
    ax1.set_xticklabels(e_labels, rotation=30, ha="right")
    ax1.set_xlabel("Average electrotonic distance between synapses")
    ax1.set_ylabel("Depolization at SIZ per synapse")
    # --- Right: stacked histogram ---
    ax2 = fig.add_subplot(gs[0,1])
    vpn_data = all_data[(all_data["Group"]=="VPNs") & all_data["Type"].notna()].copy()
    total_vpns = len(vpn_data)
    if total_vpns > 0:
        vpn_data.loc[:, "Electro_bin"] = pd.cut(vpn_data["avg_electrotonic_dist_syns"], bins=e_bins,
                                                labels=e_labels, include_lowest=True)
        count_df = vpn_data.groupby(["Electro_bin","Type"], observed=False).size().unstack(fill_value=0)
        count_df = 100*count_df/total_vpns
        bottom = np.zeros(len(count_df))
        for t in vpn_types:
            vals = count_df[t] if t in count_df.columns else np.zeros(len(count_df))
            ax2.bar(count_df.index, vals, bottom=bottom, color=type_colors.get(t,"gray"), 
                    label=t, edgecolor='black', alpha=0.9)
            bottom += vals
        ax2.set_xticks(range(len(e_labels)))
        ax2.set_xticklabels(e_labels, rotation=30, ha="right")
        ax2.set_xlabel("Electrotonic distance bins")
        ax2.set_ylabel("% of VPNs per electrotonic bin")
        ax2.set_title(f"{neuron_name} — VPN distribution across electrotonic bins")
        ax2.legend(frameon=False)

    plt.tight_layout()
    plt.show()

    # --- Levene all-vs-all ---
    valid_data = all_data[all_data["Peak_Depol_norm"].notna()].copy()
    group_bin_pairs = valid_data.groupby(["Combined_Group","Electro_bin"], observed=False).size().reset_index(name="count")
    group_bin_pairs = group_bin_pairs[group_bin_pairs["count"]>=3]

    results = []
    for (g1,b1),(g2,b2) in combinations(group_bin_pairs[["Combined_Group","Electro_bin"]].itertuples(index=False,name=None),2):
        vals1 = valid_data[(valid_data["Combined_Group"]==g1)&(valid_data["Electro_bin"]==b1)]["Peak_Depol_norm"]
        vals2 = valid_data[(valid_data["Combined_Group"]==g2)&(valid_data["Electro_bin"]==b2)]["Peak_Depol_norm"]
        if len(vals1)<3 or len(vals2)<3: continue
        stat,p = levene(vals1,vals2,center="median")
        results.append({"Group1":g1,"Bin1":b1,"n1":len(vals1),
                        "Group2":g2,"Bin2":b2,"n2":len(vals2),
                        "Levene_Stat":stat,"p_value":p})

    results_df = pd.DataFrame(results)
    if export_csv:
        if csv_path is None:
            csv_path = f"{neuron_name}_Levene_results.csv"
        results_df.to_csv(csv_path, index=False)
        print(f"✅ Levene results exported to {csv_path}")

    return all_data, results_df

#final published version of the plot_partner_close_all_combined function
def plot_partner_close_all_combined(
    neuron_name, export_csv=True, csv_path=None
):
    import os
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from itertools import combinations
    from scipy.stats import levene

    # --- Color map for VPNs ---
    type_colors = {
        "LC4": "#0000FF",
        "LPLC2": "#FFA500",
        "LPLC1": "#FF0D00",
        "LC22": "#FFF000",
        "LPLC4": "#00FF00",
        "LC6": "#FF00E9",
    }

    # --- Synapse count lookup ---
    syn_counts = {
        "DNp01": [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 22, 25, 29, 31
        ],
        "DNp03": [
            1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14,
            15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 29, 32
        ],
        "DNp01_hemi": [
            5, 11, 12, 15, 16, 18, 19, 20, 21, 22, 23, 25,
            26, 27, 28, 29, 30, 31, 32, 35, 36, 39, 47
        ],
        "DNp03_hemi": [
            6, 8, 9, 10, 11, 12, 13, 14, 15, 16,
            17, 18, 20, 21, 22, 23, 27, 30, 38
        ],
    }

    if neuron_name not in syn_counts:
        raise ValueError(f"Unknown neuron_name: {neuron_name}")

    all_syncounts = syn_counts[neuron_name]

    def try_read_csv(candidates):
        for p in candidates:
            if os.path.exists(p):
                try:
                    return pd.read_csv(p)
                except Exception as e:
                    print(f"Could not read {p}: {e}")
        return pd.DataFrame()

    # --- Load Random and Close data: unchanged from V5 ---
    def load_for_synct(sc):
        parts = []

        for Dist in range(10, 71, 10):
            fname = f"{neuron_name}_{sc}_syn_spread{Dist}_final.csv"
            df = try_read_csv([
                os.path.join(
                    "datafiles", "simulationData",
                    f"{neuron_name}_final_sims",
                    "synapse_spread", fname
                ),
                fname,
            ])

            if not df.empty:
                df["Group"] = "Random"
                df["Synct"] = sc
                parts.append(df)

        close_fname = f"{neuron_name}_close{sc}_syn_spread_final.csv"
        close_df = try_read_csv([
            os.path.join(
                "datafiles", "simulationData",
                f"{neuron_name}_final_sims",
                "CloseSims", close_fname
            ),
            close_fname,
        ])

        if not close_df.empty:
            close_df["Group"] = "Close"
            close_df["Synct"] = sc
            parts.append(close_df)

        return (
            pd.concat(parts, ignore_index=True)
            if parts else pd.DataFrame()
        )

    dfs = [load_for_synct(sc) for sc in all_syncounts]
    dfs = [df for df in dfs if not df.empty]

    # --- Load actual VPN data using the working V6 method ---
    partner_rows = []

    partner_dir = (
        f"datafiles/simulationData/{neuron_name}_final_sims/"
        f"rand_vs_partner/{neuron_name}_partner"
    )
    quintList = iterQuintsV4(partner_dir)

    for k, quint in enumerate(quintList):
        print(f"Processing group {k + 1} (_partner)")
        time_df, vSIZ_np, vsoma_np, synSite = quint

        if synSite.empty:
            continue

        idx_30 = np.argmin(np.abs(time_df - 30))
        peakV_siz = np.max(vSIZ_np[0, idx_30:])
        Vrest_siz = vSIZ_np[0, idx_30]

        partner_rows.append({
            "Group": synSite["type"].iloc[0],
            "Synct": synSite.shape[0],
            "Peak_Depol_SIZ": peakV_siz - Vrest_siz,
        })

    partner_df = pd.DataFrame(partner_rows)

    if not partner_df.empty:
        dfs.append(partner_df)

    if not dfs:
        print("No data found.")
        return None

    all_data = pd.concat(dfs, ignore_index=True)

    if (
        "Peak_Depol_SIZ" not in all_data.columns
        or "Synct" not in all_data.columns
    ):
        raise ValueError("Missing Peak_Depol_SIZ or Synct columns.")

    # --- Normalize ---
    all_data["Peak_Depol_norm"] = (
        all_data["Peak_Depol_SIZ"] / all_data["Synct"]
    )

    # --- Define groups ---
    base_groups = ["Close", "Random"]
    vpn_groups = [
        g for g in sorted(all_data["Group"].dropna().unique())
        if g not in base_groups
    ]
    groups = base_groups + vpn_groups

    valid = all_data[
        all_data["Peak_Depol_norm"]
        .replace([np.inf, -np.inf], np.nan)
        .notna()
    ].copy()

    # --- Boxplot of normalized depolarization ---
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, grp in enumerate(groups):
        sub = valid[valid["Group"] == grp]
        if sub.empty:
            continue

        color = type_colors.get(
            grp, "gray" if grp == "Close" else "black"
        )

        ax.boxplot(
            sub["Peak_Depol_norm"],
            positions=[i],
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            boxprops=dict(facecolor="none", edgecolor="black"),
            medianprops=dict(color="black"),
            whiskerprops=dict(color="black"),
            capprops=dict(color="black"),
        )

        jitter = (np.random.rand(len(sub)) - 0.5) * 0.3
        ax.scatter(
            np.full(len(sub), i) + jitter,
            sub["Peak_Depol_norm"],
            s=30,
            color=color,
            alpha=0.9,
        )

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.set_ylabel("Depolarization per synapse (norm.)")
    ax.set_title(f"{neuron_name}: Close, Random, and Individual VPNs")
    plt.tight_layout()
    plt.show()

    # --- FIGURE 2: Peak depolarization vs synapse number ---
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    fig.suptitle(
        f"{neuron_name}: Peak Depolarization vs Synapse Number",
        fontsize=14,
        fontweight="bold",
    )

    rand_df = all_data[all_data["Group"] == "Random"]
    close_df = all_data[all_data["Group"] == "Close"]
    vpn_df = all_data[all_data["Group"].isin(vpn_groups)]

    # Random
    if not rand_df.empty:
        axs[0].scatter(
            rand_df["Synct"], rand_df["Peak_Depol_SIZ"],
            color="black", alpha=0.6
        )
    axs[0].set_title("Random")
    axs[0].set_xlabel("Synapse Number")
    axs[0].set_ylabel("Peak Depolarization (mV)")

    # Close
    if not close_df.empty:
        axs[1].scatter(
            close_df["Synct"], close_df["Peak_Depol_SIZ"],
            color="gray", alpha=0.6
        )
    axs[1].set_title("Close")
    axs[1].set_xlabel("Synapse Number")

    # VPNs
    if not vpn_df.empty:
        for vpn_type, sub in vpn_df.groupby("Group"):
            axs[2].scatter(
                sub["Synct"], sub["Peak_Depol_SIZ"],
                color=type_colors.get(vpn_type, "black"),
                label=vpn_type, alpha=0.8
            )
        axs[2].legend(title="VPN Type", fontsize=8)

    axs[2].set_title("Actual VPNs")
    axs[2].set_xlabel("Synapse Number")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    # --- FIGURE 3: Overlay all groups ---
    fig, ax = plt.subplots(figsize=(8, 6))

    if not rand_df.empty:
        ax.scatter(
            rand_df["Synct"], rand_df["Peak_Depol_SIZ"],
            color="black", alpha=0.6, label="Random"
        )

    if not close_df.empty:
        ax.scatter(
            close_df["Synct"], close_df["Peak_Depol_SIZ"],
            color="gray", alpha=0.6, label="Close"
        )

    if not vpn_df.empty:
        for vpn_type, sub in vpn_df.groupby("Group"):
            ax.scatter(
                sub["Synct"], sub["Peak_Depol_SIZ"],
                color=type_colors.get(vpn_type, "black"),
                label=vpn_type, alpha=0.8
            )

    ax.set_xlabel("Synapse Number")
    ax.set_ylabel("Peak Depolarization (mV)")
    ax.set_title(f"{neuron_name}: Peak Depol Overlay")
    ax.legend(title="Group/VPN Type", fontsize=8)
    plt.tight_layout()
    plt.show()

    # --- Pairwise median-centered Levene tests ---
    levene_results = []

    for g1, g2 in combinations(groups, 2):
        vals1 = valid.loc[
            valid["Group"] == g1, "Peak_Depol_norm"
        ].astype(float)

        vals2 = valid.loc[
            valid["Group"] == g2, "Peak_Depol_norm"
        ].astype(float)

        n1 = len(vals1)
        n2 = len(vals2)

        if n1 >= 3 and n2 >= 3:
            levene_stat, levene_p = levene(
                vals1, vals2, center="median"
            )
            tested = True
            exclusion_reason = ""
        else:
            levene_stat = np.nan
            levene_p = np.nan
            tested = False
            exclusion_reason = (
                "At least one group had fewer than 3 valid observations"
            )

        levene_results.append({
            "Group1": g1,
            "Group2": g2,
            "N_Group1": n1,
            "N_Group2": n2,
            "Levene_Center": "median",
            "Levene_Stat": levene_stat,
            "Levene_p": levene_p,
            "Significant_p_lt_0.05": (
                levene_p < 0.05 if pd.notna(levene_p) else False
            ),
            "Multiple_Test_Correction": "None",
            "Tested": tested,
            "Exclusion_Reason": exclusion_reason,
        })

    levene_df = pd.DataFrame(levene_results)

    # --- Export Levene results ---
    if export_csv:
        if csv_path is None:
            csv_path = f"{neuron_name}_Levene_results.csv"

        csv_path = str(csv_path)
        if not csv_path.lower().endswith(".csv"):
            csv_path += ".csv"

        csv_path = os.path.abspath(csv_path)
        levene_df.to_csv(csv_path, index=False)

        print(f"\nExported Levene results to:\n{csv_path}")

    return all_data, levene_df

def print_significant_levene(csv_file):
    """
    Reads a CSV of Levene test results and prints:
      1. VPN vs Random/Close within the same bin (only significant)
      2. VPN vs VPN comparisons within the same bin (including ns)

    VPNs: ["LC4", "LC22", "LC6", "LPLC1", "LPLC2", "LPLC4"]

    Star notation:
        *    p < 0.05
        **   p < 0.01
        ***  p < 0.001
        **** p < 0.0001
        ns   p >= 0.05

    Expected CSV columns: ['Group1', 'Bin1', 'Group2', 'Bin2', 'Levene_Stat', 'p_value']
    """
    import pandas as pd

    df = pd.read_csv(csv_file)

    required_cols = {'Group1','Bin1','Group2','Bin2','p_value'}
    if not required_cols.issubset(df.columns):
        print("❌ CSV missing required columns")
        return

    VPNs = ["LC4", "LC22", "LC6", "LPLC1", "LPLC2", "LPLC4"]

    # Convert p-values to stars (including ns)
    def p_to_stars(p):
        if p < 0.0001:
            return "****"
        elif p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    df['stars'] = df['p_value'].apply(p_to_stars)

    # --- VPN vs Random/Close within the same bin (significant only) ---
    vpn_vs_rc = df[
        (
            (df['Group1'].isin(VPNs) & (df['Group2'] == 'Random+Close')) |
            (df['Group2'].isin(VPNs) & (df['Group1'] == 'Random+Close'))
        ) &
        (df['Bin1'] == df['Bin2'])
    ]
    vpn_vs_rc_sig = vpn_vs_rc[vpn_vs_rc['stars'] != "ns"]

    if not vpn_vs_rc_sig.empty:
        print("\n=== VPN vs Random/Close within same bin (significant) ===")
        for _, row in vpn_vs_rc_sig.iterrows():
            print(f"{row['Group1']} ({row['Bin1']}) vs {row['Group2']} ({row['Bin2']}): "
                  f"{row['stars']} (p={row['p_value']:.4g})")
    else:
        print("\nNo significant VPN vs Random/Close comparisons found within the same bin.")

    # --- VPN vs VPN within same bin (include ns) ---
    vpn_vs_vpn = df[
        (df['Group1'].isin(VPNs)) &
        (df['Group2'].isin(VPNs)) &
        (df['Bin1'] == df['Bin2'])
    ]

    if not vpn_vs_vpn.empty:
        print("\n=== VPN vs VPN within same bin ===")
        for _, row in vpn_vs_vpn.iterrows():
            print(f"{row['Group1']} ({row['Bin1']}) vs {row['Group2']} ({row['Bin2']}): "
                  f"{row['stars']} (p={row['p_value']:.4g})")
    else:
        print("\nNo VPN vs VPN comparisons found within the same bin.")

######
#Misc at the end of the file
#Analysis of depolarizations and integration of multiple synapses comparing actual VPN activations with the individual linear sum of its synapses at SIZ

def Peak_depol_per_VPN(neuron_name):
    """
    For a given neuron, iterate over partner quint files and compute:
    - peak depol (vSIZ max after rest time minus baseline Vrest)
    - average depol per synapse
    Returns a sorted DataFrame.
    """

    dir_path = f"datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_partner"

    if not os.path.isdir(dir_path):
        print(f"⚠️ Directory not found: {dir_path}")
        return pd.DataFrame()

    data_list = []

    for quint in iterQuintsV10(dir_path):
        time_df, vSIZ_np, vsoma_np, synSite_df, stim_summary_df = quint

        # --- safety checks ---
        if vSIZ_np.size == 0:
            continue
        if not isinstance(synSite_df, pd.DataFrame) or synSite_df.empty:
            continue

        # Find index where time is closest to 48 ms
        try:
            time_array = time_df.to_numpy().flatten()
            rest_idx = (numpy.abs(time_array - 48)).argmin()
        except Exception as e:
            print(f"⚠️ Failed to find rest index: {e}")
            continue

        # Resting value at that time point
        try:
            Vrest_siz = vSIZ_np[0, rest_idx]
        except IndexError:
            print(f"⚠️ Rest index {rest_idx} out of range — skipping trial")
            continue

        # Peak depol after rest time
        try:
            peakV_siz = numpy.max(vSIZ_np[:, rest_idx:])
        except Exception as e:
            print(f"⚠️ Failed to calculate peak: {e}")
            continue

        synct = len(synSite_df)

        cell_type = synSite_df["type"].iloc[0] if "type" in synSite_df.columns else None
        if cell_type in ["LC4", "LC6", "LC22", "LPLC2", "LPLC1", "LPLC4"]:
            pre_id = synSite_df["pre"].iloc[0] if "pre" in synSite_df.columns else None
            depol = peakV_siz - Vrest_siz
            avg_syn_depol = depol / synct if synct > 0 else numpy.nan

            data_list.append({
                "pre": pre_id,
                "type": cell_type,
                "synct": synct,
                "depol_together": depol,
                "avg_syn_depol": avg_syn_depol
            })

    df = pd.DataFrame(data_list)

    if df.empty:
        print(f"⚠️ No valid data for {neuron_name}")
        return df

    # Sort by type and synapse count
    return df.sort_values(by=["type", "synct"]).reset_index(drop=True)

def individual_syn_depol_per_VPN(neuron_name):
    # Determine file paths dynamically
    base_path = f"datafiles/simulationData/{neuron_name}_final_sims/singlesynactivation/"
    if not os.path.isdir(base_path):
        print(f"⚠️ Directory not found: {base_path}")
        return pd.DataFrame()

    try:
        synSite_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_synapse_sites.csv"), dtype={'pre': str})
        synCurr_siz_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_synaptic_currents_siz.csv"), header=None)
        time_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_time.csv"), header=None)
    except Exception as e:
        print(f"⚠️ Failed to load files for {neuron_name}: {e}")
        return pd.DataFrame()

    # Find index where time is closest to 48 ms
    try:
        time_array = time_df.to_numpy().flatten()
        rest_idx = (np.abs(time_array - 48)).argmin()
    except Exception as e:
        print(f"⚠️ Failed to find rest index: {e}")
        return pd.DataFrame()

    result_list = []

    for index, row in synCurr_siz_df.iterrows():
        # Peak depol after rest
        try:
            peak_value = row.iloc[rest_idx:].max()
            rest_value = row.iloc[rest_idx]
        except IndexError:
            print(f"⚠️ Rest index {rest_idx} out of range for row {index} — skipping")
            continue

        # Depolarization = peak after rest minus rest
        result = peak_value - rest_value
        result_list.append(result)

    synSite_df['summed_syn_depol'] = result_list

    sum_syn_depol = synSite_df.groupby('pre')['summed_syn_depol'].sum().reset_index()
    sum_syn_depol['syn_count'] = synSite_df.groupby('pre').size().reset_index(name='syn_count')['syn_count']

    new_df = pd.DataFrame(sum_syn_depol)
    return new_df

def combine_indiv_syn_depol_and_group_act(neuron_name):
    df = Peak_depol_per_VPN(neuron_name)
    df2 = individual_syn_depol_per_VPN(neuron_name)
    df['pre'] = df['pre'].astype(str)
    merged_df = pd.merge(df, df2, on='pre', how='inner')
    print(df)
    print(df2)
    unique_types = merged_df['type'].unique()
    print(unique_types)
    fig, axes = plt.subplots(nrows=len(unique_types), ncols=1, figsize=(8, 6))

    # Plot for each unique type
    for i, type_val in enumerate(unique_types):
        ax = axes[i] if len(unique_types) > 1 else axes  # Use the appropriate axis for the subplot
        group = merged_df[merged_df['type'] == type_val]
        ax.plot(group['syn_count'], group['summed_syn_depol'], label=type_val, color = "red")
        ax.plot(group['syn_count'], group['depol_together'], label=type_val, color = "blue")
        ax.set_xlabel('syn_count')
        ax.set_ylabel('EPSP Amplitude')
        ax.set_title(f'Plot for Type: {type_val}')
        ax.legend()

    # Adjust layout
    plt.tight_layout()

    type_colors = {'LC4': '#0000FF',
        'LPLC2': '#FFA500',
        'LPLC1': '#FF0D00',
        'LC22': '#FFF000',
        'LPLC4': '#00FF00',
        'LC6': '#FF00E9'}  # Add more colors as needed

    # Create the scatter plot
    plt.figure(figsize=(8, 6))
    for type_val in unique_types:
        group = merged_df[merged_df['type'] == type_val]
        plt.scatter(group['syn_count'], group['summed_syn_depol'] - group['depol_together'], label=type_val, color=type_colors.get(type_val, 'black'))

    # Set labels and title
    plt.xlabel('syn_count')
    plt.ylabel('Difference (summed_syn_depol - depol_together)')
    plt.title('Difference between summed_syn_depol and depol_together vs. syn_count')

    # Add legend
    plt.legend()

    # Create the scatter plot
    plt.figure(figsize=(8, 6))
    for type_val in unique_types:
        group = merged_df[merged_df['type'] == type_val]
        plt.scatter(group['syn_count'], group['avg_syn_depol'], label=type_val, color=type_colors.get(type_val, 'black'))

    # Set labels and title
    plt.xlabel('syn_count')
    plt.ylabel('Average syn_depol (from actual VPN activation)')
    plt.title('Average syn_depol vs. syn_count')

    # Add legend
    plt.legend()
    # Show plot
    plt.show()

def combined_close_syn_depol(neuron_name, n_select=10):
    """
    For a given neuron, calculate depolarizations for close synapse activations.
    Collects depol when activated together (close sims) and compares to the
    linear sum of individual synapses, pulled directly from single-site data.
    """

    if neuron_name == "DNp01":
        syn_counts = [2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 22, 25, 31, 32]
    elif neuron_name == "DNp03":
        syn_counts = [2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17,
                      18, 19, 20, 21, 22, 23, 24, 29, 32]
    elif neuron_name == "DNp01_hemi":
        syn_counts = [5, 11, 12, 15, 16, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 36, 39, 47]
    elif neuron_name == "DNp03_hemi":
        syn_counts = [6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 27, 30, 38]

    # === Load single synapse activation data once ===
    base_path = f"datafiles/simulationData/{neuron_name}_final_sims/singlesynactivation/"
    try:
        single_sites_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_synapse_sites.csv"))
        single_curr_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_synaptic_currents_siz.csv"), header=None)
        time_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_time.csv"), header=None)
    except Exception as e:
        print(f"⚠️ Failed to load single synapse activation data: {e}")
        return pd.DataFrame()

    try:
        time_array = time_df.to_numpy().flatten()
        rest_idx = (np.abs(time_array - 48)).argmin()
    except Exception as e:
        print(f"⚠️ Failed to find rest index in single synapse data: {e}")
        return pd.DataFrame()

    # Compute single-synapse depols (add to site table)
    single_sites_df['summed_syn_depol'] = [
        row.iloc[rest_idx:].max() - row.iloc[rest_idx] for _, row in single_curr_df.iterrows()
    ]

    type_colors = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    all_results = []

    # === Loop over synapse counts ===
    for i, x in enumerate(syn_counts, start=1):
        dir_path = f"datafiles/simulationData/{neuron_name}_final_sims/CloseSims/{neuron_name}_close{x}record150"
        print(f"Processing {i}/{len(syn_counts)} → syn_count = {x}")
        if not os.path.isdir(dir_path):
            continue

        trial_paths = list(iterQuintsV9(dir_path))
        if not trial_paths:
            continue

        # Quantize: pick ~10 evenly spaced trials
        n_trials = len(trial_paths)
        step = max(1, n_trials // n_select)
        selected_trials = trial_paths[::step][:n_select]

        for quint in selected_trials:
            time_df_trial, vSIZ_np, vsoma_np, synSite_df, stim_summary_df = quint
            if vSIZ_np.size == 0 or synSite_df.empty:
                continue

            # depol together
            time_array_trial = time_df_trial.to_numpy().flatten()
            rest_idx_trial = (np.abs(time_array_trial - 48)).argmin()
            Vrest_siz = vSIZ_np[0, rest_idx_trial]
            peakV_siz = np.max(vSIZ_np[:, rest_idx_trial:])
            depol_together = peakV_siz - Vrest_siz

            # === Linear sum from single-site data (no merging) ===
            summed_syn_depol = 0.0
            for _, syn in synSite_df.iterrows():
                match = single_sites_df[
                    (single_sites_df['pre_x'] == syn['pre_x']) &
                    (single_sites_df['pre_y'] == syn['pre_y']) &
                    (single_sites_df['pre_z'] == syn['pre_z']) &
                    (single_sites_df['post_x'] == syn['post_x']) &
                    (single_sites_df['post_y'] == syn['post_y']) &
                    (single_sites_df['post_z'] == syn['post_z'])
                ]
                if not match.empty:
                    summed_syn_depol += match['summed_syn_depol'].iloc[0]

            # type assignment
            unique_types = synSite_df['type'].unique()
            syn_type = unique_types[0] if len(unique_types) == 1 else "Mixed"

            all_results.append({
                'depol_together': depol_together,
                'summed_syn_depol': summed_syn_depol,
                'syn_count': x,
                'type': syn_type
            })

    if not all_results:
        print("⚠️ No valid data across synapse counts")
        return pd.DataFrame()

    final_df = pd.DataFrame(all_results)
    print(final_df)

    # === Plotting ===
    plt.figure(figsize=(10, 7))
    seen_labels = set()
    for _, row in final_df.iterrows():
        label = row['type']
        color = type_colors.get(label, 'black')
        diff_val = row['summed_syn_depol'] - row['depol_together']
        if label not in seen_labels:
            plt.scatter(row['syn_count'], diff_val,
                        label=label, color=color, alpha=0.7)
            seen_labels.add(label)
        else:
            plt.scatter(row['syn_count'], diff_val,
                        color=color, alpha=0.7)

    plt.xlabel('Synapse Count')
    plt.ylabel('Difference (summed_syn_depol - depol_together)')
    plt.title(f'{neuron_name} — Close Synapse Depolarizations ({n_select} trials per count)')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return final_df

def combined_random_syn_depol(neuron_name):
    """
    For a given neuron, calculate depolarizations for random synapse activations.
    Collects depol when activated together (random sims) and compares to the
    linear sum of individual synapses, pulled directly from single-site data.
    """

    # === Load single synapse activation data once ===
    base_path = f"datafiles/simulationData/{neuron_name}_final_sims/singlesynactivation/"
    try:
        single_sites_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_synapse_sites.csv"))
        single_curr_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_synaptic_currents_siz.csv"), header=None)
        time_df = pd.read_csv(os.path.join(base_path, "SINGLE_EXP2_time.csv"), header=None)
    except Exception as e:
        print(f"⚠️ Failed to load single synapse activation data: {e}")
        return pd.DataFrame()

    try:
        time_array = time_df.to_numpy().flatten()
        rest_idx = (np.abs(time_array - 48)).argmin()
    except Exception as e:
        print(f"⚠️ Failed to find rest index in single synapse data: {e}")
        return pd.DataFrame()

    # Compute single-synapse depols (add to site table)
    single_sites_df['summed_syn_depol'] = [
        row.iloc[rest_idx:].max() - row.iloc[rest_idx] for _, row in single_curr_df.iterrows()
    ]

    type_colors = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    all_results = []

    # === One folder for random trials ===
    dir_path = f"datafiles/simulationData/{neuron_name}_final_sims/rand_vs_partner/{neuron_name}_rand"
    if not os.path.isdir(dir_path):
        print(f"⚠️ Folder not found: {dir_path}")
        return pd.DataFrame()

    trial_paths = list(iterQuintsV11(dir_path))
    if not trial_paths:
        print("⚠️ No trials found in folder")
        return pd.DataFrame()

    # === Loop over trials with counter ===
    for i, quint in enumerate(trial_paths, start=1):
        print(f"Processing trial {i}/{len(trial_paths)}")

        time_df_trial, vSIZ_np, vsoma_np, synSite_df, stim_summary_df = quint
        if vSIZ_np.size == 0 or synSite_df.empty:
            continue

        # depol together
        time_array_trial = time_df_trial.to_numpy().flatten()
        rest_idx_trial = (np.abs(time_array_trial - 48)).argmin()
        Vrest_siz = vSIZ_np[0, rest_idx_trial]
        peakV_siz = np.max(vSIZ_np[:, rest_idx_trial:])
        depol_together = peakV_siz - Vrest_siz

        # === Linear sum from single-site data (no merging) ===
        summed_syn_depol = 0.0
        for _, syn in synSite_df.iterrows():
            match = single_sites_df[
                (single_sites_df['pre_x'] == syn['pre_x']) &
                (single_sites_df['pre_y'] == syn['pre_y']) &
                (single_sites_df['pre_z'] == syn['pre_z']) &
                (single_sites_df['post_x'] == syn['post_x']) &
                (single_sites_df['post_y'] == syn['post_y']) &
                (single_sites_df['post_z'] == syn['post_z'])
            ]
            if not match.empty:
                summed_syn_depol += match['summed_syn_depol'].iloc[0]

        # type assignment
        unique_types = synSite_df['type'].unique()
        syn_type = unique_types[0] if len(unique_types) == 1 else "Mixed"

        all_results.append({
            'depol_together': depol_together,
            'summed_syn_depol': summed_syn_depol,
            'syn_count': len(synSite_df),   # now directly count synapses
            'type': syn_type
        })

    if not all_results:
        print("⚠️ No valid results from trials")
        return pd.DataFrame()

    final_df = pd.DataFrame(all_results)
    print(final_df)

    # === Plotting ===
    plt.figure(figsize=(10, 7))
    seen_labels = set()
    for _, row in final_df.iterrows():
        label = row['type']
        color = type_colors.get(label, 'black')
        diff_val = row['summed_syn_depol'] - row['depol_together']
        if label not in seen_labels:
            plt.scatter(row['syn_count'], diff_val,
                        label=label, color=color, alpha=0.7)
            seen_labels.add(label)
        else:
            plt.scatter(row['syn_count'], diff_val,
                        color=color, alpha=0.7)

    plt.xlabel('Synapse Count')
    plt.ylabel('Difference (summed_syn_depol - depol_together)')
    plt.title(f'{neuron_name} — Random Synapse Depolarizations ({len(trial_paths)} trials)')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return final_df

def compare_random_vs_actual_syn_depol(neuron_name):
    """
    Combines random and actual VPN group depolarizations on the same plots.
    Plots:
      (1) Difference between summed_syn_depol and depol_together vs syn_count
      (2) Average synaptic depol vs syn_count
    """

    # === Load actual VPN activation data ===
    df_actual = Peak_depol_per_VPN(neuron_name)
    df_indiv = individual_syn_depol_per_VPN(neuron_name)

    if df_actual.empty or df_indiv.empty:
        print("⚠️ Missing VPN activation data.")
        return

    df_actual['pre'] = df_actual['pre'].astype(str)
    merged_df = pd.merge(df_actual, df_indiv, on='pre', how='inner')
    merged_df['source'] = 'Actual'

    # === Load random activation data (reuse your logic) ===
    rand_df = combined_random_syn_depol(neuron_name)
    if rand_df.empty:
        print("⚠️ Missing random activation data.")
        return

    rand_df['avg_syn_depol'] = rand_df['summed_syn_depol'] / rand_df['syn_count']
    rand_df['source'] = 'Random'

    # === Combine both into one dataframe for plotting ===
    combined_df = pd.concat([merged_df, rand_df], ignore_index=True)

    type_colors = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    # === Figure 1: Difference (summed_syn_depol - depol_together) ===
    plt.figure(figsize=(8, 6))
    for src, style in zip(['Actual', 'Random'], ['o', 's']):
        subset = combined_df[combined_df['source'] == src]
        for type_val in subset['type'].unique():
            group = subset[subset['type'] == type_val]
            plt.scatter(
                group['syn_count'],
                group['summed_syn_depol'] - group['depol_together'],
                label=f"{type_val} ({src})",
                color=type_colors.get(type_val, 'black'),
                alpha=0.7,
                marker=style
            )
    plt.xlabel('Synapse Count')
    plt.ylabel('Difference (summed_syn_depol - depol_together)')
    plt.title(f'{neuron_name} — Random vs Actual VPN Synapse Depolarizations')
    plt.legend()
    plt.tight_layout()

    # === Figure 2: Average synaptic depol ===
    plt.figure(figsize=(8, 6))
    for src, style in zip(['Actual', 'Random'], ['o', 's']):
        subset = combined_df[combined_df['source'] == src]
        for type_val in subset['type'].unique():
            group = subset[subset['type'] == type_val]
            plt.scatter(
                group['syn_count'],
                group['avg_syn_depol'],
                label=f"{type_val} ({src})",
                color=type_colors.get(type_val, 'black'),
                alpha=0.7,
                marker=style
            )
    plt.xlabel('Synapse Count')
    plt.ylabel('Average Synaptic Depol (mV)')
    plt.title(f'{neuron_name} — Average Synaptic Depolarization')
    plt.legend()
    plt.tight_layout()
    plt.show()

def compare_close_vs_actual_syn_depol(neuron_name, n_select=10):
    """
    Combines close-group and actual VPN depolarizations on the same plots.
    Plots:
      (1) Difference between summed_syn_depol and depol_together vs syn_count
      (2) Average synaptic depol vs syn_count
    """

    # === Load actual VPN activation data ===
    df_actual = Peak_depol_per_VPN(neuron_name)
    df_indiv = individual_syn_depol_per_VPN(neuron_name)

    if df_actual.empty or df_indiv.empty:
        print("⚠️ Missing VPN activation data.")
        return

    df_actual['pre'] = df_actual['pre'].astype(str)
    merged_df = pd.merge(df_actual, df_indiv, on='pre', how='inner')
    merged_df['source'] = 'Actual'

    # === Load close synapse activation data ===
    close_df = combined_close_syn_depol(neuron_name, n_select=n_select)
    if close_df.empty:
        print("⚠️ Missing close-group activation data.")
        return

    close_df['avg_syn_depol'] = close_df['summed_syn_depol'] / close_df['syn_count']
    close_df['source'] = 'Close'

    # === Combine both datasets ===
    combined_df = pd.concat([merged_df, close_df], ignore_index=True)

    type_colors = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    # === FIGURE 1: Difference plot ===
    plt.figure(figsize=(8, 6))
    for src, style in zip(['Actual', 'Close'], ['o', 'D']):
        subset = combined_df[combined_df['source'] == src]
        for type_val in subset['type'].unique():
            group = subset[subset['type'] == type_val]
            plt.scatter(
                group['syn_count'],
                group['summed_syn_depol'] - group['depol_together'],
                label=f"{type_val} ({src})",
                color=type_colors.get(type_val, 'black'),
                alpha=0.7,
                marker=style
            )
    plt.xlabel('Synapse Count')
    plt.ylabel('Difference (summed_syn_depol - depol_together)')
    plt.title(f'{neuron_name} — Close vs Actual VPN Synapse Depolarizations')
    plt.legend()
    plt.tight_layout()

    # === FIGURE 2: Average synaptic depol ===
    plt.figure(figsize=(8, 6))
    for src, style in zip(['Actual', 'Close'], ['o', 'D']):
        subset = combined_df[combined_df['source'] == src]
        for type_val in subset['type'].unique():
            group = subset[subset['type'] == type_val]
            plt.scatter(
                group['syn_count'],
                group['avg_syn_depol'],
                label=f"{type_val} ({src})",
                color=type_colors.get(type_val, 'black'),
                alpha=0.7,
                marker=style
            )
    plt.xlabel('Synapse Count')
    plt.ylabel('Average Synaptic Depol (mV)')
    plt.title(f'{neuron_name} — Average Synaptic Depolarization (Close vs Actual)')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return combined_df

def plot_actual_VPN_syn_depol(neuron_name):
    """
    Plots depolarization properties of actual VPN activations:
      (1) Linear sum of individual synaptic depols vs synapse count
      (2) Coactivated synaptic depol vs synapse count
      (3) Overlay of both on the same axes

    Uses:
      Peak_depol_per_VPN(neuron_name)
      individual_syn_depol_per_VPN(neuron_name)
    """

    # === Load data ===
    df_group = Peak_depol_per_VPN(neuron_name)
    df_indiv = individual_syn_depol_per_VPN(neuron_name)

    if df_group.empty or df_indiv.empty:
        print("⚠️ Missing VPN activation data for:", neuron_name)
        return

    df_group['pre'] = df_group['pre'].astype(str)
    merged_df = pd.merge(df_group, df_indiv, on='pre', how='inner')

    # === Prepare plot aesthetics ===
    type_colors = {
        'LC4': '#0000FF', 'LPLC2': '#FFA500', 'LPLC1': '#FF0D00',
        'LC22': '#FFF000', 'LPLC4': '#00FF00', 'LC6': '#FF00E9'
    }

    unique_types = merged_df['type'].unique()

    # === Figure 1: Linear sum only ===
    plt.figure(figsize=(8, 6))
    for t in unique_types:
        group = merged_df[merged_df['type'] == t]
        plt.scatter(group['syn_count'], group['summed_syn_depol'],
                    label=t, color=type_colors.get(t, 'black'), alpha=0.8)
    plt.xlabel('Synapse Count')
    plt.ylabel('Linear Sum of Individual Syn Depols (mV)')
    plt.title(f'{neuron_name} — Linear Sum of Synaptic Depolarizations')
    plt.legend()
    plt.tight_layout()

    # === Figure 2: Coactivation depol only ===
    plt.figure(figsize=(8, 6))
    for t in unique_types:
        group = merged_df[merged_df['type'] == t]
        plt.scatter(group['syn_count'], group['depol_together'],
                    label=t, color=type_colors.get(t, 'black'), alpha=0.8)
    plt.xlabel('Synapse Count')
    plt.ylabel('Depol from Coactivated Synapses (mV)')
    plt.title(f'{neuron_name} — Coactivation Depolarization at SIZ')
    plt.legend()
    plt.tight_layout()

    # === Figure 3: Overlay both ===
    plt.figure(figsize=(8, 6))
    for t in unique_types:
        group = merged_df[merged_df['type'] == t]
        plt.scatter(group['syn_count'], group['summed_syn_depol'],
                    label=f'{t} (Linear Sum)',
                    color=type_colors.get(t, 'black'), alpha=0.6, marker='o')
        plt.scatter(group['syn_count'], group['depol_together'],
                    label=f'{t} (Coactivated)',
                    color=type_colors.get(t, 'black'), alpha=0.9, marker='x')
    plt.xlabel('Synapse Count')
    plt.ylabel('Depolarization (mV)')
    plt.title(f'{neuron_name} — Linear Sum vs Coactivation Comparison')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return merged_df

######


def main():
    print(os.getcwd())

    ##### 

    #Running receptive field analysis for COM in particular hemispheres
    #Figure 3 Supplemental 1
    # VPN = "LPLC2"
    # COM_receptive_field_split_and_visualization(VPN)

    #####
    #Define and intialize neuron model
    neuron_name = "DNp01" # "DNp01", "DNp01_hemi", "DNp02", "DNp03", "DNp03_hemi", "DNp04", "DNp06"

    cell, allSections_py, allSections_nrn, somaSection, erev, axonList, tetherList, dendList, sizSection = initializeModel(neuron_name)#sizSection, erev, axonList = initializeModel(neuron_name)

    if neuron_name == "DNp01":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp01_morphData/synMap_DNp01_8_19_2025.csv', dtype={'pre': str})
    elif neuron_name == "DNp01_hemi":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp01_morphData/synMap_DNp01_hemibrain.csv', dtype={'pre': str})  
    elif neuron_name == "DNp02":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp02_morphData/synMap_DNp02_8_19_2025.csv', dtype={'pre': str}) 
    elif neuron_name == "DNp03":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp03_morphData/synMap_DNp03_8_19_2025.csv', dtype={'pre': str})    
    elif neuron_name == "DNp03_hemi":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp03_morphData/synMap_DNp03_hemibrain.csv', dtype={'pre': str})  
    elif neuron_name == "DNp04":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp04_morphData/synMap_DNp04_8_19_2025.csv', dtype={'pre': str})    
    elif neuron_name == "DNp06":
        synSite_df = pd.read_csv('datafiles/morphologyData/DNp06_morphData/synMap_DNp06_8_19_2025.csv', dtype={'pre': str})

    synSite_df = removeAxonalSynapses(synSite_df, dendList, axonList=axonList, tetherList=tetherList, somaSection=somaSection, sizSection=sizSection)
    

    ########
    #Running retinotopy analysis, Figure 4 D-E and Figure 4 supplementals

    #Define VPN population to look at for hemibrain models follow the naming convention below, FAFB models use the same names without the _hemi suffix

    # VPN = "LPLC1_hemi"
    # # #Running retinotopy analysis, Figure 4 D-E and Figure 4 supplemental and hemibrain models
    # dist_df = calculate_distance(neuron_name, VPN, synSite_df)
    # neuron_pair_dist =  calculate_average_distance_between_synapses_per_neuron_pair(neuron_name, VPN, synSite_df)
    # nearest_neighbor_df, avg_distance_per_neuron_df = nearest_neighbor(neuron_name, VPN, sizSection, synSite_df)
    # nn_and_COM_col_coded, all_syns_dist_df, neuron_syn_pair_dists = NN_COM_merge(VPN, nearest_neighbor_df, dist_df, neuron_pair_dist)
    # plot_retinotopy(neuron_syn_pair_dists, VPN, neuron_name)

    #Makes all figures needed for receptive field mapping and electrotonic analysis
    # syn_spread_attenuation(neuron_name, VPN="LPLC2")

    #############
    # # #Figure 8 and supplemental

    # plotAllSynsAndhighlightsynsFig8(neuron_name)
    # plot_peak_depol_at_SIZ_vs_dist_to_SIZ(neuron_name, sizSection)

    #############

    #Figure 9 and supplemental

    # plot_partner_VPN_activations(neuron_name, MODE= 'VPN/')
    # plot_rand_vs_partner_SIZ(neuron_name, MODE='VPN/') 
    # plot_rand_vs_partner_SOMA(neuron_name, MODE= 'VPN/')

    #Preprint version
    # Plot_DN_VPN_syns(neuron_name, location = 'SIZ', csv_filename=f"{neuron_name}_VPN_population_depol_data.csv")
    #Looks at all populations
    DNs = ("DNp01", "DNp03")
    # DNs = ("DNp01_hemi", "DNp03_hemi")
    # shunting_rate_vs_density_with_exp_curve(neuron_list = DNs, dataset_type="hemibrain")  
    # shunting_rate_vs_density(neuron_list = DNs, gap=0, dataset_type="FAFB") 

    #New version plotting simulation data by neighboring VPNs based on receptive field selections
    # DNp01
    # df = Plot_DN_VPN_neighborhood_grid(
    # neuron_name  = neuron_name,
    # VPN_list     = ["LC4", "LPLC2"],
    # location     = "SIZ",
    # csv_filename = f"{neuron_name}_LC4_LPLC2_neighborhood.csv",  max_neighbors=40)

    #DNp03
    # df = Plot_DN_VPN_neighborhood_grid(
    # neuron_name  = neuron_name,
    # VPN_list     = ["LC4", "LPLC1", "LPLC4"],
    # location     = "SIZ",
    # csv_filename = f"{neuron_name}_LC4_LPLC1_LPLC4_neighborhood.csv",  max_neighbors=15)

    #Calculate density of synapses per VPN population
    #calculate surface area of branches containing synapses from a given type, and export those branches to an SWC file to visualize
    #Only provide on gap at a time with a swc filename to save the branches used to calculate surface area for that gap, otherwise it will overwrite
    # gap = (0, 1) # , 1, 2, 3, 5, 10, 15
    # for x in gap:
    #     target_type = "LC22"
    #     print(f"gap for {target_type} {x}")
    #     area, sections = calculate_branch_surface_area_and_export_swc(
    #         synSite_df,
    #         dendList,
    #         target_type,
    #         swc_filename=f"{neuron_name}_{target_type}_selected_branches.swc", #f"{neuron_name}_{target_type}_selected_branches.swc"
    #         gap_allowance=x 
    #     )
    # Looking for on average how many sections in between synapses to account for total surface area
    # plot_nearest_neighbor_section_gaps(synSite_df, target_type="LC4", show_plot=True)

    #This function below takes in the data from the above function, the data was manually copied from the terminal output from the above function
    # plot_density_of_VPN_synsapses(csv_filename="VPN-DNs_synapse_density.csv")

    #New analysis including VPN retinotopic activations without fitting the shunting constant
    # res_df = shunting_rate_vs_density_with_neighborhood(
    #     neuron_list  = ["DNp01", "DNp03"],
    #     VPN_config   = {
    #         "DNp01": ["LC4", "LPLC2"],
    #         "DNp03": ["LC4", "LPLC1", "LPLC4"]
    #     },
    #     gap          = 0,
    #     dataset_type = "FAFB",
    #     max_neighbors= 15)

    # combined_res_df = shunting_rate_vs_density_with_neighborhood_fit(
    # neuron_list  = ["DNp01", "DNp03"],
    # VPN_config   = {
    #     "DNp01": ["LC4", "LPLC2"],
    #     "DNp03": ["LC4", "LPLC1", "LPLC4"]
    # },
    # gap          = 0,           # which gap to show in Figure 1 right subplot
    # dataset_type = "FAFB",      # or "hemibrain"
    # max_neighbors= 40)

    neuron_list = ["DNp01", "DNp03"]

    VPN_config = {
        "DNp01": ["LC4_hemi", "LPLC2_hemi"],
        "DNp03": ["LC4_hemi", "LPLC1_hemi", "LC22_hemi", "LPLC4_hemi"]
    }

    # VPN_FAFB_config = {
    #         "DNp01": ["LC4", "LPLC2"],
    #         "DNp03": ["LC4", "LPLC1", "LC22", "LPLC4"]
    #     }

    # nd_res_df = shunting_rate_vs_density_neighborhood_only(
    #     neuron_list=neuron_list,
    #     VPN_config=VPN_config,
    #     gap=0,
    #     dataset_type="FAFB",
    #     max_neighbors=60,
    # )

    #############
    #Figure 10

    #Number of VPNs per synapse count
    # plot_syn_spread_vs_dist_to_siz_partner_vs_rand(neuron_name, plots='Partner', sizSection=sizSection) # must run before the next functions
    # plot_syn_spread_vs_dist_to_siz_partner_vs_rand(neuron_name, plots='Rand', sizSection=sizSection) # must run before the next functions

    # VPN_per_syn_count(dir_path=f'datafiles/simulationData/{neuron_name}_final_sims', neuron_name=neuron_name) #used to determine which synapse counts were run

    #############
    #Generate final csv's for figure 10, simulations by synapse spread 
    # syn_list_DNp01 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 22, 25, 29, 31] 
    # syn_list_DNp03 = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 29, 32] 
    # syn_list_DNp01_hemi = [5, 11, 12, 15, 16, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 35, 36, 39, 47] 
    # syn_list_DNp03_hemi = [6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 27, 30, 38] 

    # #DNp01
    # plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=sizSection, Dist=None, syncount=None, close_analysis = False, partner_analysis= True, trials = 150)
    # for x in syn_list_DNp01: 
    #     plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=sizSection, Dist=None, syncount=x, close_analysis = True, partner_analysis = False, trials = 150)
    # for y in [10, 20, 30, 40, 50, 60, 70]:
    #         plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_physical(neuron_name, sizSection=sizSection, Dist=y, syncount=x) #old do not need, this is physical distance
    #         plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=sizSection, Dist=y, syncount=x, close_analysis = True, partner_analysis = False, trials = None)

    # #DNp03
    # plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=sizSection, Dist=None, syncount=None, close_analysis = False, partner_analysis= True, trials = 150)
    # for x in syn_list_DNp03:
    #     plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=sizSection, Dist=None, syncount=x, close_analysis = True, partner_analysis = False, trials = 150)
    #     for y in [10, 20, 30, 40, 50, 60, 70]:
    #         # plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_physical(neuron_name, erev, sizSection=sizSection, Dist=y, syncount=x) #old do not need, this is physical distance
    #         plot_syn_spread_vs_dist_to_siz_by_average_syn_spread_electrotonic(neuron_name, sizSection=sizSection, Dist=y, syncount=x, close_analysis = False, partner_analysis = False, trials = None)
            
    #final plots        
    # plot_partner_close_all_electro_bins(neuron_name, n_e_bins=7, export_csv=True, csv_path=f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_Levene_results.csv") # shows all synapse, and adds levenes test for variance
    # plot_partner_close_all_combined(neuron_name, export_csv=True, csv_path=None) #this is the final version for the paper

    #stats
    # run_levene_test(f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_Levene_results.csv")
    # print_significant_levene(f"datafiles/simulationData/{neuron_name}_final_sims/{neuron_name}_Levene_results.csv")


    #### Additional code removed from main analysis: for review purposes only

    #running data analysis for nearest neighboring syns
    #Random ingeger seeds are [266,	935, 4, 598, 466, 752, 963, 582, 680, 530]
    # VPN = "LPLC1"
    # rndm_inrts = [266, 935, 4, 598, 466, 752, 963, 582, 680, 530]
    # for x in rndm_inrts:
    #     print("random seed: ", x)
    #     synMap_df, shuffled_df = shuffle_data(neuron_name, VPN, synSite_df, int_seed= x)
    #     nearest_neighbor_df_shuffled, avg_distance_per_neuron_df_shuffled = nearest_neighbor(neuron_name, VPN, sizSection, shuffled_df)
    #     nearest_neighbor_analysis(neuron_name, nearest_neighbor_df_shuffled) 

    # #Run this once, to get the original data values below for each VPN-DN pair.
    # # nearest_neighbor_df, avg_distance_per_neuron_df = nearest_neighbor(neuron_name, VPN, sizSection, synSite_df)
    # # nearest_neighbor_analysis(neuron_name, nearest_neighbor_df) #Run this once, to get the original data values below for each VPN-DN pair.

   
    # bar_whisker_plot_shuffled_vs_original_data(dataset= "FAFBv783_updated")

    # #FAFB shuffled data for NN analysis, preprint data using codex version 630
    # # LC4_shuffled_DNp01  = [4.42, 2.31, 1.35, 1.73, 1.15, 0.96, 2.31, 3.46, 2.12, 1.73] # Original 44.56
    # # LPLC2_shuffled_DNp01 = [0.66,2.49, 1.83, 1, 1.33, 0.83, 0.89, 1, 1.16, 2.16] # Original 36.88
    # # LC4_shuffled_DNp02  = [2.8, 1.78, 2.8, 1.65, 3.56, 3.05, 0.89, 2.04, 2.67, 3.31] # Original 51.02 
    # # LC4_shuffled_DNp03  = [1.01, 1.01, 3.03, 3.03, 4.04, 3.03, 6.06, 2.02, 3.03, 4.04] # Original 49.8
    # # LPLC2_shuffled_DNp03  = [1.96, 1.76, 2.75, 0.78, 1.76, 4.12, 1.76, 1.96, 1.76, 1.37] # Original 43.75
    # # LC22_shuffled_DNp03  = [13.18, 8.18, 3.85, 7.27, 1.82, 9.09, 4.09, 3.64, 5.91, 6.82] # Original  57.27
    # # LPLC1_shuffled_DNp03  = [2.19, 3.01, 2.1, 1.37, 2.28, 2.56, 1.64, 1, 2.19, 1.74] # Original 40
    # # LPLC4_shuffled_DNp03  = [1.77, 2.28, 0.84, 1.94, 2.02, 1.69, 1.94, 2.7, 1.94, 2.02] # Original 58.52
    # # LPLC2_shuffled_DNp04 = [0.59, 1.3, 1.65, 1.18, 1.54, 0.71, 1.06, 1.14, 1.54, 1.42] # Original 54.13
    # # LC4_shuffled_DNp04  = [2.21, 1.64, 3.02, 1.7, 2.21, 2.08, 2.39, 1.89, 2.39, 2.21] # Original 46.93
    # # LC4_shuffled_DNp06  = [1.01, 1.01, 3.03, 3.03, 4.04, 3.03, 6.06, 2.02, 3.03, 4.04] # Original 65.66
    # # LPLC2_shuffled_DNp06  = [1.09, 1.09, 1.4, 1.24, 0.78, 2.33, 1.86, 1.4, 1.09, 1.24] # Original 44.72
    # # LC6_shuffled_DNp06  = [9.92, 1.92, 3.85, 9.62, 15.38	, 5.77, 15.38, 5.77, 5.77, 5.77] #Original 53.85
    # # LPLC1_shuffled_DNp06  = [3.9, 2.31, 0.87, 1.15, 2.45, 2.31, 2.16, 1.3, 2.45, 1.15] #Original 45.02
    # # original_data = [44.56, 36.88, 51.02, 49.8, 43.75, 57.27, 40, 58.52, 54.13, 46.93, 65.66, 44.72, 53.85, 45.02]
    # # combined_data = [LC4_shuffled_DNp01, LPLC2_shuffled_DNp01, LC4_shuffled_DNp02, LC4_shuffled_DNp03, LPLC2_shuffled_DNp03,LC22_shuffled_DNp03,
    # #                 LPLC1_shuffled_DNp03,  LPLC4_shuffled_DNp03, LC4_shuffled_DNp04, LPLC2_shuffled_DNp04, LC4_shuffled_DNp06, LPLC2_shuffled_DNp06,
    # #                 LC6_shuffled_DNp06, LPLC1_shuffled_DNp06]

    # #FAFB shuffled data for NN analysis, finalized elife, data using codex version 783, from codex
    # LC4_shuffled_DNp01  = [2.14, 2.14, 2.67, 1.34, 1.34, 2.41, 2.41, 2.14, 2.41, 2.94] # Original 17.38 
    # LPLC2_shuffled_DNp01 = [1.75, 1.97, 2.84, 1.31, 0.66, 1.09, 0.44, 1.09, 1.75, 1.31] # Original 16.81
    # LC4_shuffled_DNp02  = [2.10, 2.76, 4.19, 1.52, 2.29, 1.33, 2.67, 2.48, 2.76, 2.10] # Original 18.10
    # LC4_shuffled_DNp03  = [1.64, 4.10, 1.64, 2.73, 1.91, 1.37, 2.19,1.91,1.64, 4.10] # Original 19.13
    # LPLC2_shuffled_DNp03  = [0, 0, 0, 0, 9.09, 0, 0, 9.09, 0, 0 ] # Original 9.09
    # LC22_shuffled_DNp03  = [10.97, 7.74, 5.16, 5.16, 5.81, 8.39, 7.10, 1.94, 9.68, 4,52] # Original  38.71
    # LPLC1_shuffled_DNp03  = [2.36, 1.42, 2.48, 2.48, 1.30, 1.53, 1.53, 1.77, 2.36, 1.89] # Original 14.39
    # LPLC4_shuffled_DNp03  = [1.62, 1.19, 2.92, 1.51, 1.51, 1.30, 2.16, 2.16, 1.88, 1.84] # Original 30.78
    # LC4_shuffled_DNp04  = [2.54, 2.44, 1.22, 3.20, 2.16, 2.63, 2.26, 2.87, 1.50, 1.97] # Original 24.72
    # LPLC2_shuffled_DNp04 = [1.98, 1.15, 0.82, 1.48, 0.82, 0.82, 0.99, 1.81, 1.32, 2.14] # Original 15.49
    # LC4_shuffled_DNp06  = [1.72, 6.90, 1.72, 0, 6.90, 5.17, 1.72, 3.45, 5.17, 1.72] # Original 29.31
    # LPLC2_shuffled_DNp06  = [0.41, 1.65, 0.82, 1.24, 1.65, 1.85, 0.82, 0.82, 1.65, 0.41] # Original 18.76
    # LC6_shuffled_DNp06  = [4.44, 6.67, 6.67, 2.22, 4.44, 8.89, 6.67, 0, 0, 11.11] #Original 33.33
    # LPLC1_shuffled_DNp06  = [3.14, 1.77, 1.38, 2.16, 2.55, 2.55, 1.18, 1.77, 2.16, 1.57] #Original 16.70
    # original_data = [17.38, 16.81, 18.10, 19.13, 9.09, 38.71, 14.39, 30.78, 24.72, 15.49, 29.31, 18.76, 33.33, 16.70]
    # combined_data = [LC4_shuffled_DNp01, LPLC2_shuffled_DNp01, LC4_shuffled_DNp02, LC4_shuffled_DNp03, LPLC2_shuffled_DNp03,LC22_shuffled_DNp03,
    #                 LPLC1_shuffled_DNp03,  LPLC4_shuffled_DNp03, LC4_shuffled_DNp04, LPLC2_shuffled_DNp04, LC4_shuffled_DNp06, LPLC2_shuffled_DNp06,
    #                 LC6_shuffled_DNp06, LPLC1_shuffled_DNp06]

    # # #FAFB shuffled data for NN analysis, finalized elife, data using codex version 783, from newly updated synapses 8-19-2025
    # LC4_shuffled_DNp01  = [2.14, 1.78, 2.64, 1.36, 2.36, 1.43, 2.50, 1.64, 1.21, 0.86] # Original 14.92 
    # LPLC2_shuffled_DNp01 = [1.46, 0.98, 0.65, 0.65, 0.98, 2.28, 1.95, 0.81, 0.65, 1.63] # Original 9.92
    # LC4_shuffled_DNp02  = [2.99, 1.96, 1.77, 1.96, 1.87, 2.05, 1.87, 1.40, 1.96, 2.05, ] # Original 13.63
    # LC4_shuffled_DNp03  = [3.09, 1.12, 2.25, 0.42, 2.25, 1.12, 1.69, 1.69, 2.11, 3.37] # Original 18.40
    # LPLC2_shuffled_DNp03  = [0, 0, 0, 0, 7.14, 0, 0, 0, 0, 0 ] # Original 28.57
    # LC22_shuffled_DNp03  = [4.74, 6.84, 3.68, 9.47, 7.89, 3.68, 5.79, 5.26, 4.21, 4.74] # Original  29.47
    # LPLC1_shuffled_DNp03  = [2.52, 2.29, 2.21, 2.13, 2.29, 2.21, 2.05, 1.58, 1.73, 1.18] # Original 12.37
    # LPLC4_shuffled_DNp03  = [1.99, 1.18, 1.72, 2.72, 1.27, 1.63, 1.63, 1.72, 2.81, 1.54] # Original 23.03
    # LC4_shuffled_DNp04  = [2.19, 1.64, 1.81, 2.54, 2.02, 2.37, 1.39, 1.46, 2.16, 2.09] # Original 18.45
    # LPLC2_shuffled_DNp04 = [0.38, 1.38, 1.51, 1.13, 2.14, 0.88, 0.25, 1.26, 0.75] # Original 10.68
    # LC4_shuffled_DNp06  = [4.46, 2.68, 2.68, 3.57, 0, 0, 2.68, 0.89, 3.57, 2.68] # Original 34.82
    # LPLC2_shuffled_DNp06  = [0.18, 0, 1.76, 0.70, 0.53, 0.88, 0.88, 1.23, 0.88, 1.58] # Original 17.62
    # LC6_shuffled_DNp06  = [4.55, 2.27, 4.55, 13.64, 0, 4.55, 4.55, 0, 0, 0] #Original 20.45
    # LPLC1_shuffled_DNp06  = [1.24, 1.71, 0.78, 1.55, 2.02, 1.71, 1.55, 2.17, 2.33, 1.40] #Original 11.32
    # original_data = [14.92, 9.92, 13.63, 18.40, 28.57, 29.47, 12.37, 23.03, 18.45, 10.68, 34.82, 17.62, 20.45, 11.32]
    # combined_data = [LC4_shuffled_DNp01, LPLC2_shuffled_DNp01, LC4_shuffled_DNp02, LC4_shuffled_DNp03, LPLC2_shuffled_DNp03,LC22_shuffled_DNp03,
    #                 LPLC1_shuffled_DNp03,  LPLC4_shuffled_DNp03, LC4_shuffled_DNp04, LPLC2_shuffled_DNp04, LC4_shuffled_DNp06, LPLC2_shuffled_DNp06,
    #                 LC6_shuffled_DNp06, LPLC1_shuffled_DNp06]


    # # # #hemibrain shuffled data for NN analysis
    # # # # LC4_shuffled_DNp01_hemibrain  = [1.63, 2.11, 1.36, 1.45, 1.89, 1.32, 1.93, 1.10, 1.23, 1.41] # Original 13.97
    # # # # LPLC2_shuffled_DNp01_hemibrain = [1.22, 1.56, 1.29, 1.29, 1.70, 1.70, 1.22, 0.95, 1.29, 1.50] # Original 9.45
    # # # # original_data = [13.97, 9.45]
    # # # # combined_data = [LC4_shuffled_DNp01_hemibrain, LPLC2_shuffled_DNp01_hemibrain]

    # # # #Plotting of figure 5C
    # # # # Using FAFB data
    # labels = [
    #     "LC4-DNp01", "LPLC2-DNp01", "LC4-DNp02", "LC4-DNp03", "LPLC2-DNp03", "LC22-DNp03",
    #     "LPLC1-DNp03", "LPLC4-DNp03", "LC4-DNp04", "LPLC2-DNp04", "LC4-DNp06", "LPLC2-DNp06",
    #     "LC6-DNp06", "LPLC1-DNp06"
    # ]
    # # # # #Using hemibrain data
    # # # # # labels = ["LC4-DNp01", "LPLC2-DNp01"]
    # plot_vpn_dn_comparison_bargraph_with_stats(original_data, combined_data, labels)

    #######

    #Looking at the integration of synapses and how it compares to linear summation 

    # combine_indiv_syn_depol_and_group_act(neuron_name) #for actual VPNs
    # combined_close_syn_depol(neuron_name, n_select=10)
    # combined_random_syn_depol(neuron_name)# for random synapses

    #compare random vs actual VPNs
    # compare_random_vs_actual_syn_depol(neuron_name)

    #compare close vs actual VPNs
    # compare_close_vs_actual_syn_depol(neuron_name, n_select=20)

    # plot_actual_VPN_syn_depol(neuron_name)

    ######


main()


