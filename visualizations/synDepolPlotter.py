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
from scipy.stats import ttest_ind, shapiro, mannwhitneyu, linregress
from scipy import stats
from scipy.stats import ttest_1samp, levene
from scikit_posthocs import posthoc_dunn
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

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
    #if false, the section is broken up until into additional segments until the largest segment is no greater in length than 1/10 lambda

    #NOTE TO SELF, NOT IDEAL BC ONLY INCREASES, DOESNT DECREASE, WRITE CODE TO DECREMENT SEGMENT NUMBER IF APPLICABLE
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

    #DNp03 axon len - 151 | SIZ - axon[2]
    #DNp04 axon len - 77 | "SIZ" - axon[73]
    #DNp01 axon len - 33 | "SIZ" - axon[32]
    #TODO fill out remaining info

    if neuron_name == "DNp01":
        sizIndex = 0#32
    elif neuron_name == "DNp02":
        sizIndex = 2#1#128
    elif neuron_name == "DNp03":
        sizIndex = 2#1#2
    elif neuron_name == "DNp04":
        sizIndex = 4#12#73
    elif neuron_name == "DNp06":
        sizIndex = 4#5#43

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


    allSections_py = createAxon(axonEnd, allSections_py, neuron_name)

    if neuron_name == "DNp01":
        change_Ra(34)#33.2617)
        change_gLeak(3.4e-4, erev=-66.4298 - 0.2)#4.4e-9)    
        change_memCap(1)#1.4261)
        erev = -66.4298 - 0.2
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
        # equivCylHeight = 360.346/10
        # equivCylDiam = 1.199*2
        equivCylHeight = 184.656
        equivCylDiam = 1.34*2
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
    #TODO: pass figure as input/return as output so you can overlay traces

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
    '''
    ##########################
    ### VPN Pop Traces Fig ###
    ##########################

    synTypeSet = set(synTypeList)
    uniqueSynTypes = list(synTypeSet)

    gridTracesPop = gridspec.GridSpec(1, 2)
    for synType in uniqueSynTypes:
        if synType == "LC4":
            colorCode = 'g'
        elif synType == "LPLC1":
            colorCode = 'b'
        elif synType == "LPLC4":
            colorCode = 'r'
        elif synType == "LPLC2":
            colorCode = 'm'
        elif synType == "LC22":
            colorCode = 'c'
        else:
            colorCode = 'k'
        figTracesPop = plt.figure()
        axTracesPopSynSeg = figTracesPop.add_subplot(gridTracesPop[0, 0]) # row 0, col 0
        axTracesPopSIZ = figTracesPop.add_subplot(gridTracesPop[0, 1])
        for trial in range(len(syn_tVecList)):
            if synTypeList[trial] == synType:
                vSyn = syn_vSynList[trial]#.to_python()
                vSyn_np = numpy.array(vSyn)
                vSIZ = syn_vSIZList[trial]#.to_python()
                vSIZ_np = numpy.array(vSIZ)
                tVec = syn_tVecList[trial]#.to_python()
                tVec_np = numpy.array(tVec)
                axTracesPopSynSeg.plot(tVec_np, vSyn_np, colorCode)#, linewidth=0.6)
                axTracesPopSIZ.plot(tVec_np, vSIZ_np, colorCode)#, linewidth=0.6)
            else:
                pass

        axTracesPopSynSeg.set_title("synSeg recording for {} syn activations".format(synType), fontsize = 16)
        axTracesPopSIZ.set_title("SIZ recording for {} syn activations".format(synType), fontsize = 16)

        axTracesPopSynSeg.set_ylabel('Voltage (mV)')
        axTracesPopSynSeg.set_xlabel('Time (ms)')
        axTracesPopSIZ.set_ylabel('Voltage (mV)')
        axTracesPopSIZ.set_xlabel('Time (ms)')

        figTracesPop.tight_layout()
        '''

def plotMorphology(allSections_py, somaSection, axonList, tetherList, dendList, sizSection):#, sizSection):
    morph_shape_window = h.PlotShape(h.SectionList(allSections_py))           # Create a shape plot
    morph_shape_window.exec_menu('Show Diam')    # Show diameter of sections in shape plot
    morph_shape_window.color_all(9)

    #colorR for axon sections | colorB for test work | colorG for soma | colorK for presumed SIZ
    colorR = h.SectionList()
    colorB = h.SectionList()
    colorG = h.SectionList()
    colorK = h.SectionList()
    colorV = h.SectionList()

    colorG.append(somaSection)
    colorV.append(sizSection)
    #colorK.append(sizSection)
    
    for sec in axonList:
        colorR.append(sec)
    for sec in tetherList:
        colorK.append(sec)
    for sec in dendList:
        colorB.append(sec)

        # if len(sec.children()) < 1:
        #     colorB.append(sec)
        #    axonEnd = sec

    morph_shape_window.color_list(colorR, 2)
    morph_shape_window.color_list(colorG, 4)  
    morph_shape_window.color_list(colorK, 1)
    morph_shape_window.color_list(colorB, 3)
    # morph_shape_window.color_list(colorV, 7)
    # shape_window.color_list(colorB, 3)

    return morph_shape_window

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
    return syn_shape_window

def str2sec(sec_as_str):
    for sec in h.allsec():
        if sec.name() == sec_as_str:
            sec_as_sec = sec
            break
    return sec_as_sec

def selectSynapses(subsetInfo, synMap_df):
    #TODO: save synapse subset df alongside trace?
    synSubset_df = pd.DataFrame()
    
    if subsetInfo[0] == "RAND":
        #TODO: if rand is nan, REDO GEN
        #random.seed(a=1)
        randN = subsetInfo[1]
        randSelectedSynIdxs = []
        randSelectedSyns = pd.DataFrame()
        for randSynCt in range(randN):
            RANDVAR = random.randint(0, synMap_df.shape[0]-1)
            print(RANDVAR)
            randSelectedSyns = pd.concat([randSelectedSyns, synMap_df.iloc[[RANDVAR]]])
        synSubset_df = randSelectedSyns

    elif subsetInfo[0] == "CLOSEST":
        closestN = subsetInfo[1]
        closestMethod = "PATH" #"EUCLID"
        STARTSYN = synMap_df.iloc[subsetInfo[2]]
        #print(STARTSYN)

        if numpy.isnan(STARTSYN.loc['pre']):
                print("ERROR: RAND SYN IS NAN")
                raise IndexError
        else:
            for sec in h.allsec():
                if sec.name() == STARTSYN.loc["mappedSection"]:
                    STARTSYN_SEC = sec
                    break
        
        SYNDIST_DATA = []

        for index, row in synMap_df.iterrows():

            if numpy.isnan(row.loc['pre']):
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
    
    elif subsetInfo[0] == "BRANCH":
        #retrieveSynBranch(synMap_df, somaSection):
        for sec in h.allsec():
            #print(sec.name())
            if sec.name() == "dend[849]":
            #if sec.name() == "dend[416]":
                STARTSEC = sec
                break
        
        colorB = h.SectionList()
        tSec = somaSection
        listOfKidsToCheck = []

        #print(h.psection(sec=STARTSEC))
        tSec = STARTSEC
        c_sref = h.SectionRef(sec=STARTSEC)
        # print("PARENT:", c_sref.parent)
        # print("CHILDREN", STARTSEC.children())
        while True:

            while_sref = h.SectionRef(tSec)
            parent = while_sref.parent

            print("PARENT:", while_sref.parent)
            # print("CHILDREN", tSec.children()) 
            
            if "soma" in parent.name():
                MYSEC = sec
                break
            else:
                tSec = parent
                colorB.append(tSec)

        SUBSET_DF = pd.DataFrame()
        for sec in colorB:
            DF_TO_APPEND = synMap_df.loc[synMap_df['mappedSection'] == sec.name()]
            if DF_TO_APPEND.empty:
                pass
            else:
                SUBSET_DF = pd.concat([SUBSET_DF, synMap_df.loc[synMap_df['mappedSection'] == sec.name()]])
        # return colorB, SUBSET_DF


    return synSubset_df

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

def returnQuads(dir_path):
    dirList = os.listdir(dir_path)
    numQuads =  int(len([entry for entry in dirList if os.path.isfile(os.path.join(dir_path, entry))])/4)

    dirList.sort()
    quadList = []

    for quad_set in range(numQuads):
        print(quad_set, "/250 IN RETURNQUADS")
        synSite_df = dirList[0+(quad_set*4)]
        synVoltSIZ_df = dirList[1+(quad_set*4)]
        synVoltSoma_df = dirList[2+(quad_set*4)]
        time_df = dirList[3+(quad_set*4)]
        
        time_df = pd.read_csv(dir_path + '/' + time_df, header=None)
        synVoltSIZ_df = pd.read_csv(dir_path + '/' + synVoltSIZ_df, header=None)
        synVoltSoma_df = pd.read_csv(dir_path + '/' + synVoltSoma_df, header=None)
        synSite_df = pd.read_csv(dir_path + '/' + synSite_df)
        
        quadList.append([time_df, synVoltSIZ_df, synVoltSoma_df, synSite_df])
    
    return quadList

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

def boxplot_2d(x,y, ax, whis=1.5, color = None):
    xlimits = [np.percentile(x, q) for q in (25, 50, 75)]
    ylimits = [np.percentile(y, q) for q in (25, 50, 75)]

    ##the box
    box = Rectangle(
        (xlimits[0],ylimits[0]),
        (xlimits[2]-xlimits[0]),
        (ylimits[2]-ylimits[0]),
        ec = color,
        zorder=0
    )
    ax.add_patch(box)

    ##the x median
    vline = Line2D(
        [xlimits[1],xlimits[1]],[ylimits[0],ylimits[2]],
        color=color,
        zorder=1
    )
    ax.add_line(vline)

    ##the y median
    hline = Line2D(
        [xlimits[0],xlimits[2]],[ylimits[1],ylimits[1]],
        color=color,
        zorder=1
    )
    ax.add_line(hline)

    ##the central point
    ax.plot([xlimits[1]],[ylimits[1]], color=color, marker='o')

    ##the x-whisker
    ##defined as in matplotlib boxplot:
    ##As a float, determines the reach of the whiskers to the beyond the
    ##first and third quartiles. In other words, where IQR is the
    ##interquartile range (Q3-Q1), the upper whisker will extend to
    ##last datum less than Q3 + whis*IQR). Similarly, the lower whisker
    ####will extend to the first datum greater than Q1 - whis*IQR. Beyond
    ##the whiskers, data are considered outliers and are plotted as
    ##individual points. Set this to an unreasonably high value to force
    ##the whiskers to show the min and max values. Alternatively, set this
    ##to an ascending sequence of percentile (e.g., [5, 95]) to set the
    ##whiskers at specific percentiles of the data. Finally, whis can
    ##be the string 'range' to force the whiskers to the min and max of
    ##the data.
    iqr = xlimits[2]-xlimits[0]

    ##left
    left = np.min(x[x > xlimits[0]-whis*iqr])
    whisker_line = Line2D(
        [left, xlimits[0]], [ylimits[1],ylimits[1]],
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_line)
    whisker_bar = Line2D(
        [left, left], [ylimits[0],ylimits[2]],
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_bar)

    ##right
    right = np.max(x[x < xlimits[2]+whis*iqr])
    whisker_line = Line2D(
        [right, xlimits[2]], [ylimits[1],ylimits[1]],
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_line)
    whisker_bar = Line2D(
        [right, right], [ylimits[0],ylimits[2]],
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_bar)

    ##the y-whisker
    iqr = ylimits[2]-ylimits[0]

    ##bottom
    bottom = np.min(y[y > ylimits[0]-whis*iqr])
    whisker_line = Line2D(
        [xlimits[1],xlimits[1]], [bottom, ylimits[0]], 
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_line)
    whisker_bar = Line2D(
        [xlimits[0],xlimits[2]], [bottom, bottom], 
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_bar)

    ##top
    top = np.max(y[y < ylimits[2]+whis*iqr])
    whisker_line = Line2D(
        [xlimits[1],xlimits[1]], [top, ylimits[2]], 
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_line)
    whisker_bar = Line2D(
        [xlimits[0],xlimits[2]], [top, top], 
        color = color,
        zorder = 1
    )
    ax.add_line(whisker_bar)

    ##outliers
    mask = (x<left)|(x>right)|(y<bottom)|(y>top)
    ax.scatter(
        x[mask],y[mask],
        facecolors=color, edgecolors='k'
    )


#########
#Plotting of single VPN activated synapses 

def plotallsynsexp2syn(neuron_name):
    if neuron_name == "DNp01":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synapse_sites.csv', dtype={'pre': str})
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_siz.csv', header=None)
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_soma.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_time.csv', header=None)

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
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synapse_sites.csv', dtype={'pre': str})
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_siz.csv', header=None)
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_soma.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_time.csv', header=None)
        
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

#Figure 9
def plotAllSynsAndhighlightsynsFig9(neuron_name, erev):
    # erev= -66.6544
    erev = -66.4298 - 0.2
    if neuron_name == "DNp01":
        synSite_df = pd.read_csv("datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synapse_sites.csv", dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_time.csv', header=None)

    elif neuron_name == "DNp03":
        synSite_df = pd.read_csv("datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synapse_sites.csv", dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_time.csv', header=None)
        
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

    SIZMAX = max(synCurr_siz_df.max(axis=1))
    
    #NEED THIS ONE
    TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=None, gridPos=[2,1], axList=None)
    # print(synSite_df.head())
    TESTWINDOW = plotSynapseLocations(synSite_df, colorKey=1, sizeArg=6)
    # clock.sleep(1200)
    if neuron_name == "DNp01":
        # RSYN_IDX_1 = 385
        # RSYN_IDX_2 = 5119
        # RSYN_IDX_3 = 856
        RSYN_IDX_1 = 501
        RSYN_IDX_2 = 39
        RSYN_IDX_3 = 481
        RSYN_IDX_4 = 1053
        RSYN_IDX_5 = 670
        RSYN_IDX_6 = 746
    elif neuron_name == "DNp02":
        RSYN_IDX_1 = 3855
        RSYN_IDX_2 = 3834
        RSYN_IDX_3 = 2852
    elif neuron_name == "DNp03":
        RSYN_IDX_1 = 5
        RSYN_IDX_2 = 150
        #RSYN_IDX_3 = 400
        RSYN_IDX_4 = 530
        RSYN_IDX_5 = 535
        #RSYN_IDX_6 = 555
        RSYN_IDX_7 = 920
        RSYN_IDX_8 = 1000
        #RSYN_IDX_9 = 1231
        RSYN_IDX_10 = 1828
        RSYN_IDX_11 = 1832
        #RSYN_IDX_12 = 1842
        RSYN_IDX_13 = 1948
        RSYN_IDX_14 = 1985
        #RSYN_IDX_15 = 2200
    elif neuron_name == "DNp04":
        RSYN_IDX_1 = 292
        RSYN_IDX_2 = 153
        RSYN_IDX_3 = 775
    elif neuron_name == "DNp06":
        RSYN_IDX_1 = 2695
        RSYN_IDX_2 = 15276
        RSYN_IDX_3 = 525
    else:
        RSYN_IDX_1 = random.randint(0, synSite_df.shape[0])
        RSYN_IDX_2 = random.randint(0, synSite_df.shape[0])
        RSYN_IDX_3 = random.randint(0, synSite_df.shape[0])

    if neuron_name == "DNp01":
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_1], colorKey=3, sizeArg=4)#syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_2], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=4)
        #TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_3], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=4)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_4], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_5], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)
        #TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_6], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=4)

        TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_siz_df.iloc[RSYN_IDX_4], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_siz_df.iloc[RSYN_IDX_5], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_6], synCurr_siz_df.iloc[RSYN_IDX_6], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_siz_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_siz_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_3], synCurr_siz_df.iloc[RSYN_IDX_3], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df, synCurr_synpt_df, '#899499', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_synpt_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_synpt_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_3], synCurr_synpt_df.iloc[RSYN_IDX_3], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_synpt_df.iloc[RSYN_IDX_4], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_synpt_df.iloc[RSYN_IDX_5], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_6], synCurr_synpt_df.iloc[RSYN_IDX_6], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF_AXLIST[1].set_xlim([47.5, 60])
        TF_AXLIST[1].set_ylim([erev-0.001, SIZMAX+0.02])
        TF_AXLIST[0].set_xlim([47.5, 60])
        TF_AXLIST[0].set_ylim([erev-0.001, SYNPTMAX+0.02])
        TF_AXLIST[0].set_title('Recording at Synapse Site')
        TF_AXLIST[1].set_title('Recording at SIZ')

        TF_AXLIST[0].set_xlabel('Time (ms)')
        TF_AXLIST[1].set_xlabel('Time (ms)')
        TF_AXLIST[0].set_ylabel('Voltage (mV)')
        TF_AXLIST[1].set_ylabel('Voltage (mV)')


    elif neuron_name =="DNp03":
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_1], colorKey=3, sizeArg=6)#syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_2], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=6)
        #TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_3], colorKey=3, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_4], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_5], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        #TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_6], colorKey=8, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_7], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_8], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        #TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_9], colorKey=2, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_10], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_11], colorKey=5, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_13], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_14], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        #TESTWINDOW = plotSynapseLocations(synSite_df.iloc[RSYN_IDX_12], colorKey=4, syn_shape_window=TESTWINDOW, sizeArg=6)
        TF, TF_AXLIST = plotTraces(time_df, synCurr_siz_df, '#899499', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_siz_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_siz_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_3], synCurr_siz_df.iloc[RSYN_IDX_3], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_siz_df.iloc[RSYN_IDX_4], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_siz_df.iloc[RSYN_IDX_5], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_6], synCurr_siz_df.iloc[RSYN_IDX_6], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_7], synCurr_siz_df.iloc[RSYN_IDX_7], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_8], synCurr_siz_df.iloc[RSYN_IDX_8], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_9], synCurr_siz_df.iloc[RSYN_IDX_9], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_siz_df.iloc[RSYN_IDX_10], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_siz_df.iloc[RSYN_IDX_11], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_siz_df.iloc[RSYN_IDX_13], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_siz_df.iloc[RSYN_IDX_14], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_12], synCurr_siz_df.iloc[RSYN_IDX_12], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 1], axList=TF_AXLIST)
        
        TF, TF_AXLIST = plotTraces(time_df, synCurr_synpt_df, '#899499', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_1], synCurr_synpt_df.iloc[RSYN_IDX_1], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_2], synCurr_synpt_df.iloc[RSYN_IDX_2], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_3], synCurr_synpt_df.iloc[RSYN_IDX_3], colorCode='#0000FF', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_4], synCurr_synpt_df.iloc[RSYN_IDX_4], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_5], synCurr_synpt_df.iloc[RSYN_IDX_5], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_6], synCurr_synpt_df.iloc[RSYN_IDX_6], colorCode='#FFF000', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_7], synCurr_synpt_df.iloc[RSYN_IDX_7], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_8], synCurr_synpt_df.iloc[RSYN_IDX_8], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_9], synCurr_synpt_df.iloc[RSYN_IDX_9], colorCode='#00FF00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_synpt_df.iloc[RSYN_IDX_10], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_synpt_df.iloc[RSYN_IDX_11], colorCode='#FFA500', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_10], synCurr_synpt_df.iloc[RSYN_IDX_13], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_11], synCurr_synpt_df.iloc[RSYN_IDX_14], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)
        #TF, TF_AXLIST = plotTraces(time_df.iloc[RSYN_IDX_12], synCurr_synpt_df.iloc[RSYN_IDX_12], colorCode='#FF0D00', FIGURE=TF, gridPos=[2, 0], axList=TF_AXLIST)

        TF_AXLIST[1].set_xlim([47.5, 60])
        TF_AXLIST[1].set_ylim([erev+5.7, SIZMAX+0.04])
        TF_AXLIST[0].set_xlim([47.5, 60])
        TF_AXLIST[0].set_ylim([erev+5.5, SYNPTMAX+0.04])
        TF_AXLIST[0].set_title('Recording at Synapse Site')
        TF_AXLIST[1].set_title('Recording at SIZ')

        TF_AXLIST[0].set_xlabel('Time (ms)')
        TF_AXLIST[1].set_xlabel('Time (ms)')
        TF_AXLIST[0].set_ylabel('Voltage (mV)')
        TF_AXLIST[1].set_ylabel('Voltage (mV)')
    # TESTWINDOW.point_mark(sizSection(0.5), 2, "O", 10)
    # TEST_XORG = somaSection.x3d(1)
    # TEST_YORG = somaSection.y3d(1)
    # TEST_ZORG = somaSection.z3d(1)
    # TESTWINDOW.rotate(TEST_XORG, TEST_YORG, TEST_ZORG, 0.5, 0.5, 0.5)

    #F4BB44 #899499

    # TF.savefig('datafiles/figures/allSyn_628/'+neuron_name+'_allSyn/'+neuron_name+'_allTraces_'+str(RSYN_IDX_1)+'_'+str(RSYN_IDX_2)+'_'+str(RSYN_IDX_3)+'_7623.svgz', bbox_inches='tight')
    # TESTWINDOW.printfile('datafiles/figures/morphVPN_actSingle_forAMS_NBioDM_'+neuron_name+'_synLocations_'+str(RSYN_IDX_1)+'_'+str(RSYN_IDX_2)+'_'+str(RSYN_IDX_3)+'.cdr')

    fig = plt.figure()
    peakDepol_synpt = synCurr_synpt_df.max(axis=1)
    peakDepol_siz = synCurr_siz_df.max(axis=1)
    plt.plot(peakDepol_synpt.index, peakDepol_synpt - erev, 'go')
    plt.plot(1150, numpy.mean(peakDepol_synpt) - erev, 'ro', markersize=10)
    print('\n\navg depol is', numpy.mean(peakDepol_synpt) - erev-5.5, 'mV')
    plt.plot(peakDepol_siz.index, peakDepol_siz - erev-5., 'ko')
    plt.ylabel('Magnitude of depolarization (mV)')
    plt.show()

#Inset of Figure 9
def plot_peak_depol_at_SIZ_vs_dist_to_SIZ(neuron_name, sizSection):
    if neuron_name == "DNp01":
        synSite_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synapse_sites.csv', dtype={'pre': str})
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df= pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_time.csv', header=None)

    if neuron_name == "DNp03":
        synSite_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synapse_sites.csv', dtype={'pre': str})
        synCurr_soma__df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_soma.csv', header=None)
        synCurr_siz_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_synpt.csv', header=None)
        time_df = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_time.csv', header=None)
    gs_fig = gridspec.GridSpec(2, 2)
    TEST_IND = plt.figure()
    ax1_fig = TEST_IND.add_subplot(gs_fig[0, :])
    ax2_fig = TEST_IND.add_subplot(gs_fig[1, :])

    dist_to_SIZ =[]
    peak_depol_at_SIZ = []
    for index, row in synSite_df.iterrows():
        matching_rows = synSite_df[(synSite_df['mappedSection'] == row['mappedSection']) & (synSite_df['mappedSegRangeVar'] == row['mappedSegRangeVar'])]
        index = matching_rows.index[0]
        rest_SIZ = np.array(synCurr_siz_df.iloc[index][2002])
        rest_syn = np.array(synCurr_synpt_df.iloc[index][2002])
        synSec = str2sec(row['mappedSection'])
        distToSIZ = h.distance(sizSection(0.05), synSec(row['mappedSegRangeVar']))
        dist_to_SIZ.append(distToSIZ)
        peak_depol_syn = numpy.max(numpy.array(synCurr_synpt_df.iloc[index]))- rest_syn
        peak_depol_SIZ = numpy.max(numpy.array(synCurr_siz_df.iloc[index]))- rest_SIZ
        peak_depol_at_SIZ.append(peak_depol_SIZ)
        ax1_fig.plot(distToSIZ, peak_depol_syn, 'ko')
        ax2_fig.plot(distToSIZ, peak_depol_SIZ, 'ko')
    ax1_fig.set_xlabel('Distance from Synapse to SIZ (um)')
    ax2_fig.set_xlabel('Distance from Synapse to SIZ (um)')
    ax1_fig.set_ylabel('Depolarization (mV)')
    ax2_fig.set_ylabel('Depolarization (mV)')
    ax1_fig.set_title('Individual Synapse Depolarizations at Synapse Location ({})'.format(neuron_name))
    ax2_fig.set_title('Individual Synapse Depolarizations at SIZ ({})'.format(neuron_name))

    coeff = numpy.polyfit(dist_to_SIZ, peak_depol_at_SIZ, 1)
    line = numpy.polyval(coeff, dist_to_SIZ)
    ax2_fig.plot(dist_to_SIZ, line, '-', color='red')
    slope, intercept, r_value, p_value, std_err = linregress(dist_to_SIZ, peak_depol_at_SIZ)

    # Calculate R^2
    R_squared = r_value ** 2

    print("R^2 value:", R_squared)
    print("p value:", p_value)

    TEST_IND.tight_layout()
    plt.show()

#Figure 10
def plot_partner_VPN_activations(neuron_name, MODE= 'all/'):
    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    quintList = returnQuints(dir_path)
    fig = plt.figure()

    LC4_0 = True
    LC6_0 = True
    LC22_0 = True
    LPLC1_0 = True
    LPLC2_0 = True
    LPLC4_0 = True
    Rand_0 = True

    for quint_ct in range(len(quintList)):
        quint = quintList[quint_ct]
        
        vSIZ_np = numpy.array(quint[1])
        peakV_siz = numpy.max(vSIZ_np)
        Vrest_siz = vSIZ_np[0, 2002]
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

    plt.gcf().get_axes()[0].spines['top'].set_visible(False)
    plt.gcf().get_axes()[0].spines['right'].set_visible(False)
    if neuron_name == "DNp01":
        plt.xlim([0, 26])
    elif neuron_name == "DNp03":
        plt.xlim([0, 40])
    plt.ylabel("Peak Depolarization at SIZ (mV)")
    plt.xlabel("Number of Synapses")
    plt.title("Recording at SIZ")
    plt.suptitle("Voltage Response of {} to Activations of Partnered VPN Synapses".format(neuron_name))
    plt.tight_layout()
    plt.legend()
    plt.show()

def plot_rand_vs_partner_SIZ(neuron_name, MODE='all/'):
    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    
    quintList = returnQuints(dir_path)
    fig = plt.figure()

    dataList = []

    LC4_0 = True
    LC6_0 = True
    LC22_0 = True
    LPLC1_0 = True
    LPLC2_0 = True
    LPLC4_0 = True

    Rand_0 = True

    for quint_ct in range(len(quintList)):
        quint = quintList[quint_ct]
        vSIZ_np = numpy.array(quint[1])
        peakV_siz = numpy.max(vSIZ_np)
        Vrest_siz = vSIZ_np[0, 2002]
        synct = quint[4].shape
        synct = synct[0]
        #TW = plotSynapseLocations(quint[4], colorKey=2, sizeArg=4, syn_shape_window=None)
        JIT = numpy.random.default_rng().normal(0, 0.1, 1)
        # fig = plotTraces(time_df, vSIZ, colorCode='k', FIGURE=fig)
        
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

    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_rand'
    if neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_rand'
    quintList = returnQuints(dir_path)
    
    for quint_ct in range(len(quintList)):
        quint = quintList[quint_ct]

        vSIZ_np = numpy.array(quint[1])
        peakV_siz = numpy.max(vSIZ_np)
        Vrest_siz = vSIZ_np[0, 2002]
        synct = quint[4].shape
        synct = synct[0]
        if synct >= 40:
            break

        JIT = numpy.random.default_rng().normal(0, 0.1, 1)
        if Rand_0:
            plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'ko', label='Random')#, markersize=2)
            Rand_0 = False
        else:
            plt.plot(synct+JIT+0.25, peakV_siz-Vrest_siz, 'ko')

    plt.gcf().get_axes()[0].spines['top'].set_visible(False)
    plt.gcf().get_axes()[0].spines['right'].set_visible(False)
    if neuron_name == "DNp01":
        plt.xlim([0, 25])
    elif neuron_name == "DNp03":
        plt.xlim([0, 45])
    plt.ylabel("Peak Depolarization at SIZ (mV)")
    plt.xlabel("Number of Synapses")
    plt.title("Recording at SIZ")
    plt.suptitle("Voltage Response of {} to Activations of Partnered vs Random Synapses".format(neuron_name))
    plt.tight_layout()
    plt.legend()
    plt.show()

def plot_rand_vs_partner_SOMA(neuron_name, MODE='all/'):
    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_partner'
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
        Vrest_soma = vSoma_np[0, 2002]
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
        Vrest_soma = vSoma_np[0, 2002]
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

#Figure 11
def Plot_DN_VPN_syns(neuron_name, erev, location = None):
    dir_path_lin30VPN = 'datafiles/simulationData/'+neuron_name+'_incRandSynCt_VPN'

    quintList_lin30VPN = returnQuints(dir_path_lin30VPN)


    lin30VPN_fig = plt.figure()
    for quint_ct in range(len(quintList_lin30VPN)):
        print('plotting trial', quint_ct, '/', len(quintList_lin30VPN))
        quint = quintList_lin30VPN[quint_ct]
        
        vSIZ_np = numpy.array(quint[1])
        peakV_siz = numpy.max(vSIZ_np)
        vSoma_np = numpy.array(quint[2])
        peakV_soma = numpy.max(vSoma_np)
        Vrest_siz = vSIZ_np[0, 2002]
        Vrest_soma = vSoma_np[0, 2002]
        synct = quint[4].shape
        if quint[4].loc[0]['type'] == 'na':
            colorCode = '#000000'
        elif quint[4].loc[0]['type'] == 'LC4':
            colorCode = '#0000FF'
        elif quint[4].loc[0]['type'] == 'LPLC2':
            colorCode = '#FFA500'
        elif quint[4].loc[0]['type'] == 'LPLC1':
            colorCode = '#00FF00'
        elif quint[4].loc[0]['type'] == 'LC22':
            colorCode = '#FFF000'
        elif quint[4].loc[0]['type'] == 'LPLC4':
            colorCode = '#FF0D00'
        elif quint[4].loc[0]['type'] == 'LC6':
            colorCode = '#FF00E9'
        else:
            print("ERROR:"+quint[4].loc[0]['type'])
            break



        synct = synct[0]
        if location =="Soma":
            plt.plot(synct, peakV_soma-Vrest_soma, colorCode, marker='o')
        elif location == "SIZ":
            plt.plot(synct, peakV_siz-Vrest_siz, colorCode, marker='o')
        elif location == "Both":
             plt.plot(synct, peakV_soma-Vrest_soma, colorCode, marker='o')
             plt.plot(synct, peakV_siz-Vrest_siz, colorCode, marker='s')

    plt.xlabel('Number of synapses activated')
    plt.ylabel('Voltage deflection from resting membrane potential (mV)')
    plt.title('Depolarization of {} as a Function of Number of Synapses Activated'.format(neuron_name))
    plt.show()

#This is a work around, in case of memory issues with the above, but requires the manual creation of the folders with correct data seperated by LC/LPLC types
def Plot_DNp03_VPN_syns_by_type(neuron_name, location = None, type = None):
    if type == "LC":
        dir_path_lin30VPN = 'datafiles/simulationData/'+neuron_name+'_incRandSynCt650_VPNnoJitter_FIX2'
    elif type == "LPLC":
        dir_path_lin30VPN = 'datafiles/simulationData/'+neuron_name+'_incRandSynCt650_VPNnoJitter_FIX3'

    quintList_lin30VPN = returnQuints(dir_path_lin30VPN)


    lin30VPN_fig = plt.figure()
    for quint_ct in range(len(quintList_lin30VPN)):
        print('plotting trial', quint_ct, '/', len(quintList_lin30VPN))
        quint = quintList_lin30VPN[quint_ct]
        
        vSIZ_np = numpy.array(quint[1])
        peakV_siz = numpy.max(vSIZ_np)
        vSoma_np = numpy.array(quint[2])
        peakV_soma = numpy.max(vSoma_np)
        Vrest_siz = vSIZ_np[0, 2002]
        Vrest_soma = vSoma_np[0, 2002]
        synct = quint[4].shape
        if quint[4].loc[0]['type'] == 'na':
            colorCode = '#000000'
        elif quint[4].loc[0]['type'] == 'LC4':
            colorCode = '#0000FF'
        elif quint[4].loc[0]['type'] == 'LPLC2':
            colorCode = '#FFA500'
        elif quint[4].loc[0]['type'] == 'LPLC1':
            colorCode = '#00FF00'
        elif quint[4].loc[0]['type'] == 'LC22':
            colorCode = '#FFF000'
        elif quint[4].loc[0]['type'] == 'LPLC4':
            colorCode = '#FF0D00'
        elif quint[4].loc[0]['type'] == 'LC6':
            colorCode = '#FF00E9'
        else:
            print("ERROR:"+quint[4].loc[0]['type'])
            break



        synct = synct[0]
        if location =="Soma":
            plt.plot(synct, peakV_soma-Vrest_soma, colorCode, marker='o')
        elif location == "SIZ":
            plt.plot(synct, peakV_siz-Vrest_siz, colorCode, marker='o')
        elif location == "Both":
             plt.plot(synct, peakV_soma-Vrest_soma, colorCode, marker='o')
             plt.plot(synct, peakV_siz-Vrest_siz, colorCode, marker='s')
    plt.xlim(0, 1200)
    plt.ylim(0, 40)
    plt.xlabel('Number of synapses activated')
    plt.ylabel('Voltage deflection from resting membrane potential (mV)')
    plt.title('Depolarization of {} as a Function of Number of Synapses Activated'.format(neuron_name))
    plt.show()

#Figure 12
#Must run this first before running the below, and run partner and random seperately.
def plot_syn_spread_vs_dist_to_siz_partner_vs_rand(neuron_name, erev, sizSection=None, plots = None):
    if plots == "Rand":
        dir_path = 'datafiles/simulationData/'+neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+'_rand'
        quintList = returnQuints(dir_path)

    elif plots =="Partner":
        dir_path = 'datafiles/simulationData/'+neuron_name+'_final_sims/rand_vs_partner/'+neuron_name+'_partner'
        quintList = returnQuints(dir_path)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))
    synct_list = []
    peakV_siz_list= []
    peakV_soma_list =[]
    avg_dist_list = []
    avg_syn_spread_list = []
    synTypeList = []
    synpreList = []
    k = 0

    for quint_ct, quint in enumerate(quintList):
            k+=1
            print(k)
            synSegRangeVarList = []
            synSecList = []
            synSite = quint[4]
            vSIZ_np = numpy.array(quint[1])
            peakV_siz = numpy.max(vSIZ_np)
            Vrest_siz = vSIZ_np[0, 2002]
            vsoma_np = numpy.array(quint[2])
            peakV_soma = numpy.max(vsoma_np)
            Vrest_soma = vsoma_np[0, 2002]

            synct = quint[4].shape[0]
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
    dir_path = 'datafiles/simulationData'
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
        quintList = returnQuints(dir_path_rand)
    elif plots =="Close":
        dir_path_close = 'datafiles/simulationData/'+neuron_name+f'_close{synnum}record{trials}'
        quintList = returnQuints(dir_path_close)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))
    synct_list = []
    peakV_siz_list= []
    avg_dist_list = []
    avg_syn_spread_list = []
    
    synTypeList = []
    synpreList = []
    k = 0

    for quint_ct, quint in enumerate(quintList):
            synSegRangeVarList = []
            synSecList = []
            k+=1
            print(k)
            synSite = quint[4]
            vSIZ_np = numpy.array(quint[1])
            peakV_siz = numpy.max(vSIZ_np)
            Vrest_siz = vSIZ_np[0, 2002]
            synct = quint[4].shape[0]
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

    plt.savefig(f'{plots}_syn_spread_vs_dist_to_siz_vs_SIZ_depol_{trials}_trials_{synnum}_synapses.png')
    plt.show()

def plot_syn_spread_vs_dist_to_siz_from_csv_rand_vs_partner(neuron_name, synnum = None, trials = None):
    dir_path = 'datafiles/simulationData'
    rand_df = pd.read_csv(dir_path + '/'+neuron_name+'_'+f'Rand_syn_spread_vs_dist_to_siz_data_{trials}_trials_{synnum}_synapses.csv')
    close_df = pd.read_csv(dir_path + '/'+neuron_name+'_'+f'Close_syn_spread_vs_dist_to_siz_data_{trials}_trials_{synnum}_synapses.csv')
    partner_df = pd.read_csv(dir_path + '/'+neuron_name+'_'+f'Partner_syn_spread_vs_dist_to_siz_data.csv')
    if synnum != None:
        subset = synnum
        partner_df = partner_df.loc[partner_df['Synct'] == subset]
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15))
    print(rand_df)
    print(close_df)
    print(partner_df)

    subset_close = close_df
    subset_rand = rand_df
    subset_partner = partner_df
    ax1.plot(subset_rand['Avg_syn_spread'], subset_rand['Peak_Depol_SIZ'], 's', color='black', label = "Rand")
    ax1.plot(subset_close['Avg_syn_spread'], subset_close['Peak_Depol_SIZ'], 'o', color='red', label = "Close")
    ax1.plot(subset_partner['Avg_syn_spread'], subset_partner['Peak_Depol_SIZ'], '^', color='blue', label = "Partner")
    # ax2.plot(subset_rand['avg_dist_to_SIZ'], subset_rand['Avg_syn_spread'], 's', color='black', label = "Rand")
    # ax2.plot(subset_close['avg_dist_to_SIZ'], subset_close['Avg_syn_spread'], 'o', color='red',label = "Close")
    # ax2.plot(subset_partner['avg_dist_to_SIZ'], subset_partner['Avg_syn_spread'], '^', color='blue',label = "Partner")
    ax3.plot(subset_rand['avg_dist_to_SIZ'], subset_rand['Peak_Depol_SIZ'], 's', color='black', label = "Rand")
    ax3.plot(subset_close['avg_dist_to_SIZ'], subset_close['Peak_Depol_SIZ'], 'o', color='red', label = "Close")
    ax3.plot(subset_partner['avg_dist_to_SIZ'], subset_partner['Peak_Depol_SIZ'], '^', color='blue', label = "Partner")
    
    ax1.set_xlabel("Average synapse spread")
    ax1.set_ylabel("Peak depolarization at SIZ")
    ax1.legend()

    ax2.set_xlabel("Average synapse spread")
    ax2.set_ylabel("Peak depolarization at SIZ")
    ax2.legend()

    coeff_rand = numpy.polyfit(rand_df['avg_dist_to_SIZ'], rand_df['Peak_Depol_SIZ'], 1)
    line_rand = numpy.polyval(coeff_rand, rand_df['avg_dist_to_SIZ'])

    coeff_close = numpy.polyfit(close_df['avg_dist_to_SIZ'], close_df['Peak_Depol_SIZ'], 1)
    line_close = numpy.polyval(coeff_close, close_df['avg_dist_to_SIZ'])

    ax3.set_xlabel("Average distance of synapses to SIZ")
    ax3.set_ylabel("Peak depolarization at SIZ")
    ax3.legend()

    # Assuming rand_df and close_df have the same length
    x_values_rand = rand_df['avg_dist_to_SIZ']
    y_values_rand = rand_df['Peak_Depol_SIZ']
    x_values_rand_syn_spread = rand_df['Avg_syn_spread']

    x_values_close = close_df['avg_dist_to_SIZ']
    y_values_close = close_df['Peak_Depol_SIZ']
    x_values_close_syn_spread = close_df['Avg_syn_spread']

    x_values_partner = partner_df['avg_dist_to_SIZ']
    y_values_partner = partner_df['Peak_Depol_SIZ']
    x_values_partner_syn_spread = partner_df['Avg_syn_spread']

    # Fitting a linear line to the data
    coeff_rand = numpy.polyfit(x_values_rand, y_values_rand, 1)
    line_rand = numpy.polyval(coeff_rand, x_values_rand)

    coeff_close = numpy.polyfit(x_values_close, y_values_close, 1)
    line_close = numpy.polyval(coeff_close, x_values_close)

    coeff_partner = numpy.polyfit(x_values_partner, y_values_partner, 1)
    line_partner = numpy.polyval(coeff_partner, x_values_partner)

    # Extrapolate the lines
    x_values_extrapolate = numpy.linspace(min(min(x_values_rand), min(x_values_close), min(x_values_partner)),
                                    max(max(x_values_rand), max(x_values_close), max(x_values_partner)), 50)

    line_rand_extrapolated = numpy.polyval(coeff_rand, x_values_extrapolate)
    line_close_extrapolated = numpy.polyval(coeff_close, x_values_extrapolate)
    line_partner_extrapolated = numpy.polyval(coeff_partner, x_values_extrapolate)
    

    ax3.plot(x_values_extrapolate, line_rand_extrapolated, '-', color='black', label = "Rand")
    ax3.plot(x_values_extrapolate, line_close_extrapolated, '-', color='red', label = "Close")
    ax3.plot(x_values_extrapolate, line_partner_extrapolated, '-', color='blue', label = "Partner")

    plt.figure(figsize=(8,6))
    plt.boxplot([rand_df['Peak_Depol_SIZ'], close_df['Peak_Depol_SIZ'], partner_df['Peak_Depol_SIZ']], labels=['Random Data', 'Close Data', "Single VPN"])
    plt.title('Box Plot of Peak_Depol_SIZ')
    plt.ylabel('Peak_Depol_SIZ')

    plt.figure(figsize=(8,6))
    plt.boxplot([rand_df['Avg_syn_spread'], close_df['Avg_syn_spread'], partner_df['Avg_syn_spread']], labels=['Random Data', 'Close Data', "Single VPN"])
    plt.title('Box Plot of Synapse spread')
    plt.ylabel('Avg_syn_spread')

    #doing the box plot
    boxplot_2d(x_values_rand_syn_spread,y_values_rand,ax=ax2, whis=1, color = 'black')
    boxplot_2d(x_values_close_syn_spread,y_values_close,ax=ax2, whis=1, color = 'red')
    boxplot_2d(x_values_partner_syn_spread,y_values_partner,ax=ax2, whis=1, color = 'blue')

    statistic_depol, p_value_depol = levene(rand_df['Peak_Depol_SIZ'], close_df['Peak_Depol_SIZ'])

    # Perform Levene's test for Avg_syn_spread
    statistic_spread, p_value_spread = levene(rand_df['Avg_syn_spread'], close_df['Avg_syn_spread'])

    statistic_dist_to_SIZ, p_value_dist = levene(rand_df['avg_dist_to_SIZ'], partner_df['avg_dist_to_SIZ'])

    print("Levene's test for Peak_Depol_SIZ (rand vs close):")
    print("Statistic:", statistic_depol)
    print("p-value:", p_value_depol)

    print("\nLevene's test for Avg_syn_spread (rand vs close):")
    print("Statistic:", statistic_spread)
    print("p-value:", p_value_spread)

    print("\nLevene's test for SIZ Distance (rand vs partner):")
    print("Statistic:", statistic_dist_to_SIZ)
    print("p-value:", p_value_dist)

    range_dist_rand = rand_df['avg_dist_to_SIZ'].max() - rand_df['avg_dist_to_SIZ'].min()
    range_dist_close = close_df['avg_dist_to_SIZ'].max() - close_df['avg_dist_to_SIZ'].min()

    # Print minimum and maximum values along with ranges
    print("\nRange for avg_dist_to_SIZ - Random Group:")
    print("Minimum:", rand_df['avg_dist_to_SIZ'].min())
    print("Maximum:", rand_df['avg_dist_to_SIZ'].max())
    print("Range:", range_dist_rand)

    print("\nRange for avg_dist_to_SIZ - Close Group:")
    print("Minimum:", close_df['avg_dist_to_SIZ'].min())
    print("Maximum:", close_df['avg_dist_to_SIZ'].max())
    print("Range:", range_dist_close)

    # Perform linear regression for rand
    slope_rand, intercept_rand, r_value_rand, p_rand, stderr_rand = linregress(x_values_rand, y_values_rand)

    # Perform linear regression for close
    slope_close, intercept_close, r_value_close, p_close, stderr_close = linregress(x_values_close, y_values_close)

    # Perform linear regression for VPNs
    slope_partner, intercept_partner, r_value_partner, p_partner, stderr_partner = linregress(x_values_partner, y_values_partner)

    print(f"R^2 value rand {r_value_rand**2}")
    print(f"slope p value rand {p_rand}")
    print(f"R^2 value close {r_value_close**2}")
    print(f"slope p value close {p_close}")
    print(f"R^2 value partner {r_value_partner**2}")
    print(f"slope p value partner{p_partner}")
    plt.show()

    close_normality_p_value_peak_depol = shapiro(close_df['Peak_Depol_SIZ'])[1]
    rand_normality_p_value_peak_depol = shapiro(rand_df['Peak_Depol_SIZ'])[1]
    if len(partner_df['Peak_Depol_SIZ']) >= 3:
        partner_normality_p_value_peak_depol = shapiro(partner_df['Peak_Depol_SIZ'])[1]
    else:
        pass

    print("Shapiro wilk test for normality p values peak depol p values")
    print(f'close normality p value{close_normality_p_value_peak_depol}')
    print(f'random normality p value{rand_normality_p_value_peak_depol}')
    if len(partner_df['Peak_Depol_SIZ']) >= 3:
        print(f'partner normality p value{partner_normality_p_value_peak_depol}')
    else:
        pass

    close_normality_p_value_peak_depol = shapiro(close_df['Avg_syn_spread'])[1]
    rand_normality_p_value_peak_depol = shapiro(rand_df['Avg_syn_spread'])[1]
    if len(partner_df['Avg_syn_spread']) >= 3:
        partner_normality_p_value_peak_depol = shapiro(partner_df['Avg_syn_spread'])[1]
    else:
        pass

    print("Shapiro wilk test for normality p values avg syn spread")
    print(f'close normality p value{close_normality_p_value_peak_depol}')
    print(f'random normality p value{rand_normality_p_value_peak_depol}')
    if len(partner_df['Avg_syn_spread']) >= 3:
        print(f'partner normality p value{partner_normality_p_value_peak_depol}')
    else:
        pass

    peak_depol = [close_df['Peak_Depol_SIZ'], partner_df['Peak_Depol_SIZ'], rand_df['Peak_Depol_SIZ']]
    syn_spread = [close_df['Avg_syn_spread'], partner_df['Avg_syn_spread'], rand_df['Avg_syn_spread']]

    statistic, p_value_peak_depol = stats.kruskal(peak_depol[0], peak_depol[1], peak_depol[2])
    statistic2, p_value_syn_spread = stats.kruskal(syn_spread[0], syn_spread[1], syn_spread[2])

    # Step 2: Post Hoc Analysis (Dunn's Test)
    if p_value_peak_depol < 0.05:
        # If Kruskal-Wallis test is significant, perform post hoc analysis
        print("Peak depolarizations bonferroni adjusted p values")
        posthoc_res = posthoc_dunn(peak_depol, p_adjust='bonferroni')  # Adjust p-values for multiple comparisons
        print(posthoc_res)
    else:
        print("No significant differences between groups according to Kruskal-Wallis test.")

        # Step 2: Post Hoc Analysis (Dunn's Test)
    if p_value_syn_spread < 0.05:
        # If Kruskal-Wallis test is significant, perform post hoc analysis
        print("Synapse spread bonferroni adjusted p values")
        posthoc_res = posthoc_dunn(syn_spread, p_adjust='bonferroni')  # Adjust p-values for multiple comparisons
        print(posthoc_res)
    else:
        print("No significant differences between groups according to Kruskal-Wallis test.")

######
#Analysis of depolarizations and integration of multiple synapses comparing actual VPN activations with the individual linear sum of its synapses at SIZ

def Peak_depol_per_VPN(neuron_name):
    if neuron_name == "DNp01":
        dir_path = 'datafiles/simulationData/DNp01_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    elif neuron_name == "DNp03":
        dir_path = 'datafiles/simulationData/DNp03_final_sims/rand_vs_partner/'+neuron_name+'_partner'
    
    quintList = returnQuints(dir_path)
    data_list = []
    for quint_ct in range(len(quintList)):
        quint = quintList[quint_ct]
        vSIZ_np = numpy.array(quint[1])
        peakV_siz = numpy.max(vSIZ_np)
        Vrest_siz = vSIZ_np[0, 2002]
        synct = quint[4].shape
        synct = synct[0]
        if quint[4]['type'].iloc[0] in ["LC4", "LC6", "LC22", "LPLC2", "LPLC1", "LPLC4"]:
            id = (quint[4]['pre']).iloc[0]
            type = quint[4]['type'].iloc[0]
            depol = peakV_siz - Vrest_siz
            avg_syn_depol = depol / synct
            data_list.append({'pre': id, 'type': type, 'synct': synct, 'depol_together': depol, 'avg_syn_depol': avg_syn_depol})

    df = pd.DataFrame(data_list)

    # Sort the DataFrame by 'type' and then by 'synct'
    df_sorted = df.sort_values(by=['type', 'synct'])

    # Print the sorted DataFrame
    return df_sorted

def individual_syn_depol_per_VPN(neuron_name):
    if neuron_name == "DNp01":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synapse_sites.csv', dtype={'pre': str})
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/ssinglesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_siz.csv', header=None)
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_soma.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp01_final_sims/singlesynactivation/SINGLE_EXP2_20240214-181219_time.csv', header=None)
        result_list = []
    
    elif neuron_name == "DNp03":
        synSite_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synapse_sites.csv', dtype={'pre': str})
        synCurr_soma_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_soma.csv', header=None)
        synCurr_siz_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_siz.csv', header=None)
        synCurr_synpt_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_synaptic_currents_synpt.csv', header=None)
        time_df_exp2 = pd.read_csv('datafiles/simulationData/DNp03_final_sims/singlesynactivation/SINGLE_EXP2_20240214-123501_time.csv', header=None)
        result_list = []

    # Iterate through each row
    for index, row in synCurr_siz_df_exp2.iterrows():
        # Find the maximum value in the row
        max_value = row.max()

        # Get the value from the specified position 
        column_value = row.iloc[2002]  
        
        # Subtract the max value from the value at the specified position for the current row
        result = column_value - max_value
        
        # Append the result to the list
        result_list.append(result * -1)

    # Append the result list as a new column to the synSite_df_exp2 DataFrame
    synSite_df_exp2['summed_syn_depol'] = result_list

    sum_syn_depol = synSite_df_exp2.groupby('pre')['summed_syn_depol'].sum().reset_index()
    sum_syn_depol['syn_count'] = synSite_df_exp2.groupby('pre').size().reset_index(name='syn_count')['syn_count']

    # Creating a new DataFrame with the sum of 'syn_depol' and 'syn_count' for each 'pre' value
    new_df = pd.DataFrame(sum_syn_depol)

    return new_df

def combine_indiv_syn_depol_and_group_act(neuron_name):
    df = Peak_depol_per_VPN(neuron_name)
    df2 = individual_syn_depol_per_VPN(neuron_name)
    df['pre'] = df['pre'].astype(str)
    merged_df = pd.merge(df, df2, on='pre', how='inner')
    unique_types = merged_df['type'].unique()
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

    type_colors = {'LC4': 'blue', 'LPLC2': 'orange', 'LPLC1': 'red', 'LPLC4': 'green', 'LC22': 'yellow'}  # Add more colors as needed

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

######
######


#############
    #Retinotopy set up and analysis Figure 3 Supplemental 1, Figure 4D-E, Figure 4 Supplemental 1, Figure 5 

def calculate_distance(neuron_name, VPN, synMap_df):
    synSecList = []
    synTypeList = []
    synpreList = []
    synSegRangeVarList = []

    if VPN not in ["LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4"]:
        raise ValueError("Invalid VPN value")

    mask = synMap_df['type'] == VPN
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
    if VPN not in ["LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4"]:
        raise ValueError("Invalid VPN value")

    mask = synMap_df['type'] == VPN
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
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC6":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC6_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC22":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC1":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC2":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC4":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
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

def bar_whisker_plot_shuffled_vs_original_data():
    #Data  shown here was retrieved from the NN analysis above which generates csv's with the following values. 
    data1 = [42.88, 51.02, 49.8, 49.8, 43.75, 57.27, 40, 58.52, 54.13, 46.93, 65.66, 44.72, 53.85, 45.02]
    data2 = [2.154, 2.455, 3.03, 3.03, 1.998, 6.385, 2.008, 1.914, 2.174, 1.213, 3.03, 1.352, 7.915, 2.005]

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
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC6":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC6_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LC22":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LC22_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC1":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC1_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC2":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC2_receptivefield_COM.xlsx", dtype={'updated_ids': str})
    elif VPN == "LPLC4":
        VPN_COM = pd.read_excel("datafiles/Receptive_field_files/LPLC4_receptivefield_COM.xlsx", dtype={'updated_ids': str})
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
    'Dorsal': [33, np.nan, 32, np.nan, 30, np.nan, 32, np.nan, 60, np.nan, 30, np.nan],
    'Ventral': [21, np.nan, 30, np.nan, 31, np.nan, 35, np.nan, 47, np.nan, 24, np.nan],
    'Anterior': [np.nan, 27, np.nan, 36, np.nan, 38, np.nan, 27, np.nan, 46, np.nan, 32],
    'Posterior': [np.nan, 27, np.nan, 26, np.nan, 23, np.nan, 40, np.nan, 61, np.nan, 22]
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

#############

def main():
    print(os.getcwd())

    ##### 

    #Running Statistics comparing nearest neighbors in original dataset vs shuffled. 

    #Figure 5B
    #bar_whisker_plot_shuffled_vs_original_data()

    # LPLC2_shuffled_DNp01 = [0.66,2.49, 1.83, 1, 1.33, 0.83, 0.89, 1, 1.16, 2.16] # Original 36.88
    # LC4_shuffled_DNp01  = [4.42, 2.31, 1.35, 1.73, 1.15, 0.96, 2.31, 3.46, 2.12, 1.73] # Original 42.88
    # LC4_shuffled_DNp02  = [2.8, 1.78, 2.8, 1.65, 3.56, 3.05, 0.89, 2.04, 2.67, 3.31] # Original 51.02 
    # LC4_shuffled_DNp03  = [1.01, 1.01, 3.03, 3.03, 4.04, 3.03, 6.06, 2.02, 3.03, 4.04] # Original 49.8
    # LPLC2_shuffled_DNp03  = [1.96, 1.76, 2.75, 0.78, 1.76, 4.12, 1.76, 1.96, 1.76, 1.37] # Original 43.75
    # LC22_shuffled_DNp03  = [13.18, 8.18, 3.85, 7.27, 1.82, 9.09, 4.09, 3.64, 5.91, 6.82] # Original  57.27
    # LPLC1_shuffled_DNp03  = [2.19, 3.01, 2.1, 1.37, 2.28, 2.56, 1.64, 1, 2.19, 1.74] # Original 40
    # LPLC4_shuffled_DNp03  = [1.77, 2.28, 0.84, 1.94, 2.02, 1.69, 1.94, 2.7, 1.94, 2.02] # Original 58.52
    # LPLC2_shuffled_DNp04 = [0.59, 1.3, 1.65, 1.18, 1.54, 0.71, 1.06, 1.14, 1.54, 1.42] # Original 54.13
    # LC4_shuffled_DNp04  = [2.21, 1.64, 3.02, 1.7, 2.21, 2.08, 2.39, 1.89, 2.39, 2.21] # Original 46.93
    # LC4_shuffled_DNp06  = [1.01, 1.01, 3.03, 3.03, 4.04, 3.03, 6.06, 2.02, 3.03, 4.04] # Original 65.66
    # LPLC2_shuffled_DNp06  = [1.09, 1.09, 1.4, 1.24, 0.78, 2.33, 1.86, 1.4, 1.09, 1.24] # Original 44.72
    # LC6_shuffled_DNp06  = [9.92, 1.92, 3.85, 9.62, 15.38	, 5.77, 15.38, 5.77, 5.77, 5.77] #Original 53.85
    # LPLC1_shuffled_DNp06  = [3.9, 2.31, 0.87, 1.15, 2.45, 2.31, 2.16, 1.3, 2.45, 1.15] #Original 45.02

    #Statistics in Figure 5C
    # NN_stats(original_data= 45.02, shuffled_data=LPLC1_shuffled_DNp06)

    ##### 
    #Running receptive field analysis for COM in particular hemispheres
    #Figure 3 Supplemental 1
    # VPN = "LPLC2"
    # COM_receptive_field_split_and_visualization(VPN)

    ####
    #Define and intialize neuron model
    neuron_name = "DNp01"

    cell, allSections_py, allSections_nrn, somaSection, erev, axonList, tetherList, dendList, sizSection = initializeModel(neuron_name)#sizSection, erev, axonList = initializeModel(neuron_name)

    if neuron_name == "DNp01":
        synSite_df = pd.read_csv('morphologyData/DNp01_morphData/synMap_DNp01_021324.csv', dtype={'pre': str})
    elif neuron_name == "DNp02":
        synSite_df = pd.read_csv('morphologyData/DNp02_morphData/synMap_DNp02_7623.csv', dtype={'pre': str}) 
    elif neuron_name == "DNp03":
        synSite_df = pd.read_csv('morphologyData/DNp03_morphData/synMap_DNp03_021324.csv', dtype={'pre': str})    
    elif neuron_name == "DNp04":
        synSite_df = pd.read_csv('morphologyData/DNp04_morphData/synMap_DNp04_7623.csv', dtype={'pre': str})    
    elif neuron_name == "DNp06":
        synSite_df = pd.read_csv('morphologyData/DNp06_morphData/synMap_DNp06_7623.csv', dtype={'pre': str})    
    synSite_df = removeAxonalSynapses(synSite_df, dendList, axonList=axonList, tetherList=tetherList, somaSection=somaSection, sizSection=sizSection)
    

    ########

    #Running retinotopy analysis, Figure 4 D-E and Figure 4 supplemental
    #Define VPN population to look at

    #VPN = "LPLC2"
    # dist_df = calculate_distance(neuron_name, VPN, synSite_df)
    # neuron_pair_dist =  calculate_average_distance_between_synapses_per_neuron_pair(neuron_name, VPN, synSite_df)
    # nearest_neighbor_df, avg_distance_per_neuron_df = nearest_neighbor(neuron_name, VPN, sizSection, synSite_df)
    # nn_and_COM_col_coded, all_syns_dist_df, neuron_syn_pair_dists = NN_COM_merge(VPN, nearest_neighbor_df, dist_df, neuron_pair_dist)
    # plot_retinotopy(neuron_syn_pair_dists, VPN, neuron_name)

    ########

    # Printing peak depolarization per VPN
    #Peak_depol_per_VPN(neuron_name, MODE='VPN/')

    # Printing single synapse depolarization at SIZ per VPN

    #plot_peak_depol_at_SIZ_vs_dist_to_SIZ(neuron_name, sizSection)

    # Combing the two dataframes
    #Looking how the integration of synapses and compares to linear summation (integration is sublinear)

    #combine_indiv_syn_depol_and_group_act(neuron_name)

    #individual_syn_depol_per_VPN(neuron_name)

    #############
    #Figure 9
    #plotAllSynsAndhighlightsynsFig9(neuron_name, erev)

    #Figure 10
    #plot_partner_VPN_activations(neuron_name, MODE= 'VPN/')
    #plot_rand_vs_partner_SIZ(neuron_name, MODE='VPN/')
    #plot_rand_vs_partner_SOMA(neuron_name, MODE= 'VPN/')

    #Figure 11
    #Plot_DN_VPN_syns(neuron_name, erev, location = 'SIZ')
    #Plot_DNp03_VPN_syns_by_type(neuron_name, location = "SIZ", type = 'LC')
    #Plot_DNp03_VPN_syns_by_type(neuron_name, location = "SIZ", type = 'LPLC')

    #Figure 12
    #plot_syn_spread_vs_dist_to_siz(neuron_name, erev, MODE='all/', sizSection=sizSection, plots='Rand', synnum =None, trials = None)
    #plot_syn_spread_vs_dist_to_siz_partner_vs_rand(neuron_name, erev, plots='Rand', sizSection=sizSection)
    #plot_syn_spread_vs_dist_to_siz_partner_vs_rand(neuron_name, erev, plots='Partner', sizSection=sizSection)
    #plot_syn_spread_vs_dist_to_siz_from_csv_rand_vs_partner(neuron_name, synnum = 15, trials = 50)


main()


