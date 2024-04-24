from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from neuron import h, gui
from neuron.units import ms, mV
import time as clock
import math
import scipy.optimize as opt
from tkinter import Tk
import tkinter.filedialog as fd



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

def loadEphysData(filename):
    #TODO: confirm this works as expected
    with open(filename) as f:
        #need to create np array of times and voltages

        #count num of data points to preallocate array size
        arraySize = sum(1 for line in open(filename))
        timeArr = np.zeros([arraySize])
        voltageArr = np.zeros([arraySize])


        #HJ ephys data as recorded and converted to .dat from MATLAB has 3 space chars b/w time pt data and voltage data
        #TODO PLEASE ADJUST THIS TO FIT TO FORMAT OF EPHYS DATA BEING READ | CONFIRM BELOW FOR LOOP IS APPROPRIATE
        for i, line in enumerate(f):
            line = line.strip()
            line = line.replace("   ", "")
            line_split = line.split(' ')
            #print(line_split)
            #TODO double check formatting for ephys data to make sure no leading or trailing whitespace
            time = line_split[0]
            voltage = line_split[1]
            timeArr[i] = time
            voltageArr[i] = voltage
        # print(timeArr)
        # print(voltageArr)

    return timeArr, voltageArr

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

def runSim(sectionList_py, electrodeSec, somaSection, exp_tData=None, exp_vData=None, current=-0.04, erev=-72.5, continueRun=1200, injDur=1000):
    #sectionList_py is passed in case you want to set up a recording electrode/vector at a particular indexed section
    #example call: vec.record(sectionList_py[section_index](0.5)._ref_v)


    stimobj = h.IClamp(somaSection(0.5))
    stimobj.delay = 100
    stimobj.dur = injDur
    stimobj.amp = current
    stimobj.i = 0
    ampInjvect = h.Vector().record(stimobj._ref_i)

    vInjVec = h.Vector()
    tInjVec = h.Vector()
    vInjVec.record(somaSection(0.5)._ref_v)
    tInjVec.record(h._ref_t)

    eInjVec = h.Vector()
    eInjVec.record(electrodeSec(0.5)._ref_v)


    h.finitialize(erev * mV)
    h.continuerun(continueRun * ms)
    aInj = ampInjvect.to_python()
    vInj = vInjVec.to_python()
    tInj = tInjVec.to_python()
    vInj_np = np.array(vInj)
    tInj_np = np.array(tInj)
    eInjVec_np = np.array(eInjVec)
    aInjVec_np = np.array(aInj)

    #############################################################
    ### plotting code for exp data vs sim data for given trial ###
    #############################################################

    #only plots if exp time and voltage data vectors are passed to the runSim function
    if exp_tData is not None and exp_vData is not None:
        fig = plt.plot(stimobj, vInj_np)
        plt.plot(exp_tData*1000, exp_vData)
        plt.plot(tInj_np+14500, vInj_np)
        plt.xlabel('time (ms)')
        plt.ylabel('mV')
        plt.show()
        h.stoprun
        #TODO: is h.stoprun necessary?
        # plt.close('all')

    return tInj_np, vInj_np, eInjVec_np, aInjVec_np

def initializeModel(morph_file, neuron_name, hasElectrode):

    cell = instantiate_swc(morph_file)

    allSections_nrn = h.SectionList()
    for sec in h.allsec():
        allSections_nrn.append(sec=sec)
    
    # Create a Python list from this SectionList
    # Select sections from the list by their index

    allSections_py = [sec for sec in allSections_nrn]    

    #colorR for axon sections | colorB for test work | colorG for soma | colorK for presumed SIZ
    colorR = h.SectionList()
    colorB = h.SectionList()
    colorG = h.SectionList()
    colorK = h.SectionList()
    colorV = h.SectionList()

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
            colorG.append(somaSection)
        elif "axon" in sec.name():
            axonList.append(sec)
        elif "dend_11" in sec.name():
            tetherList.append(sec)
        elif "dend_6" in sec.name():
            axonEnd = sec
        else:
            dendList.append(sec)

    i = 0
    for sec in axonList:
        if i == sizIndex:
            sizSection = sec
        i += 1
        
    # i = 0
    # for sec in axonList:
    #     if len(sec.children()) < 1:
    #         # print(sec)
    #         # colorV.append(sec)
    #         axonEnd = sec
    #         # break
    #     i += 1
    colorV.append(axonEnd)

    allSections_py = createAxon(axonEnd, allSections_py, neuron_name)
    if hasElectrode:
        allSections_py, electrodeSec = createElectrode(somaSection, allSections_py, neuron_name)
    else:
        electrodeSec = None

    shape_window = h.PlotShape(h.SectionList(allSections_py))           # Create a shape plot
    shape_window.exec_menu('Show Diam')    # Show diameter of sections in shape plot
    shape_window.color_all(9)

    shape_window.color_list(axonList, 2)
    shape_window.color_list(colorG, 4)  
    shape_window.color_list(tetherList, 1)
    shape_window.color_list(dendList, 3)
    shape_window.color_list(colorV, 7)


    if neuron_name == "DNp01":
        change_Ra(10.3294)#19.7102)#34)#33.2617)
        change_gLeak(0.0002982717, erev=-72.5)#0.000584725, erev=-72.5)#3.4e-4, erev=-72.5)#4.4e-9)    
        change_memCap(1, electrodeSec=None)#2.0934, electrodeSec=None)#1)#1.4261)
        erev = -72.5
    elif neuron_name == "DNp02": #or neuron_name == "DNp06":
        change_Ra(91)
        change_gLeak(0.0002415, erev=-70.8)  
        change_memCap(1)
        erev = -70.8#-70.812 TOOK FIRST 5000 DATAPOINTS FROM EACH HJ CURRENT STEP AND TOOK MEAN
    elif neuron_name == "DNp03":
        # 30.964169523649325,0.00013449623136848623,0.7553324245115953
        change_Ra(ra=30.964)#41.0265)
        change_gLeak(gleak=0.000134496, erev=-58.75)#0.000239761, erev=-58.75)#0.00034, erev=-58.75)
        change_memCap(memcap=0.75533)#1.3765)
        erev = -58.75
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
    
    return cell, allSections_py, allSections_nrn, somaSection, sizSection, axonEnd, erev, axonList, tetherList, dendList, electrodeSec#, shape_window

def createElectrode(somaSection, pySectionList, neuron_name=None):
    electrodeSec = h.Section()
    electrodeSec.L = 10
    electrodeSec.diam = 1
    electrodeSec.connect(somaSection, 0)

    pySectionList.append(electrodeSec)

    # equivCylAxon.connect(axonEnd, 1)
    # print(equivCylAxon(0.5).area())

    return pySectionList, electrodeSec

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

    lineDir = np.array([lineEnd_X, lineEnd_Y, lineEnd_Z]) - np.array([lineStart_X, lineStart_Y, lineStart_Z])
    line_direction_norm = lineDir / np.linalg.norm(lineDir)

    # TODO: ADD OTHERS
    if neuron_name == "DNp01":
        # equivCylHeight = 360.346
        # equivCylDiam = 1.199*2
        equivCylHeight = 241.69#184.656
        equivCylDiam = 3.32*2#1.34*2
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

    new_point = np.array([lineEnd_X, lineEnd_Y, lineEnd_Z]) + equivCylHeight * line_direction_norm
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

def nsegDiscretization(sectionListToDiscretize):
    #TODO: SEE NOTE
    #NOTE TO SELF, NOT IDEAL BC ONLY INCREASES, DOESNT DECREASE, WRITE CODE TO DECREMENT SEGMENT NUMBER IF APPLICABLE
    #TO SAVE ON COMPUTATIONAL COMPLEXITY
    for sec in sectionListToDiscretize:
        #secMorphDict = sec.psection()['morphology']
        secDensityMechDict = sec.psection()['density_mechs']
        secLambda = math.sqrt( ( (1 / secDensityMechDict['pas']['g'][0]) * sec.diam) / (4*sec.Ra) )
        ###print("\nlength of section", sec, "is", (sec.L/sec.nseg)/secLambda, "lambda | sec.L:", sec.L, "| sec.nseg:", sec.nseg, "| lambda:", secLambda)
        if (sec.L/sec.nseg)/secLambda > 0.1:
            ###print('had to fix section', sec)
            ###print("\nlength of section", sec, "is", (sec.L/sec.nseg)/secLambda, "lambda")
            numSeg_log2 = math.log2(sec.L/secLambda / 0.1)
            numSeg = math.ceil(2**numSeg_log2)
            if numSeg % 2 == 0:
                numSeg += 1
            sec.nseg = numSeg
            ###print("\nlength of section", sec, "is now", (sec.L/sec.nseg)/secLambda, "lambda")
            ###print("FIXED by doing", sec.nseg, "segments")

    return

def extractData(exp_timeData, exp_voltageData, sim_timeData, sim_voltageData):

    EXP_INJ_START = 10000/1000
    EXP_INJ_END = 11000/1000

    exp_durInj_idx = np.where(np.logical_and(exp_timeData >= EXP_INJ_START, exp_timeData <= (EXP_INJ_END+(200/1000))))
    exp_durInj_timeData = exp_timeData[exp_durInj_idx]
    exp_durInj_voltageData = exp_voltageData[exp_durInj_idx]

    SIM_INJ_START = 100
    SIM_INJ_END = 1200

    sim_durInj_idx = np.where(np.logical_and(sim_timeData >= SIM_INJ_START, sim_timeData <= (SIM_INJ_END+0.0)))
    sim_durInj_timeData = sim_timeData[sim_durInj_idx]
    sim_durInj_voltageData = sim_voltageData[sim_durInj_idx]
    sim_durInj_timeData = sim_durInj_timeData

    RS_sim_durInj_timeData = sim_durInj_timeData[0::2]
    RS_sim_durInj_voltageData = sim_durInj_voltageData[0::2]

    print(RS_sim_durInj_timeData.shape)

    return [((exp_durInj_timeData-10)*1000)[0:1000], exp_durInj_voltageData[0:1000]], [(RS_sim_durInj_timeData-100)[0:1000], RS_sim_durInj_voltageData[0:1000]],[((exp_durInj_timeData-9.9)*1000)[-4000:-3000], exp_durInj_voltageData[-4000:-3000]], [(RS_sim_durInj_timeData)[-2000:-1000], RS_sim_durInj_voltageData[-2000:-1000]]

def func_decay(x, a, tau_d, c):
     return a * np.exp(-x / tau_d) + c

def func_recovery(x, a, tau_r, c):
     return -a * np.exp(-x / tau_r) + c

def calcTau_decay(time, voltage, t, v):
    optimizedParameters, pcov = opt.curve_fit(func_decay, time, voltage)
    plt.plot(time, func_decay(time, *optimizedParameters), 'r-', label="fit")
    plt.plot(time, voltage, "b-", label="Data")
    plt.plot(t, v, 'g-')
    plt.axvline(x=time[0]+0.25, color='k')
    plt.xlim([time[0], time[0]+10])
    plt.xlabel('Time (ms)')
    plt.ylabel('Votlage (mV)')
    plt.legend(['curve fit to simulated trace', 'simulated trace', 'experimental data','cutoff (0.25 ms)'])
    plt.title("DNp03 decay")
    plt.show()
    tau_d = optimizedParameters[1]
    return tau_d

def calcTau_recovery(time, voltage, t, v):
    optimizedParameters, pcov = opt.curve_fit(func_recovery, time-time[0], voltage)
    plt.plot(time, func_recovery(time-time[0], *optimizedParameters), 'r-', label="fit")
    plt.plot(time, voltage, "b-", label="Data")
    plt.plot(t, v, 'g-')
    plt.axvline(x=time[0]+0.25, color='k')
    plt.xlim([time[0], time[0]+10])
    plt.xlabel('Time (ms)')
    plt.ylabel('Votlage (mV)')
    plt.legend(['curve fit to simulated trace', 'simulated trace', 'experimental data','cutoff (0.25 ms)'])
    plt.title('DNp03 recovery')
    plt.show()
    tau_r = optimizedParameters[1]
    return tau_r

def calculateRMSE_justDecay(exp_timeData, exp_voltageData, sim_timeData, sim_voltageData):

    # EXP_INJ_START = 500/1000
    # EXP_INJ_END = 1000/1000

    EXP_INJ_START = 10000/1000
    EXP_INJ_END = 11000/1000

    exp_durInj_idx = np.where(np.logical_and(exp_timeData >= EXP_INJ_START, exp_timeData <= (EXP_INJ_END+(200/1000))))
    exp_durInj_timeData = exp_timeData[exp_durInj_idx]
    exp_durInj_voltageData = exp_voltageData[exp_durInj_idx]

    SIM_INJ_START = 100
    SIM_INJ_END = 500

    sim_durInj_idx = np.where(np.logical_and(sim_timeData >= SIM_INJ_START, sim_timeData <= (SIM_INJ_END+0.0)))
    
    sim_durInj_timeData = sim_timeData[sim_durInj_idx]
    sim_durInj_voltageData = sim_voltageData[sim_durInj_idx]
    sim_durInj_timeData = sim_durInj_timeData

    # print(exp_durInj_voltageData.shape)
    # print(sim_durInj_voltageData.shape)

    RS_sim_durInj_timeData = sim_durInj_timeData[0::2]
    RS_sim_durInj_voltageData = sim_durInj_voltageData[0::2]

    # NEW = plt.figure()
    # fig, ax = plt.subplots()
    # ax1 = plt.subplot(3, 1, 1)
    # ax2 = plt.subplot(3, 1, 2)
    # ax3 = plt.subplot(3, 1, 3)
    # ax1.plot(exp_durInj_timeData/1000, exp_durInj_voltageData)

    # ax2.plot(sim_durInj_timeData/1000-0.5, sim_durInj_voltageData)

    # ax3.plot(exp_timeData, exp_voltageData, 'g-')
    # ax3.plot(sim_timeData/1000, sim_voltageData, 'b-')

    lenExpData = len(exp_durInj_voltageData)

    #exp_voltageData_fitCurve = func(exp_durInj_timeData/1000, *exp_popt)
    
    RMSE_val = np.sqrt(np.mean((exp_durInj_voltageData[(0+3):150] - RS_sim_durInj_voltageData[(0+3):150]) ** 2))#0:lenExpData-1]) ** 2))
    RMSE_val_100 = np.sqrt(np.mean((exp_durInj_voltageData[(0+3):100] - RS_sim_durInj_voltageData[(0+3):100]) ** 2))
    # print(RS_sim_durInj_timeData[3])
    # print(RS_sim_durInj_timeData[150])
    # print(RS_sim_durInj_timeData[1000])


    plt.figure()
    plt.plot(((exp_durInj_timeData-9.9)*1000)[3:100], exp_durInj_voltageData[3:100], 'k')
    plt.plot((RS_sim_durInj_timeData)[3:100], RS_sim_durInj_voltageData[3:100], 'g')
    plt.plot(((exp_durInj_timeData-9.9)*1000)[100:150], exp_durInj_voltageData[100:150], 'k--')
    plt.plot((RS_sim_durInj_timeData)[100:150], RS_sim_durInj_voltageData[100:150], 'g--')
    # plt.plot((exp_durInj_timeData-10)*1000, exp_durInj_voltageData, 'b-')
    # plt.plot(sim_durInj_timeData/1000-0.5, sim_durInj_voltageData, 'r')

    # plt.plot(((exp_durInj_timeData-10)*1000)[0:150], exp_durInj_voltageData[0:150], 'k')
    # plt.plot((RS_sim_durInj_timeData-100)[0:1000], RS_sim_durInj_voltageData[0:1000], 'r')
    # plt.plot((RS_sim_durInj_timeData-100)[0:150], RS_sim_durInj_voltageData[0:150], 'g')
    # plt.axvline(x=((exp_durInj_timeData-10)*1000)[0]+0.25, color='b')
    # plt.axvline(x=(RS_sim_durInj_timeData-100)[5], color='r')
    # plt.title('just decay')
    # # print(exp_durInj_timeData[30])
    # # plt.plot(sim_durInj_timeData/1000-0.5, RS_sim_durInj_voltageData[0:lenExpData-1], 'g')

    plt.show()
    # plt.close('all')

    return RMSE_val, RMSE_val_100

def calculateRMSE_justRecovery(exp_timeData, exp_voltageData, sim_timeData, sim_voltageData):
    # plt.figure()
    EXP_INJ_START = 10000/1000
    EXP_INJ_END = 11000/1000

    exp_durInj_idx = np.where(np.logical_and(exp_timeData >= EXP_INJ_START, exp_timeData <= (EXP_INJ_END+(200/1000))))
    exp_durInj_timeData = exp_timeData[exp_durInj_idx]
    exp_durInj_voltageData = exp_voltageData[exp_durInj_idx]

    SIM_INJ_START = 100
    SIM_INJ_END = 1200

    sim_durInj_idx = np.where(np.logical_and(sim_timeData >= SIM_INJ_START, sim_timeData <= (SIM_INJ_END+000)))
    sim_durInj_timeData = sim_timeData[sim_durInj_idx]
    sim_durInj_voltageData = sim_voltageData[sim_durInj_idx]
    sim_durInj_timeData = sim_durInj_timeData

    RS_sim_durInj_timeData = sim_durInj_timeData[0::2]
    RS_sim_durInj_voltageData = sim_durInj_voltageData[0::2]

    # fig, ax = plt.subplots()
    # ax1 = plt.subplot(3, 1, 1)
    # ax2 = plt.subplot(3, 1, 2)
    # ax3 = plt.subplot(3, 1, 3)
    # ax1.plot(exp_durInj_timeData/1000, exp_durInj_voltageData)

    # ax2.plot(sim_durInj_timeData/1000-0.5, sim_durInj_voltageData)

    # ax3.plot(exp_timeData, exp_voltageData, 'g-')
    # ax3.plot(sim_timeData/1000, sim_voltageData, 'b-')

    lenExpData = len(exp_durInj_voltageData)

    #exp_voltageData_fitCurve = func(exp_durInj_timeData/1000, *exp_popt)
    
    # RMSE_val = np.sqrt(np.mean((exp_durInj_voltageData[(-4000+5):-3000] - RS_sim_durInj_voltageData[(-2000+5):-1000]) ** 2))#0:lenExpData-1]) ** 2))
    RMSE_val = np.sqrt(np.mean((exp_durInj_voltageData[(-4000+3):-3850] - RS_sim_durInj_voltageData[(-2000+3):-1850]) ** 2))#0:lenExpData-1]) ** 2))
    RMSE_val_100 = np.sqrt(np.mean((exp_durInj_voltageData[(-4000+3):-3900] - RS_sim_durInj_voltageData[(-2000+3):-1900]) ** 2))#0:lenExpData-1]) ** 2))
    # print(exp_durInj_timeData[-6000], exp_durInj_timeData[-5999], RS_sim_durInj_timeData[-4000], RS_sim_durInj_timeData[-3999])

    # print(RS_sim_durInj_timeData[-2000+3])
    # print(RS_sim_durInj_timeData[-1850])
    # print(RS_sim_durInj_timeData[-1000])
    # print(exp_durInj_timeData[-4000+3])
    # print(exp_durInj_timeData[-3000])
    # print(exp_durInj_voltageData.shape)
    # print(sim_durInj_voltageData.shape)

    plt.figure()
    plt.plot(((exp_durInj_timeData-9.9)*1000)[-4000:-3900], exp_durInj_voltageData[-4000:-3900], 'k')
    plt.plot(((exp_durInj_timeData-9.9)*1000)[-3900:-3850], exp_durInj_voltageData[-3900:-3850], 'k--')
    plt.plot((RS_sim_durInj_timeData)[-2000:-1900], RS_sim_durInj_voltageData[-2000:-1900], 'g')
    plt.plot((RS_sim_durInj_timeData)[-1900:-1850], RS_sim_durInj_voltageData[-1900:-1850], 'g--')

    # plt.plot(((exp_durInj_timeData-9.9)*1000)[-4000:-3000], exp_durInj_voltageData[-4000:-3000], 'k')
    # plt.plot((RS_sim_durInj_timeData+700)[-2000:-1000], RS_sim_durInj_voltageData[-2000:-1000], 'g')
    # plt.axvline(x=((exp_durInj_timeData-9.9)*1000)[-4000]+0.25, color='b')
    # plt.axvline(x=(RS_sim_durInj_timeData+700)[-1995], color='r')
    # plt.title('just recovery')

    plt.show()

    return RMSE_val, RMSE_val_100


def plotMorphColorCode_wSIZ(allSections_py, somaSection, axonList, tetherList, dendList, sizSection, neuron_name):
    morph_shape_window = h.Shape()#h.SectionList(allSections_py))           # Create a shape plot
    # morph_shape_window.exec_menu('Show Diam')    # Show diameter of sections in shape plot
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


    if sizSection is not None:
        if neuron_name == "DNp01":
            synTemp = h.AlphaSynapse(sizSection(0.0551186))
  
    return morph_shape_window



def main():
    hasElectrode=True
    neuron_name = "DNp01"
    Tk().withdraw()
    fd_title = "Select morphology file to use for synapse mapping"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"morphologyData", title=fd_title)

    # morph_file = 'datafiles/morphologyData/' + neuron_name + '_morphData/' + neuron_name + '_um_model.swc'
    cell, allSections_py, allSections_nrn, somaSection, sizSection, axonEnd, erev, axonList, tetherList, dendList, electrodeSec = initializeModel(morph_file, neuron_name, hasElectrode)

    '''
    # electrodeArea1 = electrodeSec(0).area()
    # electrodeArea = electrodeSec(0.5).area()

    # totalCap = 2e-12


    # secSA_nTB = 0
    # for pt in range(electrodeSec.n3d()):
    #     # print(sec, pt, "/", sec.n3d())
    #     #print(pt, sec.n3d())
    #     if pt+1 == electrodeSec.n3d():
    #         break
    #     # print("pair:", pt, pt+1)
    #     d1 = electrodeSec.diam3d(pt)
    #     r1 = d1/2
    #     r1 = r1
    #     #print(r1)
    #     d2 = electrodeSec.diam3d(pt+1)
    #     r2 = d2/2
    #     r2 = r2
    #     #TODO: SEC.L has to be distance
    #     x1 = electrodeSec.x3d(pt)
    #     x2 = electrodeSec.x3d(pt+1)
    #     y1 = electrodeSec.y3d(pt)
    #     y2 = electrodeSec.y3d(pt+1)
    #     z1 = electrodeSec.z3d(pt)
    #     z2 = electrodeSec.z3d(pt+1)
    #     DIST = math.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))
    #     #print(DIST)
    #     DIST = DIST
    #     SA_minus_TopBottom_p1 = math.pi * (r1+r2)
    #     SA_minus_TopBottom_p2 = math.sqrt(((r1-r2)*(r1-r2))+(DIST*DIST))
    #     SA_minus_TopBottom = SA_minus_TopBottom_p1*SA_minus_TopBottom_p2
    #     secSA_nTB += SA_minus_TopBottom
    # #print(secSA_nTB, sec(0.5).area())
    # # SA_byNeuron += electrodeSec(0.5).area()
    # # totSA += secSA_nTB

    # print(electrodeArea, secSA_nTB)
    '''




    msw = plotMorphColorCode_wSIZ(allSections_py, somaSection, axonList, tetherList, dendList, sizSection, neuron_name)

    timestr = clock.strftime("%Y%m%d-%H%M%S")
    if neuron_name == "DNp01":
        DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray = loadEphysData('ephysData/DNp01_hp_withoutBiasCurrent_avg_fly1thru6_BRIDGE_FIXED.dat')
        DNp01_fly1thru6_60pA_hpol_noHold_timeArray = DNp01_fly1thru6_60pA_hpol_noHold_timeArray - 0.00025
        expData = [[DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray]]

    elif neuron_name == "DNp03":
        #DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray = loadEphysData('ephysData/DNp03_7423_hp_withoutBiasCurrent_avg_fly4567.dat')
        DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray = loadEphysData('ephysData/DNp03_7423_hp_withoutBiasCurrent_avg_fly457.dat')
        expData = [[DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray]]

    for trialSet in expData:
        exp_erev = round(np.mean(trialSet[1][0:5000]), 2)
        exp_erev = -66.4
        trialSet.append(exp_erev)
    

    if neuron_name == "DNp01": 
        erev = -66.4298 - 0.2
        raVal = 212
        gleakVal = 1/2300
        cmVal = 0.7

    elif neuron_name == "DNp03":
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -61.15

    elec_raVal = 235.6#150                     
    elec_gleakVal = 0
    elec_cmVal = 6.4 # Membrane capacitance in micro Farads / cm^2
    #electode geom props: l = 10um | d = 1um

    sealCon_8GOhm = 0.0003978
    sealCon_2GOhm = 0.0016
    elec_gleakVal = sealCon_8GOhm


    change_Ra(ra=raVal, electrodeSec=electrodeSec, electrodeVal = elec_raVal)
    change_gLeak(gleak=gleakVal, erev=erev, electrodeSec=electrodeSec, electrodeVal = elec_gleakVal)
    change_memCap(memcap=cmVal, electrodeSec=electrodeSec, electrodeVal = elec_cmVal)


    if neuron_name == "DNp01":
        current = -0.06
    elif neuron_name == "DNp03":
        current = -0.002
    else:
        print("set current value")
        quit(0)

    #records from the soma section to fit
    # tInj_np_hpol_DecayRecov, vInj_np_hpol_DecayRecov = runSim(allSections_py, somaSection, exp_tData=None, exp_vData=None, current=-0.002, erev = expData[0][2], continueRun=1200)#-60)#current=-0.0333, erev = -58.75)
    
    #records from the electode section for fitting.
    tInj_np_hpol_justDecay, vInj_np_hpol_justDecay, eInj_np_hpol_justDecay, aInjVec_np_hpol_current = runSim(allSections_py, electrodeSec, somaSection, exp_tData=None, exp_vData=None, current=current, erev = -59.55  , continueRun=1200)
    



    ERROR_decay, ERR_decay_100 = calculateRMSE_justDecay(expData[0][0], expData[0][1], tInj_np_hpol_justDecay, vInj_np_hpol_justDecay)
    if (np.isinf(ERROR_decay)): 
            return float('inf'),
    ERROR_recov, ERR_recov_100 = calculateRMSE_justRecovery(expData[0][0], expData[0][1], tInj_np_hpol_justDecay, vInj_np_hpol_justDecay)
    if (np.isinf(ERROR_recov)): 
            return float('inf'),

    print(ERROR_decay, ERROR_recov, ERROR_decay+ERROR_recov)
    print(ERR_decay_100, ERR_recov_100, ERR_decay_100+ERR_recov_100)

    errVal = ERROR_decay+ERROR_recov
    print(errVal)

    ERROR_decay_elec, ERR_decay_100_elec = calculateRMSE_justDecay(expData[0][0], expData[0][1], tInj_np_hpol_justDecay, eInj_np_hpol_justDecay)
    if (np.isinf(ERROR_decay_elec)): 
            return float('inf'),
    ERROR_recov_elec, ERR_recov_100_elec = calculateRMSE_justRecovery(expData[0][0], expData[0][1], tInj_np_hpol_justDecay, eInj_np_hpol_justDecay)
    if (np.isinf(ERROR_recov_elec)): 
            return float('inf'),

    print(ERROR_decay_elec, ERROR_recov_elec, ERROR_decay_elec+ERROR_recov_elec)
    print(ERR_decay_100_elec, ERR_recov_100_elec, ERR_decay_100_elec+ERR_recov_100_elec)

    errVal_elec = ERROR_decay_elec+ERROR_recov_elec
    print(errVal_elec)
    
    gs = gridspec.GridSpec(3, 2)
    fig_DR = plt.figure()

    ax_full = fig_DR.add_subplot(gs[0, :]) # row 1, span all columns
    ax_decay = fig_DR.add_subplot(gs[1, :])
    ax_recov = fig_DR.add_subplot(gs[2, :])
    if neuron_name == "DNp01":
        plt.plot(tInj_np_hpol_justDecay, aInjVec_np_hpol_current)
        ax_full.plot(DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray)
        ax_full.plot((tInj_np_hpol_justDecay/1000)+9.9, eInj_np_hpol_justDecay)
        ax_full.plot((tInj_np_hpol_justDecay/1000)+9.9, vInj_np_hpol_justDecay)
        ax_decay.plot(DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray)
        ax_decay.plot((tInj_np_hpol_justDecay/1000)+9.9, eInj_np_hpol_justDecay)
        ax_decay.plot((tInj_np_hpol_justDecay/1000)+9.9, vInj_np_hpol_justDecay)
        ax_recov.plot(DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray)
        ax_recov.plot((tInj_np_hpol_justDecay/1000)+9.9, eInj_np_hpol_justDecay)
        ax_recov.plot((tInj_np_hpol_justDecay/1000)+9.9, vInj_np_hpol_justDecay)
    elif neuron_name == "DNp03":
        plt.plot(tInj_np_hpol_justDecay, aInjVec_np_hpol_current)
        ax_full.plot(DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray)
        ax_full.plot((tInj_np_hpol_justDecay/1000)+9.9, eInj_np_hpol_justDecay)
        ax_decay.plot(DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray)
        ax_decay.plot((tInj_np_hpol_justDecay/1000)+9.9, eInj_np_hpol_justDecay)
        ax_recov.plot(DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray)
        ax_recov.plot((tInj_np_hpol_justDecay/1000)+9.9, eInj_np_hpol_justDecay)

    TEST = ((tInj_np_hpol_justDecay/1000))[150]

    ax_full.spines['top'].set_visible(False)
    ax_full.spines['right'].set_visible(False)
    ax_decay.spines['top'].set_visible(False)
    ax_decay.spines['right'].set_visible(False)
    ax_recov.spines['top'].set_visible(False)
    ax_recov.spines['right'].set_visible(False)

    ax_full.set_title(f'Simulated vs Experimental {neuron_name} Trace (Hyperpolarization)')

    ax_full.set_xlim(9.75, 11.25)
    ax_decay.set_xlim(9.975, 10.025)
    ax_recov.set_xlim(10.975, 11.025)

    ax_full.set_ylabel('Voltage (mV)')
    ax_full.set_xlabel('Time (ms)')

    plt.suptitle('Optimization of {} (ra={}, gleak={}, cm={}, electrode={}, elec_erev=0, elec_gleak={}, elec_ra={}, elec_cm={}, err={})'.format(neuron_name, round(raVal, 3), round(gleakVal, 6), round(cmVal, 3), hasElectrode, elec_gleakVal, elec_raVal, elec_cmVal, errVal))
    plt.tight_layout()

   # plt.savefig('DNp01_finalPassiveProp.svgz', format='svgz')
    plt.show()

main()