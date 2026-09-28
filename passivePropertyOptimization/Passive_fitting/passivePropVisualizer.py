from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from neuron import h, gui
from neuron.units import ms, mV
import time as clock
import math
import scipy.optimize as opt
from tkinter import Tk
import tkinter.filedialog as fd
import seaborn as sns
import pandas as pd
from scipy.optimize import curve_fit



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
    with open(filename) as f:
        #need to create np array of times and voltages

        #count num of data points to preallocate array size
        arraySize = sum(1 for line in open(filename))
        timeArr = np.zeros([arraySize])
        voltageArr = np.zeros([arraySize])
        STDArr = np.zeros([arraySize])

        for i, line in enumerate(f):
            line = line.strip()
            line = line.replace("   ", "")
            line_split = line.split(' ')
            time = line_split[0]
            voltage = line_split[1]
            STD = line_split[2]
            timeArr[i] = time
            voltageArr[i] = voltage
            STDArr[i] = STD

    return timeArr, voltageArr, STDArr

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

    # i = 0
    # for sec in axonList:
    #     if i == sizIndex:
    #         sizSection = sec
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
        change_Ra(10.3294)
        change_gLeak(0.0002982717, erev=-72.5)
        change_memCap(1, electrodeSec=None)
        erev = -72.5
    elif neuron_name == "DNp02": 
        change_Ra(91)
        change_gLeak(0.0002415, erev=-70.8)  
        change_memCap(1)
        erev = -70.8
    elif neuron_name == "DNp03":
        change_Ra(ra=30.964)
        change_gLeak(gleak=0.000134496, erev=-58.75)
        change_memCap(memcap=0.75533)#
        erev = -58.75
    elif neuron_name == "DNp06":
        change_Ra(91)
        change_gLeak(0.0002415, erev=-60)  
        change_memCap(1)
        erev = -60
    else:
        change_Ra()
        change_gLeak()
        change_memCap()
        erev=-72.5
    
    return cell, allSections_py, allSections_nrn, somaSection, sizSection, axonEnd, erev, axonList, tetherList, dendList, electrodeSec#, shape_window

def createElectrode(somaSection, pySectionList, neuron_name=None):
    electrodeSec = h.Section()
    electrodeSec.L = 10
    electrodeSec.diam = 1
    electrodeSec.connect(somaSection, 0)

    pySectionList.append(electrodeSec)

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
        equivCylHeight = 241.69
        equivCylDiam = 3.32
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
    # d_lambda = 0.1  # fraction of lambda
    # frequency = 100  # Hz
    for sec in sectionListToDiscretize:
        secDensityMechDict = sec.psection()['density_mechs']
        secLambda = math.sqrt( ( (1 / secDensityMechDict['pas']['g'][0]) * sec.diam) / (4*sec.Ra) )
        if (sec.L/sec.nseg)/secLambda > 0.1:
            numSeg_log2 = math.log2(sec.L/secLambda / 0.1)
            numSeg = math.ceil(2**numSeg_log2)
            if numSeg % 2 == 0:
                numSeg += 1
            sec.nseg = numSeg
    # for sec in sectionListToDiscretize:
    #     secDensityMechDict = sec.psection()['density_mechs']
    #     c_m = secDensityMechDict.get('pas', {}).get('e', [0.8])[0] if 'cm' in dir(sec) else 0.8
    #     c_m = sec.cm  # Use the section's cm directly
    #     r_a = sec.Ra
    #     diameter = sec.diam
    #     length = sec.L
        
    #     lambda_f = 1e5 * math.sqrt(diameter / (4 * math.pi * frequency * c_m * r_a))
    #     ncomp = int((length / (d_lambda * lambda_f) + 0.9) / 2) * 2 + 1
    #     sec.nseg = ncomp
    return

def calculateRMSE_justDecay(exp_timeData, exp_voltageData, sim_timeData, sim_voltageData):

    # First RMSE window: 248 ms to 252 ms
    start_ms = 250.2
    end_ms = 253

    # Second RMSE window: 248 ms to 257 ms (adds 5 ms)
    end_ms_extended = 257

    # --- First window ---
    exp_mask = (exp_timeData >= start_ms) & (exp_timeData <= end_ms)
    sim_mask = (sim_timeData >= start_ms) & (sim_timeData <= end_ms)

    exp_voltage_window = exp_voltageData[exp_mask]
    sim_voltage_window = sim_voltageData[sim_mask][::2]  # downsample

    min_len = min(len(exp_voltage_window), len(sim_voltage_window))
    RMSE_val = np.sqrt(np.mean((exp_voltage_window[:min_len] - sim_voltage_window[:min_len]) ** 2))

    # --- Second window ---
    exp_mask_ext = (exp_timeData >= start_ms) & (exp_timeData <= end_ms_extended)
    sim_mask_ext = (sim_timeData >= start_ms) & (sim_timeData <= end_ms_extended)

    exp_voltage_window_ext = exp_voltageData[exp_mask_ext]
    sim_voltage_window_ext = sim_voltageData[sim_mask_ext][::2]  # downsample

    min_len_ext = min(len(exp_voltage_window_ext), len(sim_voltage_window_ext))
    RMSE_val_extended = np.sqrt(np.mean((exp_voltage_window_ext[:min_len_ext] - sim_voltage_window_ext[:min_len_ext]) ** 2))

    # Plotting
    plt.figure()
    plt.plot(exp_timeData[exp_mask], exp_voltage_window, 'k', label=f'Experimental ({start_ms}-{end_ms} ms)')
    plt.plot(sim_timeData[sim_mask][::2][:min_len], sim_voltage_window[:min_len], 'g', label=f'Simulation ({start_ms}-{end_ms} ms)')
    plt.plot(exp_timeData[exp_mask_ext], exp_voltage_window_ext, 'k--', label='Experimental (248–257)')
    plt.plot(sim_timeData[sim_mask_ext][::2][:min_len_ext], sim_voltage_window_ext[:min_len_ext], 'g--', label=f'Simulation {start_ms}-{end_ms_extended} ms)')
    plt.xlabel("Time (ms)")
    plt.ylabel("Voltage")
    plt.title(f"Voltage Comparison: RMSE for {start_ms}- {end_ms} ms and {start_ms}-{end_ms_extended} ms")
    plt.legend()
    # plt.show()

    return RMSE_val, RMSE_val_extended

def calculateRMSE_justRecovery(exp_timeData, exp_voltageData, sim_timeData, sim_voltageData):

    # Define recovery window
    start_ms = 300.2
    end_ms = 303
    end_ms_extended = 305  # extended by 5 ms

    # --- First window: 297–300 ms ---
    exp_mask = (exp_timeData >= start_ms) & (exp_timeData <= end_ms)
    sim_mask = (sim_timeData >= start_ms) & (sim_timeData <= end_ms)

    exp_voltage_window = exp_voltageData[exp_mask]
    sim_voltage_window = sim_voltageData[sim_mask][::2]  # downsample

    min_len = min(len(exp_voltage_window), len(sim_voltage_window))
    RMSE_val = np.sqrt(np.mean((exp_voltage_window[:min_len] - sim_voltage_window[:min_len]) ** 2))

    # --- Second window: 297–305 ms ---
    exp_mask_ext = (exp_timeData >= start_ms) & (exp_timeData <= end_ms_extended)
    sim_mask_ext = (sim_timeData >= start_ms) & (sim_timeData <= end_ms_extended)

    exp_voltage_window_ext = exp_voltageData[exp_mask_ext]
    sim_voltage_window_ext = sim_voltageData[sim_mask_ext][::2]  # downsample

    min_len_ext = min(len(exp_voltage_window_ext), len(sim_voltage_window_ext))
    RMSE_val_extended = np.sqrt(np.mean((exp_voltage_window_ext[:min_len_ext] - sim_voltage_window_ext[:min_len_ext]) ** 2))

    # Plotting
    plt.figure()
    plt.plot(exp_timeData[exp_mask], exp_voltage_window, 'k', label='Exp (297–300)')
    plt.plot(sim_timeData[sim_mask][::2][:min_len], sim_voltage_window[:min_len], 'g', label='Sim (297–300)')
    plt.plot(exp_timeData[exp_mask_ext], exp_voltage_window_ext, 'k--', label='Exp (297–305)')
    plt.plot(sim_timeData[sim_mask_ext][::2][:min_len_ext], sim_voltage_window_ext[:min_len_ext], 'g--', label='Sim (297–305)')

    plt.xlabel("Time (ms)")
    plt.ylabel("Voltage")
    plt.title("Recovery Phase RMSE")
    plt.legend()
    # plt.show()

    return RMSE_val, RMSE_val_extended

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
    colorK.append(sizSection)
    
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

def apply_liquid_junction_correction(expData, ljp_correction_mV):
    """
    Applies liquid junction potential correction to the experimental voltage data.

    Parameters:
        expData (list): A list containing [timeArray, voltageArray]
        ljp_correction_mV (float): The LJP correction value in millivolts (to subtract)

    Returns:
        corrected_time (np.array): Unchanged time array
        corrected_voltage (np.array): Voltage array after LJP correction
    """
    time_array, voltage_array = expData

    # Apply LJP correction
    corrected_voltage = voltage_array - ljp_correction_mV

    return time_array, corrected_voltage

def normalized_rmse(rmse, signal):
    return 100 * rmse / np.ptp(signal)  # percent RMSE relative to signal range

def analyze_experimental_properties(neuron_name):
    # Select path based on neuron
    if neuron_name == "DNp01":
        path = "ephysData/DNp01_ephysData/DNp01_flies1_2-4-8_avg_traces.csv"
        current = -0.06 #nA
    elif neuron_name == "DNp03":
        path = "ephysData/DNp03_ephysData/DNp03_flies1-8_avg_traces.csv"
        current = -0.002 #nA
    else:
        raise ValueError(f"Unknown neuron_name: {neuron_name}")

    # Load data
    df = pd.read_csv(path)
    time = df["Time"].values
    fly_columns = [col for col in df.columns if col != "Time"]

    # --- Injection parameters (set according to your protocol) ---
    delay = 250       # ms, current step onset
    injDur = 50      # ms, duration of injection
    stim_amp = current  # nA, injected current amplitude (adjust to match your stimobj.amp)

    results = {"Fly": [], "RestVm": [], "Rin": [], "Tau": [], "Ctotal": []}

    # Analyze each fly trace
    for fly in fly_columns:
        exp_vData = df[fly].values

        # Baseline before current injection
        baseline_window = (time > (delay - 250)) & (time < delay)
        baseline = exp_vData[baseline_window].mean()

        # Steady state during injection (midpoint ± 12.5 ms)
        mid_start = delay + (injDur / 2) - 12.5
        mid_end   = delay + (injDur / 2) + 12.5
        ss_window = (time > mid_start) & (time < mid_end)
        steady_state = exp_vData[ss_window].mean()

        deltaV = steady_state - baseline  # mV
        deltaI = stim_amp                 # nA
        Rin = deltaV / deltaI if deltaI != 0 else np.nan  # MΩ

        # Fit exponential to tau during injection phase
        rise_window = (time > 250) & (time < 260)
        t_rise = time[rise_window] - delay   # subtract delay to align with injection onset
        v_rise = exp_vData[rise_window]

        try:
            popt, _ = curve_fit(exp_decay, t_rise, v_rise, p0=[deltaV, 20, baseline])
            tau = popt[1]  # ms
        except RuntimeError:
            tau = np.nan

        # Capacitance (pF) from tau / Rin
        C_total = tau / Rin * 1000 if (not np.isnan(tau) and Rin != 0) else np.nan

        # Store results
        results["Fly"].append(fly)
        results["RestVm"].append(baseline-13) # adjust for LJP of 13 mV
        results["Rin"].append(Rin)
        results["Tau"].append(tau)
        results["Ctotal"].append(C_total)

        print(f"{neuron_name} {fly}: Vm={baseline:.2f} mV, Rin={Rin:.2f} MΩ, tau={tau:.2f} ms, C={C_total:.2f} pF")

    # --- Convert to DataFrame ---
    res_df = pd.DataFrame(results)

    # --- Plot scatter + box for each property ---
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    props = [("RestVm", "Resting Vm (mV)"),
             ("Rin", "Input Resistance (MΩ)"),
             ("Tau", "Tau (ms)"),
             ("Ctotal", "Capacitance (pF)")]

    for ax, (prop, label) in zip(axes, props):
        data = res_df[prop].values
        # Add horizontal jitter
        x_jitter = 1 + 0.3 * (np.random.rand(len(data)) - 0.5)  # ±0.025 jitter
        ax.scatter(x_jitter, data, c='k', edgecolors='none', s=50, alpha=0.7)
        ax.boxplot(data, positions=[1], widths=0.3, showfliers=False)
        ax.set_title(label)
        ax.set_xticks([])
        ax.grid(True, linestyle="--", alpha=0.5)

    fig.suptitle(f"Experimental properties: {neuron_name}", fontsize=14)
    plt.tight_layout()
    plt.show()

    return res_df

def plot_with_std(ax, t, v, std, color='black', label=None):
                if std is not None and len(std) == len(v):
                    ax.fill_between(t, v - std, v + std, color='gray', alpha=0.4, edgecolor='none')
                ax.plot(t, v, color=color, label=label)
                
def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C

def runSim(sectionList_py, electrodeSec, somaSection, 
           exp_tData=None, exp_vData=None, STD=None, 
           current=-0.04, erev=-72.5, continueRun=1200, injDur=1000, delay=100):

    stimobj = h.IClamp(somaSection(0.5))
    stimobj.delay = delay
    stimobj.dur = injDur
    stimobj.amp = current

    ampInjvect = h.Vector().record(stimobj._ref_i)
    vInjVec = h.Vector().record(somaSection(0.5)._ref_v)
    tInjVec = h.Vector().record(h._ref_t)
    eInjVec = h.Vector().record(electrodeSec(0.5)._ref_v)

    h.finitialize(erev)
    h.continuerun(continueRun)

    # Convert to numpy arrays
    aInjVec_np = np.array(ampInjvect)
    vInj_np = np.array(vInjVec)
    tInj_np = np.array(tInjVec)
    eInjVec_np = np.array(eInjVec)

    # --- Simulation measurements ---
    baseline_window = (tInj_np > (delay - 100)) & (tInj_np < delay)
    baseline = vInj_np[baseline_window].mean()

    mid_start = delay + (injDur / 2) - 12.5
    mid_end   = delay + (injDur / 2) + 12.5
    ss_window = (tInj_np > mid_start) & (tInj_np < mid_end)
    steady_state = vInj_np[ss_window].mean()

    deltaV = steady_state - baseline  # mV
    deltaI = stimobj.amp              # nA
    Rin = deltaV / deltaI             # MΩ

    # Estimate tau
    rise_window = (tInj_np > 250) & (tInj_np < 260)
    t_rise = tInj_np[rise_window] - delay
    v_rise = vInj_np[rise_window]

    try:
        popt, _ = curve_fit(exp_decay, t_rise, v_rise, p0=[deltaV, 20, baseline])
        tau = popt[1]
    except RuntimeError:
        tau = np.nan

    C_total = tau / Rin * 1000 if not np.isnan(tau) else np.nan
    print(f"Simulation: Rin_MOhm: {Rin:.3f}, tau_ms: {tau:.3f}, C_total_pF: {C_total:.3f}")

    #############################################################
    ### plotting + RMSE calculation if exp data is provided #####
    #############################################################
    if exp_tData is not None and exp_vData is not None:
        rmse_decay, rmse_decay_ext = calculateRMSE_justDecay(exp_tData, exp_vData, tInj_np, eInjVec_np)
        rmse_recovery, rmse_recovery_ext = calculateRMSE_justRecovery(exp_tData, exp_vData, tInj_np, eInjVec_np)

        fig = plt.figure(figsize=(15, 7))
        gs = gridspec.GridSpec(3, 3, height_ratios=[0.4, 0.5, 0.1])

        def plot_with_std(ax, t, v, std, color='black', label=None):
            if std is not None and len(std) == len(v):
                ax.fill_between(t, v - std, v + std, color='gray', alpha=0.4, edgecolor='none')
            ax.plot(t, v, color=color, label=label)
        # --- Downsample simulation to match experimental points if needed ---
        factor = len(tInj_np) // len(exp_vData)
        if factor > 1:
            tInj_np = tInj_np[::factor]
            vInj_np = vInj_np[::factor]
            eInjVec_np = eInjVec_np[::factor]
            aInjVec_np = aInjVec_np[::factor]

        # Trim if off by one sample
        if len(tInj_np) > len(exp_vData):
            tInj_np   = tInj_np[:len(exp_vData)]
            vInj_np   = vInj_np[:len(exp_vData)]
            eInjVec_np = eInjVec_np[:len(exp_vData)]
            aInjVec_np = aInjVec_np[:len(exp_vData)]
        elif len(exp_vData) > len(tInj_np):
            exp_vData = exp_vData[:len(tInj_np)]
            if STD is not None:
                STD = STD[:len(tInj_np)]
        baseline_window = (tInj_np > (delay - 100)) & (tInj_np < delay)
        baseline = exp_vData[baseline_window].mean()

        mid_start = delay + (injDur / 2) - 12.5
        mid_end   = delay + (injDur / 2) + 12.5
        ss_window = (tInj_np > mid_start) & (tInj_np < mid_end)
        steady_state = exp_vData[ss_window].mean()

        deltaV = steady_state - baseline  # mV
        deltaI = stimobj.amp              # nA
        Rin = deltaV / deltaI

        rise_window = (tInj_np > 250) & (tInj_np < 260)
        t_rise = tInj_np[rise_window] - delay
        v_rise = exp_vData[rise_window]

        try:
            popt, _ = curve_fit(exp_decay, t_rise, v_rise, p0=[deltaV, 20, baseline])
            tau = popt[1]
        except RuntimeError:
            tau = np.nan

        C_total = tau / Rin * 1000 if not np.isnan(tau) else np.nan
        print(f"Experimental: Rin_MOhm: {Rin:.3f}, tau_ms: {tau:.3f}, C_total_pF: {C_total:.3f}")

        # Full trace
        ax1 = fig.add_subplot(gs[0:2, 0])
        plot_with_std(ax1, exp_tData, exp_vData, STD, color='black', label='Experimental')
        ax1.plot(tInj_np, eInjVec_np, color='red', label='Simulated')
        ax1.set_title('Full Voltage Trace')
        ax1.set_ylabel('mV')
        ax1.legend()
        ax1.tick_params(labelbottom=False)
        ax1.set_xlim(220, 320)
        # ax1.set_ylim(-72, -71.2)

        # Current injection
        ax2 = fig.add_subplot(gs[2, 0])
        ax2.plot(tInj_np, aInjVec_np, color='blue')
        ax2.set_title('Injected Current')
        ax2.set_xlabel('Time (ms)')
        ax2.set_ylabel('nA')
        ax2.set_xlim(220, 320)

        # Decay region
        ax3 = fig.add_subplot(gs[:, 1])
        decay_start, decay_end, decay_ext_end = 250.2, 253, 257

        # Main decay
        mask_exp = (exp_tData >= decay_start) & (exp_tData <= decay_end)
        mask_sim = (tInj_np >= decay_start) & (tInj_np <= decay_end)
        plot_with_std(ax3, exp_tData[mask_exp], exp_vData[mask_exp],
                      STD[mask_exp] if STD is not None else None, color='black')
        ax3.plot(tInj_np[mask_sim], eInjVec_np[mask_sim], color='red')

        mask_exp_ext = (exp_tData > decay_end) & (exp_tData <= decay_ext_end)
        mask_sim_ext = (tInj_np > decay_end) & (tInj_np <= decay_ext_end)
        # Extended decay
        mask_exp_ext = (exp_tData > decay_end) & (exp_tData <= decay_ext_end)
        if STD is not None and len(STD) == len(exp_vData):
            ax3.fill_between(
                exp_tData[mask_exp_ext],
                exp_vData[mask_exp_ext] - STD[mask_exp_ext],
                exp_vData[mask_exp_ext] + STD[mask_exp_ext],
                color='gray', alpha=0.4, edgecolor='none'
            )
        # Plot lines
        ax3.plot(exp_tData[mask_exp_ext], exp_vData[mask_exp_ext], color='black')
        ax3.plot(tInj_np[mask_sim_ext], eInjVec_np[mask_sim_ext], color='red')

        ax3.set_xlim(decay_start, decay_ext_end)
        ax3.set_title(f'Decay Region\nRMSE = {rmse_decay:.3f}, Ext = {rmse_decay_ext:.3f}')
        ax3.set_xlabel('Time (ms)')
        ax3.set_ylabel('mV')

        # Recovery region
        ax4 = fig.add_subplot(gs[:, 2])
        recovery_start, recovery_end, recovery_ext_end = 300.2, 303, 305

        # Main recovery
        mask_exp = (exp_tData >= recovery_start) & (exp_tData <= recovery_end)
        mask_sim = (tInj_np >= recovery_start) & (tInj_np <= recovery_end)
        plot_with_std(ax4, exp_tData[mask_exp], exp_vData[mask_exp],
                      STD[mask_exp] if STD is not None else None, color='black')
        ax4.plot(tInj_np[mask_sim], eInjVec_np[mask_sim], color='red')

        # Extended recovery
        mask_exp_ext = (exp_tData > recovery_end) & (exp_tData <= recovery_ext_end)
        mask_sim_ext = (tInj_np > recovery_end) & (tInj_np <= recovery_ext_end)
        if STD is not None and len(STD) == len(exp_vData):
            ax4.fill_between(
                exp_tData[mask_exp_ext],
                exp_vData[mask_exp_ext] - STD[mask_exp_ext],
                exp_vData[mask_exp_ext] + STD[mask_exp_ext],
                color='gray', alpha=0.4, edgecolor='none'
            )

        # Plot lines
        ax4.plot(exp_tData[mask_exp_ext], exp_vData[mask_exp_ext], color='black')
        ax4.plot(tInj_np[mask_sim_ext], eInjVec_np[mask_sim_ext], color='red')

        ax4.set_xlim(recovery_start, recovery_ext_end)
        ax4.set_title(f'Recovery Region\nRMSE = {rmse_recovery:.3f}, Ext = {rmse_recovery_ext:.3f}')
        ax4.set_xlabel('Time (ms)')
        ax4.set_ylabel('mV')

        plt.tight_layout()
        plt.show()

    return tInj_np, vInj_np, eInjVec_np, aInjVec_np


def main():
    hasElectrode=True
    neuron_name = "DNp03" # 

    # analyze_experimental_properties(neuron_name) # Analyze experimental properties for the given neuron

    fly_num = "avg" #DNp03 fly 3 and 5, DNp01 fly 4 and 6, 
    Tk().withdraw()
    fd_title = "Select morphology file to use for passive property fitting"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"morphologyData", title=fd_title)

    # morph_file = 'datafiles/morphologyData/' + neuron_name + '_morphData/' + neuron_name + '_um_model.swc'
    cell, allSections_py, allSections_nrn, somaSection, sizSection, axonEnd, erev, axonList, tetherList, dendList, electrodeSec = initializeModel(morph_file, neuron_name, hasElectrode)


    msw = plotMorphColorCode_wSIZ(allSections_py, somaSection, axonList, tetherList, dendList, sizSection, neuron_name)

    timestr = clock.strftime("%Y%m%d-%H%M%S")
    if neuron_name in ["DNp01", "DNp01_hemi"]:
        if fly_num == "4":
            DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray, STD = loadEphysData('ephysData/DNp01_ephysData/DNp01_avg_trace_fly4.dat') # higher rmp
        elif fly_num == "6":
            DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray, STD = loadEphysData('ephysData/DNp01_ephysData/DNp01_avg_trace_fly6.dat')# lower rmp
        elif fly_num == "avg":
            DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray, STD = loadEphysData('ephysData/DNp01_ephysData/DNp01_avg_trace_flies1_2-4-8.dat')
        expData = [DNp01_fly1thru6_60pA_hpol_noHold_timeArray, DNp01_fly1thru6_60pA_hpol_noHold_voltageArray]
        ljp = 13  # from gouwen and wilson 2009
        corrected_time, corrected_voltage = apply_liquid_junction_correction(expData, ljp)
        expData = [corrected_time, corrected_voltage]

    elif neuron_name in ["DNp03", "DNp03_hemi"]:
        if fly_num == "2":
            DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray, STD = loadEphysData('ephysData/DNp03_ephysData/DNp03_avg_trace_fly2.dat')# higher rmp
        elif fly_num == "6": 
            DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray, STD = loadEphysData('ephysData/DNp03_ephysData/DNp03_avg_trace_fly6.dat') # lower rmp
        elif fly_num == "avg":
            DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray, STD = loadEphysData('ephysData/DNp03_ephysData/DNp03_avg_trace_flies1-8.dat')
        expData = [DNp03_fly4567_2pA_hpol_noHold_timeArray, DNp03_fly4567_2pA_hpol_noHold_voltageArray]
        ljp = 13  # from gouwen and wilson 2009
        corrected_time, corrected_voltage = apply_liquid_junction_correction(expData, ljp)
        expData = [corrected_time, corrected_voltage]


    #new passive properties, after liquid junction potential correction, avg trace 
    if neuron_name == "DNp01" and fly_num == "avg": 
        #avg trace, flies 1,2,4-8
        initial = -76.75
        erev = -77.5
        raVal = 250
        gleakVal = 1/2675
        cmVal = 0.8
    elif neuron_name == "DNp01" and fly_num == "5":
        initial = -66.3
        erev = -66.94
        raVal = 300
        gleakVal =  1/3275
        cmVal = 0.79
    elif neuron_name == "DNp01" and fly_num == "3":
        initial = -87.3
        erev = -88.05
        raVal = 250
        gleakVal = 1/2820
        cmVal = 0.50
    elif neuron_name == "DNp01_hemi" and fly_num == "avg": 
        initial = -76.75
        erev = -77.5
        raVal = 190
        gleakVal = 1/1900
        cmVal = 0.85
    elif neuron_name == "DNp01_hemi" and fly_num == "5": 
        initial = -66.3
        erev = -66.95
        raVal = 220
        gleakVal =  1/2375
        cmVal = 0.85
    elif neuron_name == "DNp01_hemi" and fly_num == "3": 
        initial = -87.2
        erev = -88.05
        raVal = 210
        gleakVal =  1/1875
        cmVal = 0.55
    elif neuron_name == "DNp03" and fly_num == "avg":
        #avg trace, flies 1-8
        initial = -71.35
        raVal = 55
        gleakVal = 1/2725
        cmVal = 1.1
        erev = -73.71
    elif neuron_name == "DNp03" and fly_num == "2":
        initial = -68.4
        raVal = 50
        gleakVal = 1/3150
        cmVal = 0.8
        erev = -70.70
    elif neuron_name == "DNp03" and fly_num == "6":
        initial = -73.15
        raVal = 30
        gleakVal = 1/2700
        cmVal = 2
        erev = -74.8
    elif neuron_name == "DNp03_hemi" and fly_num == "avg":
        #avg trace, flies 1-8
        initial = -71.35
        raVal = 320
        gleakVal = 1/3755
        cmVal = 0.7
        erev = -73.73
    elif neuron_name == "DNp03_hemi" and fly_num == "2":
        initial = -67.75
        raVal = 350
        gleakVal = 1/3550
        cmVal = 0.7
        erev = -70.75
    elif neuron_name == "DNp03_hemi" and fly_num == "6":
        initial = -73.25
        raVal = 120
        gleakVal = 1/4500
        cmVal = 1.1
        erev = -74.85

    elec_raVal = 235.6                     
    elec_gleakVal = 0
    elec_cmVal = 6.4 # Membrane capacitance in micro Farads / cm^2
    #electode geom props: l = 10um | d = 1um

    sealCon_8GOhm = 0.0003978
    elec_gleakVal = sealCon_8GOhm


    change_Ra(ra=raVal, electrodeSec=electrodeSec, electrodeVal = elec_raVal)
    change_gLeak(gleak=gleakVal, erev=erev, electrodeSec=electrodeSec, electrodeVal = elec_gleakVal)
    change_memCap(memcap=cmVal, electrodeSec=electrodeSec, electrodeVal = elec_cmVal)
    nsegDiscretization(allSections_py)

    if neuron_name in ["DNp01", "DNp01_hemi"]:
        current = -0.06 #60 pA 
    elif neuron_name in ["DNp03", "DNp03_hemi"]:
        current = -0.002 # 2pA
    else:
        print("set current value")
        quit(0)
    
    #records from the electode section for fitting, Figure 6:
    tInj_np_hpol_justDecay, vInj_np_hpol_justDecay, eInj_np_hpol_justDecay, aInjVec_np_hpol_current = runSim(allSections_py, electrodeSec, somaSection, exp_tData=expData[0], exp_vData=expData[1], STD= STD, current=current, erev = initial, continueRun=550, injDur=50,  delay = 250)
    
main()