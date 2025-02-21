%By Alexander Vasserman

%%
clear all
close all
clear
%% load data fly1-6 hpol
folder = 'C:\Users\antho\VPN-DN-synapse-normalization\Matlab code';
filename = 'fly01_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat';
File = fullfile(folder, filename);
load(File)
fly1_hpol_responsePlots = [response_plot01, response_plot02, response_plot03,...
    response_plot04, response_plot05, response_plot06,...
    response_plot07, response_plot08, response_plot09,...
    response_plot10, response_plot11, response_plot12,...
    response_plot13, response_plot14, response_plot15,...
    response_plot16, response_plot17, response_plot18,...
    response_plot19, response_plot20];

t_hpol = 0:numel(response_plot01)-1;
t_hpol = t_hpol./20000;
t_hpol = t_hpol';


folder = 'C:\Users\antho\VPN-DN-synapse-normalization\Matlab code';
filename = 'fly02_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat';File = fullfile(folder, filename);
load(File)
fly2_hpol_responsePlots = [response_plot01, response_plot02, response_plot03,...
    response_plot04, response_plot05, response_plot06,...
    response_plot07, response_plot08, response_plot09,...
    response_plot10, response_plot11, response_plot12,...
    response_plot13, response_plot14, response_plot15,...
    response_plot16, response_plot17, response_plot18,...
    response_plot19, response_plot20];

folder = 'C:\Users\antho\VPN-DN-synapse-normalization\Matlab code';
filename = 'fly03_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat';
File = fullfile(folder, filename);
load(File)
fly3_hpol_responsePlots = [response_plot01, response_plot02, response_plot03,...
    response_plot04, response_plot05, response_plot06,...
    response_plot07, response_plot08, response_plot09,...
    response_plot10, response_plot11, response_plot12,...
    response_plot13, response_plot14, response_plot15,...
    response_plot16, response_plot17, response_plot18,...
    response_plot19, response_plot20];

% folder = 'C:\Users\antho\VPN-DN-synapse-normalization\Matlab code';
% filename = 'fly04_DNp01_hyperpolarization current injection 60pA_300ms_duration_with 5sec interval.mat';
% File = fullfile(folder, filename);
% load(File)
% fly4_hpol_responsePlots = [response_plot01, response_plot02, response_plot03,...
%     response_plot04, response_plot05, response_plot06,...
%     response_plot07, response_plot08, response_plot09,...
%     response_plot10, response_plot11, response_plot12,...
%     response_plot13, response_plot14, response_plot15,...
%     response_plot16, response_plot17, response_plot18,...
%     response_plot19, response_plot20];
% 
% folder = 'C:\Users\Sasha\Research\EphysDataPrep\DNp01 current injection 07-28-2023 upload\60 pA hyperpolarization or depolarization applied\5th';
% filename = 'fly05_DNp01_hyperpolarization current injection 60pA_300ms_duration_with 5sec interval.mat';
% File = fullfile(folder, filename);
% load(File)
% fly5_hpol_responsePlots = [response_plot01, response_plot02, response_plot03,...
%     response_plot04, response_plot05, response_plot06,...
%     response_plot07, response_plot08, response_plot09,...
%     response_plot10, response_plot11, response_plot12,...
%     response_plot13, response_plot14, response_plot15,...
%     response_plot16, response_plot17, response_plot18,...
%     response_plot19, response_plot20];
% 
% folder = 'C:\Users\Sasha\Research\EphysDataPrep\DNp01 current injection 08-03-2023 upload\60 pA hyperpolarization or depolarization applied\6th';
% filename = 'fly06_DNp01_hyperpolarization current injection 60pA_300ms_duration_with 5sec interval.mat';
% File = fullfile(folder, filename);
% load(File)
% fly6_hpol_responsePlots = [response_plot01, response_plot02, response_plot03,...
%     response_plot04, response_plot05, response_plot06,...
%     response_plot07, response_plot08, response_plot09,...
%     response_plot10, response_plot11, response_plot12,...
%     response_plot13, response_plot14, response_plot15,...
%     response_plot16, response_plot17, response_plot18,...
%     response_plot19, response_plot20];

t_hpol_300 = 0:numel(response_plot01)-1;
t_hpol_300 = t_hpol_300./20000;
t_hpol_300 = t_hpol_300';


%% make average traces
fly1_hpol_avg = returnAvgTrace(fly1_hpol_responsePlots);
fly2_hpol_avg = returnAvgTrace(fly2_hpol_responsePlots);
fly3_hpol_avg = returnAvgTrace(fly3_hpol_responsePlots);
% fly4_hpol_avg = returnAvgTrace(fly4_hpol_responsePlots);
% fly5_hpol_avg = returnAvgTrace(fly5_hpol_responsePlots);
% fly6_hpol_avg = returnAvgTrace(fly6_hpol_responsePlots);

% fly_45_72823_300ms = returnAvgTrace([fly4_72823_hpol_avg, fly5_72823_hpol_avg]);

%% quick plot of traces: note the different injection durations for some
%since all we care about for fit is preserve membrane decay/rise properties
%we artificially extend the shorter traces such that the decay and rise
%portions align (so they can be averaged), since those are consistent
figure(1);

time = t_hpol;
time_300ms = t_hpol_300;
colorGrey = [0.7, 0.7, 0.7];

subplot(3, 1, 1)
%note the subtraction of 7.5 seconds to facilitate alignment with the
%shorter trials so the injection duration begins at the same timepoint
plot(time-7.5, fly1_hpol_avg, 'r')
hold on
plot(time-7.5, fly2_hpol_avg, 'r')
plot(time-7.5, fly3_hpol_avg, 'r')
xlabel('sec')
ylabel('mV')
title("Fly 1-3 (longer inj)")

subplot(3, 1, 2)

% plot(time_300ms, fly4_hpol_avg, 'b')
% hold on
% plot(time_300ms, fly5_hpol_avg, 'b')
% plot(time_300ms, fly6_hpol_avg, 'b')
ylabel('mV')
xlabel('sec')
title("Fly 4-6 (shorter inj)")

subplot(3, 1, 3)
plot(time-7.5, fly1_hpol_avg, 'r')
hold on
plot(time-7.5, fly2_hpol_avg, 'g')
plot(time-7.5, fly3_hpol_avg, 'b')
% plot(time_300ms, fly4_hpol_avg, 'b')
% plot(time_300ms, fly5_hpol_avg, 'b')
% plot(time_300ms, fly6_hpol_avg, 'b')
ylabel('mV')
xlabel('sec')
title("Fly 1-6")

sgtitle("Quick plot of traces")

%% plot individual flies variability
figure(2);
plot(time-7.5, fly1_hpol_responsePlots,'Color', colorGrey)
hold on
plot(time-7.5, fly1_hpol_avg, 'r')
xlabel('sec')
ylabel('mV')
xlim([2 4])
title("Fly 1 (longer inj)")

figure(3);
plot(time-7.5, fly2_hpol_responsePlots, 'Color', colorGrey)
hold on
plot(time-7.5, fly2_hpol_avg, 'g')
xlabel('sec')
ylabel('mV')
xlim([2 4])
title("Fly 2 (longer inj)")

figure(4);
plot(time-7.5, fly3_hpol_responsePlots, 'Color', colorGrey)
hold on
plot(time-7.5, fly3_hpol_avg, 'b')
xlabel('sec')
ylabel('mV')
xlim([2 4])
title("Fly 3 (longer inj)")

%% fix bridge error
%bridge error is manually removed using indices determined by visual
%inspection

%(NOTE: ALL DONE WITHOUT BIAS AKA HOLDING CURRENT)

%both sets of flies (1-3 and 4-6) had bridge error, but due to diff inj
%durations, bridge error was corrected separately for avg traces of each
%set

%lots of commented plot code, left for posterity, can be removed

injStart_short = find(time_300ms == 2.3);
injEnd_short = find(time_300ms == 3);
injStart_long = find(time == 9.8);
injEnd_long = find(time == 11.2);

%fly 1-3 bridge error correction

avg_123_unbias = returnAvgTrace([fly1_hpol_avg, fly2_hpol_avg, fly3_hpol_avg]);
%separate into before and after drop 
PART1 = avg_123_unbias(1:200003);
PART2 = avg_123_unbias(200004:420000);
%calculate the voltage drop due to bridge error from avg traces
BRIDGE_DROP = avg_123_unbias(200003) - avg_123_unbias(200001);
%account for bridge drop + stich before and after sections of trace
avg_123_unbias_bridgeFix = [PART1(1:end-2); PART1(end-2); PART1(end-2); PART1(end-2); PART1(end-2); PART1(end-2); PART2(4:19998)+abs(BRIDGE_DROP); PART2(19998)+abs(BRIDGE_DROP); PART2(19998)+abs(BRIDGE_DROP); PART2(19998)+abs(BRIDGE_DROP); PART2(20002:end)];


%for 4-6, the traces need to be extended as well so the decays and rises
%can align

%first, calculate during injection indices for flies 4-6

% durInj_fly4_t = time_300ms(injStart_short:injEnd_short);
% durInj_fly5_t = time_300ms(injStart_short:injEnd_short);
% durInj_fly6_t = time_300ms(injStart_short:injEnd_short);

durInj_fly1_v = fly1_hpol_avg(injStart_long:injEnd_long);
durInj_fly2_v = fly2_hpol_avg(injStart_long:injEnd_long);
durInj_fly3_v = fly3_hpol_avg(injStart_long:injEnd_long);
% durInj_fly4_v = fly4_hpol_avg(injStart_short:injEnd_short);
% durInj_fly5_v = fly5_hpol_avg(injStart_short:injEnd_short);
% durInj_fly6_v = fly6_hpol_avg(injStart_short:injEnd_short);

%extend traces using custom function (see end of script)
% [extend_fly4_t, extend_fly4_v] = extend300msTrace(durInj_fly4_v, durInj_fly4_t);
% [extend_fly5_t, extend_fly5_v] = extend300msTrace(durInj_fly5_v, durInj_fly5_t);
% [extend_fly6_t, extend_fly6_v] = extend300msTrace(durInj_fly6_v, durInj_fly6_t);

%now correct for bridge error in 4-6 using same protocol as for 1-3
% avg_456_unbias = returnAvgTrace([extend_fly4_v, extend_fly5_v, extend_fly6_v]);
% P1_S_END = find(extend_fly4_t+7.5==time(200003));
% PART1_SHORT = avg_456_unbias(1:P1_S_END);
% PART2_SHORT = avg_456_unbias(P1_S_END+1:28000);
% BRIDGE_DROP_456 = avg_456_unbias(P1_S_END) - avg_456_unbias(P1_S_END-2);
%pretty sure the commented line below was my initial, incorrect, attempt as
%restiching the traces: i cannot remember why i did not delete it, so it
%might be useful for some unknown future purpose, but it can likely be
%removed
% avg_456_unbias_bridgeFix = [PART1_SHORT(1:end-2); PART1_SHORT(end-2); PART1_SHORT(end-2); PART2_SHORT(1:19998)+abs(BRIDGE_DROP_456); PART2_SHORT(19998)+abs(BRIDGE_DROP); PART2_SHORT(19998)+abs(BRIDGE_DROP_456);...
%     PART2_SHORT(19998)+abs(BRIDGE_DROP_456); PART2_SHORT(19998)+abs(BRIDGE_DROP_456); PART2_SHORT(20002:end)];

% avg_456_unbias_bridgeFix = [PART1_SHORT(1:end-2); PART1_SHORT(end-2);...
%     PART1_SHORT(end-2); PART1_SHORT(end-2); PART1_SHORT(end-2);...
%     PART2_SHORT(3:19998)+abs(BRIDGE_DROP_456);...
%     PART2_SHORT(19998)+abs(BRIDGE_DROP_456);...
%     PART2_SHORT(19998)+abs(BRIDGE_DROP_456);...
%     PART2_SHORT(19998)+abs(BRIDGE_DROP_456);...
%     PART2_SHORT(19998)+abs(BRIDGE_DROP_456);...
%     PART2_SHORT(19998)+abs(BRIDGE_DROP_456); PART2_SHORT(20003:end)];


%create master avg trace of flies 1-6 using bridge fixed 1-3 avg trace and
%extended, bridge fixed 4-6 avg trace 
% avg_1thru6_BF = returnAvgTrace([avg_123_unbias_bridgeFix(injStart_long:injEnd_long), avg_456_unbias_bridgeFix]);

%% plot bridge fixed traces

figure(5);
subplot(3, 3, 1)

% colorGrey = [0.7, 0.7, 0.7];
% plot(time, avg_123_unbias_bridgeFix, 'r')
% hold on
% plot(extend_fly6_t+7.5, avg_456_unbias_bridgeFix, 'b')
% plot(extend_fly6_t+7.5, avg_1thru6_BF, 'k')

plot(time, fly1_hpol_avg, 'color', colorGrey);
xline(10+0.25/1000, 'g')
plot(time, fly2_hpol_avg, 'color', colorGrey);
plot(time, fly3_hpol_avg,'color', colorGrey);
% plot(extend_fly4_t+7.5, extend_fly4_v, 'color', colorGrey)
% plot(extend_fly5_t+7.5, extend_fly5_v, 'color', colorGrey)
% plot(extend_fly6_t+7.5, extend_fly6_v, 'color', colorGrey)

legend('fly1-3 bridge fixed', 'fly4-6 bridge fixed', 'fly1-6 bridge fix', 'indiv traces 1-6 (no bridge fix)', 'decay start used in python sim')
xlim([9.995, 10.005])
xlabel('sec')
ylabel('mV')
title('decay phase (zoomed) | bridge fixed 1-3 and 4-6 with indiv (non-bridge fixed) traces')

subplot(3, 3, 2)
plot(time, avg_123_unbias_bridgeFix, 'r')
hold on
plot(time, avg_123_unbias, 'k');
xlim([9.999, 10.01])
xlabel('sec')
ylabel('mV')
legend("bridge fix", "no bridge fix")
title("fly 1-3 bridge fix vs no bridge fix")

% subplot(3,3,3)
% plot(extend_fly6_t+7.5, avg_456_unbias_bridgeFix, 'b')
% hold on
% plot(extend_fly4_t+7.5, avg_456_unbias, 'k');
% xlim([9.999, 10.01])
% xlabel('sec')
% ylabel('mV')
% legend("bridge fix", "no bridge fix")
% title("fly 4-6 bridge fix vs no bridge fix")

% plot(extend_fly4_t+7.5, avg_456_unbias, 'r', 'LineWidth', 1, 'LineWidth', 2)
% % plot(extend_fly6_t+7.5, avg_456_unbias, 'b')
% plot(extend_fly6_t(20002:end)+7.5, avg_456_unbias_bridgeFix(20002:end), 'g')
% plot(time_300ms+7.5, fly4_hpol_avg, 'b')
% plot(time_300ms+7.5, fly5_hpol_avg, 'b')
% plot(time_300ms+7.5, fly6_hpol_avg, 'b')

subplot(3, 3, 4)
colorGrey = [0.7, 0.7, 0.7];
plot(time, avg_123_unbias_bridgeFix, 'r')
hold on
% plot(extend_fly6_t+7.5, avg_456_unbias_bridgeFix, 'b')
plot(extend_fly6_t+7.5, avg_1thru6_BF, 'k')

plot(time, fly1_hpol_avg, 'color', colorGrey);
plot(time, fly2_hpol_avg, 'color', colorGrey);
plot(time, fly3_hpol_avg,'color', colorGrey);
% plot(extend_fly4_t+7.5, extend_fly4_v, 'color', colorGrey)
% plot(extend_fly5_t+7.5, extend_fly5_v, 'color', colorGrey)
% plot(extend_fly6_t+7.5, extend_fly6_v, 'color', colorGrey)

legend('fly1-3 bridge fixed', 'fly4-6 bridge fixed', 'fly1-6 bridge fix', 'indiv traces 1-6 (no bridge fix)', 'Location', 'Southeast')
xlim([10.995, 11.005])
xlabel('sec')
ylabel('mV')
title('rise phase (zoomed) | bridge fixed 1-3 and 4-6 with indiv (non-bridge fixed) traces')

subplot(3,3,5)
plot(time, avg_123_unbias_bridgeFix, 'r')
hold on
plot(time, avg_123_unbias, 'k');
xlim([10.995, 11.005])
xlabel('sec')
ylabel('mV')
legend("bridge fix", "no bridge fix")
title("fly 1-3 bridge fix vs no bridge fix")
% plot(time, avg_123_unbias, 'r', 'LineWidth', 1, 'LineWidth', 2)
% plot(time, avg_123_unbias_bridgeFix, 'k');
% plot(extend_fly6_t(24003:end)+7.5, avg_456_unbias_bridgeFix(24003:end), 'g')
% plot(extend_fly5_t+7.5, avg_456_unbias, 'm')
% plot(time, fly1_hpol_avg, 'r');
% plot(time, fly2_hpol_avg, 'r');
% plot(time, fly3_hpol_avg, 'r');
% plot(extend_fly4_t+7.5, extend_fly4_v, 'b')
% plot(extend_fly5_t+7.5, extend_fly5_v, 'b')
% plot(extend_fly6_t+7.5, extend_fly6_v, 'b')

% subplot(3,3,6)
% plot(extend_fly6_t+7.5, avg_456_unbias_bridgeFix, 'b')
% hold on
% plot(extend_fly4_t+7.5, avg_456_unbias, 'k');
% xlim([10.995, 11.005])
% xlabel('sec')
% ylabel('mV')
% legend("bridge fix", "no bridge fix")
% title("fly 4-6 bridge fix vs no bridge fix")

subplot(3, 1, 3)

plot(time, avg_123_unbias_bridgeFix, 'r')
% hold on
% plot(extend_fly4_t+7.5, avg_456_unbias_bridgeFix, 'b');
% plot(extend_fly6_t+7.5, avg_1thru6_BF, 'k')

xlim([9.9, 11.1])
xlabel('sec')
ylabel('mV')
legend("1-3 bridge fix", "4-6 bridge fix", "1-6 bridge fix")
title("all bridge fix traces")
sgtitle("fly 1-3, fly4-6, fly1-6 bridge fixed (without bias current) average traces at various zoom levels")

%uncomment line to write data to .dat file (NOTE: line currently written
%for bridgefix 1-6, but 1-3 is what was used for our actual sims
% writeEphysData(extend_fly6_t+7.5, avg_1thru6_BF, 'DNp01_hp_withoutBiasCurrent_avg_fly1thru6_BRIDGE_FIXED.dat')

%% calculate indices used for curve fitting

%if you want to calculate taus, use the indices to index out the data 
%(TOFIT_123, TOFIT_456, TOFIT_1thru6 or the equivlanet indices for any
%other indiv or avg traces) and then use curve fitter to fit it to the 
%correct exponential form voltage = a*exp(-t/b)+c where b is the tau and t
%is time (make sure to center the time vector you fit to so that point 1 is
%at 0)

STARTIDX_123 = find(time >= 10);
STARTIDX_123 = STARTIDX_123(1)
ENDIDX_123 = find(time >= 10.02);
ENDIDX_123 = ENDIDX_123(1)

STARTIDX_456 = find(extend_fly4_t+7.5 >= 10);
STARTIDX_456 = STARTIDX_456(1)
ENDIDX_456 = find(extend_fly4_t+7.5 >= 10.02);
ENDIDX_456 = ENDIDX_456(1)

STARTIDX_1thru6 = find(extend_fly6_t+7.5 >= 10);
STARTIDX_1thru6 = STARTIDX_1thru6(1)
ENDIDX_1thru6 = find(extend_fly6_t+7.5 >= 10.02);
ENDIDX_1thru6 = ENDIDX_1thru6(1)

TOFIT_123 = avg_123_unbias_bridgeFix(STARTIDX_123:ENDIDX_123);
TOFIT_456 = avg_456_unbias_bridgeFix(STARTIDX_456:ENDIDX_456);
TOFIT_1thru6 = avg_1thru6_BF(STARTIDX_1thru6:ENDIDX_1thru6);

%% junk plot code

% plot(extend_fly4_t+7.5, extend_fly4_v, 'color', colorGrey)
% plot(extend_fly5_t+7.5, extend_fly5_v, 'color', colorGrey)
% plot(extend_fly6_t+7.5, extend_fly6_v, 'color', colorGrey)
% plot(extend_fly4_t+7.5, avg_456_unbias, 'r', 'LineWidth', 1, 'LineWidth', 2)
% plot(extend_fly4_t+7.5, avg_456_unbias_bridgeFix, 'k');


% plot(time, avg_123_unbias, 'm', 'LineWidth', 1)
% plot(extend_fly6_t+7.5, avg_456_unbias_bridgeFix, 'b')
% plot(time, fly1_hpol_avg, 'color', colorGrey);
% % hold on
% plot(time, fly2_hpol_avg, 'color', colorGrey);
% plot(time, fly3_hpol_avg,'color', colorGrey);
% plot(time, avg_123_unbias, 'r', 'LineWidth', 1, 'LineWidth', 2)
% plot(time, avg_123_unbias_bridgeFix, 'k');
% plot(extend_fly6_t(24003:end)+7.5, avg_456_unbias_bridgeFix(24003:end), 'g')
% plot(extend_fly5_t+7.5, avg_456_unbias, 'm')
% plot(time, fly1_hpol_avg, 'r');
% plot(time, fly2_hpol_avg, 'r');
% plot(time, fly3_hpol_avg, 'r');
% plot(extend_fly4_t+7.5, extend_fly4_v, 'b')
% plot(extend_fly5_t+7.5, extend_fly5_v, 'b')
% plot(extend_fly6_t+7.5, extend_fly6_v, 'b')

% plot(time, avg_123_unbias_bridgeFix, 'r');
% plot(extend_fly6_t+7.5, avg_1thru6_BF, 'k')
% 
% legend("avg 1-3", "fly1", "fly2", "fly3")

% FULLAVG = returnAvgTrace([injHpol_fly1_wHold_center-fly1_hpol_wHold_offset,...
%     injHpol_fly2_wHold_center-fly2_hpol_wHold_offset,...
%     injHpol_fly3_wHold_center-fly3_hpol_wHold_offset,...
%     injHpol_fly4_wHold_center-fly4_hpol_wHold_offset,...
%     injHpol_fly6_wHold_center-fly6_hpol_wHold_offset,...
%     injHpol_fly7_wHold_center-fly7_hpol_wHold_offset]);
% 
% plot(injTime_center, FULLAVG)

% avg_fly4_hpol = returnAvgTrace(fly4_hpol_responsePlots);
% figure(1);
% subplot(2, 1, 1)
% plot(t, avg_fly4_hpol)
% xlim([9.95, 10.05])
% subplot(2, 1, 2)
% plot(t, avg_fly4_hpol)
% xlim([10, 10.05])
% sgtitle("fly 4 no bias average trace at various zooms")

%% junk code for plotting voltage deflection

% fly1_noBias_ind = [];
% fly2_noBias_ind = [];
% fly3_noBias_ind = [];
% 
% for i = 1:20
%     fly1_noBias_ind(i) = calcAmp(fly1_hpol_responsePlots(:, i));
%     fly2_noBias_ind(i) = calcAmp(fly2_hpol_responsePlots(:, i));
%     fly3_noBias_ind(i) = calcAmp(fly3_hpol_responsePlots(:, i));
% end
% 
% fly1_deflection_noBias = calcAmp(fly1_hpol_avg);
% fly2_deflection_noBias = calcAmp(fly2_hpol_avg);
% fly3_deflection_noBias = calcAmp(fly3_hpol_avg);
% 
% 
% figure(3);
% boxplot([fly1_noBias_ind',...
%     fly2_noBias_ind', fly3_noBias_ind'],...
%     'Labels', {'f1 no bias', 'f2 with bias', 'f3 no bias'})
% ylabel('deflection from resting (mV)')
% hold on
% plot(1, fly1_deflection_noBias, 'kd', 'LineWidth', 2)
% plot(2, fly2_deflection_noBias, 'kd', 'LineWidth', 2)
% plot(3, fly3_deflection_noBias, 'kd', 'LineWidth', 2)

% writeEphysData(t, avg_123_unbias_copy_FIXMAN, 'DNp01_7423_hp_withoutBiasCurrent_avg_fly123_BRIDGE_FIXED.dat')

%% functions
function deflection = calcAmp(trace)
restAmp = mean(trace(1:100000));
durInjAmp = mean(trace(200001:220001));
deflection = restAmp-durInjAmp;
end

function writeEphysData(xData, yData, filename)
dataMat = [xData, yData];
formatSpec = '%2.5f %4.2f\n';
fileID = fopen(filename, 'w');
for i = 1:numel(xData)
    temp = dataMat(i, :);
    fprintf(fileID, formatSpec, temp);
end
end

function avgTrace = returnAvgTrace(traceList)
sumTraces = sum(traceList, 2);
avgTrace = sumTraces./size(traceList, 2);
end

function [extendedTime, extendedTrace] = extend300msTrace(shortTrace, time)
idxStart = find(time == 2.55);
idxEnd = find(time == 2.75);
extensionVoltVal = mean(shortTrace(idxStart:idxEnd));
extensionTime = (1:14000).*(1/20000);
extendedTime = [time(1:idxStart); extensionTime'+2.55; time(idxStart+1:end)+extensionTime(end)];
close all
% plot(1:5001, extendedTime(1:5001), 'r');
% hold on
% plot(5002:19002, extendedTime(5002:19002), 'b');
% plot(19003:28001, extendedTime(19003:end), 'g');
extendedTrace = [shortTrace(1:idxStart); ones(14000, 1).*extensionVoltVal; shortTrace(idxStart+1:end)];
% plot(extendedTime, extendedTrace)
end
