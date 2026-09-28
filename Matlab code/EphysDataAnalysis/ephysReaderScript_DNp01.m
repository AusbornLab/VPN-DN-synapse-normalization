%By Alexander Vasserman/Anthony Moreno-Sanchez

%%
clear all
close all
clear
%% load data fly1-5 hpol
folder = pwd;

% [fly1_hpol_responsePlots, t_hpol] = loadFlyData_old(folder, 'fly01_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
% [fly2_hpol_responsePlots, ~] = loadFlyData_old(folder, 'fly02_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
% [fly3_hpol_responsePlots, ~] = loadFlyData_old(folder, 'fly03_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
% [fly4_hpol_responsePlots, ~] = loadFlyData(folder, 'fly04_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
% [fly5_hpol_responsePlots, ~] = loadFlyData(folder, 'fly05_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
[fly1_raw, ~] = loadFlyData_old(folder, 'fly01_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
[fly2_raw, ~] = loadFlyData_old(folder, 'fly02_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');

% Apply bridge correction to each trace in raw data (assuming each column is a trace)
fly1_raw_corrected = zeros(size(fly1_raw));
for i = 1:size(fly1_raw, 2)
    fly1_raw_corrected(:, i) = applyBridgeFix(fly1_raw(:, i));
end

fly2_raw_corrected = zeros(size(fly2_raw));
for i = 1:size(fly2_raw, 2)
    fly2_raw_corrected(:, i) = applyBridgeFix(fly2_raw(:, i));
end

% Reshape corrected data
[fly1_hpol_responsePlots, t_hpol] = reshapeFlyResponseToRows(fly1_raw_corrected);
[fly2_hpol_responsePlots, ~] = reshapeFlyResponseToRows(fly2_raw_corrected);

% [fly3_raw, ~] = loadFlyData_old(folder, 'fly03_DNp01_hyperpolarization current injection 60pA_with 20sec interval.mat');
% [fly3_hpol_responsePlots, ~] = reshapeFlyResponseToRows(fly3_raw);

[fly4_hpol_responsePlots, fly5_hpol_responsePlots, fly6_hpol_responsePlots, fly7_hpol_responsePlots, fly8_hpol_responsePlots, t_hpol] = loadFlyData(folder, 'combined_fly_traces_DNp01.mat');

colors = {
    [1, 0, 0];        % red
    [0, 1, 1];        % cyan
    [0, 1, 0];        % green
    [0, 0, 1];        % blue
    [0.54, 0, 0.54];  % violet
    [1, 1, 0];        % yellow
    [0, 0, 0];        % black
    [1, 0, 1];        % magenta
};
%% make average traces
[fly1_hpol_avg, STD1] = returnAvgTrace(fly1_hpol_responsePlots);
[fly2_hpol_avg, STD2] = returnAvgTrace(fly2_hpol_responsePlots);
% [fly3_hpol_avg, STD3] = returnAvgTrace(fly3_hpol_responsePlots);
[fly4_hpol_avg, STD4] = returnAvgTrace(fly4_hpol_responsePlots);
[fly5_hpol_avg, STD5] = returnAvgTrace(fly5_hpol_responsePlots);
[fly6_hpol_avg, STD6] = returnAvgTrace(fly6_hpol_responsePlots);
[fly7_hpol_avg, STD7] = returnAvgTrace(fly7_hpol_responsePlots);
[fly8_hpol_avg, STD8] = returnAvgTrace(fly8_hpol_responsePlots);

time = t_hpol;
colorGrey = [0.7, 0.7, 0.7];

% %% Plotting of average traces per fly
figure(1);
plot(time, fly1_hpol_avg, 'Color', colors{1}, 'DisplayName', 'Fly 1');
plot(time, fly2_hpol_avg, 'Color', colors{2}, 'DisplayName', 'Fly 2');
% plot(time, fly3_hpol_avg, 'Color', colors{3}, 'DisplayName', 'Fly 3');
plot(time, fly4_hpol_avg, 'Color', colors{4}, 'DisplayName', 'Fly 4');
plot(time, fly5_hpol_avg, 'Color', colors{5}, 'DisplayName', 'Fly 5');
plot(time, fly6_hpol_avg, 'Color', colors{6}, 'DisplayName', 'Fly 6');
plot(time, fly7_hpol_avg, 'Color', colors{7}, 'DisplayName', 'Fly 7');
plot(time, fly8_hpol_avg, 'Color', colors{8}, 'DisplayName', 'Fly 8');
xlim([249 254])
title('Flies 1-8');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
legend('Fly 1', 'Fly 2', 'Fly 4', 'Fly 5', 'Fly 6','Fly 7','Fly 8');
hold off;

avgTraces = {fly1_hpol_avg; fly2_hpol_avg; fly4_hpol_avg; fly5_hpol_avg; fly6_hpol_avg; fly7_hpol_avg; fly8_hpol_avg};

tauValues = fitTauAndPlotBox(avgTraces, t_hpol, colors);
hold on;

% Plotting with STD per fly
xmin = 200;
xmax = 400;
plotWithStd(time, fly1_hpol_avg, STD1, colors{1}, "Fly 1", xmin, xmax);
plotWithStd(time, fly2_hpol_avg, STD2, colors{2}, "Fly 2", xmin, xmax);
% plotWithStd(time, fly3_hpol_avg, STD3, colors{3}, "Fly 3", xmin, xmax);
plotWithStd(time, fly4_hpol_avg, STD4, colors{4}, "Fly 4", xmin, xmax);
plotWithStd(time, fly5_hpol_avg, STD5, colors{5}, "Fly 5", xmin, xmax);
plotWithStd(time, fly6_hpol_avg, STD6, colors{6}, "Fly 6", xmin, xmax);
plotWithStd(time, fly7_hpol_avg, STD7, colors{7}, "Fly 7", xmin, xmax);
plotWithStd(time, fly8_hpol_avg, STD8, colors{8}, "Fly 8", xmin, xmax);

% Plotting all trials of a single fly with their average
plot_individual_traces_and_average(time, fly1_hpol_avg, fly1_hpol_responsePlots, 'Fly 1 trials', colors{1})
plot_individual_traces_and_average(time, fly2_hpol_avg, fly2_hpol_responsePlots, 'Fly 2 trials', colors{2})
% plot_individual_traces_and_average(time, fly3_hpol_avg, fly3_hpol_responsePlots, 'Fly 3 trials', colors{3})
plot_individual_traces_and_average(time, fly4_hpol_avg, fly4_hpol_responsePlots, 'Fly 4 trials', colors{4})
plot_individual_traces_and_average(time, fly5_hpol_avg, fly5_hpol_responsePlots, 'Fly 5 trials', colors{5})
plot_individual_traces_and_average(time, fly6_hpol_avg, fly6_hpol_responsePlots, 'Fly 6 trials', colors{6})
plot_individual_traces_and_average(time, fly7_hpol_avg, fly7_hpol_responsePlots, 'Fly 7 trials', colors{7})
plot_individual_traces_and_average(time, fly8_hpol_avg, fly8_hpol_responsePlots, 'Fly 8 trials', colors{8})

%% Baseline subtract the average traces and plot them
datapoints = 200:600;

initial_avg_per_trace = [
    mean(fly1_hpol_avg(datapoints)), mean(fly2_hpol_avg(datapoints)), mean(fly4_hpol_avg(datapoints)), mean(fly5_hpol_avg(datapoints)), mean(fly6_hpol_avg(datapoints)), mean(fly7_hpol_avg(datapoints)), mean(fly8_hpol_avg(datapoints))
];

% Compute overall initial average
overall_avg_initial = mean(initial_avg_per_trace);

% Normalized traces
fly1_hpol_norm = fly1_hpol_avg - initial_avg_per_trace(1) + overall_avg_initial;
fly2_hpol_norm = fly2_hpol_avg - initial_avg_per_trace(2) + overall_avg_initial;
% fly3_hpol_norm = fly3_hpol_avg - initial_avg_per_trace(3) + overall_avg_initial;
fly4_hpol_norm = fly4_hpol_avg - initial_avg_per_trace(3) + overall_avg_initial;
fly5_hpol_norm = fly5_hpol_avg - initial_avg_per_trace(4) + overall_avg_initial;
fly6_hpol_norm = fly6_hpol_avg - initial_avg_per_trace(5) + overall_avg_initial;
fly7_hpol_norm = fly7_hpol_avg - initial_avg_per_trace(6) + overall_avg_initial;
fly8_hpol_norm = fly8_hpol_avg - initial_avg_per_trace(7) + overall_avg_initial;

% Plot normalized traces
figure;
subplot(2,1,1);
hold on;
plot(time, fly1_hpol_avg, 'Color', colors{1}, 'DisplayName', 'Fly 1');
plot(time, fly2_hpol_avg, 'Color', colors{2}, 'DisplayName', 'Fly 2');
% plot(time, fly3_hpol_avg, 'Color', colors{3}, 'DisplayName', 'Fly 3');
plot(time, fly4_hpol_avg, 'Color', colors{4}, 'DisplayName', 'Fly 4');
plot(time, fly5_hpol_avg, 'Color', colors{5}, 'DisplayName', 'Fly 5');
plot(time, fly6_hpol_avg, 'Color', colors{6}, 'DisplayName', 'Fly 6');
plot(time, fly7_hpol_avg, 'Color', colors{7}, 'DisplayName', 'Fly 7');
plot(time, fly8_hpol_avg, 'Color', colors{8}, 'DisplayName', 'Fly 8');
xlim([249 254])
title('Original Data');
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
legend('Fly 1', 'Fly 2', 'Fly 4', 'Fly 5', 'Fly 6','Fly 7','Fly 8');
hold off;

subplot(2,1,2);
hold on;
plot(time, fly1_hpol_norm, 'Color', colors{1}, 'DisplayName', 'Fly 1');
plot(time, fly2_hpol_norm, 'Color', colors{2}, 'DisplayName', 'Fly 2');
% plot(time, fly3_hpol_norm, 'Color', colors{3}, 'DisplayName', 'Fly 3');
plot(time, fly4_hpol_norm, 'Color', colors{4}, 'DisplayName', 'Fly 4');
plot(time, fly5_hpol_norm, 'Color', colors{5}, 'DisplayName', 'Fly 5');
plot(time, fly6_hpol_norm, 'Color', colors{6}, 'DisplayName', 'Fly 6');
plot(time, fly7_hpol_norm, 'Color', colors{7}, 'DisplayName', 'Fly 7');
plot(time, fly8_hpol_norm, 'Color', colors{8}, 'DisplayName', 'Fly 8');
xlim([249 254])
title('Normalized Data');
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
legend('Fly 1', 'Fly 2', 'Fly 4', 'Fly 5', 'Fly 6','Fly 7','Fly 8');


fly_nobaseline = [fly1_hpol_avg; fly2_hpol_avg; fly4_hpol_avg; fly5_hpol_avg; fly6_hpol_avg; fly7_hpol_avg; fly8_hpol_avg];
fly_baseline = [fly1_hpol_norm; fly2_hpol_norm; fly4_hpol_norm; fly5_hpol_norm; fly6_hpol_norm; fly7_hpol_norm; fly8_hpol_norm];

[fly_nobaseline_avg_trace, STD_nobaseline] = returnAvgTrace(fly_nobaseline);
[fly_baseline_avg_trace, STD_baseline] = returnAvgTrace(fly_baseline);

plotWithStd(time, fly_nobaseline_avg_trace, STD_nobaseline, 'r', "Flies 1-8 not baseline subtracted", xmin, xmax);
plotWithStd(time, fly_baseline_avg_trace, STD_baseline, 'r', "Flies 1-8 baseline subtracted", xmin, xmax);

% fly3 was excluded from analysis due to low input resistance, flies 4 and
% 6 were fit considering upper and lower bound of resting membrane
% potentials and time constants
% writeEphysData_with_STD(time, fly_baseline_avg_trace, STD_baseline, 'DNp01_avg_trace_flies1_2-4-8.dat')
% writeEphysData_with_STD(time, fly4_hpol_avg, STD4, 'DNp01_avg_trace_fly4.dat')
% writeEphysData_with_STD(time, fly6_hpol_avg, STD6, 'DNp01_avg_trace_fly6.dat')
writeEphysData_toCSV(time, fly_nobaseline, 'DNp01_flies1_2-4-8_avg_traces.csv', [1 2 4 5 6 7 8]); 


flyTraces = {
    fly1_hpol_responsePlots
    fly2_hpol_responsePlots
    fly4_hpol_responsePlots
    fly5_hpol_responsePlots
    fly6_hpol_responsePlots
    fly7_hpol_responsePlots
    fly8_hpol_responsePlots
};

[results, summary] = analyze_experimental_properties( ...
    'DNp01', t_hpol, flyTraces, [1 2 4 5 6 7 8]);

%% functions
function writeEphysData_with_STD(xData, yData, STD, filename)
    % Ensure all inputs are column vectors
    xData = xData(:);
    yData = yData(:);
    STD   = STD(:);

    % Check for consistent vector lengths
    if numel(xData) ~= numel(yData) || numel(xData) ~= numel(STD)
        error('xData, yData, and STD must have the same number of elements.');
    end

    % Combine data into matrix
    dataMat = [xData, yData, STD];

    % Updated format specifier for three columns
    formatSpec = '%2.5f %4.5f %4.5f\n';

    % Write to file
    fileID = fopen(filename, 'w');
    for i = 1:numel(xData)
        fprintf(fileID, formatSpec, dataMat(i, :));
    end
    fclose(fileID);
end

function writeEphysData_toCSV(timeVec, traceMatrix, filename, flyNums)
    % Ensure timeVec is a column vector
    timeVec = timeVec(:);

    % Ensure traceMatrix rows match time length
    if size(traceMatrix, 2) ~= numel(timeVec) && size(traceMatrix, 1) ~= numel(timeVec)
        error('timeVec length must match one dimension of traceMatrix.');
    end

    % If traces are in rows, transpose so each trace is a column
    if size(traceMatrix, 1) ~= numel(timeVec)
        traceMatrix = traceMatrix.'; 
    end

    % Check flyNums matches number of traces
    if nargin < 4
        flyNums = 1:size(traceMatrix,2); % default: sequential
    elseif numel(flyNums) ~= size(traceMatrix,2)
        error('Number of flyNums must match number of traces (columns).');
    end

    % Combine time and data into one matrix
    dataMat = [timeVec, traceMatrix];

    % Create variable names: Time, FlyX for each fly number
    varNames = ["Time", "Fly" + string(flyNums)];

    % Convert to table for cleaner CSV export
    T = array2table(dataMat, 'VariableNames', varNames);

    % Write to CSV
    writetable(T, filename);

    fprintf('Data written to %s\n', filename);
end

function [avgTrace, S] = returnAvgTrace(traceList)
sumTraces = sum(traceList, 1);
avgTrace = sumTraces./size(traceList, 1);
S = std(traceList,0,1);
end

% old reading function
function [responsePlots, t_hpol] = loadFlyData_old(folder, filename)
    % Constructs full file path and loads the .mat file
    File = fullfile(folder, filename);
    data = load(File);
    
    % Extracts all response plot variables dynamically
    fieldNames = fieldnames(data);
    responsePlots = [];

    for i = 1:numel(fieldNames)
        if contains(fieldNames{i}, 'response_plot')
            responsePlots = [responsePlots, data.(fieldNames{i})]; % Concatenate response plots
        end
    end
    
    % Construct time vector if at least one response plot exists
    if ~isempty(responsePlots)
        t_hpol = (0:numel(responsePlots(:,1))-1) ./ 20000;
        t_hpol = t_hpol';
    else
        t_hpol = [];
    end
end

function [reshapedData, t_hpol] = reshapeFlyResponseToRows(dataLong)
% Reshapes long format fly response data into [trials x 11000] with 550 ms segments
% around a current injection event from 10s to 11s.
%
% INPUT:
%   dataLong: [numSamples x numTrials] data sampled at 20 kHz (i.e., 0.05 ms/sample)
%
% OUTPUT:
%   reshapedData: [numTrials x 11000] (each row is one trial of 550 ms)
%   t_hpol:       [1 x 11000] time vector in milliseconds, relative to start of slice

    [numSamples, numTrials] = size(dataLong);
    Fs = 20000;  % Sampling rate in Hz
    msPerSample = 1e3 / Fs;  % 0.05 ms

    % Time window around injection
    injStart_ms = 10000;   % 10 s
    injEnd_ms = 11000;     % 11 s

    % Define segments in ms
    preInj_ms = 250;
    injStartSlice_ms = 25;
    injEndSlice_ms = 25;
    postInj_ms = 250;

    % Convert to sample indices
    preInj_samples = round(preInj_ms / msPerSample);               % 5000
    injStartSlice_samples = round(injStartSlice_ms / msPerSample); % 500
    injEndSlice_samples = round(injEndSlice_ms / msPerSample);     % 500
    postInj_samples = round(postInj_ms / msPerSample);             % 5000

    % Injection sample indices
    injStart_idx = round(injStart_ms / msPerSample); % 200000
    injEnd_idx = round(injEnd_ms / msPerSample);     % 220000

    % Slices
    idx_preInj = (injStart_idx - preInj_samples):(injStart_idx - 1);
    idx_injStart = injStart_idx:(injStart_idx + injStartSlice_samples - 1);
    idx_injEnd = (injEnd_idx - injEndSlice_samples):(injEnd_idx - 1);
    idx_postInj = injEnd_idx:(injEnd_idx + postInj_samples - 1);

    % Final concatenated index
    allIndices = [idx_preInj, idx_injStart, idx_injEnd, idx_postInj]; % Total: 11000

    if max(allIndices) > numSamples
        error('Not enough samples after injection for slicing.');
    end

    % Reshape
    reshapedData = zeros(numTrials, length(allIndices));
    for i = 1:numTrials
        reshapedData(i, :) = dataLong(allIndices, i);
    end

    % Time vector in ms (relative to beginning of segment)
    t_hpol = (0:length(allIndices)-1) * msPerSample;
end

% old function
% function [responseData1, responseData2, responseData3, t_hpol] = loadFlyData(folder, filename)
%     % Constructs full file path and loads the .mat file
%     File = fullfile(folder, filename);
%     data = load(File);
%     
%     % Extracts all response variables dynamically
%     fieldNames = fieldnames(data);
%     
%     % Regular expression to match trial patterns
%     trialPattern = 'response_trail_(\d+)_\d+';
%     matchedFields = fieldNames(~cellfun(@isempty, regexp(fieldNames, trialPattern, 'tokens')));
%     
%     % Sort fields manually for consistent order
%     matchedFields = sort(matchedFields);
%     
%     % Initialize structures for separated data
%     responseData1 = [];
%     responseData2 = [];
%     responseData3 = [];
%     
%     % Separate data by trials
%     for i = 1:numel(matchedFields)
%         trialNumber = regexp(matchedFields{i}, trialPattern, 'tokens');
%         trialNumber = str2double(trialNumber{1}{1});
%         
%         dataMatrix = data.(matchedFields{i});
%         
%         if trialNumber == 1
%             responseData1 = [responseData1, dataMatrix]; % Concatenate along columns
%         elseif trialNumber == 2
%             responseData2 = [responseData2, dataMatrix];
%         elseif trialNumber == 3
%             responseData3 = [responseData3, dataMatrix];
%         end
%     end
%     
%     % Construct time vector where each point represents 1 ms
%     if ~isempty(responseData1)
%         t_hpol = (0:size(responseData1, 1) - 1)/2; % Time in milliseconds
%         t_hpol = t_hpol';
%     else
%         t_hpol = [];
%     end
% end

function [responseData1, responseData2, responseData3, responseData4, responseData5, t_hpol] = loadFlyData(folder, filename)
    % Constructs full file path and loads the .mat file
    File = fullfile(folder, filename);
    data = load(File);
    
    % Extract field names that match the pattern *_all_traces
    fieldNames = fieldnames(data);
    tracePattern = 'Fly\d+_all_traces';
    matchedFields = fieldNames(~cellfun(@isempty, regexp(fieldNames, tracePattern, 'once')));
    
    % Sort matched fields for consistent order (e.g., Fly01 before Fly03)
    matchedFields = sort(matchedFields);
    
    % Initialize output variables
    responseData1 = [];
    responseData2 = [];
    responseData3 = [];
    responseData4 = [];
    responseData5 = [];
    
    % Assign matched data to outputs
    for i = 1:numel(matchedFields)
        thisData = data.(matchedFields{i});
        
        switch i
            case 1
                responseData1 = thisData;
            case 2
                responseData2 = thisData;
            case 3
                responseData3 = thisData;
            case 4
                responseData4 = thisData;
            case 5
                responseData5 = thisData;
        end
    end

    % Generate time vector (based on number of columns = time samples)
    if ~isempty(responseData1)
        numTimePoints = size(responseData1, 2);
        t_hpol = (0:numTimePoints - 1) * 0.05;  % in milliseconds (20 kHz)
        t_hpol = t_hpol';
    else
        t_hpol = [];
    end
end

function tauValues = fitTauAndPlotBox(avgTraces, t_hpol, colors)
    % fitTauAndPlotBox - Fit exponential decay to multiple traces and plot boxplot of time constants
    %
    % Inputs:
    %   avgTraces - cell array of average traces (ordered for flies 1,2,4,5,6,7,8)
    %   t_hpol    - time vector (N x 1), in ms
    %   colors    - 1x8 cell array of color codes (full list, including unused color{3})
    %
    % Output:
    %   tauValues - array of fitted tau values

    tauValues = [];
    flyIndex = [];

    % These are the actual fly labels (i.e., indices into colors)
    flyLabels = [1, 2, 4, 5, 6, 7, 8];

    % Fit range (ms)
    fitStart = 250.2; 
    fitEnd = 260;

    figure(100); clf; hold on;
    title('Fit Range Data');
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');

    for i = 1:numel(avgTraces)
        flyNum = flyLabels(i);  % Actual fly label (for display + color indexing)
        trace = avgTraces{i};

        if isrow(trace), trace = trace'; end

        idx = t_hpol >= fitStart & t_hpol <= fitEnd;
        t_fit = t_hpol(idx);
        y_fit = trace(idx);

        if numel(t_fit) < 5 || range(y_fit) < 0.1
            warning('Fly %d: insufficient data variation for fitting.', flyNum);
            tauValues(end+1) = NaN;
            flyIndex(end+1) = flyNum;
            continue;
        end

        t_fit_shifted = t_fit - t_fit(1);
        plot(t_fit, y_fit, 'Color', colors{flyNum}, ...
             'DisplayName', sprintf('Fly %d (raw)', flyNum));

        startA = max(y_fit) - min(y_fit);
        startC = min(y_fit);
        decay_ratio = (y_fit(end) - startC) / startA;

        if decay_ratio <= 0
            startTau = 1;
        else
            decay_idx = find(y_fit <= startC + startA * exp(-1), 1);
            startTau = t_fit_shifted(decay_idx);
            if isempty(startTau), startTau = 1; end
        end

        expModel = fittype('A*exp(-x/tau) + C', ...
            'independent', 'x', ...
            'coefficients', {'A', 'tau', 'C'});

        try
            fitResult = fit(t_fit_shifted, y_fit, expModel, ...
                'StartPoint', [startA, startTau, startC], ...
                'Lower', [0, 0.001, -Inf], ...
                'Upper', [Inf, Inf, Inf]);

            tau = fitResult.tau;
            y_model = fitResult(t_fit_shifted);

            plot(t_fit, y_model, '--', 'Color', colors{flyNum}, ...
                'DisplayName', sprintf('Fly %d (fit)', flyNum));
        catch
            warning('Fit failed for Fly %d.', flyNum);
            tau = NaN;
        end

        tauValues(end+1) = tau;
        flyIndex(end+1) = flyNum;
    end

    legend show;

    % --- Boxplot with scatter ---
    figure('Name', 'Tau Boxplot'); clf;
    boxplot(tauValues, 'Colors', [0 0 0], 'Symbol', '');
    hold on;

    for i = 1:numel(tauValues)
        flyNum = flyIndex(i);
        scatter(1, tauValues(i), 60, colors{flyNum}, ...
            'filled', 'MarkerEdgeColor', 'k');
    end

    ylabel('\tau (ms)', 'FontSize', 12);
    title('Fitted Time Constants (\tau)', 'FontSize', 14);
    set(gca, 'XTickLabel', {'All Flies'}, 'FontSize', 12);
    grid on;
    hold off
end

function plot_individual_traces_and_average(time, avg_trace, all_traces, label_str, color_str)
    figure;
    % Plot all individual traces with light color
    hold on;
    plot(time, all_traces, 'Color', [0.8 0.8 0.8]); % Light gray for individual traces

    % Plot the average trace in specified color
    plot(time, avg_trace, 'Color', color_str, 'LineWidth', 2);

    % Formatting
    title(label_str, 'Interpreter', 'none');
    xlabel('Time');
    ylabel('Response');
    grid on;
    box on;
    hold off;
end

function plotWithStd(time, avg_trace, std_trace, color, titleText, xmin, xmax)
    % Function to plot average traces with shaded standard deviation bounds
    % Inputs:
    %   time       - Time vector
    %   avg_trace  - Average trace data
    %   std_trace  - Standard deviation data
    %   color      - Line color for the average trace
    %   titleText  - Title of the plot

    % Define color for shading
    colorGrey = [0.5, 0.5, 0.5];

    % Compute upper and lower bounds
    upper_bound = avg_trace + std_trace;
    lower_bound = avg_trace - std_trace;

    % Ensure vectors have the same orientation
    if isrow(time)
        time = time(:);
    end
    if isrow(upper_bound)
        upper_bound = upper_bound(:);
    end
    if isrow(lower_bound)
        lower_bound = lower_bound(:);
    end

    % Shift time for alignment
    time_shifted = time 

    % Create a new figure
    figure;
    
    % Fill area between upper and lower bounds
    fill([time_shifted; flipud(time_shifted)], ...
         [upper_bound; flipud(lower_bound)], ...
         colorGrey, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    hold on;

    % Plot the average trace
    plot(time_shifted, avg_trace, 'Color', color, 'LineWidth', 1.5);

    % Labels and formatting
    xlabel('Time (s)');
    ylabel('Membrane potential (mV)');
    xlim([xmin xmax]);
    title(titleText);
    hold off;
end

function corrected_trace = applyBridgeFix(trace)
    % Hardcoded drop index to match original loop
    drop_idx = 200003;

    % Define segments based on your logic
    PART1 = trace(1:drop_idx);
    PART2 = trace(drop_idx+1:end);

    % Bridge drop = difference between 200003 and 200001
    BRIDGE_DROP = trace(drop_idx) - trace(drop_idx - 2);

    % Apply correction exactly like original code
    corrected_trace = [ ...
        PART1(1:end-2); ...
        repmat(PART1(end-2), 5, 1); ...
        PART2(4:19998) + abs(BRIDGE_DROP); ...
        repmat(PART2(19998) + abs(BRIDGE_DROP), 3, 1); ...
        PART2(20002:end) ...
    ];
    
end

function [results, summary] = analyze_experimental_properties(neuron_name, time, flyTraces, flyNums)
    % flyTraces: cell array containing trials-by-time matrices.
    % time: milliseconds; voltage: mV.
    % results: measurements for each injection and each average trace.
    % summary: numeric per-fly means, SDs, and average-trace properties.

    switch neuron_name
        case 'DNp01'
            current = -0.06;   % nA
        case 'DNp03'
            current = -0.002;  % nA
        otherwise
            error('Unknown neuron: %s', neuron_name);
    end

    time = time(:);
    nFlies = numel(flyTraces);

    assert(nFlies > 0 && nFlies == numel(flyNums), ...
        'Provide one fly number per trace matrix.');

    % Same analysis windows as the Python function.
    baselineWindow = time > 0 & time < 250;
    steadyWindow = time > 262.5 & time < 287.5;
    fitWindow = time > 250 & time < 260;
    tFit = time(fitWindow) - 250;

    % For resistance from ONLY the final injection segment, use:
    % steadyWindow = time > 275 & time < 300;

    expModel = fittype('A*exp(-x/tau)+C', ...
        'independent', 'x', ...
        'coefficients', {'A', 'tau', 'C'});

    results = table();
    summaryValues = nan(nFlies, 7);
    flyMeans = zeros(nFlies, numel(time));

    for f = 1:nFlies
        traces = flyTraces{f};

        assert(~isempty(traces) && size(traces, 2) == numel(time), ...
            'Each matrix must contain trials in rows and time in columns.');

        nInjections = size(traces, 1);
        flyMeans(f, :) = mean(traces, 1);

        % Analyze individual injections, followed by the average trace.
        traces = [traces; flyMeans(f, :)];
        labels = [compose("Injection %d", (1:nInjections)'); ...
                  "Average trace"];

        measurements = nan(nInjections + 1, 4);

        for k = 1:size(traces, 1)
            voltage = traces(k, :)';

            baseline = mean(voltage(baselineWindow));
            steadyState = mean(voltage(steadyWindow));
            deltaV = steadyState - baseline;

            Rin = deltaV / current;       % MOhm
            restingVm = baseline - 13;    % LJP corrected, mV
            tau = NaN;
            capacitance = NaN;

            try
                fitted = fit(tFit, voltage(fitWindow), expModel, ...
                    'StartPoint', [max(-deltaV, 0.001), 1, steadyState], ...
                    'Lower', [0, 0.001, -Inf], ...
                    'Upper', [Inf, Inf, Inf]);

                tau = fitted.tau;         % ms
            catch ME
                warning('Fly %d, %s: tau fit failed: %s', ...
                    flyNums(f), labels(k), ME.message);
            end

            if isfinite(tau) && isfinite(Rin) && Rin > 0
                capacitance = 1000 * tau / Rin;  % pF
            end

            measurements(k, :) = [restingVm, Rin, tau, capacitance];
        end

        flyResults = array2table(measurements, ...
            'VariableNames', {'RestVm_LJP_mV', 'Rin_MOhm', ...
                              'Tau_ms', 'Ctotal_pF'});

        flyResults = addvars(flyResults, ...
            repmat(string(flyNums(f)), nInjections + 1, 1), labels, ...
            'Before', 1, 'NewVariableNames', {'Fly', 'Trace'});

        fprintf('\n%s — Fly %d\n', neuron_name, flyNums(f));
        disp(flyResults);
        results = [results; flyResults];

        % Summarize injections ONLY: exclude the appended average trace.
        restingValues = measurements(1:nInjections, 1);
        resistanceValues = measurements(1:nInjections, 2);

        restingValues = restingValues(isfinite(restingValues));
        resistanceValues = resistanceValues(isfinite(resistanceValues));

        meanRest = mean(restingValues);
        meanRin = mean(resistanceValues);

        sdRest = NaN;
        sdRin = NaN;

        if numel(restingValues) > 1
            sdRest = std(restingValues, 0);       % Sample SD
        end
        if numel(resistanceValues) > 1
            sdRin = std(resistanceValues, 0);     % Sample SD
        end

        if numel(restingValues) < nInjections || ...
                numel(resistanceValues) < nInjections
            warning('Fly %d: nonfinite values excluded from mean/SD.', ...
                flyNums(f));
        end

        tauAverageTrace = measurements(end, 3);
        capacitanceAverage = NaN;

        if isfinite(tauAverageTrace) && isfinite(meanRin) && meanRin > 0
            capacitanceAverage = 1000 * tauAverageTrace / meanRin;
        end

        summaryValues(f, :) = [ ...
            nInjections, meanRin, sdRin, meanRest, sdRest, ...
            tauAverageTrace, capacitanceAverage];
    end

    % Numeric summary for subsequent analysis.
    summary = array2table(summaryValues, ...
        'VariableNames', { ...
            'N_injections', ...
            'Rin_mean_MOhm', 'Rin_SD_MOhm', ...
            'RestVm_mean_mV', 'RestVm_SD_mV', ...
            'Tau_average_trace_ms', 'Ctotal_pF'});

    summary = addvars(summary, string(flyNums(:)), ...
        'Before', 1, 'NewVariableNames', 'Fly');

    % Across-fly averages: do not average the within-fly SDs.
    averageValues = nan(1, 7);
    meanColumns = [2, 4, 6, 7];
    averageValues(meanColumns) = ...
        mean(summaryValues(:, meanColumns), 1, 'omitnan');

    averageRow = array2table(averageValues, ...
        'VariableNames', summary.Properties.VariableNames(2:end));
    averageRow = addvars(averageRow, "Average", ...
        'Before', 1, 'NewVariableNames', 'Fly');

    summary = [summary; averageRow];

    % Formatted console table.
    injectionText = [string(summaryValues(:, 1)); "N/A"];

    rinText = [compose("%.2f ± %.2f", ...
        summaryValues(:, 2), summaryValues(:, 3)); ...
        compose("%.3f", averageValues(2))];

    restText = [compose("%.2f ± %.2f", ...
        summaryValues(:, 4), summaryValues(:, 5)); ...
        compose("%.3f", averageValues(4))];

    displayTable = table( ...
        summary.Fly, injectionText, rinText, restText, ...
        summary.Tau_average_trace_ms, summary.Ctotal_pF, ...
        'VariableNames', { ...
            'Fly', 'Current_injections', ...
            'Rin_MOhm_mean_SD', 'RestVm_LJP_mV_mean_SD', ...
            'Tau_ms', 'Capacitance_pF'});

    fprintf('\n%s — PER-FLY SUMMARY (mean ± SD across injections)\n', ...
        neuron_name);
    disp(displayTable);
end