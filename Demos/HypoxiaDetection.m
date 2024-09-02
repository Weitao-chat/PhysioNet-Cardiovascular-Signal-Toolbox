%	OVERVIEW:
%       This code analyzes a segment of data collected in the 
%       fetal sheep which contains ECG signals 
%
%   OUTPUT:
%       ECG Metrics exported to .csv files
%
%   DEPENDENCIES & LIBRARIES:
%       https://github.com/Weitao-chat/Advancing-Fetal-Surveillance-with-Physiological-Sensing-Detecting-Hypoxia-in-Fetal-Sheep
%   APA REFERENCE: 
%       Weitao Tang. (2024). Weitao-chat/Advancing-Fetal-Surveillance-with-Physiological-Sensing-Detecting-Hypoxia-in-Fetal-Sheep: V1.0.1 (V1.0.1). 
%       Zenodo. https://doi.org/10.5281/zenodo.12540027
%	REPO:       
%       https://github.com/Weitao-chat/Advancing-Fetal-Surveillance-with-Physiological-Sensing-Detecting-Hypoxia-in-Fetal-Sheep
%   ORIGINAL SOURCE AND AUTHORS:     
%       Weitao Tang  
%	COPYRIGHT (C) 2024
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

run(['..' filesep 'startup.m'])

% Remove old files generated by this code
OldFolder = [pwd,filesep, 'OutputData', filesep, 'ResultsFetalSheep'];
if exist(OldFolder, 'dir')
    rmdir(OldFolder, 's');
    fprintf('Old Demo Folder deleted \n');
end


HRVparams = InitializeHRVparams('Hypoxia_Detection'); % include the project name

HRVparams.PeakDetect.THRES = 0.3;
HRVparams.PeakDetect.REF_PERIOD = 0.15;
HRVparams.PeakDetect.windows = 5;

HRVparams.poincare.on = 0; % Poincare analysis off for this code
HRVparams.MSE.on = 0; % MSE analysis off for this demo
HRVparams.HRT.on = 0; % HRT analysis off for this code


% 1. Load Raw Fetal Sheep Data (ECG Waveform)
filePath = fullfile(pwd, '..', 'Demos', 'TestData', 'Weitao_issue_data', '20112_Baseline.mat');
ECG_data = load(filePath);

% Get the field names (assuming there's only one field)
fieldNames = fieldnames(ECG_data);

% Access the data dynamically
ecg_segment = ECG_data.(fieldNames{1});

jqrs_ann = run_qrsdet_by_seg(ecg_segment, HRVparams);
GainQrsDetect = 200000000; % Default value for gain (adu/physical unit)

sqrs_ann = run_sqrs(ecg_segment * GainQrsDetect, HRVparams, 0);
wqrs_ann = wqrsm_fast(ecg_segment * GainQrsDetect, HRVparams.Fs);

disp(length(jqrs_ann)); % check the detection length of jqrs
disp(length(sqrs_ann)); % check the detection length of sqrs 
disp(length(wqrs_ann)); % check the detection length of wqrs 

plot(ecg_segment);
title('ECG Signal');
xlabel('Samples');
ylabel('Amplitude');

[SQIjw, StartSQIwindows_jw] = bsqi(jqrs_ann(:), wqrs_ann(:), HRVparams);


ConvertRawDataToRRIntervals(ecg_segment,HRVparams,'3')

% Set the file ID and segment to be analyzed
FileID = '20112'; % Example
Segment = 'Baseline'; % Example

% Define and create a subfolder for saving the results
outputFolder = fullfile(pwd, '..', 'Demos', 'OutputData', 'ResultsHD', [FileID '_' Segment]);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Update the writedata path in HRVparams
HRVparams.writedata = outputFolder;

% 2. Analyze data using HRV PhysioNet Cardiovascular Signal Toolbox
[results, resFilename] = Main_HRV_Analysis(ecg_segment, [], 'ECGWaveform', HRVparams, [FileID '_' Segment]);

fprintf('Results for %s_%s saved in %s\n', FileID, Segment, outputFolder);

