%% Test REST (go to 'REST' folder)
%% set path and initialize bcilab
bcilab_path = which('bcilab.m');
if isempty(bcilab_path)
    current_path = pwd;
    addpath(['dependencies' filesep 'BCILAB']); bcilab
    cd(current_path);
    addpath(genpath('./'));
end

%% refresh workspace
% clear all will break bcilab and require it to restart as it uses global
% variables 
close all
delete(timerfind)
clear

%% define opts structure
% whether to customize pipeline 
opts.customize_pipeline = true;

% (optional) define config file name
opts.config = 'Config_ORICA_SleepHeadband';

% (optional) channel location file
load(['data' filesep 'chanlocs' filesep 'SleepHeadband_8_Stream.mat']); 
opts.chanlocs = chanlocs;

% point to headModel
% opts.headModel = ['data' filesep 'head_models' filesep 'quick20HeadModel'];

% (optional) path to calibration data and select time window
opts.calibration_data = ['data' filesep 'sleepBand_sample.set'];
opts.calibration_window = [0,10]; % sec

% use playback data
opts.playback = 1;

%% start REST
REST(opts)
