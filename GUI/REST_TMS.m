function varargout = REST_TMS(varargin)
%REST_TMS MATLAB code file for REST_TMS.fig
%      REST_TMS, by itself, creates a new REST_TMS or raises the existing
%      singleton*.
%
%      H = REST_TMS returns the handle to a new REST_TMS or the handle to
%      the existing singleton*.
%
%      REST_TMS('Property','Value',...) creates a new REST_TMS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to REST_TMS_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      REST_TMS('CALLBACK') and REST_TMS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in REST_TMS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help REST_TMS

% Last Modified by GUIDE v2.5 16-Jun-2017 17:57:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @REST_TMS_OpeningFcn, ...
    'gui_OutputFcn',  @REST_TMS_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before REST_TMS is made visible.
function REST_TMS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Set parameters
handles.nplots = 8;         % number of sources or channels to plot
handles.isSourceMode = 0;   % flag for display mode: channel 0, source 1
handles.isSwitched = 0;
handles.curIC = 1;          % IC index for source localization
handles.lock = [];
handles.color_lock = [0.5 1 0.5];
handles.reject = [];
handles.color_reject = [1 0.5 0.5];

% Parse varagin
[handles,calibData,customize_pipeline] = startup_check_inputs(handles,varargin);

% Intialize ORICA
handles = startup_initializeORICA(handles,calibData,customize_pipeline);

% Create ORICA timer
oricaTimer = timer('Period',.1,'ExecutionMode','fixedSpacing','TimerFcn',{@onl_filtered_ORICA,handles.streamName},'StartDelay',0.1,'Tag','oricaTimer','Name','oricaTimer');%,'BusyMode','queue');

% Start EEG stream
% [~, handles.bufferName] = vis_stream_ORICA('figurehandles',handles.figure1,'axishandles',handles.axisEEG,'streamname',handles.streamName);
[~, handles.bufferName] = vis_stream_ORICA('figurehandles',handles.figure1,'streamname',handles.streamName);
% eegTimer = timerfind('Name','eegTimer');
% set(eegTimer,'UserData',{hObject,1})

% Gather and disperse pipeline function names
funs = get_pipeline_functions();
funnames = cell(length(funs)-2,1);
for it = 2:length(funs)
    temp = arg_report('properties',funs{it});
    funnames{it-1} = temp.name;
    if iscell(funnames{it-1})
        funnames{it-1} = funnames{it-1}{1}; end
end
% set(handles.popupmenuEEG,'String',['Raw Data'; funnames; 'ICA Cleaned'])
% buffer = evalin('base',handles.bufferName);
% buffer.funs = funs;
% assignin('base',handles.bufferName,buffer);

% Find if channels have been removed
funsstr = cellfun(@func2str,funs,'uniformoutput',false');
if any(strcmp(funsstr,'flt_selchans'))
    pipeline = evalin('base','pipeline'); %#ok<NASGU>
    numparts = find(flipud(strcmp(funsstr,'flt_selchans')))-1;
    eval(['removed = pipeline' repmat('.parts{2}',1,numparts) '.parts{4};']);
    handles.rmchan_index = ismember({handles.chanlocs.labels},removed);
    % adjust chanlocs
    handles.urchanlocs = handles.chanlocs;
    handles.chanlocs(handles.rmchan_index) = [];
    handles.nic = length(handles.chanlocs);
    handles.ics = 1:handles.nic;
    % adjust headModel
    if isfield(handles,'headModel')
        handles.urheadModel = handles.headModel;
        handles.headModel.dropChannels(handles.rmchan_index); % !!! had to change the headModel contructor
        handles.K(handles.rmchan_index,:) = [];
        %     LFM = load(handles.headModel.leadFieldFile);
        %     LFM.K(handles.rmchan_index,:) = [];
        %     save(handles.headModel.leadFieldFile,'-struct','LFM')
    end
end

if handles.isSourceMode
    % Populate scalp maps
    for it = 1:handles.nplots
        set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)])); clf;
        [~,Zi,~,Xi,Yi,intx,inty] = topoplotFast_LRBF(rand(size(handles.chanlocs)), handles.chanlocs);
        set(handles.(['panelIC' int2str(it)]),'Title',['IC' int2str(handles.ics(it))])
    end
else
    % Populate channel maps and scalp maps
    for it = 1:handles.nplots
        set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]))
        topoplotCH([],handles.chanlocs(it), 'style', 'blank', 'electrodes','on', ...
            'electcolor',[1 0 0],'emarkersize',20);
        set(handles.(['panelIC' int2str(it)]),'Title',handles.chanlocs(it).labels);
    end
end

% Create scalp map timer
topoTimer = timer('Period',round(1/handles.nplots*1000)/1000,'ExecutionMode','fixedRate','TimerFcn',{@vis_topo,hObject},'StartDelay',0.2,'Tag','topoTimer','Name','topoTimer');

% Create data timer (starts as power spectrum)
% infoTimer = timer('Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');
% infoTimer = timer('Period',round(1/handles.nplots*1000)/1000,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');
infoTimer = timer('Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');

% Set panel and button colors
handles.color_bg = get(handles.figure1,'Color');
names = fieldnames(handles);
ind = find(any([strncmpi(names,'panel',5),strncmpi(names,'toggle',6),strncmpi(names,'push',4),strncmpi(names,'popup',5)],2));
for it = 1:length(ind)
    set(handles.(names{ind(it)}),'BackgroundColor',handles.color_bg)
end

% Save timers
% handles.pauseTimers = [eegTimer,topoTimer,infoTimer];
handles.pauseTimers = [topoTimer,infoTimer];


% Choose default command line output for REST_TMS
handles.output = hObject;

% Update handles structure
handles.intialized = true;
guidata(hObject, handles);

% Start timers
start(oricaTimer);
% start(eegTimer);
start(topoTimer);
start(infoTimer);

% UIWAIT makes REST_TMS wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function [funs] = get_pipeline_functions(p)
if ~exist('p','var')
    p = evalin('base','pipeline'); end
funs = {};
if p.subnodes
    % recursive call for deeper parts of pipeline
    for k=p.subnodes
        [funs] = get_pipeline_functions(p.parts{k}); end
end
% build the outputs
funs = [funs;{p.head}];



function [handles, calibData, customize_pipeline] = startup_check_inputs(handles,in)
% !!! need to add appropriate errors and explanations

% set display mode
if handles.isSourceMode
    set(handles.uibuttonMode,'selectedobject',handles.radiobuttonIC);
else
    set(handles.uibuttonMode,'selectedobject',handles.radiobuttonCH);
end

% check channel locations
if isfield(in{1},'chanlocs')
    handles.chanlocs = in{1}.chanlocs; end

% check calibration data
if isfield(in{1},'calibration_data')
    if isstruct(in{1}.calibration_data)
        calibData = in{1}.calibration_data;
    elseif ischar(in{1}.calibration_data)
        calibData = pop_loadset(in{1}.calibration_data);
    end
    if ~isfield(handles,'chanlocs') && isfield(calibData,'chanlocs') && ~isempty(calibData.chanlocs)
        handles.chanlocs = calibData.chanlocs; end
else
    calibData = [];
end

% check if playback is requested
if isfield(in{1},'playback') && in{1}.playback
    playbackStream = play_eegset_lsl(calibData,'REST_playback_data','REST_playback_markers',[],true); end

% shorten calibration data if requested
if isfield(in{1},'calibration_window')
    if isscalar(in{1}.calibration_window)
        calibData = pop_select(calibData,'time',[0 in{1}.calibration_window]);
    else
        calibData = pop_select(calibData,'time',in{1}.calibration_window);
    end
end

% check if localization is possible and adjust GUI accordingly
handles = startup_check_localization(handles,in{1});

if ~isfield(handles,'chanlocs') && isfield(handles,'headModel')
    for it = 1:size(handles.headModel.channelSpace,1)
        handles.chanlocs(it).labels = handles.headModel.channelLabel{it};
        handles.chanlocs(it).X = handles.headModel.channelSpace(it,1);
        handles.chanlocs(it).Y = handles.headModel.channelSpace(it,2);
        handles.chanlocs(it).Z = handles.headModel.channelSpace(it,3);
    end
    handles.chanlocs = convertlocs(handles.chanlocs);
end

% if still no chanlocs, error
if ~isfield(handles,'chanlocs')
    error('REST: No channel location information provided!')
end

% check whether to open pipeline arg_guipanel
if isfield(in{1},'customize_pipeline')
    customize_pipeline = in{1}.customize_pipeline;
else
    customize_pipeline = false;
end

% check if config file is defined
if isfield(in{1},'config')
    handles.config = in{1}.config;
else
    handles.config = 'Config_ORICA';
end



function handles = startup_check_localization(handles,in) % !!! combine headmodel in localization
% if no headModel provided, remove localization button
if ~isfield(in,'headModel') || isempty(in.headModel)
    set(handles.pushbuttonLocalize,'string','Localize (off)')
    
    % if the provided headModel is a string, load the headModel
elseif isa(in.headModel,'char')
    handles.headModel = headModel.loadFromFile(in.headModel);
    temp = load(handles.headModel.surfaces);
    handles.nVertices = size(temp.surfData(3).vertices,1);
    
    % if there is not an accompanying *_SSPEB.mat file containing the
    % cropped lead-field matrix, laplacian matrix, and valid vertex indices,
    % then calulate them. (takes a while).
    if ~exist([in.headModel '_SSPEB.mat'],'file')
        [~,K,L,rmIndices] = getSourceSpace4PEB(handles.headModel);
        hmInd = setdiff(1:handles.nVertices,rmIndices);
        save([in.headModel '_SSPEB.mat'],'K','L','hmInd')
        
        handles.K = K;
        handles.L = L;
        handles.hmInd = hmInd;
    else
        temp = load([in.headModel '_SSPEB.mat']);
        handles.K = temp.K;
        handles.L = temp.L;
        handles.hmInd = temp.hmInd;
    end
    
    if handles.isSourceMode
        set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']); end

    % if the provided headModel is an object, use it
elseif isa(in.headModel,'headModel')
    handles.headModel = in.headModel;
    
    % find the number of vertices in the model
    if ischar(handles.headModel.surfaces)
        temp = load(handles.headModel.surfaces);
        handles.nVertices = size(temp.surfData(3).vertices,1);
    else
        handles.nVertices = size(handles.headModel.surfaces(3).vertices,1);
    end
    
    % if K, L, and hmInd are not provided in opts, then calulate them.
    if ~isfield(in,'K') || ~isfield(in,'L') || ~isfield(in,'hmInd')
        [~,handles.K,handles.L,rmIndices] = ...
            getSourceSpace4PEB(handles.headModel);
        handles.hmInd = setdiff(1:handles.nVertices,rmIndices);
    else
        handles.K = in.K;
        handles.L = in.L;
        handles.hmInd = in.hmInd;
    end
    if handles.isSourceMode
        set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']); end

end



function handles = startup_initializeORICA(handles,calibData,customize_pipeline)

% load LSL
if ~exist('env_translatepath','file')
    % standalone case
    lib = lsl_loadlib();
else
    % if we're within BCILAB we want to make sure that the library is also found if the toolbox is compiled
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
end

% find streams
streams = lsl_resolve_all(lib,3);
streamnames = cellfun(@(s)s.name(),streams ,'UniformOutput',false);
if isempty(streamnames)
    error('There is no stream visible on the network.');
elseif length(streamnames) == 1
    assignin('base','streamname',streamnames); % save stream name in base workspace
else
    % if more than 2 (EEG) streams, pop up a GUI to select one.
    selStream(streamnames); % result 'streamname' is saved in base workspace
    streamnames = evalin('base','streamname');
    streamnames = streamnames{1};
end
run_readlsl_ORICA('MatlabStream',streamnames,'DataStreamQuery', ['name=''' streamnames '''']);
if ~isvarname(streamnames), streamnames = streamnames(~ismember(streamnames,['-' ' '])); end
handles.streamName = streamnames;
opts.lsl.StreamName = streamnames;

% create learning rate buffer
bufflen = 60; % seconds
handles.srate = getfield(onl_peek(opts.lsl.StreamName,1,'samples'),'srate');
handles.nbchan = getfield(onl_peek(opts.lsl.StreamName,1,'samples'),'nbchan');
assignin('base','learning_rate',nan(1,bufflen*handles.srate));

% set stream info in GUI
set(handles.textStreamName,'string',streamnames);
set(handles.textNumAll,'string',num2str(handles.nbchan));
set(handles.textSrate,'string',[num2str(handles.srate) ' Hz']);

% look for pre-existing config file for pipeline
REST_path = fileparts(fileparts(which('REST')));
opts.BCILAB_PipelineConfigFile = ...
    [REST_path filesep 'data' filesep 'config' filesep handles.config '.mat']; % make sure this file doesn't have 'signal' entry

% define the pipeline configuration
tic
try
    fltPipCfg = exp_eval(io_load(opts.BCILAB_PipelineConfigFile));
    if customize_pipeline
        fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
            'Parameters',[{'signal',onl_peek(opts.lsl.StreamName,1,'samples')} fltPipCfg], ...
            'PanelOnly',false);
    end
catch
    disp('-- no existing pipeline or fail loading pipeline--');
    fltPipCfg = {};
end

% open pipeline configuration gui if no settings found or if user requested
if isempty(fltPipCfg)
    fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
        'Parameters',[{'signal',onl_peek(opts.lsl.StreamName,1,'samples')} fltPipCfg], ...
        'PanelOnly',false);
end

if isfield(fltPipCfg,'pselchans')
    if isfield(calibData.etc,'badChLabels')
        fltPipCfg.pselchans.channels = calibData.etc.badChLabels;
    end
end

% save the configuration %!!! maybe disable this?
if ~isempty(fltPipCfg)
    if isfield(fltPipCfg,'signal')
        fltPipCfg = rmfield(fltPipCfg,'signal'); end
    save(env_translatepath(opts.BCILAB_PipelineConfigFile), ...
        '-struct','fltPipCfg');
end

% grab calib data from online stream if there is none
if isempty(calibData)
    disp('Collecting calibration data from online stream... please wait 10 seconds...');
    pause(10-toc); % uh oh!
    calibData = onl_peek(opts.lsl.StreamName,10,'seconds');
end

% check for bad channels
calibData = warmStartWithBadChRemoved(calibData);

% set number of EEG channels in GUI
if isfield(fltPipCfg,'pselchans') && ~isempty(fltPipCfg.pselchans.channels)
    set(handles.textNumEEG,'string',num2str(calibData.nbchan - length(fltPipCfg.pselchans.channels)));
else
    set(handles.textNumEEG,'string',num2str(calibData.nbchan));
end

% run pipline on calibration data
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);
assignin('base','pipeline',pipeline);



function vis_topo(varargin)

% get handles
handles = guidata(varargin{3});

% update index
it = mod(get(varargin{1},'TasksExecuted')-1,handles.nplots)+1;
hstr = ['axesIC' int2str(it)];
hand = get(handles.(hstr),'children');

if handles.isSourceMode     % display IC scalp maps
    % get the updated stream buffer
    W = evalin('base','pipeline.state.icaweights');
    sphere = evalin('base','pipeline.state.icasphere');
    
    try
        Winv = inv(W*sphere);
        
        %     [map, cmin, cmax] = topoplotUpdate(Winv(:,handles.ics(it)), handles.chanlocs,'electrodes','off','gridscale',32);
        map = zeros(handles.topoNPixel,1);
        map(~handles.topoNaNMask) = handles.topoMat*Winv(:,handles.ics(it));
        maxabs = max(abs(map));
        cmin = -maxabs;
        cmax =  maxabs;
        map(handles.topoNaNMask) = NaN;
        map = reshape(map,sqrt(handles.topoNPixel),[]);
        
        surfInd = strcmp(get(hand,'type'),'surface');
        set(hand(surfInd),'CData',map);
        set(handles.(hstr),'CLim',[cmin cmax]);
    end
    
    % update name and buttons
    lock = any(handles.lock==handles.ics(it));
    reject = any(handles.reject==handles.ics(it));
    set(handles.(['panelIC' int2str(it)]),'Title',['IC' int2str(handles.ics(it))])
    set(handles.(['togglebuttonLock' int2str(it)]),'Value',lock,...
        'BackgroundColor',handles.color_lock*lock + handles.color_bg*(1-lock))
    set(handles.(['togglebuttonReject' int2str(it)]),'Value',reject,...
        'BackgroundColor',handles.color_reject*reject + handles.color_bg*(1-reject))
end


function infoPSD(varargin)

% plot PSD of selected IC
try  %#ok<*TRYNC>
    handles = guidata(varargin{3});
    srate = evalin('base',[handles.bufferName '.srate']);
    data = evalin('base',[handles.bufferName '.data{end}']); % FIXME: make this more robust
    
    % parameters for computing PSD
    dataLen = 5 * srate;        % seconds
    windowLen = 2.5 * srate;      % seconds
    noverlap = 2 * srate;    % seconds
    nfft = 2^nextpow2(srate);
    
    % select data
    if all(data(:,end)==0)
        mark=length(data);
        while true
            mark = find(data(1,1:mark)~=0,1,'last');
            if all(data(:,mark)~=0)
                break; end
        end
        data = data(:,max(1,mark-dataLen+1):mark);
    else
        data = data(:,max(1,end-dataLen+1):end);
    end
    
    data = bsxfun(@minus,data,mean(data,2));
    
    if handles.isSourceMode
        W = evalin('base','pipeline.state.icaweights');
        sphere = evalin('base','pipeline.state.icasphere');
        data = W*sphere*data;
    end
    
    % FIXME: use online PSD estimation
    % compute fft
    [pxx,f] = pwelch(data',round(windowLen),round(noverlap),nfft,srate);
    maxFreq = 40;
    ind = f <= maxFreq;
    
    % update all plots
    for it = 1:handles.nplots
        hstr = ['axes' int2str(it)];
        plot(handles.(hstr),f(ind),db(pxx(ind,it))) % !!!?
        grid(handles.(hstr),'on');
        axis(handles.(hstr),'tight')
        set(handles.(hstr),'XTick',[0 10:10:f(end)],'box','off')
        %     xlabel(handles.(hstr),'Frequency (Hz)')
        %     ylabel(handles.(hstr),'Power/Frequency (dB/Hz)')
    end
    
end



function genICSelectGUI(hObject,handles)
if handles.isSourceMode
    temp = get(handles.figure1,'Position');
    fhandle = figure('toolbar','none','Menubar','none','Name','Select IC','position',[1 1 temp(3:4)],'Resize','on','Colormap',colormap('jet'),'DeleteFcn',{@closeFigIC,hObject});
    W = evalin('base','pipeline.state.icaweights');
    sphere = evalin('base','pipeline.state.icasphere');
    Winv = inv(W*sphere);
    
    rowcols(2) = ceil(sqrt(handles.nic));
    rowcols(1) = ceil(handles.nic/rowcols(2));
    scaleMatTopo = [1 0 0 0;0 1 0 0;0 0 1 0;0 .2 0 .8];
    buttonGap = .1;
    scaleMatReject = [1 0 0 0;0 1 0 0;0 0 .5-buttonGap/2 0;.5+buttonGap/2 0 0 .2];
    for it = 1:handles.nic
        h(it) = subaxis(rowcols(1),rowcols(2),it,'MR',.025,'ML',.025,'MT',.025,'MB',.025,'SH',0,'SV',0.02);
        tempPos = get(h(it),'Position');
        set(h(it),'position',get(h(it),'position')*scaleMatTopo)
        topoplotFast_LRBF(Winv(:,it),handles.chanlocs);
        title(['IC' int2str(it)])
        
        lock = any(handles.lock==it);
        reject = any(handles.reject==it);
        buttonLock(it) = uicontrol('Style', 'togglebutton', 'String', 'Lock',...
            'Units','normalize','Position', tempPos.*[1 1 .5-buttonGap/2 .2],...
            'Callback', {@lockIC,it,hObject},'Value',lock,...
            'BackgroundColor',handles.color_lock*lock + handles.color_bg*(1-lock));
        buttonReject(it) = uicontrol('Style', 'togglebutton', 'String', 'Reject',...
            'Units','normalize','Position', tempPos*scaleMatReject,...
            'Callback', {@rejectIC,it,hObject},'Value',reject,...
            'BackgroundColor',handles.color_reject*reject + handles.color_bg*(1-reject));
    end
    handles.figIC.buttonLock = buttonLock;
    handles.figIC.buttonReject = buttonReject;
    handles.figIC.handle = fhandle;
end
guidata(hObject,handles);


function closeFigIC(varargin)
hObject = varargin{3};
% load handles
handles = guidata(hObject);
if isfield(handles,'figIC')
    handles = rmfield(handles,'figIC'); end
guidata(hObject,handles);



function lockIC(varargin)
ic = varargin{3};
button = varargin{1};
% load handles
if numel(varargin)>3
    hObject = varargin{4};
else
    hObject = get(button,'parent');
end
handles = guidata(hObject);
if get(button,'Value') % turned lock on
    handles.lock = sort([handles.lock ic]);
    set(button,'BackgroundColor',[0.5 1 0.5])
    if isfield(handles,'figIC')
        set(handles.figIC.buttonLock(ic),'Value',1,'BackgroundColor',handles.color_lock); end
else % turned lock off
    handles.lock(handles.lock==ic) = [];
    set(button,'BackgroundColor',handles.color_bg)
    if isfield(handles,'figIC')
        set(handles.figIC.buttonLock(ic),'value',0,'BackgroundColor',handles.color_bg); end
end
% save handles
guidata(hObject,handles);
% update ics to plot
updateICs(hObject,false)



function rejectIC(varargin)
ic = varargin{3};
button = varargin{1};
% load handles
if numel(varargin)>3
    hObject = varargin{4};
else
    hObject = get(button,'parent');
end
handles = guidata(hObject);
if get(button,'Value') % turned reject on
    handles.reject = sort([handles.reject ic]);
    set(button,'BackgroundColor',[1 0.5 0.5])
    if isfield(handles,'figIC')
        set(handles.figIC.buttonReject(ic),'Value',1,'BackgroundColor',handles.color_reject); end
else % turned reject off
    handles.reject(handles.reject==ic) = [];
    set(button,'BackgroundColor',handles.color_bg)
    if isfield(handles,'figIC')
        set(handles.figIC.buttonReject(ic),'value',0,'BackgroundColor',handles.color_bg); end
end
% save handles
guidata(hObject,handles);



function updateICs(hObject,flag_sort)
% load handles
handles = guidata(hObject);
if flag_sort
    % load info
    S = evalin('base','pipeline.state.icasphere');
    W = evalin('base','pipeline.state.icaweights');
    V = evalin('base','pipeline.state.Var');
    other = setdiff(1:handles.nic,[handles.lock]);
    % sort by residual variance
    meanvar = mean(pinv(W*S).^2).*V';
    [~,ind_lock] = sort(meanvar(handles.lock),'descend');
    [~,ind_other] = sort(meanvar(other),'descend');
    handles.ics = [handles.lock(ind_lock) other(ind_other)];
else
    handles.ics = [handles.lock setdiff(1:handles.nic,[handles.lock])];
end
guidata(hObject,handles);



% --- Outputs from this function are returned to the command line.
function varargout = REST_TMS_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'intialized') && handles.intialized
    if isfield(handles,'figIC')
        try
            close(handles.figIC.handle); end, end
    if isfield(handles,'figLoc')
        try
            close(handles.figIC.handle); end, end
    
    timerNames = {'eegTimer','oricaTimer','topoTimer','infoTimer','locTimer',[handles.streamName '_timer']};
    % warning off MATLAB:timer:deleterunning
    for it = 1:length(timerNames)
        temp = timerfind('name',timerNames{it});
        if isempty(temp)
            continue; end
        stop(temp)
        delete(temp)
    end
end




% --- Executes on button press in pushbuttonSelect.
function pushbuttonSelect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'figIC')
    figure(handles.figIC.handle)
    return
end

if strcmpi(get(handles.pushbuttonPause,'string'),'Resume')
    genICSelectGUI(hObject,handles)
else
    stop(handles.pauseTimers)
    genICSelectGUI(hObject,handles)
    start(handles.pauseTimers)
end


% --- Executes on button press in pushbuttonPause.
function pushbuttonPause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(get(handles.pushbuttonPause,'string'),'Pause')
    set(handles.pushbuttonPause,'string','Resume');
    stop(handles.pauseTimers)
else
    set(handles.pushbuttonPause,'string','Pause');
    start(handles.pauseTimers);
end


% --- Executes on button press in pushbuttonLocalize.
function pushbuttonLocalize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLocalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.isSourceMode
    % if no existing headModel, return
    if ~isfield(handles,'headModel') || isempty(handles.headModel)
        disp('No existing headModel')
        return
    end
    
    % if it already exists then just focus on it
    if isfield(handles,'figLoc')
        figure(handles.figLoc.handle.hFigure)
        return
    end
    
    % otherwise pause if necessary and then create the window
    method =  questdlg('Choose method for source Localization','Localization', ...
        'LORETA','Dipole Fit','LORETA');
    
    flag_resume = false;
    if strcmpi(get(handles.pushbuttonPause,'string'),'Pause')
        stop(handles.pauseTimers);
        flag_resume = true;
    end
    
    switch method
        case 'LORETA'
            figLoc_gen_LORETA(hObject,handles)
        case 'Dipole Fit'
            figLoc_gen_dipolefit(hObject,handles)
    end
    
    if flag_resume %#ok<*ALIGN>
        handles = guidata(hObject);
        start(handles.pauseTimers); 
    end
    
else
    disp('Source localization is only supported in the Source Mode.')
end


function figLoc_gen_LORETA(hObject,handles)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% run dynamicLoreta once to generate state and initial localization
[U,S,V] = svd(handles.K/handles.L,'econ');
Ut = U';
s2 = diag(S).^2;
iLV = handles.L\V;
options = struct('maxTol',1e-3,'maxIter',100,'gridSize',100,'verbose',false,'history',true,'useGPU',false,'initNoiseFactor',0.001);
[J,sigma2,tau2] = dynamicLoreta(Winv(:,handles.curIC), Ut, s2,iLV,[],[], options);

% create headModel plot
if length(handles.hmInd) == length(J)
    Jest = zeros(handles.nVertices,1);
    Jest(handles.hmInd,:) = J;
else
    Jest = zeros(handles.nVertices,3);
    Jest(handles.hmInd,:) = reshape(J,[],3);
end
fhandle = handles.headModel.plotOnModel(Jest(:),Winv(:,handles.curIC),sprintf('IC %d Localization (LORETA)',handles.curIC));
set(fhandle.hFigure,'DeleteFcn',{@closeFigLoc,hObject},'name',['IC' num2str(handles.curIC)]);
colorbar('hide');

% create timer
locTimer = timer('Period',3,'StartDelay',3,'ExecutionMode','fixedRate','TimerFcn',{@figLoc_update_LORETA,hObject},'Tag','locTimer','Name','locTimer');

% save headModel plot, timer, and dynamicLoreta parameters to handles
handles.pauseTimers = [handles.pauseTimers,locTimer];
handles.figLoc.handle = fhandle;
handles.figLoc.IC = handles.curIC;
handles.figLoc.Ut = Ut;
handles.figLoc.s2 = s2;
handles.figLoc.iLV = iLV;
handles.figLoc.tau2 = tau2;
handles.figLoc.sigma2 = sigma2;
handles.figLoc.options = options;
temp = load(handles.headModel.surfaces);
handles.figLoc.scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,temp.surfData(1).vertices,32);
if range(handles.figLoc.scalp)<.5
    handles.figLoc.scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,temp.surfData(1).vertices,32/1000);
end
guidata(hObject,handles);



function figLoc_update_LORETA(varargin)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% parse inputs
handles = guidata(varargin{3});

% run dynamicLoreta
Ut = handles.figLoc.Ut;
s2 = handles.figLoc.s2;
iLV = handles.figLoc.iLV;
tau2 = handles.figLoc.tau2;
sigma2 = handles.figLoc.sigma2;
options = handles.figLoc.options;
[J,sigma2,tau2] = dynamicLoreta(Winv(:,handles.figLoc.IC), Ut, s2,iLV,sigma2,tau2, options);

% update figure and related object values
if length(handles.hmInd) == length(J)
    Jest = zeros(handles.nVertices,1);
    Jest(handles.hmInd) = J;
    handles.figLoc.handle.sourceMagnitud = Jest;
    set(handles.figLoc.handle.hCortex,'FaceVertexCData',handles.figLoc.handle.sourceMagnitud)
else
    Jest = zeros(handles.nVertices,3);
    Jest(handles.hmInd,:) = reshape(J,[],3);
    handles.figLoc.handle.sourceOrientation = Jest;
    handles.figLoc.handle.sourceMagnitud = squeeze(sqrt(sum(Jest.^2,2)));
    set(handles.figLoc.handle.hVector,'udata',Jest(:,1),'vdata',Jest(:,2),'wdata',Jest(:,3))
    set(handles.figLoc.handle.hCortex,'FaceVertexCData',handles.figLoc.handle.sourceMagnitud)
end
scalp_val = handles.figLoc.scalp*Winv(:,handles.figLoc.IC);
set(handles.figLoc.handle.hScalp,'FaceVertexCData',scalp_val)

% upade color limits
maxabs = max(abs(scalp_val));
handles.figLoc.handle.clim.scalp = [-maxabs maxabs];
maxabs = max(abs(handles.figLoc.handle.sourceMagnitud));
handles.figLoc.handle.clim.source = [-maxabs maxabs];


% save dynamicLoreta parameters to handles
handles.figLoc.tau2 = tau2;
handles.figLoc.sigma2 = sigma2;
guidata(varargin{3},handles);



function figLoc_gen_dipolefit(hObject,handles)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% load surfaces
sd = load(handles.headModel.surfaces);
vertices = sd.surfData(3).vertices(handles.hmInd,:);

% calculate scalp potential transform and guess units
scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,sd.surfData(1).vertices,32);
if range(scalp)<.5
    scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,sd.surfData(1).vertices,32/1000);
end

% create figure and dipole plot
fhandle = handles.headModel.plotDipoles([0 0 0],[0 0 0]);
set(fhandle.hFigure,'DeleteFcn',{@closeFigLoc,hObject},'name',['IC' num2str(handles.curIC)]);
hold(fhandle.hAxes,'on');

% dipole fit
[dipoles, ~, moments, rv, state] = dipole_fit(Winv(:,handles.curIC),handles.K,vertices);

% set scalp potentials
set(fhandle.hScalp,'FaceVertexCData',scalp*Winv(:,handles.curIC),'facecolor','interp')

% plot dipole
arrowlen = max(range(sd.surfData(3).vertices))/10;
if length(moments)/3==length(vertices)
    dipoles.moment = dipoles.moment/norm(dipoles.moment)*arrowlen;
    handles.figLoc.fixed_dip = false;
    arrows_dip = quiver3(dipoles.location(1),dipoles.location(2),dipoles.location(3), ...
        dipoles.moment(1),dipoles.moment(2),dipoles.moment(3),'ko','filled');
else
    handles.figLoc.fixed_dip = true;
    normals = geometricTools.getSurfaceNormals(sd.surfData(end).vertices,sd.surfData(end).faces,false);
    arrows_dip = quiver3(dipoles.location(1),dipoles.location(2),dipoles.location(3), ...
        normals(dipoles.L,1)*arrowlen,normals(dipoles.L,2)*arrowlen,normals(dipoles.L,3)*arrowlen,'ko','filled');
end
colormap(bipolar(512, 0.99))
maxabs = max(abs(vec(scalp*Winv(:,handles.curIC))));
caxis(fhandle.hAxes,[-maxabs maxabs]);

% create residual variance text
handles.figLoc.axisRV = axes('parent',fhandle.hFigure,'position',[.05 .9 .3 .1],'hittest','off');
axis(handles.figLoc.axisRV,'off')
handles.figLoc.textRV = text(0,0.5,sprintf('Residual Variance: %04.1f%%',rv*100), ...
    'parent',handles.figLoc.axisRV,'fontweight','bold','fontsize',16,'hittest','off');

% create timer
locTimer = timer('Period',3,'StartDelay',3,'ExecutionMode','fixedRate', ...
    'TimerFcn',{@figLoc_update_dipolefit,hObject},'Tag','locTimer','Name','locTimer');

% save headModel plot, timer, and bfpf state to handles
handles.figLoc.state = state;
handles.pauseTimers = [handles.pauseTimers,locTimer];
handles.figLoc.handle = fhandle;
handles.figLoc.IC = handles.curIC;
handles.figLoc.arrows_dip = arrows_dip;
handles.figLoc.scalp = scalp;
handles.figLoc.arrowlen = arrowlen;
if handles.figLoc.fixed_dip
    handles.figLoc.normals = normals; end
guidata(hObject,handles);



function figLoc_update_dipolefit(varargin)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% parse inputs
handles = guidata(varargin{3});
if handles.figLoc.fixed_dip
    normals = handles.figLoc.normals; end

% load surfaces
temp = load(handles.headModel.surfaces);
vertices = temp.surfData(3).vertices(handles.hmInd,:);

% update dipole
[dipoles, ~, ~, rv, handles.figLoc.state] = ...
    dipole_fit(Winv(:,handles.figLoc.IC),handles.K,vertices,handles.figLoc.state);

% update dipole plot
if ~handles.figLoc.fixed_dip
    dipoles.moment = dipoles.moment/norm(dipoles.moment)*handles.figLoc.arrowlen;
    set(handles.figLoc.arrows_dip, ...
        'XData',dipoles.location(1), ...
        'YData',dipoles.location(2), ...
        'ZData',dipoles.location(3), ...
        'UData',dipoles.moment(1), ...
        'VData',dipoles.moment(2), ...
        'WData',dipoles.moment(3));
else
    set(handles.figLoc.arrows_dip, ...
        'XData',dipoles.location(1), ...
        'YData',dipoles.location(2), ...
        'ZData',dipoles.location(3), ...
        'UData',normals(dipoles.L,1)*handles.figLoc.arrowlen, ...
        'VData',normals(dipoles.L,2)*handles.figLoc.arrowlen, ...
        'WData',normals(dipoles.L,3)*handles.figLoc.arrowlen);
end

% update scalp plot
set(handles.figLoc.handle.hScalp,'FaceVertexCData', ...
    handles.figLoc.scalp*Winv(:,handles.figLoc.IC))
maxabs = max(abs(vec(handles.figLoc.scalp*Winv(:,handles.figLoc.IC))));
caxis(handles.figLoc.handle.hAxes,[-maxabs maxabs]);

% update RV text
set(handles.figLoc.textRV,'string',sprintf('Residual Variance: %04.1f%%',rv*100));



function closeFigLoc(varargin)
hObject = varargin{3};

% load handles
handles = guidata(hObject);

% delete figure handle from handles
if isfield(handles,'figLoc')
    handles = rmfield(handles,'figLoc'); end

% delete timer and remove from pauseTimers
locTimerInd = strcmp(get(handles.pauseTimers,'Name'),'locTimer');
stop(handles.pauseTimers(locTimerInd));
delete(handles.pauseTimers(locTimerInd));
handles.pauseTimers(locTimerInd) = [];

% save handles
guidata(hObject,handles);



% --- Executes on button press in togglebuttonLock2.
function togglebuttonLock2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock2


% --- Executes on button press in togglebuttonReject2.
function togglebuttonReject2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject2


% --- Executes on button press in togglebuttonReject1.
function togglebuttonReject1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject1


% --- Executes on button press in togglebuttonLock1.
function togglebuttonLock1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock1


% --- Executes on selection change in popupmenuEEG.
function popupmenuEEG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuEEG


% --- Executes during object creation, after setting all properties.
function popupmenuEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuInfo.
function popupmenuInfo_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInfo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInfo


% --- Executes during object creation, after setting all properties.
function popupmenuInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSortICs.
function pushbuttonSortICs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSortICs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    updateICs(hObject,true); end




% --- Executes on mouse press over axes background.
function axesIC1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(1);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(2);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(3);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(4);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(5);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(6);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(7);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isSourceMode
    handles.curIC = handles.ics(8);
    set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);
end
% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in togglebuttonReject8.
function togglebuttonReject8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject8


% --- Executes on button press in togglebuttonLock8.
function togglebuttonLock8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock8




% --- Executes on button press in togglebuttonReject7.
function togglebuttonReject7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject7


% --- Executes on button press in togglebuttonLock7.
function togglebuttonLock7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock7




% --- Executes on button press in togglebuttonReject6.
function togglebuttonReject6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject6


% --- Executes on button press in togglebuttonLock6.
function togglebuttonLock6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock6




% --- Executes on button press in togglebuttonReject5.
function togglebuttonReject5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject5


% --- Executes on button press in togglebuttonLock5.
function togglebuttonLock5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock5




% --- Executes on button press in togglebuttonReject4.
function togglebuttonReject4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject4


% --- Executes on button press in togglebuttonLock4.
function togglebuttonLock4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock4




% --- Executes on button press in togglebuttonReject3.
function togglebuttonReject3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject3


% --- Executes on button press in togglebuttonLock3.
function togglebuttonLock3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock3


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in togglebuttonLock2.
function togglebutton17_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock2


% --- Executes on button press in togglebuttonReject2.
function togglebutton18_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject2


% --- Executes on button press in togglebuttonReject3.
function togglebutton22_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject3


% --- Executes on button press in togglebuttonLock3.
function togglebutton21_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock3


% --- Executes on button press in togglebuttonLock4.
function togglebutton26_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonLock4


% --- Executes on button press in togglebuttonReject4.
function togglebutton25_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonReject4


% --- Executes on button press in togglebuttonSwitch.
function togglebuttonSwitch_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSwitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSwitch


% --- Executes on button press in radiobuttonCH.
function radiobuttonCH_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.isSourceMode = 0;

% Populate channel maps
for it = 1:handles.nplots
    set(handles.(['axesIC' int2str(it)]),'NextPlot','add');
    cla(handles.(['axesIC' int2str(it)]));
    set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]));
    topoplotCH([],handles.chanlocs(it), 'style', 'blank', 'electrodes','on', ...
        'electcolor',[1 0 0],'emarkersize',20);
end
for it = 1:handles.nplots
    set(handles.(['panelIC' int2str(it)]),'Title',handles.chanlocs(it).labels);
end

set(handles.pushbuttonLocalize,'string','Localize');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in radiobuttonIC.
function radiobuttonIC_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.isSourceMode = 1;

% Populate scalp maps
for it = 1:handles.nplots
    set(handles.(['axesIC' int2str(it)]),'NextPlot','add');
    cla(handles.(['axesIC' int2str(it)]));
    set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]));
    [~,Zi,~,Xi,Yi,intx,inty] = topoplotFast_LRBF(rand(size(handles.chanlocs)), handles.chanlocs);
    set(handles.(['panelIC' int2str(it)]),'Title',['IC' int2str(handles.ics(it))])
end

% Generate scalp map interpolation matrix (jerry rigged)
nChan = length(handles.chanlocs);
in = eye(nChan);
out = zeros(32^2,nChan);
for it = 1:nChan
    op = rbfcreate(double([inty;intx]),in(:,it)','RBFFunction', 'linear');
    out(:,it) = rbfinterp(double([Xi(:),Yi(:)]'), op);
end
handles.topoMat = out/in;
handles.topoNaNMask = isnan(Zi);
handles.topoNPixel = size(out,1);
handles.topoMat(handles.topoNaNMask,:) = [];

set(handles.pushbuttonLocalize,'string',['Localize (IC ' int2str(handles.curIC) ')']);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbuttonViewStream.
function pushbuttonViewStream_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonViewStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

