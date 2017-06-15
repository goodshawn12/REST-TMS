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

% Last Modified by GUIDE v2.5 13-Jun-2017 17:53:51

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
handles.ntopo = 8;
handles.curIC = 1; % cursor for PSD plot
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

% Populate scalp maps
for it = 1:handles.ntopo
    set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]))
    [~,Zi,~,Xi,Yi,intx,inty] = topoplotFast_LRBF(rand(size(handles.chanlocs)), handles.chanlocs);
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

% Create scalp map timer
topoTimer = timer('Period',round(1/handles.ntopo*1000)/1000,'ExecutionMode','fixedRate','TimerFcn',{@vis_topo,hObject},'StartDelay',0.2,'Tag','topoTimer','Name','topoTimer');

% Create data timer (starts as power spectrum)
% infoTimer = timer('Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');
infoTimer = timer('Period',round(1/handles.ntopo*1000)/1000,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');

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
    set(handles.pushbuttonLocalize,'HitTest','off','visible','off')
    
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
assignin('base','learning_rate',nan(1,bufflen*handles.srate));

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

% run pipline on calibration data
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);
assignin('base','pipeline',pipeline);



function vis_topo(varargin)
% get the updated stream buffer
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');

% get handles
handles = guidata(varargin{3});

% update topo plot
it = mod(get(varargin{1},'TasksExecuted')-1,handles.ntopo)+1;
% if it==1
%     updateICs(varargin{3}); end
hstr = ['axesIC' int2str(it)];
hand = get(handles.(hstr),'children');

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


function infoPSD(varargin)

% plot PSD of selected IC
try  %#ok<*TRYNC>
    secs2samp = 5; % seconds
    
    W = evalin('base','pipeline.state.icaweights');
    % if isempty(W), W = evalin('base','Wn'); end
    sphere = evalin('base','pipeline.state.icasphere');
    handles = guidata(varargin{3});
        
    srate = evalin('base',[handles.bufferName '.srate']);
    data = evalin('base',[handles.bufferName '.data{end}']); % !!! make this more robust
    if all(data(:,end)==0)
        mark=length(data);
        while true
            mark = find(data(1,1:mark)~=0,1,'last');
            if all(data(:,mark)~=0)
                break; end
        end
        data = data(:,max(1,mark-srate*secs2samp+1):mark);
    else
        data = data(:,max(1,end-srate*secs2samp+1):end);
    end
    
    data = bsxfun(@minus,data,mean(data,2));
    data = W*sphere*data;
    
    % find index to be updated
    it = mod(get(varargin{1},'TasksExecuted')-1,handles.ntopo)+1;
    hstr = ['axes' int2str(it)];

    data = data(it,:);    
    [data,f] = pwelch(data,[],[],[],srate);
    
    maxFreq = 40;
    ind = f <= maxFreq;
    plot(handles.(hstr),f(ind),db(data(ind))) % !!!?
    grid(handles.(hstr),'on');
    axis(handles.(hstr),'tight')
    set(handles.(hstr),'XTick',[0 10:10:f(end)],'box','off')
    %     xlabel(handles.(hstr),'Frequency (Hz)')
    %     ylabel(handles.(hstr),'Power/Frequency (dB/Hz)')

end



function genICSelectGUI(hObject,handles)
temp = get(handles.figure1,'Position');
fhandle = figure('toolbar','none','Menubar','none','Name','IC Select','position',[1 1 temp(3:4)],'Resize','on','Colormap',colormap('jet'),'DeleteFcn',{@closeFigIC,hObject});
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


% --- Executes on mouse press over axes background.
function axesIC1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonSelect.
function pushbuttonSelect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonPause.
function pushbuttonPause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonLocalize.
function pushbuttonLocalize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLocalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axesIC2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on mouse press over axes background.
function axesIC8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on mouse press over axes background.
function axesIC7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on mouse press over axes background.
function axesIC6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on mouse press over axes background.
function axesIC5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on mouse press over axes background.
function axesIC4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on mouse press over axes background.
function axesIC3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

% Hint: get(hObject,'Value') returns toggle state of radiobuttonCH


% --- Executes on button press in radiobuttonIC.
function radiobuttonIC_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonIC


% --- Executes on button press in pushbuttonViewStream.
function pushbuttonViewStream_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonViewStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
