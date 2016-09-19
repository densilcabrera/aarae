% DO NOT EDIT THIS INITIALIZATION FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!
function varargout = audio_recorder(varargin)
% AUDIO_RECORDER MATLAB code for audio_recorder.fig
%      AUDIO_RECORDER, by itself, creates a new AUDIO_RECORDER or raises the existing
%      singleton*.
%
%      H = AUDIO_RECORDER returns the handle to a new AUDIO_RECORDER or the handle to
%      the existing singleton*.
%
%      AUDIO_RECORDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUDIO_RECORDER.M with the given input arguments.
%
%      AUDIO_RECORDER('Property','Value',...) creates a new AUDIO_RECORDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before audio_recorder_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to audio_recorder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help audio_recorder

% Last Modified by GUIDE v2.5 25-Jul-2015 18:45:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @audio_recorder_OpeningFcn, ...
                   'gui_OutputFcn',  @audio_recorder_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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






% Note: for AARAE Release 7 (July 2015), channel mapping was introduced
% instead of specifying the number of contiguous channels from 1. However,
% the relevant function names and variable names have been retained from
% the previous versions.




% --- Executes just before audio_recorder is made visible.
function audio_recorder_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to audio_recorder (see VARARGIN)


% This next couple of lines checks if the GUI is being called from the main
% window, otherwise it doesn't run.
dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'main_stage1'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.main_stage1 = varargin{mainGuiInput+1};
    
    fontsize
    
    % Call the 'desktop'
    hMain = getappdata(0,'hMain');
    handles.signaldata = getappdata(hMain,'testsignal');
    handles.recording = [];
    handles.position = get(handles.audio_recorder,'Position');
    handles.OUT_axes_position = get(handles.OUT_axes,'Position');
    UserData.state = false;
    set(handles.stop_btn,'UserData',UserData);
    mainHandles = guidata(handles.main_stage1);
    inputdevinfo = dspAudioDeviceInfo('inputs');
    inputnames = {inputdevinfo.name}';
    inputnames = regexprep(inputnames,'\s\(Windows DirectSound\)','');
    inputnames = regexprep(inputnames,'\s\(ASIO\)','');
    inputnames = regexprep(inputnames,'\s\(Core Audio\)','');
    set(handles.inputdev_popup,'String',inputnames)
    outputdevinfo = dspAudioDeviceInfo('outputs');
    outputnames = {outputdevinfo.name}';
    outputnames = regexprep(outputnames,'\s\(Windows DirectSound\)','');
    outputnames = regexprep(outputnames,'\s\(ASIO\)','');
    outputnames = regexprep(outputnames,'\s\(Core Audio\)','');
    set(handles.outputdev_popup,'String',outputnames)
    if ~isempty(handles.signaldata) && ismatrix(handles.signaldata.audio) && ~strcmp(handles.signaldata.datatype,'syscal')% If there's a signal loaded in the 'desktop'...
        % Allow visibility of playback option along with the specs of
        % the playback signal
        mainHandles = guidata(handles.main_stage1);
        selectednode = mainHandles.mytree.getSelectedNodes;
        set(handles.pb_enable,'Visible','on','Value',1);
        set(handles.IN_name,'String',['rec_' selectednode(1).getName.char])
        handles.outputdata = handles.signaldata;
        if mainHandles.alternate==1 && isfield(handles.outputdata,'audio2')
            % swap audio with audio2
            audio = handles.outputdata.audio2;
            handles.outputdata.audio2 = handles.outputdata.audio;
            handles.outputdata.audio = audio;
        end
        handles.dur = length(handles.outputdata.audio)/handles.outputdata.fs;
        handles.t = linspace(0,handles.dur,length(handles.outputdata.audio));
        output_settings{1} = ['Playback audio loaded: ' selectednode(1).getName.char];
        output_settings{2} = ['Number of audio channels: ' num2str(size(handles.outputdata.audio,2))];
        output_settings{3} = ['Sampling frequency = ',num2str(handles.outputdata.fs),' samples/s'];
        %output_settings{4} = ['Bit depth = ',num2str(handles.outputdata.nbits)];
        output_settings{4} = ['Duration = ',num2str(handles.dur),' s'];
        set(handles.IN_numchsout,'String','1')
        % currently channel mapping of output is only allowed for mono
        % signals
        if size(handles.outputdata.audio,2) == 1
            set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','on')
        else
            set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','off')
        end
        pixels = get_axes_width(handles.OUT_axes);
        [t, line] = reduce_to_width(handles.t', handles.outputdata.audio, pixels, [-inf inf]);
        plot(handles.OUT_axes,t,line)
        set(handles.OUT_axes,'tag','OUT_axes')
        set(handles.output_settings,'String',output_settings);
        set(handles.inputdev_popup,'Value',getappdata(hMain,'audio_recorder_input'));
        set(handles.outputdev_popup,'Value',getappdata(hMain,'audio_recorder_output'));
        set(handles.IN_numchs,'String',num2str(getappdata(hMain,'audio_recorder_numchs')));
        set(handles.text1,'String','Add time');
        set(handles.IN_duration,'String',num2str(getappdata(hMain,'audio_recorder_duration')));
        set(handles.IN_fs,'Enable','off');
        set(handles.IN_fs,'String','-');
        %set(handles.IN_nbits,'Enable','off');
        %set(handles.IN_nbits,'String','-');
        set(handles.IN_qdur,'String',num2str(getappdata(hMain,'audio_recorder_qdur')));
        set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
        handles.numchs = str2num(get(handles.IN_numchs,'String'));
        handles.addtime = str2double(get(handles.IN_duration,'String'));
        set(handles.SilenceRequestCheckBox,'Value',getappdata(hMain,'audio_recorder_silencerequest'));
        handles.silenceplease.audio = mainHandles.silenceplease.audio;
        handles.silenceplease.fs = mainHandles.silenceplease.fs;
        handles.thankyou.audio = mainHandles.thankyou.audio;
        handles.thankyou.fs = mainHandles.thankyou.fs;
        handles.fs = handles.outputdata.fs;
        %handles.nbits = handles.outputdata.nbits;
        handles.qdur = str2double(get(handles.IN_qdur,'String'));
        handles.buffer = str2double(get(handles.IN_buffer,'String'));
        xlim(handles.IN_axes,[0 round(handles.dur+handles.addtime)])
        xlim(handles.OUT_axes,[0 round(handles.dur+handles.addtime)])
    else
        % If there's no signal loaded in the desktop just allocate memory
        % space for the signal to be recorded
        set(handles.pb_enable,'Visible','off','Value',0);
        set(handles.output_panel,'Visible','off');
        set(handles.text1,'String','Duration');
        set(handles.IN_name,'String','recording')
        set(handles.inputdev_popup,'Value',getappdata(hMain,'audio_recorder_input'));
        set(handles.IN_numchs,'String',num2str(getappdata(hMain,'audio_recorder_numchs')));
        set(handles.IN_duration,'String',num2str(getappdata(hMain,'audio_recorder_duration')));
        set(handles.IN_fs,'Enable','on');
        set(handles.IN_fs,'String',num2str(getappdata(hMain,'audio_recorder_fs')));
        %set(handles.IN_nbits,'Enable','on');
        %set(handles.IN_nbits,'String',num2str(getappdata(hMain,'audio_recorder_nbits')));
        set(handles.IN_qdur,'String',num2str(getappdata(hMain,'audio_recorder_qdur')));
        set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
        set(handles.OUT_axes,'Visible','off');
        set(handles.audio_recorder,'Position',handles.position-[0 0 0 handles.OUT_axes_position(4)]);
        set(handles.SilenceRequestCheckBox,'Value',getappdata(hMain,'audio_recorder_silencerequest'));
        handles.silenceplease.audio = mainHandles.silenceplease.audio;
        handles.silenceplease.fs = mainHandles.silenceplease.fs;
        handles.thankyou.audio = mainHandles.thankyou.audio;
        handles.thankyou.fs = mainHandles.thankyou.fs;
        handles.numchs = str2num(get(handles.IN_numchs,'String'));
        handles.dur = str2double(get(handles.IN_duration,'String'));
        handles.fs = str2double(get(handles.IN_fs,'String'));
        %handles.nbits = str2double(get(handles.IN_nbits,'String'));
        handles.qdur = str2double(get(handles.IN_qdur,'String'));
        handles.buffer = str2double(get(handles.IN_buffer,'String'));
        xlim(handles.IN_axes,[0 round(handles.dur)])
        xlim(handles.OUT_axes,[0 round(handles.dur)])
    end
    if isfield(mainHandles,'syscalstats') && ~isempty(mainHandles.syscalstats)
        handles.savenewsyscal = 0;
        handles.syscalstats = mainHandles.syscalstats;
        if isfield(handles.syscalstats,'latency')
            if get(handles.pb_enable,'Value') == 1
                set(handles.delay_chk,'Enable','on','Value',1)
            else
                set(handles.delay_chk,'Enable','off','Value',0)
            end
            set(handles.delaytext,'String',[num2str(handles.syscalstats.latency) ' samples'])
        end
        if isfield(handles.syscalstats,'cal')
            set(handles.cal_chk,'Enable','on','Value',1)
            set(handles.caltext,'String',[num2str(handles.syscalstats.cal) ' dB'])
        end
        if isfield(handles.syscalstats,'audio')
            if get(handles.pb_enable,'Value') == 1
                set(handles.invfilter_chk,'Enable','on','Value',1)
            else
                set(handles.invfilter_chk,'Enable','off','Value',0)
            end
            set(handles.invftext,'String','Available')
        end
    else
        handles.savenewsyscal = 1;
        handles.syscalstats = struct([]);
    end
end

% Update handles structure
guidata(hObject, handles);

if dontOpen
   disp('-----------------------------------------------------');
   disp('This function is part of the AARAE framework, it is') 
   disp('not a standalone function. To call this function,')
   disp('click on the appropriate calling button on the main');
   disp('Window. E.g.:');
   disp('   Record');
   disp('-----------------------------------------------------');
else
   uiwait(hObject);
end

% UIWAIT makes audio_recorder wait for user response (see UIRESUME)
% uiwait(handles.audio_recorder);


% --- Outputs from this function are returned to the command line.
function varargout = audio_recorder_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMain = getappdata(0,'hMain');
if handles.savenewsyscal == 1 && ~isempty(handles.syscalstats)
    setappdata(hMain,'savenewsyscalstats',1);
    setappdata(hMain,'syscalstats',handles.syscalstats)
else
    setappdata(hMain,'savenewsyscalstats',0)
end

% Get default command line output from handles structure
varargout{1} = handles.recording;
delete(hObject);


% --- Executes on button press in record_btn.
function record_btn_Callback(hObject, ~, handles) %#ok : called with record button
% hObject    handle to record_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call handles from main window
%mainHandles = guidata(handles.main_stage1);
set([handles.cancel_btn handles.load_btn handles.preview_btn handles.syscal_btn],'Enable','off')


% ASK FOR SILENCE
if get(handles.SilenceRequestCheckBox,'Value') == 1
    %SilenceRequestCheckBox
    set(hObject,'Enable','off');
    pause on
    pause(0.000001)
    pause off
    % Set playback object
    if isfield(handles,'outputdata')
        fs = handles.outputdata.fs;
    else
        fs = handles.fs;
    end
    outputs = cellstr(get(handles.outputdev_popup,'String'));
    outputdevname = outputs{get(handles.outputdev_popup,'Value')};
    
        % Set playback audio
    handles.hsr1 = dsp.SignalSource;
    if handles.silenceplease.fs ~= fs
        gcd_fs = gcd(handles.silenceplease.fs,fs); % greatest common denominator
        thankyou = resample(handles.silenceplease.audio,fs/gcd_fs,handles.silenceplease.fs/gcd_fs);
    else
        thankyou = handles.silenceplease.audio;
    end
    thankyou = thankyou(:,1,1,1,1,1); % make sure it is 1 channel (it should be anyway)
    if length(str2num(get(handles.IN_numchsout,'String'))) > 1
            thankyou = repmat(thankyou,1,length(str2num(get(handles.IN_numchsout,'String'))));
    end
    handles.hap = dsp.AudioPlayer('DeviceName',outputdevname,...
        'SampleRate',fs,...
        'QueueDuration',handles.qdur,...
        'BufferSizeSource','Property',...
        'BufferSize',handles.buffer);
    
    %if get(handles.invfilter_chk,'Value') == 1, quietplease = filter(handles.syscalstats.audio2,1,quietplease); end
    handles.hsr1.Signal = [thankyou;zeros(floor(handles.hap.QueueDuration*handles.fs),1)];
    handles.hsr1.SamplesPerFrame = 2048;
    guidata(hObject,handles)
    ncycles = ceil(length(handles.hsr1.Signal)/handles.hsr1.SamplesPerFrame);
    set(hObject,'BackgroundColor','red');
    set(handles.stop_btn,'Visible','on');
    % Initialize playback routine
    pause on
        %UserData = get(handles.stop_btn,'UserData');
        h = helpdlg('Silence request. Press OK to start recording.','AARAE info');
        for i = 1:ncycles
           UserData = get(handles.stop_btn,'UserData');
           if UserData.state == false
               step(handles.hap,step(handles.hsr1));
           else
               break
           end
           %waitbar(i/ncycles,h)
           pause(0.0000001)
        end
        %if i == ncycles, delete(h); end
    pause off
    % Release playback and audio data objects
    release(handles.hap)
    release(handles.hsr1)
    uiwait(h)
end


if get(handles.pb_enable,'Value') == 1
    % Simultaneous playback and record routine
    set(hObject,'Enable','off');
    pause on
    pause(0.000001)
    pause off
    % Set playback object
    outputs = cellstr(get(handles.outputdev_popup,'String'));
    outputdevname = outputs{get(handles.outputdev_popup,'Value')};
    handles.hap = dsp.AudioPlayer('DeviceName',outputdevname,...
        'SampleRate',handles.outputdata.fs,...
        'QueueDuration',handles.qdur,...
        'BufferSizeSource','Property',...
        'BufferSize',handles.buffer);
    % Set record object
    inputs = cellstr(get(handles.inputdev_popup,'String'));
    inputdevname = inputs{get(handles.inputdev_popup,'Value')};
    handles.har = dsp.AudioRecorder('DeviceName',inputdevname,...
        'SampleRate',handles.outputdata.fs,...
        'QueueDuration',handles.qdur,...
        'OutputDataType','double',...
        'ChannelMappingSource','Property',...
        'ChannelMapping',handles.numchs,...
        'BufferSizeSource','Property',...
        'BufferSize',handles.buffer);
    % Set playback audio
    handles.hsr1 = dsp.SignalSource;
    if size(handles.outputdata.audio,2) == 1 && length(str2num(get(handles.IN_numchsout,'String'))) > 1
        if get(handles.sim_chk,'Value') == 1
            playbackaudio = repmat(handles.outputdata.audio,1,length(str2num(get(handles.IN_numchsout,'String'))));
            addtimeok = false;
        else
            playbackaudio = zeros((length(handles.outputdata.audio)+floor(handles.addtime*handles.fs))*length(str2num(get(handles.IN_numchsout,'String'))),length(str2num(get(handles.IN_numchsout,'String'))));
            j = 1;
            playbackaudio(1:length(handles.outputdata.audio),j) = handles.outputdata.audio;
            for j = length(str2num(get(handles.IN_numchsout,'String')))-1
                playbackaudio((length(handles.outputdata.audio)+floor(handles.addtime*handles.fs))*j+1:(length(handles.outputdata.audio)+floor(handles.addtime*handles.fs))*j+length(handles.outputdata.audio),j+1) = handles.outputdata.audio;
            end
            addtimeok = true;
        end
    else
        playbackaudio = handles.outputdata.audio;
        addtimeok = false;
    end
    %if get(handles.invfilter_chk,'Value') == 1, playbackaudio = filter(handles.syscalstats.audio2,1,playbackaudio); end
    if ~addtimeok
        handles.hsr1.Signal = [playbackaudio;zeros(floor((handles.addtime+handles.hap.QueueDuration)*handles.fs),size(playbackaudio,2))];
    else
        handles.hsr1.Signal = [playbackaudio;zeros(floor(handles.hap.QueueDuration*handles.fs),size(playbackaudio,2))];
    end
    handles.hsr1.SamplesPerFrame = handles.har.SamplesPerFrame;
    guidata(hObject,handles)
    handles.rec = [];
    ncycles = ceil(length(handles.hsr1.Signal)/handles.har.SamplesPerFrame);
    audio = zeros(ncycles*handles.har.SamplesPerFrame,length(handles.numchs));
    set(hObject,'BackgroundColor','red');
    set(handles.stop_btn,'Visible','on');
    % Initialize playback/record routine
    pause on
    try
        UserData = get(handles.stop_btn,'UserData');
        h = waitbar(0,'Recording in progress...','Name','AARAE info');
        handles.recordtime = datestr(now);
        for i = 1:ncycles
           UserData = get(handles.stop_btn,'UserData');
           if UserData.state == false
               audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:) = step(handles.har);
               step(handles.hap,step(handles.hsr1));
           else
               break
           end
           waitbar(i/ncycles,h)
           pause(0.0000001)
        end
        if i == ncycles, delete(h); end
    catch sthgwrong
        UserData.state = true;
        handles.rec = [];
        syswarning = sthgwrong.message;
        warndlg(syswarning,'AARAE info')
    end
    pause off
    % Check recording and adjust for QueueDuration latency
    handles.rec = audio;
    if ~isempty(handles.rec)
        handles.rec = handles.rec(handles.hap.QueueDuration*handles.fs:end,:);
        if UserData.state == false
            if ~addtimeok
                handles.rec = handles.rec(1:(size(handles.outputdata.audio,1)+handles.addtime*handles.fs),:);
            else
                handles.rec = handles.rec(1:size(playbackaudio,1),:);
                handles.rec = reshape(handles.rec,(size(handles.outputdata.audio,1)+handles.addtime*handles.fs),length(handles.numchs),1,1,length(str2num(get(handles.IN_numchsout,'String'))));
            end
        else
            UserData.state = false;
            set(handles.stop_btn,'UserData',UserData);
        end
    % Plot recording
        time = linspace(0,size(handles.rec,1)/handles.fs,length(handles.rec));
        pixels = get_axes_width(handles.IN_axes);
        [t, linea] = reduce_to_width(time', handles.rec, pixels, [-inf inf]);
        plot(handles.IN_axes,t,linea);
        xlim(handles.IN_axes,[0 max(time)])
    end
    % Release playback, record and audio data objects
    release(handles.hap)
    release(handles.har)
    release(handles.hsr1)
else
    % Record-only routine
    dur = handles.dur*handles.fs;
    set(hObject,'Enable','off');
    pause on
    pause(0.000001)
    pause off
    % Set record object
    inputs = cellstr(get(handles.inputdev_popup,'String'));
    inputdevname = inputs{get(handles.inputdev_popup,'Value')};
    handles.har = dsp.AudioRecorder('DeviceName',inputdevname,...
        'SampleRate',handles.fs,...
        'OutputDataType','double',...
        'ChannelMappingSource','Property',...
        'ChannelMapping',handles.numchs,...
        'BufferSizeSource','Property',...
        'BufferSize',handles.buffer,...
        'QueueDuration',handles.qdur);
    guidata(hObject,handles)
    ncycles = ceil(dur/handles.har.SamplesPerFrame);
    handles.rec = [];
    audio = zeros(ncycles*handles.har.SamplesPerFrame,length(handles.numchs));
    set(hObject,'BackgroundColor','red');
    set(handles.stop_btn,'Visible','on');
    % Initialize record routine
    pause on
    try
        UserData = get(handles.stop_btn,'UserData');
        h = waitbar(0,'Recording in progress...','Name','AARAE info');
        handles.recordtime = datestr(now);
        for i = 1:ncycles
           UserData = get(handles.stop_btn,'UserData');
           if UserData.state == false
               audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:) = step(handles.har);
           else
               break
           end
           waitbar(i/ncycles,h)
           pause(0.0000001)
        end
        if i == ncycles, delete(h); end
    catch sthgwrong
        UserData.state = true;
        handles.rec = [];
        syswarning = sthgwrong.message;
        warndlg(syswarning,'AARAE info')
    end
    pause off
    handles.rec = audio;
    % Check recording and adjust for Duration
    if ~isempty(handles.rec)
        if UserData.state == false
            handles.rec = handles.rec(1:dur,:);
        else
            UserData.state = false;
            set(handles.stop_btn,'UserData',UserData);
        end
        % Plot recording
        time = linspace(0,size(handles.rec,1)/handles.fs,length(handles.rec));
        pixels = get_axes_width(handles.IN_axes);
        [t, linea] = reduce_to_width(time', handles.rec, pixels, [-inf inf]);
        plot(handles.IN_axes,t,linea);
        xlim(handles.IN_axes,[0 max(time)])
    end
    % Release record object
    release(handles.har)
end


% SAY THANKYOU
if get(handles.SilenceRequestCheckBox,'Value') == 1
    %SilenceRequestCheckBox
    set(hObject,'Enable','off');
    pause on
    pause(0.000001)
    pause off
    % Set playback object
    if isfield(handles,'outputdata')
        fs = handles.outputdata.fs;
    else
        fs = handles.fs;
    end
    outputs = cellstr(get(handles.outputdev_popup,'String'));
    outputdevname = outputs{get(handles.outputdev_popup,'Value')};
    
        % Set playback audio
    handles.hsr1 = dsp.SignalSource;
    if handles.thankyou.fs ~= fs
        gcd_fs = gcd(handles.thankyou.fs,fs); % greatest common denominator
        thankyou = resample(handles.thankyou.audio,fs/gcd_fs,handles.thankyou.fs/gcd_fs);
    else
        thankyou = handles.thankyou.audio;
    end
    thankyou = thankyou(:,1,1,1,1,1); % make sure it is 1 channel (it should be anyway)
    if length(str2num(get(handles.IN_numchsout,'String'))) > 1
            thankyou = repmat(thankyou,1,length(str2num(get(handles.IN_numchsout,'String'))));
    end
    handles.hap = dsp.AudioPlayer('DeviceName',outputdevname,...
        'SampleRate',fs,...
        'QueueDuration',handles.qdur,...
        'BufferSizeSource','Property',...
        'BufferSize',handles.buffer);
    
    %if get(handles.invfilter_chk,'Value') == 1, quietplease = filter(handles.syscalstats.audio2,1,quietplease); end
    handles.hsr1.Signal = [thankyou;zeros(floor(handles.hap.QueueDuration*handles.fs),1)];
    handles.hsr1.SamplesPerFrame = 2048;
    guidata(hObject,handles)
    ncycles = ceil(length(handles.hsr1.Signal)/handles.hsr1.SamplesPerFrame);
    set(hObject,'BackgroundColor','red');
    set(handles.stop_btn,'Visible','on');
    % Initialize playback routine
    pause on
        h = helpdlg('Thankyou!','AARAE info');
        for i = 1:ncycles
           UserData = get(handles.stop_btn,'UserData');
           if UserData.state == false
               step(handles.hap,step(handles.hsr1));
           else
               break
           end
           pause(0.0000001)
        end
        if i == ncycles, delete(h); end
    pause off
    % Release playback and audio data objects
    release(handles.hap)
    release(handles.hsr1)
end




set(handles.record_btn,'BackgroundColor',[0.94 0.94 0.94]);
set(handles.record_btn,'Enable','on');
set(handles.stop_btn,'Visible','off');
set([handles.load_btn handles.preview_btn handles.cancel_btn handles.syscal_btn],'Enable','on')

guidata(hObject,handles);


% --- Executes on button press in load_btn.
function load_btn_Callback(hObject, ~, handles) %#ok : called when sending recording to main window
% hObject    handle to load_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMain = getappdata(0,'hMain');
% Obtain handles using GUIDATA with the caller's handle

handles.testsignal = handles.rec;
if isfield(handles,'outputdata')
    handles.recording = handles.outputdata;
end

% Warnings...
if isempty(handles.rec)
    warndlg('No signal recorded!');
    setappdata(hMain,'testsignal',[]);
else
    pause on
    set(hObject,'BackgroundColor','red');
    set([handles.record_btn,handles.syscal_btn,handles.cancel_btn,handles.load_btn,handles.preview_btn],'Enable','off')
    pause(0.000001)
    pause off
    hMain = getappdata(0,'hMain');
    if get(handles.delay_chk,'Value') == 1, handles.rec = [handles.rec(handles.syscalstats.latency:end,:);zeros(handles.syscalstats.latency,size(handles.rec,2))]; end
    if get(handles.invfilter_chk,'Value') == 1, handles.rec = filter(handles.syscalstats.audio2,1,handles.rec); end
    handles.recording.audio = handles.rec;
    if get(handles.pb_enable,'Value') && isfield(handles.outputdata,'audio2')
        handles.recording.audio2 = handles.outputdata.audio2;
    elseif get(handles.pb_enable,'Value') && ~isfield(handles.outputdata,'audio2')
        handles.recording.audio2 = handles.outputdata.audio;
    else
        handles.recording.audio2 = [];
    end
    handles.recording.fs = handles.fs;
    %handles.recording.nbits = handles.nbits;
    try
        handles.recording.chanID = cellstr([repmat('Chan',size(handles.recording.audio,2),1) num2str((handles.numchs)')]);
    
    catch
        handles.recording.chanID = cellstr([repmat('Chan',size(handles.recording.audio,2),1) num2str((1:size(handles.recording.audio,2))')]);    
    end
    % Potentially add dim5ID or outchanID here (once its name is defined).
    % should it be a property or a first order field?
    if get(handles.cal_chk,'Value') == 1
        handles.recording.cal = handles.syscalstats.cal(1:size(handles.recording.audio,2));
        handles.recording.properties.units = handles.syscalstats.units;
        handles.recording.properties.units_ref = handles.syscalstats.units_ref;
        handles.recording.properties.units_type = handles.syscalstats.units_type;
    end
    handles.recording.history = cell(2,4);
    if isfield(handles,'outputdata')
        handles.recording.history{1,1} = handles.recordtime;
        handles.recording.history{1,2} = 'Recorded with playback';
        if isfield(handles.outputdata,'name')
            handles.recording.history{1,3} = handles.outputdata.name;
        end
        if isfield(handles.outputdata,'history')
            handles.recording.history{1,4} = handles.outputdata.history;
        end
        inputs = cellstr(get(handles.inputdev_popup,'String'));
        inputdevname = inputs{get(handles.inputdev_popup,'Value')};
        outputs = cellstr(get(handles.outputdev_popup,'String'));
        outputdevname = outputs{get(handles.outputdev_popup,'Value')};
    else
        handles.recording.history{1,1} = handles.recordtime;
        handles.recording.history{1,2} = 'Recorded without playback';
        handles.recording.history{1,3} = 'AARAE audio_recorder';
        inputs = cellstr(get(handles.inputdev_popup,'String'));
        inputdevname = inputs{get(handles.inputdev_popup,'Value')};
        outputdevname = 'none';
    end
    handles.recording.history{2,1} = handles.recordtime;
    handles.recording.history{2,2} = 'Recorded';
    handles.recording.history{2,3} = 'AARAE audio_recorder';
    handles.recording.history{2,4} = {'fs', handles.fs;...
        'dur', handles.dur;...
        'numchs', handles.numchs;...
        'qdur', handles.qdur;...
        'buffer', handles.buffer;...
        'output device', outputdevname;...
        'input device', inputdevname};
    name = matlab.lang.makeValidName(get(handles.IN_name,'String'));
    if isempty(name), name = 'untitled'; end
    setappdata(hMain,'signalname',name);
end
setappdata(hMain,'syscalstats',handles.syscalstats);
guidata(hObject,handles)
uiresume(handles.audio_recorder);

% --- Executes on button press in cancel_btn.
function cancel_btn_Callback(~, ~, handles) %#ok : called with the Cancel button
% hObject    handle to cancel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear the desktop if recording is canceled
hMain = getappdata(0,'hMain');
setappdata(hMain,'testsignal',[]);

uiresume(handles.audio_recorder);


function IN_numchs_Callback(hObject, ~, handles) %#ok number of channels input
% hObject    handle to IN_numchs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_numchs as text
%        str2double(get(hObject,'String')) returns contents of IN_numchs as a double

% allow multiple numbers to be entered for channel mapping
numchs = str2num(get(handles.IN_numchs, 'string')); 
hMain = getappdata(0,'hMain');
% Check user's input
numchs = numchs(~isnan(numchs));
numchs = numchs(numchs>0);
numchs = round(numchs);
numchs = unique(numchs);
if isempty(numchs)
    set(hObject,'String',num2str(handles.numchs));
    warndlg('All inputs MUST be real positive numbers!');
else
    handles.numchs = numchs;
    setappdata(hMain,'audio_recorder_numchs',numchs)
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function IN_numchs_CreateFcn(hObject, ~, ~) %#ok number of channels input box creation
% hObject    handle to IN_numchs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IN_duration_Callback(hObject, ~, handles) %#ok : duration input box
% hObject    handle to IN_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_duration as text
%        str2double(get(hObject,'String')) returns contents of IN_duration as a double

% Get duration input
dur = round(str2double(get(handles.IN_duration, 'string')));
hMain = getappdata(0,'hMain');
% Check user's input
if (isnan(dur)||dur<=0)
    if get(handles.pb_enable,'Value') == 1
        set(hObject,'String',num2str(handles.addtime));
    else
        set(hObject,'String',num2str(handles.dur));
    end
    warndlg('All inputs MUST be real positive numbers!');
else
    set(hObject,'String',num2str(dur))
    setappdata(hMain,'audio_recorder_duration',dur)
    handles.dur = dur;
    handles.addtime = dur;
    if get(handles.pb_enable,'Value') == 1
        xlim(handles.IN_axes,[0 round(handles.dur+handles.addtime)])
        xlim(handles.OUT_axes,[0 round(handles.dur+handles.addtime)])
    else
        xlim(handles.IN_axes,[0 round(handles.dur)])
        xlim(handles.OUT_axes,[0 round(handles.dur)])
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function IN_duration_CreateFcn(hObject, ~, ~) %#ok : duration input box creation
% hObject    handle to IN_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IN_fs_Callback(hObject, ~, handles) %#ok : sampling frequency input box
% hObject    handle to IN_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_fs as text
%        str2double(get(hObject,'String')) returns contents of IN_fs as a double

% Get sampling frequency input
fs = str2double(get(handles.IN_fs, 'string'));
hMain = getappdata(0,'hMain');
% Check user's input
if (isnan(fs)||fs<=0)
    set(hObject,'String',num2str(handles.fs))
    warndlg('All inputs MUST be real positive numbers!');
else
    handles.fs = fs;
    setappdata(hMain,'audio_recorder_fs',fs)
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function IN_fs_CreateFcn(hObject, ~, ~) %#ok : sampling frequency input box creation
% hObject    handle to IN_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%function IN_nbits_Callback(hObject, ~, handles) %#ok : bit depth input box
% hObject    handle to IN_nbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_nbits as text
%        str2double(get(hObject,'String')) returns contents of IN_nbits as a double

% Get bit depth input
%nbits = str2double(get(handles.IN_nbits, 'string'));
%hMain = getappdata(0,'hMain');
% Check user's input
%if (isnan(nbits)||nbits<=0)
%    set(hObject,'String',num2str(handles.nbits))
%    warndlg('All inputs MUST be real positive numbers!');
%else
%    handles.nbits = nbits;
%    setappdata(hMain,'audio_recorder_nbits',nbits)
%end
%guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
%function IN_nbits_CreateFcn(hObject, ~, ~) %#ok : bit depth input box creation
% hObject    handle to IN_nbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end


% --- Executes when user attempts to close audio_recorder.
function audio_recorder_CloseRequestFcn(hObject, ~, ~) %#ok: window close request function
% hObject    handle to audio_recorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(hObject);


function IN_name_Callback(~, ~, ~) %#ok : recording name input box
% hObject    handle to IN_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_name as text
%        str2double(get(hObject,'String')) returns contents of IN_name as a double


% --- Executes during object creation, after setting all properties.
function IN_name_CreateFcn(hObject, ~, ~) %#ok : recording name input box creation
% hObject    handle to IN_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_enable.
function pb_enable_Callback(hObject, ~, handles) %#ok : playback enable checkbox function
% hObject    handle to pb_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pb_enable
hMain = getappdata(0,'hMain');
if get(hObject,'Value') == 1
    if isfield(handles.syscalstats,'latency'), set(handles.delay_chk,'Enable','on','Value',1); end
    if isfield(handles.syscalstats,'audio'), set(handles.invfilter_chk,'Enable','on','Value',1); end
    set(handles.output_panel,'Visible','on');
    set(handles.IN_numchsout,'String','1')
    if size(handles.outputdata.audio,2) == 1
        set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','on')
    else
        set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','off')
    end
    mainHandles = guidata(handles.main_stage1);
    selectednode = mainHandles.mytree.getSelectedNodes;
    set(handles.IN_name,'String',['rec_' selectednode(1).getName.char])
    set(handles.IN_numchs,'String',num2str(getappdata(hMain,'audio_recorder_numchs')));
    set(handles.text1,'String','Add time');
    set(handles.IN_duration,'String',num2str(getappdata(hMain,'audio_recorder_duration')));
    set(handles.IN_fs,'Enable','off');
    set(handles.IN_fs,'String','-');
    %set(handles.IN_nbits,'Enable','off');
    %set(handles.IN_nbits,'String','-');
    set(handles.IN_qdur,'String',num2str(getappdata(hMain,'audio_recorder_qdur')));
    set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
    set(handles.audio_recorder,'Position',handles.position);
    set(handles.OUT_axes,'Visible','on');
    children = get(handles.OUT_axes,'Children');
    set(children,'Visible','on');
    handles.numchs = str2num(get(handles.IN_numchs,'String'));
    handles.addtime = str2double(get(handles.IN_duration,'String'));
    handles.fs = handles.outputdata.fs;
    %handles.nbits = handles.outputdata.nbits;
    xlim(handles.IN_axes,[0 round(handles.dur+handles.addtime)])
    xlim(handles.OUT_axes,[0 round(handles.dur+handles.addtime)])
else
    if isfield(handles.syscalstats,'latency'), set(handles.delay_chk,'Enable','off','Value',0); end
    if isfield(handles.syscalstats,'audio'), set(handles.invfilter_chk,'Enable','off','Value',0); end
    set(handles.output_panel,'Visible','off');
    set(handles.IN_name,'String','recording')
    set(handles.text1,'String','Duration');
    set(handles.IN_duration,'String',num2str(getappdata(hMain,'audio_recorder_duration')));
    set(handles.IN_fs,'Enable','on');
    set(handles.IN_fs,'String',num2str(getappdata(hMain,'audio_recorder_fs')));
    %set(handles.IN_nbits,'Enable','on');
    %set(handles.IN_nbits,'String',num2str(getappdata(hMain,'audio_recorder_nbits')));
    set(handles.IN_qdur,'String',num2str(getappdata(hMain,'audio_recorder_qdur')));
    set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
    set(handles.audio_recorder,'Position',handles.position-[0 0 0 handles.OUT_axes_position(4)]);
    set(handles.OUT_axes,'Visible','off');
    children = get(handles.OUT_axes,'Children');
    set(children,'Visible','off');
    handles.numchs = str2num(get(handles.IN_numchs,'String'));
    handles.dur = str2double(get(handles.IN_duration,'String'));
    handles.fs = str2double(get(handles.IN_fs,'String'));
    %handles.nbits = str2double(get(handles.IN_nbits,'String'));
    xlim(handles.IN_axes,[0 round(handles.dur)])
    xlim(handles.OUT_axes,[0 round(handles.dur)])
end
guidata(hObject,handles);


% --- Executes on button press in stop_btn.
function stop_btn_Callback(hObject, ~, handles) %#ok : stop recording button
% hObject    handle to stop_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData.state = true;
set(hObject,'UserData',UserData);
pause on
pause(0.000001)
pause off
guidata(hObject,handles)


% --- Executes on button press in syscal_btn.
function syscal_btn_Callback(hObject, ~, handles) %#ok : system calibration button
% hObject    handle to syscal_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
recavail = 0;
if strcmp(get(handles.load_btn,'Enable'),'on'), recavail = 1; end
if recavail == 1, set([handles.load_btn,handles.preview_btn],'Enable','off'); end
set([handles.record_btn,handles.syscal_btn,handles.cancel_btn],'Enable','off')
syscalstats = syscal('audio_recorder', handles.audio_recorder);
if recavail == 1, set([handles.load_btn,handles.preview_btn],'Enable','on'); end
set([handles.record_btn,handles.syscal_btn,handles.cancel_btn],'Enable','on')
if isfield(syscalstats,'audio2'), handles.syscalstats(1).audio2 = syscalstats.audio2; end
if isfield(syscalstats,'fs'), handles.syscalstats(1).fs = syscalstats.fs; end
if isfield(syscalstats,'latency') && ~isnan(syscalstats.latency)
    if isfield(handles.syscalstats,'latency')
        if handles.syscalstats.latency ~= syscalstats.latency
             replace_value = questdlg(['The system latency measurement has changed from ' num2str(handles.syscalstats.latency) ' samples to ' num2str(syscalstats.latency) ' samples. Would you like to replace this value?'],...
                                       'AARAE info',...
                                       'Yes', 'No', 'No');
             switch replace_value
                 case 'Yes'
                     handles.savenewsyscal = 1;
                     handles.syscalstats.latency = syscalstats.latency;
             end
        end
    else
        handles.savenewsyscal = 1;
        handles.syscalstats(1).latency = syscalstats.latency;
    end
    set(handles.delay_chk,'Enable','on','Value',1)
    set(handles.delaytext,'String',[num2str(handles.syscalstats.latency) ' samples'])
end
if isfield(syscalstats,'cal')% && ~isnan(syscalstats.cal)
    if isfield(handles.syscalstats,'cal')
        if handles.syscalstats.cal ~= syscalstats.cal
             replace_value = questdlg(['The gain calibration measurement has changed from ' num2str(handles.syscalstats.cal) ' dB to ' num2str(syscalstats.cal) ' dB. Would you like to replace this value?'],...
                                       'AARAE info',...
                                       'Yes', 'No', 'No');
             switch replace_value
                 case 'Yes'
                     handles.savenewsyscal = 1;
                     handles.syscalstats.cal = syscalstats.cal;
             end
        end
    else
        handles.savenewsyscal = 1;
        handles.syscalstats(1).cal = syscalstats.cal;
    end
    set(handles.cal_chk,'Enable','on','Value',1)
    set(handles.caltext,'String',[num2str(handles.syscalstats.cal) ' dB'])
end
if isfield(syscalstats,'units')
    handles.syscalstats.units = syscalstats.units;
else
    handles.syscalstats.units = 'Val';
end
if isfield(syscalstats,'units_ref')
    handles.syscalstats.units_ref = syscalstats.units_ref;
else
    handles.syscalstats.units_ref = 1;
end
if isfield(syscalstats,'units')
    handles.syscalstats.units_type = syscalstats.units_type;
else
    handles.syscalstats.units_type = 1;
end
if isfield(syscalstats,'invfilter') && ~isempty(syscalstats.invfilter)
    if isfield(handles.syscalstats,'invfilter')
        if handles.syscalstats.invfilter ~= syscalstats.invfilter
             replace_value = questdlg('The inverse filter has been modified. Would you like to replace it?',...
                                       'AARAE info',...
                                       'Yes', 'No', 'No');
             switch replace_value
                 case 'Yes'
                     handles.savenewsyscal = 1;
                     handles.syscalstats.audio = syscalstats.invfilter;
             end
        end
    else
        handles.savenewsyscal = 1;
        handles.syscalstats(1).audio = syscalstats.invfilter;
    end
    set(handles.invfilter_chk,'Enable','on','Value',1)
    set(handles.invftext,'String','Available')
end
guidata(hObject,handles)


% --- Executes on button press in delay_chk.
function delay_chk_Callback(~, ~, ~) %#ok : latency checkbox
% hObject    handle to delay_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of delay_chk


% --- Executes on button press in cal_chk.
function cal_chk_Callback(~, ~, ~) %#ok : calibration checkbox
% hObject    handle to cal_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cal_chk


% --- Executes on button press in invfilter_chk.
function invfilter_chk_Callback(~, ~, ~) %#ok : inverse filter checkbox
% hObject    handle to invfilter_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invfilter_chk


% --- Executes on button press in preview_btn.
function preview_btn_Callback(~, ~, handles) %#ok : executes on preview button click
% hObject    handle to preview_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%hMain = getappdata(0,'hMain');
% Obtain handles using GUIDATA with the caller's handle

handles.testsignal = handles.rec;

% Warnings...
if isempty(handles.rec)
    warndlg('No signal recorded!');
    setappdata(hMain,'testsignal',[]);
else
%    hMain = getappdata(0,'hMain');
    if get(handles.delay_chk,'Value') == 1, handles.rec = [handles.rec(handles.syscalstats.latency:end,:);zeros(handles.syscalstats.latency-1,size(handles.rec,2))]; end
    if get(handles.invfilter_chk,'Value') == 1, handles.rec = filter(handles.syscalstats.audio2,1,handles.rec); end
    time1 = linspace(0,size(handles.rec,1)/handles.fs,length(handles.rec));
    time2 = linspace(0,size(handles.testsignal,1)/handles.fs,length(handles.testsignal));
    figure
    subplot(2,1,1); plot(time2,handles.testsignal); title('Recorded signal')
    subplot(2,1,2); plot(time1,handles.rec); title('System calibrated recording')
end



function IN_qdur_Callback(hObject, ~, handles) %#ok : queue duration input box
% hObject    handle to IN_qdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_qdur as text
%        str2double(get(hObject,'String')) returns contents of IN_qdur as a double
hMain = getappdata(0,'hMain');
% Get bit depth input
qdur = str2double(get(hObject, 'string'));
% Check user's input
if (isnan(qdur) || qdur<=0)
    set(hObject,'String',num2str(handles.qdur))
    warndlg('All inputs MUST be real positive numbers!');
else
    handles.qdur = qdur;
    setappdata(hMain,'audio_recorder_qdur',qdur)
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function IN_qdur_CreateFcn(hObject, ~, ~) %#ok : queue dueration input box creation
% hObject    handle to IN_qdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IN_buffer_Callback(hObject, ~, handles) %#ok : buffer size input box
% hObject    handle to IN_buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_buffer as text
%        str2double(get(hObject,'String')) returns contents of IN_buffer as a double
hMain = getappdata(0,'hMain');
buffer = str2double(get(hObject,'String'));
if (isnan(buffer) || buffer<=0)
    set(hObject,'String',num2str(handles.buffer))
    warndlg('All inputs MUST be real positive numbers!');
else
    handles.buffer = buffer;
    setappdata(hMain,'audio_recorder_buffer',buffer)
    if isfield(handles,'syscalstats') && isfield(handles.syscalstats,'latency')
        handles.syscalstats = rmfield(handles.syscalstats,'latency');
        handles.syscalstats = rmfield(handles.syscalstats,'audio');
        handles.syscalstats = rmfield(handles.syscalstats,'fs');
        set(handles.delay_chk,'Enable','off','Value',0)
        set(handles.delaytext,'String',[])
    end
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function IN_buffer_CreateFcn(hObject, ~, ~) %#ok : buffer size input box creation
% hObject    handle to IN_buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in inputdev_popup.
function inputdev_popup_Callback(hObject, ~, ~) %#ok : input device menu
% hObject    handle to inputdev_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns inputdev_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from inputdev_popup
hMain = getappdata(0,'hMain');
setappdata(hMain,'audio_recorder_input',get(hObject,'Value'));


% --- Executes during object creation, after setting all properties.
function inputdev_popup_CreateFcn(hObject, ~, ~) %#ok : input device menu creation
% hObject    handle to inputdev_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in outputdev_popup.
function outputdev_popup_Callback(hObject, ~, ~)  %#ok : output device menu
% hObject    handle to outputdev_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns outputdev_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from outputdev_popup
hMain = getappdata(0,'hMain');
setappdata(hMain,'audio_recorder_output',get(hObject,'Value'));

% --- Executes during object creation, after setting all properties.
function outputdev_popup_CreateFcn(hObject, ~, ~) %#ok : output device menu creation
% hObject    handle to outputdev_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IN_numchsout_Callback(~, ~, ~) %#ok : Executed when number of output channels changes
% hObject    handle to IN_numchsout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_numchsout as text
%        str2double(get(hObject,'String')) returns contents of IN_numchsout as a double


% --- Executes during object creation, after setting all properties.
function IN_numchsout_CreateFcn(hObject, ~, ~) %#ok : Creates number of output channels input box
% hObject    handle to IN_numchsout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sim_chk.
function sim_chk_Callback(~, ~, ~) %#ok : Executed when simultaneous playback checkbox changes
% hObject    handle to sim_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sim_chk


% --- Executes on button press in SilenceRequestCheckBox.
function SilenceRequestCheckBox_Callback(hObject, ~, ~)
% hObject    handle to SilenceRequestCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SilenceRequestCheckBox
hMain = getappdata(0,'hMain');
setappdata(hMain,'audio_recorder_silencerequest',get(hObject,'Value'));

    


% --- Executes during object creation, after setting all properties.
function SilenceRequestCheckBox_CreateFcn(~, ~, ~)
% hObject    handle to SilenceRequestCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
