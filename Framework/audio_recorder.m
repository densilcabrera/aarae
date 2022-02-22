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

% Last Modified by GUIDE v2.5 10-Mar-2021 09:46:20

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

% refresh audio device list on opening the audio recorder window
audiodevreset;

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
    %*****    inputdevinfo = dspAudioDeviceInfo('inputs');
    try % for windows machines look for ASIO device...
        aDW = audioDeviceWriter('Driver', 'ASIO');
        outputnames = getAudioDevices(aDW);
        if strcmp(outputnames,{'No audio output device detected'})
            delete(aDW)
            error('No ASIO Device Detected')
        else
        end
        aDR = audioDeviceReader('Driver', 'ASIO');
        inputnames = getAudioDevices(aDR);
        if strcmp(inputnames{1,1}, {'No audio input device detected'})
            delete(aDR)
            error('No ASIO Device Detected')
        else
        end
        release(aDW)
        release(aDR)
        delete(aDW)
        delete(aDR)
        devinfo = audiodevinfo;
        handles.winputnames = {devinfo.input.Name};
        handles.woutputnames = {devinfo.output.Name};
        inputnames = [inputnames,handles.winputnames];
        outputnames = [outputnames,handles.woutputnames];
        handles.ASIO =1;
        if getappdata(hMain, 'audio_recorder_ASIO') == 0 && handles.ASIO == 1
            h = helpdlg('ASIO device detected. Adjust buffer in ASIO device control panel to match user input');
            uiwait(h)
            setappdata(hMain,'audio_recorder_ASIO',1)
        end
    catch
        devinfo = audiodevinfo;
        inputnames = {devinfo.input.Name};
        outputnames = {devinfo.output.Name};
        handles.winputnames = {};
        handles.woutputnames = {};
        handles.ASIO=0;
    end
    inputnames = regexprep(inputnames,'\s\(Windows DirectSound\)','');
    inputnames = regexprep(inputnames,'\s\(ASIO\)','');
    inputnames = regexprep(inputnames,'\s\(Core Audio\)','');
    inputnames = regexprep(inputnames,'\s\(ALSA\)','');
    set(handles.inputdev_popup,'String',inputnames)
    %*****    outputdevinfo = dspAudioDeviceInfo('outputs');
    outputnames = regexprep(outputnames,'\s\(Windows DirectSound\)','');
    outputnames = regexprep(outputnames,'\s\(ASIO\)','');
    outputnames = regexprep(outputnames,'\s\(Core Audio\)','');
    outputnames = regexprep(outputnames,'\s\(ALSA\)','');
    set(handles.outputdev_popup,'String',outputnames)
    if ~isempty(handles.signaldata) && ismatrix(handles.signaldata.audio) && ~strcmp(handles.signaldata.datatype,'syscal')% If there's a signal loaded in the 'desktop'...
        % Allow visibility of playback option along with the specs of
        % the playback signal
        mainHandles = guidata(handles.main_stage1);
        selectednode = mainHandles.mytree.getSelectedNodes;
        set(handles.pb_enable,'Visible','on','Value',1);
        set(handles.nretries_gui,'Visible','on');
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
        
        if size(handles.outputdata.audio,2) == 1
            set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','on')
            set(handles.IN_numchsout,'String',num2str(getappdata(hMain,'audio_recorder_numchsout')));
        else
            set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','off')
            numchsout = (1:size(handles.outputdata.audio,2));
            numchsout = numchsout(~isnan(numchsout));
            numchsout = numchsout(numchsout>0);
            numchsout = round(numchsout);
            numchsout = unique(numchsout);
            set(handles.IN_numchsout,'String',num2str(numchsout));
            set(handles.sim_chk,'Value', 1)
            
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
        set(handles.recorder_instructions,'String','Instructions');
        set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
        handles.numchs = str2num(get(handles.IN_numchs,'String'));
        handles.addtime = str2double(get(handles.IN_duration,'String'));
        set(handles.SilenceRequestCheckBox,'Value',getappdata(hMain,'audio_recorder_silencerequest'));
        handles.silenceplease.audio = mainHandles.silenceplease.audio;
        handles.silenceplease.fs = mainHandles.silenceplease.fs;
        handles.thankyou.audio = mainHandles.thankyou.audio;
        handles.thankyou.fs = mainHandles.thankyou.fs;
        set(handles.Playback_delay,'String',num2str(getappdata(hMain,'audio_recorder_playbackdelay')));
        handles.fs = handles.outputdata.fs;
        %handles.nbits = handles.outputdata.nbits;
        %         handles.recorder_instructions = str2double(get(handles.recorder_instructions,'String'));
        handles.buffer = str2double(get(handles.IN_buffer,'String'));
        xlim(handles.IN_axes,[0 handles.dur+handles.addtime])
        xlim(handles.OUT_axes,[0 handles.dur+handles.addtime])
    else
        % If there's no signal loaded in the desktop just allocate memory
        % space for the signal to be recorded
        set(handles.pb_enable,'Visible','off','Value',0);
        set(handles.nretries_gui,'Visible','off');
        set(handles.retries_text,'Visible','off');
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
        %         set(handles.recorder_instructions,'String',num2str(getappdata(hMain,'audio_recorder_qdur')));
        set(handles.recorder_instructions,'String','Instructions!');
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
        handles.recorder_instructions = get(handles.recorder_instructions,'String');
        handles.buffer = str2double(get(handles.IN_buffer,'String'));
        xlim(handles.IN_axes,[0 handles.dur])
        xlim(handles.OUT_axes,[0 handles.dur])
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
set(handles.stack_btn,'Visible','Off')
set(handles.retry, 'Value', 0)% default to off. too many errors using directsound but fine for learning... otherwise use ASIO
% 
% % look for NI GPIB-USB device for controlling B&K turntable
% try
%     % Find a VISA-GPIB object.
%     handles.ttobj = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::10::INSTR', 'Tag', '');
%     
%     % Create the VISA-GPIB object if it does not exist
%     % otherwise use the object that was found.
%     if isempty(handles.ttobj)
%         handles.ttobj = visa('NI', 'GPIB0::10::INSTR');
%         fclose(handles.ttobj);
%         delete(handles.ttobj);
%         clear handles.ttobj;
%     else
%         fclose(handles.ttobj);
%         handles.ttobj = handles.ttobj(1);
%     end
%     set(handles.turn_chk, 'Visible', 'on')
%     set(handles.turn_cont_chk, 'Visible', 'on')
% catch
    set(handles.turn_chk, 'Visible', 'off')
    set(handles.turn_cont_chk, 'Visible', 'off')
% end
 
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

% UIWAIT makes audio_recorder wait for user response (see UIRESUME)x
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
 if strcmp(handles.stack_btn.Visible, 'on')
       msg = 'WARNING: Current recording not in stack. Continue without adding to stack?';
       answer = questdlg(msg, 'Record audio',...
           'Yes', 'No', 'No');
       if ~strcmp(answer,'Yes')
           return
       end
 end
% Call handles from main window
%mainHandles = guidata(handles.main_stage1);
set([handles.cancel_btn handles.load_btn handles.preview_btn handles.syscal_btn],'Enable','off')
set([handles.stack_btn,handles.stack_edit_btn], 'Visible', 'off','Enable','off');
% set ASIO based on input (or output - no situation where there is output
% only and windows only allows same device for Input/Output) device selection.
% matlab portaudio implementation of ASIO is irritating. why not let the
% device name change the driver... ?
    inputs = cellstr(get(handles.inputdev_popup,'String'));
    inputdevname = inputs{get(handles.inputdev_popup,'Value')};
    if ~ispc ||any(strcmp(inputdevname, regexprep(handles.winputnames,'\s\(Windows DirectSound)','')))
        handles.ASIO = 0; 
    else
        handles.ASIO = 1;
    end

% % Create Turntable object and necessary settings. DISABLED. Crappy
% implementation, better to do outside of aarae - will include examples of 
% how to do this in new folder "extras"
% if get(handles.turn_chk, 'Value') == 1
%     try
%         % Find a VISA-GPIB object.
%         handles.ttobj = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::10::INSTR', 'Tag', '');
%         
%         % Create the VISA-GPIB object if it does not exist
%         % otherwise use the object that was found.
%         if isempty(handles.ttobj)
%             handles.ttobj = visa('NI', 'GPIB0::10::INSTR');
%         else
%             fclose(handles.ttobj);
%             handles.ttobj = handles.ttobj(1);
%         end
%         
%         %open turntable object
%         fopen(handles.ttobj)
%             
%         if get(handles.turn_cont_chk, 'Value') == 0
%             % get turn values, these are absolute references, turn table
%             % setup is done prior to running the recorder.
%             ttAccel = get(handles.ttAccel, 'String');
%             fprintf(handles.ttobj, ['acc: ', ttAccel]);% TODO check this new line
%             % pause?
%             turn_abs =  str2num(get(handles.turn_abs, 'String'));
%             nturns = size(turn_abs,2) % turn_abs instead of turn_rel
%             turn_return = 0; % could be optional...
%             turn_return_cmd = ['turn_abs ', num2str(turn_return)];            
%             handles.turn_ID = [];
%             for i = 1:nturns
%                 turn_abs_cmd{i} = ['turn_abs ', num2str(turn_abs(:,i))];
%                 handles.turn_ID{i} = [num2str(turn_abs(:,i)), ' deg'];
%             end
%             % move turn table to first position
% fprintf(handles.ttobj, char(turn_abs_cmd(:,1)));
%                 pause(1) % JONO: EDIT THIS
%                 fprintf(handles.ttobj, 'start')
%                 pause(2)
%      else
%             nturns = 1;
%             turn_cont = ['cont ', get(handles.turn_abs, 'String')];
%             % %             % set 0
%             % %             fprintf(handles.ttobj, 'set: 0');
%             % set continous turn speed
%             fprintf(handles.ttobj, turn_cont)
%             fprintf(handles.ttobj, 'start')
%         end      
%     % Pause after setting up turn table.
% pause(.1);          
%     catch
%         warndlg('Unable to connect to GPIB device')
%         nturns = 1;
%     end
% else
    nturns = 1; % no turntable. Turntable disabled for now anyway.
% end


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
    
    % Set playback audio for play/record routine
    handles.hsr1 = dsp.SignalSource('SamplesPerFrame', handles.buffer); 
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
    handles.hsr1.Signal = thankyou;
    handles.hsr1.SamplesPerFrame = handles.buffer;
    %    setup(handles.hsr1);
    
    switch handles.ASIO
        case 0
            handles.hap = audioDeviceWriter(...
                'BitDepth', '24-bit integer',...
                'SampleRate',fs,...
                'BufferSize', handles.buffer,...
                'SupportVariableSizeInput', 0,...
                'ChannelMappingSource', 'Property',...
                'ChannelMapping', str2num(get(handles.IN_numchsout, 'String')),...
                'Device', outputdevname);
        case 1
            % consider adding a check that the asio buffer is the same as
            % the user buffer.
            handles.hap = audioDeviceWriter(...
                'Driver', 'ASIO',...
                'BitDepth', '24-bit integer',...
                'SampleRate',fs,...
                'BufferSize', handles.buffer,...
                'SupportVariableSizeInput', 0,...
                'ChannelMappingSource', 'Property',...
                'ChannelMapping', str2num(get(handles.IN_numchsout, 'String')),...
                'Device', outputdevname);
    end
    %     setup(handles.hap,...
    %         zeros(handles.hsr1.SamplesPerFrame,...
    %         length(str2num(get(handles.IN_numchsout, 'String')))));
    
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
            handles.hap(handles.hsr1());
        else
            break
        end
        %waitbar(i/ncycles,h)
    end
    %if i == ncycles, delete(h); end
    % Pause to wait for audioDeviceWriter
    pause(length(handles.hsr1.Signal)/fs)
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
    switch handles.ASIO
        case 0
            handles.hap = audioDeviceWriter(...
                'BitDepth', '24-bit integer',...
                'SampleRate',handles.outputdata.fs,...
                'BufferSize', handles.buffer,...`
                'SupportVariableSizeInput', 0,...
                'ChannelMappingSource', 'Property',...
                'ChannelMapping', str2num(get(handles.IN_numchsout, 'String')),...
                'Device', outputdevname);
        case 1
            handles.hap = audioDeviceWriter(...
                'Driver', 'ASIO',...
                'BitDepth', '24-bit integer',...
                'SampleRate',handles.outputdata.fs,...
                'BufferSize', handles.buffer,...
                'SupportVariableSizeInput', 0,...
                'ChannelMappingSource', 'Property',...
                'ChannelMapping', str2num(get(handles.IN_numchsout, 'String')),...
                'Device', outputdevname);
    end
    
    % Set record object
    inputs = cellstr(get(handles.inputdev_popup,'String'));
    inputdevname = inputs{get(handles.inputdev_popup,'Value')};
    switch handles.ASIO
        case 0
            handles.har = audioDeviceReader(...
                'Device',inputdevname,...
                'SampleRate',handles.outputdata.fs,...
                'BitDepth', '24-bit integer',...
                'OutputDataType','double',...
                'ChannelMappingSource','Property',...
                'ChannelMapping',handles.numchs,...
                'SamplesPerFrame',handles.buffer);
        case 1
            handles.har = audioDeviceReader(...
                'Driver', 'ASIO',...
                'Device',inputdevname,...
                'SampleRate',handles.outputdata.fs,...
                'BitDepth', '24-bit integer',...
                'OutputDataType','double',...
                'ChannelMappingSource','Property',...
                'ChannelMapping',handles.numchs,...
                'SamplesPerFrame',handles.buffer);
    end
    %     setup(handles.har)
    
    % Set playback audio
    handles.hsr1 = dsp.SignalSource('SamplesPerFrame', handles.buffer);
    if size(handles.outputdata.audio,2) == 1 && length(str2num(get(handles.IN_numchsout,'String'))) > 1
        if isfield(handles.syscalstats, 'latency') && get(handles.delay_chk,'Value') == 1
            handles.toremove = zeros(handles.syscalstats.latency,1);
        else
            handles.toremove = [];
        end
        if get(handles.sim_chk,'Value') == 1
            playbackaudio = repmat([handles.outputdata.audio;handles.toremove],1,length(str2num(get(handles.IN_numchsout,'String'))));
            %             addtimeok = false;
        else
            playbackaudio = [handles.outputdata.audio;handles.toremove];
        end
    else
        if isfield(handles.syscalstats, 'latency') && get(handles.delay_chk,'Value') == 1
            handles.toremove = zeros(handles.syscalstats.latency, size(handles.outputdata.audio,2));
        else
            handles.toremove = [];
        end
        playbackaudio = [handles.outputdata.audio;handles.toremove];
    end
    playbackdelay = str2num(get(handles.Playback_delay,'String'));
    if get(handles.turn_chk, 'Value') == 1
        try
            % ttdelay = turn_abs/30;
            ttdelay = 2;
        catch
            ttdelay = 2;
        end
    else
        ttdelay = 0;
    end
    %if get(handles.invfilter_chk,'Value') == 1, playbackaudio = filter(handles.syscalstats.audio2,1,playbackaudio); end
    % addtime ok recorder_instructions
    handles.hsr1.Signal = [zeros(floor(ttdelay*handles.fs),size(playbackaudio,2));...
        zeros(floor(playbackdelay*handles.fs),size(playbackaudio,2));...
        playbackaudio;...
        zeros(floor((handles.addtime)*handles.fs),size(playbackaudio,2))];
    handles.hsr1.SamplesPerFrame = handles.har.SamplesPerFrame;
    setup(handles.hsr1)
    guidata(hObject,handles)
    handles.rec = [];
    ncycles = ceil(length(handles.hsr1.Signal)/handles.har.SamplesPerFrame);
    if get(handles.sim_chk,'Value') == 0
        if length(str2num(get(handles.IN_numchsout,'String'))) >1
            audio = zeros(ncycles*handles.har.SamplesPerFrame,length(handles.numchs),1,nturns,length(str2num(get(handles.IN_numchsout,'String'))));
        else
            audio = zeros(ncycles*handles.har.SamplesPerFrame,length(handles.numchs),1,nturns);
        end
    else
        audio = zeros(ncycles*handles.har.SamplesPerFrame,length(handles.numchs),1,nturns);
    end
    set(hObject,'BackgroundColor','red');
    set(handles.stop_btn,'Visible','on');
    
    
    % Initialize playback/record routine
    pause on
    try
        UserData = get(handles.stop_btn,'UserData');
        h = waitbar(0,'Recording in progress...','Name','AARAE info');
        handles.recordtime = datestr(now);
        OutCh = str2num(get(handles.IN_numchsout,'String'));
        totalUnderrun = 0;
        numUnderrun = zeros(ncycles,length(OutCh),nturns);
        numOverrun = zeros(ncycles,length(OutCh),nturns);
        totalOverrun=0;
        if get(handles.retry, 'Value') == 1
            nretries = str2num(get(handles.nretries_gui,'String'));
            if isempty(nretries)
                nretries = 0;
            end
            switch get(handles.sim_chk,'Value')
                case 0
                    dummy = zeros(handles.hsr1.SamplesPerFrame,1);
                    for tt = 1:nturns
                        for outind = 1:length(OutCh)
%                             tic
                            release(handles.hap)
                            handles.hap.ChannelMapping = OutCh(outind);
                            setup(handles.hap,dummy)
                            [~] = handles.har();
                            handles.hap(dummy);
%                             pause(1) % new line 06/06/17
%                         pause(.5)
                            for i = 1:ncycles
                                UserData = get(handles.stop_btn,'UserData');
                                if UserData.state == false
                                    [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:,1,tt,outind),numOverrun(i,outind,tt)] = handles.har();
                                    numUnderrun(i,outind, tt) = handles.hap(handles.hsr1());
                                    if numUnderrun(i,outind, tt) + numOverrun(i, outind, tt) >0 && nretries >0
                                        ermessage = ['Error - Output channel: ', num2str(OutCh(outind)), '... Trying again in ', num2str(handles.addtime), 's'];
                                        waitbar(((outind-1)*ncycles+i)/(length(OutCh)*ncycles),h, ermessage)
                                        pause(handles.addtime)
                                        release(handles.hsr1)
                                        reset(handles.hsr1)
                                        audio(:,:,:,tt,outind) = 0;
                                        break
                                    end
                                    totalUnderrun = totalUnderrun + numUnderrun(i,outind, tt);
                                    totalOverrun= totalOverrun+ numOverrun(i,outind, tt);
                                else
                                    break
                                end
                                waitbar(((outind-1)*ncycles*tt+i)/(length(OutCh)*ncycles*nturns),h,'Recording in progress...')
                            end
                            if nretries > 0 && numUnderrun(i,outind,tt) + numOverrun(i, outind,tt) >0
                                for j = 1:nretries
                                    ermessage = ['Retrying Output channel: ', num2str(OutCh(outind)), '... Retry: ', num2str(j)];
                                    waitbar(((outind-1)*ncycles+i)/(length(OutCh)*ncycles),h, ermessage)
                                    for i = 1:ncycles
                                        UserData = get(handles.stop_btn,'UserData');
                                        if UserData.state == false
                                            [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:,1,tt,outind),numOverrun(i,outind)] = handles.har();
                                            numUnderrun(i,outind,tt) = handles.hap(handles.hsr1());
                                            if numUnderrun(i,outind, tt) + numOverrun(i, outind, tt) >0 && nretries >0
                                                ermessage = ['Error - Output channel: ', num2str(OutCh(outind)), '... Retry: ', num2str(j), '... Trying again in ', num2str(handles.addtime), 's'];
                                                waitbar(((outind-1)*ncycles+i)/(length(OutCh)*ncycles),h, ermessage)
                                                pause(handles.addtime)
                                                release(handles.hsr1)
                                                reset(handles.hsr1)
                                                audio(:,:,:,tt,outind) = 0;
                                                break
                                            end
                                            totalUnderrun = totalUnderrun + numUnderrun(i,outind, tt);
                                            totalOverrun= totalOverrun+ numOverrun(i,outind, tt);
                                        else
                                            break
                                        end
                                        waitbar(((outind-1)*ncycles*tt+i)/(length(OutCh)*ncycles*nturns),h,'Recording in progress...')
                                    end
                                    if isDone(handles.hsr1) && numUnderrun(i,outind, tt) + numOverrun(i, outind) == 0
                                        handles.message{outind} = {'Error detected.'; 'Successful Routine.'; ['Outchan: ', num2str(OutCh(outind))]; ['retry #: ', num2str(j)]};
                                        % disp(handles.message{outind})
                                        break
                                    else
                                        handles.message{outind} = {'Error detected.'; 'Unuccessful Routine.'; ['Outchan: ', num2str(OutCh(outind))]; ['retry #: ', num2str(j)] };
                                        %disp(handles.message{outind})
                                        if j == nretries
                                            totalUnderrun = totalUnderrun + numUnderrun(i,outind, tt);
                                            totalOverrun= totalOverrun+ numOverrun(i,outind, tt);
                                        end
                                    end
                                    release(handles.hsr1)
                                    reset(handles.hsr1)
                                end
                            end
                            release(handles.hsr1)
                            reset(handles.hsr1)
%                             toc
                        end
                        if get(handles.turn_chk, 'Value') ==1 && tt < nturns
                            fprintf(handles.ttobj, char(turn_abs_cmd(:,tt+1))); %JONO: 
                            pause(1) % JONO: EDIT THIS
                            fprintf(handles.ttobj, 'start')
                            pause(2)
                            % pause(2) % removed because it causes
                            % playback/record error/latency - need to add
                            % time to playback/record signal then remove
                        end
                    end
                case 1
                    dummy = zeros(handles.hsr1.SamplesPerFrame, length(str2num(get(handles.IN_numchsout, 'String'))));
                    setup(handles.hap,dummy);
                    [~] = handles.har();
                    handles.hap(dummy);
                    for i = 1:ncycles
                        UserData = get(handles.stop_btn,'UserData');
                        if UserData.state == false
                            [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:),numOverrun(i)] = handles.har();
                            numUnderrun(i) = handles.hap(handles.hsr1());
                            if numUnderrun(i) + numOverrun(i) >0 && nretries >0
                                ermessage = ['Error... Trying again in ', num2str(handles.addtime), 's'];
                                waitbar((i/ncycles),h, ermessage)
                                pause(handles.addtime)
                                release(handles.hsr1)
                                reset(handles.hsr1)
                                audio(:,:) = 0;
                                break
                            end
                            totalUnderrun = totalUnderrun + numUnderrun(i);
                            totalOverrun= totalOverrun + numOverrun(i);
                        else
                            break
                        end
                        waitbar(i/ncycles,h)
                    end
                    if nretries > 0 && numUnderrun(i) + numOverrun(i) >0
                        for j = 1:nretries
                            for i = 1:ncycles
                                UserData = get(handles.stop_btn,'UserData');
                                if UserData.state == false
                                    [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:),numOverrun(i)] = handles.har();
                                    numUnderrun(i) = handles.hap(handles.hsr1());
                                    if numUnderrun(i) + numOverrun(i) >0 && nretries >0
                                        ermessage = ['Error... Retry: ', num2str(j), '... Trying again in ', num2str(handles.addtime), 's'];
                                        waitbar((i/ncycles),h, ermessage)
                                        pause(handles.addtime)
                                        release(handles.hsr1)
                                        reset(handles.hsr1)
                                        audio(:,:) = 0;
                                        break
                                    end
                                    totalUnderrun = totalUnderrun + numUnderrun(i);
                                    totalOverrun= totalOverrun+ numOverrun(i);
                                else
                                    break
                                end
                                waitbar(i/ncycles,h)
                            end
                            if isDone(handles.hsr1) && numUnderrun(i) + numOverrun(i) == 0
                                handles.message{1} = {'Error detected.'; 'Successful Routine.'; 'Outchan:  Synchronous'; ['retry #: ', num2str(j)]};
                                %disp(handles.message{1})
                                release(handles.hsr1)
                                reset(handles.hsr1)
                                break
                            else
                                handles.message{1} = {'Error detected.'; 'Unuccessful Routine.'; 'Outchan:  Synchronous'; ['retry #: ', num2str(j)] };
                                %disp(handles.message{1})
                                if j == nretries
                                    totalUnderrun = totalUnderrun + numUnderrun(i);
                                    totalOverrun= totalOverrun+ numOverrun(i);
                                end
                            end
                            release(handles.hsr1)
                            reset(handles.hsr1)
                        end
                    end
            end
        else
            switch get(handles.sim_chk,'Value')
                case 0
                    dummy = zeros(handles.hsr1.SamplesPerFrame,1);
                    for outind = 1:length(OutCh)
                        release(handles.hap)
                        handles.hap.ChannelMapping = OutCh(outind);
                        setup(handles.hap,dummy);
                        [~] = handles.har();
                        handles.hap(dummy);
                        for i = 1:ncycles
                            UserData = get(handles.stop_btn,'UserData');
                            if UserData.state == false
                                [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:,1,1,outind),numOverrun(i,outind)] = handles.har();
                                numUnderrun(i,outind) = handles.hap(handles.hsr1());
                                totalUnderrun = totalUnderrun + numUnderrun(i,outind);
                                totalOverrun= totalOverrun+ numOverrun(i,outind);
                            else
                                break
                            end
                            waitbar(((outind-1)*ncycles+i)/(length(OutCh)*ncycles),h)
                        end
                        release(handles.hsr1)
                        reset(handles.hsr1)
                    end
                case 1
                    dummy = zeros(handles.hsr1.SamplesPerFrame,str2num(get(handles.IN_numchsout, 'String')));
                    setup(handles.hap,dummy);
                    [~] = handles.har();
                    handles.hap(dummy);
                    for i = 1:ncycles
                        UserData = get(handles.stop_btn,'UserData');
                        if UserData.state == false
                            [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:),numOverrun(i)] = handles.har();
                            numUnderrun(i) = handles.hap(handles.hsr1());
                            totalUnderrun = totalUnderrun + numUnderrun(i);
                            totalOverrun= totalOverrun + numOverrun(i);
                        else
                            break
                        end
                        waitbar(i/ncycles,h)
                    end
            end
        end
        if get(handles.turn_chk, 'Value') == 1
            if get(handles.turn_cont_chk, 'Value') == 1
                fprintf(handles.ttobj, 'stop')
            else
                fprintf(handles.ttobj, turn_return_cmd)
                pause(2) % JONO: EDIT THIS
                fprintf(handles.ttobj, 'start')
                pause(1)
                hh = helpdlg('turntable returning to zero');
                uiwait(hh)
            end
            fclose(handles.ttobj);
            delete(handles.ttobj);
            clear handles.ttobj;
        end
        if i == ncycles, delete(h); end
    catch sthgwrong
        UserData.state = true;
        handles.rec = [];
        syswarning = sthgwrong.message;
        warndlg(syswarning,'AARAE info')
    end
    pause off
    % Check recording and adjust for QueueDuration latency, account for
    % buffer etc.
    handles.rec = audio(floor(ttdelay*handles.fs)+length(handles.toremove)+1:end,:,:,:,:);
    if ~isempty(handles.rec)
        %         handles.rec = handles.rec(handles.qdur*handles.fs:end,:);
        if UserData.state == false
            %                 handles.rec = handles.rec(1:(size(handles.outputdata.audio,1)+handles.addtime*handles.fs),:);
        else
            UserData.state = false;
            set(handles.stop_btn,'UserData',UserData);
        end
        % check errors
        if get(handles.retry, 'Value') == 1 && totalUnderrun || totalOverrun >0
            warndlg('Record/playback error: Try again or load signal and see error field');
            handles.Underrun = numUnderrun;
            handles.Overrun = numOverrun;
            handles.error = 1;
            % % % % %             f = figure('Name','Recording Errors', ...
            % % % % %                 'Position',[200 200 620 360]);
            % % % % %             inderror = 1:ncycles;
            % % % % %             daterror = [numOverrun,numUnderrun];
            % % % % %
            % % % % %             for c = 1:length(OutCh)
            % % % % %                 cnames{:,2*c-1} = ['Overrun Out Channel: ',num2str(OutCh(c))];
            % % % % %                 cnames{:,2*c} = ['Underrun Out Channel: ',num2str(OutCh(c))];
            % % % % %             end
            % % % % %
            % % % % %             rnames = num2str(inderror);
            % % % % %             t2 =uitable('Data',daterror,'ColumnName',cnames, 'RowName', rnames);
            % % % % %             set(t2,'ColumnWidth',{90});
            % % % % %             [~,tables] = disptables(f,[t2],{'Error'});
        else
            handles.error = 0;
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
    switch handles.ASIO
        case 0
            handles.har = audioDeviceReader(...
                'Device',inputdevname,...
                'SampleRate',handles.fs,...
                'BitDepth', '24-bit integer',...
                'OutputDataType','double',...
                'ChannelMappingSource','Property',...
                'ChannelMapping',handles.numchs,...
                'SamplesPerFrame',handles.buffer);
        case 1
            handles.har = audioDeviceReader(...
                'Driver', 'ASIO',...
                'Device',inputdevname,...
                'SampleRate',handles.fs,...
                'BitDepth', '24-bit integer',...
                'OutputDataType','double',...
                'ChannelMappingSource','Property',...
                'ChannelMapping',handles.numchs,...
                'SamplesPerFrame',handles.buffer);
    end
    %     setup(handles.har)
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
        numOverrun=[];
        totalOverrun=0;
        if get(handles.retry, 'Value') == 1
            for i = 1:ncycles
                UserData = get(handles.stop_btn,'UserData');
                if UserData.state == false
                    [audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:),numOverrun(i)] = handles.har();
                    totalOverrun = totalOverrun + numOverrun(i);
                else
                    break
                end
                waitbar(i/ncycles,h)
            end
        else
            for i = 1:ncycles
                UserData = get(handles.stop_btn,'UserData');
                if UserData.state == false
                    audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:) = handles.har();
                else
                    break
                end
                waitbar(i/ncycles,h)
            end
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
        % check errors
        if get(handles.retry, 'Value') == 1 && totalOverrun >0
            warndlg('Record/playback error: Try again or load signal and see error field');
            handles.Underrun = [];
            handles.Overrun = numOverrun;
            handles.error = 1;
            % % % % %             f = figure('Name','Recording Errors', ...
            % % % % %                 'Position',[200 200 620 360]);
            % % % % %             inderror = 1:ncycles;
            % % % % %             daterror = [numOverrun,numUnderrun];
            % % % % %
            % % % % %             for c = 1:length(OutCh)
            % % % % %                 cnames{:,2*c-1} = ['Overrun Out Channel: ',num2str(OutCh(c))];
            % % % % %                 cnames{:,2*c} = ['Underrun Out Channel: ',num2str(OutCh(c))];
            % % % % %             end
            % % % % %
            % % % % %             rnames = num2str(inderror);
            % % % % %             t2 =uitable('Data',daterror,'ColumnName',cnames, 'RowName', rnames);
            % % % % %             set(t2,'ColumnWidth',{90});
            % % % % %             [~,tables] = disptables(f,[t2],{'Error'});
        else
            handles.error = 0;
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
    handles.hsr1 = dsp.SignalSource('SamplesPerFrame', handles.buffer);
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
    
    
    %if get(handles.invfilter_chk,'Value') == 1, quietplease = filter(handles.syscalstats.audio2,1,quietplease); end
    handles.hsr1.Signal = [thankyou];
    handles.hsr1.SamplesPerFrame = handles.buffer;
    %     setup(handles.hsr1);
    
    switch handles.ASIO
        case 0
            handles.hap = audioDeviceWriter(...
                'BitDepth', '24-bit integer',...
                'SampleRate',fs,...
                'BufferSize', handles.buffer,...
                'SupportVariableSizeInput', 0,...
                'ChannelMappingSource', 'Property',...
                'ChannelMapping', str2num(get(handles.IN_numchsout, 'String')),...
                'Device', outputdevname);
        case 1
            % consider adding a check that the asio buffer is the same as
            % the user buffer.
            handles.hap = audioDeviceWriter(...
                'Driver', 'ASIO',...
                'BitDepth', '24-bit integer',...
                'SampleRate',fs,...
                'BufferSize', handles.buffer,...
                'SupportVariableSizeInput', 0,...
                'ChannelMappingSource', 'Property',...
                'ChannelMapping', str2num(get(handles.IN_numchsout, 'String')),...
                'Device', outputdevname);
    end
    %     setup(handles.hap,...
    %         zeros(handles.hsr1.SamplesPerFrame,...
    %         length(str2num(get(handles.IN_numchsout, 'String')))));
    
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
            handles.hap(handles.hsr1());
        else
            break
        end
    end
    if i == ncycles, delete(h); end
    pause(length(handles.hsr1.Signal)/fs)
    pause off
    % Release playback and audio data objects
    release(handles.hap)
    release(handles.hsr1)
end

try
    delete(handles.hsr1)
    delete(handles.hap)
    delete(handles.har)
catch
    delete(handles.har)
end


set(handles.record_btn,'BackgroundColor',[0.94 0.94 0.94]);
set(handles.record_btn,'Enable','on');
set(handles.stop_btn,'Visible','off');
set([handles.load_btn handles.preview_btn handles.cancel_btn handles.syscal_btn],'Enable','on')
if any(audio)
    set(handles.stack_btn, 'Visible', 'on','Enable','on');
%     handles.rec_size_string = ['rec size   = ', num2str(size(audio))];
    handles.rec_size_string = ['rec size:   [',...
    num2str(size(handles.rec,1),'%10.3e\n'),...
    '] x [', num2str(size(handles.rec,2)),...
    '] x [', num2str(size(handles.rec,3)),...
    '] x [', num2str(size(handles.rec,4)),...
    '] x [', num2str(size(handles.rec,5)),...
    '] x [', num2str(size(handles.rec,6)),']'];
set(handles.rec_size_text,'string',handles.rec_size_string, 'Enable','on');

end
guidata(hObject,handles);


% --- Executes on button press in load_btn.
function load_btn_Callback(hObject, ~, handles) %#ok : called when sending recording to main window
% hObject    handle to load_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMain = getappdata(0,'hMain');
% Obtain handles using GUIDATA with the caller's handle
if isfield(handles,'stack')
   if strcmp(handles.stack_btn.Visible, 'on')
       msg = 'WARNING: Current recording not in stack. Continue without adding to stack?';
       answer = questdlg(msg, 'Load audio',...
           'Yes', 'No', 'No');
       if ~strcmp(answer,'Yes')
           return
       end
   end
    handles.rec = handles.stack; % JONO added this, TODO change other handles as necessary, i.e. record time etc
end
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
    if get(handles.delay_chk,'Value') == 1 && isempty(handles.toremove), handles.rec = [handles.rec(handles.syscalstats.latency:end,:);zeros(handles.syscalstats.latency,size(handles.rec,2))]; end
    if get(handles.invfilter_chk,'Value') == 1, handles.rec = filter(handles.syscalstats.audio2,1,handles.rec); end
    handles.recording.audio = handles.rec;
    
    if get(handles.retry, 'Value') == 1 && handles.error == 1
        handles.recording.error.Underrun = handles.Underrun;
        handles.recording.error.Overrun = handles.Overrun;
    end
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
        handles.recording.chanID = cellstr([repmat('InChan',size(handles.recording.audio,2),1) num2str((handles.numchs)')]);
        
    catch
        handles.recording.chanID = cellstr([repmat('InChan',size(handles.recording.audio,2),1) num2str((1:size(handles.recording.audio,2))')]);
    end
    % Potentially add dim5ID or outchanID here (once its name is defined).
    % should it be a property or a first order field?
    % ***** outchanID in dim5 added
    if get(handles.sim_chk,'Value') == 0
        OutCh = str2num(get(handles.IN_numchsout,'String'));
        try
            handles.recording.OutchanID = cellstr([repmat('OutCh',size(handles.recording.audio,5),1) num2str((OutCh)')]);
            
        catch
            handles.recording.OutchanID = cellstr([repmat('OutCh',size(handles.recording.audio,5),1) num2str((1:size(handles.recording.audio,5))')]);
        end
    end
    if get(handles.turn_chk,'Value') ==1
        handles.recording.turnID = handles.turn_ID;
    end
    if get(handles.cal_chk,'Value') == 1
        %%%% handles.recording.cal = handles.syscalstats.cal(1:size(handles.recording.audio,2));
        handles.recording.cal = handles.syscalstats.cal(handles.numchs);
        if isfield(handles.syscalstats, 'units')
            handles.recording.properties.units = handles.syscalstats.units; %(handles.numchs);
            handles.recording.properties.units_ref = handles.syscalstats.units_ref; %(handles.numchs);
            handles.recording.properties.units_type = handles.syscalstats.units_type; %(handles.numchs);
        end
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
        playbackdelay = str2double(get(handles.Playback_delay,'String'));
        if get(handles.sim_chk,'Value')
            synchplayback = 'simultaneous';
        else
            synchplayback = 'sequential';
        end;
        outputdevname = outputs{get(handles.outputdev_popup,'Value')};
        outchans = get(handles.IN_numchsout,'String');
    else
        handles.recording.history{1,1} = handles.recordtime;
        handles.recording.history{1,2} = 'Recorded without playback';
        handles.recording.history{1,3} = 'AARAE audio_recorder';
        inputs = cellstr(get(handles.inputdev_popup,'String'));
        inputdevname = inputs{get(handles.inputdev_popup,'Value')};
        playbackdelay = '-';
        synchplayback = '-';
        outputdevname = 'none';
        outchans = '-';
    end
    handles.recording.history{2,1} = handles.recordtime;
    handles.recording.history{2,2} = 'Recorded';
    handles.recording.history{2,3} = 'AARAE audio_recorder';
    handles.recording.history{2,4} = {'fs', handles.fs;...
        'duration', handles.dur;...
        'input device', inputdevname;...
        'input channels', handles.numchs;...
        'buffer length', handles.buffer;...
        'output device', outputdevname;...
        'output channels', outchans;...
        'multi-channel output mode', synchplayback;...
        'playback delay', playbackdelay;...
        'silence request', get(handles.SilenceRequestCheckBox,'Value')};
    if isfield(handles,'comment')
        handles.recording.history = [handles.recording.history; handles.comment];
    end
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
if get(handles.pb_enable,'Value') == 1
    dur = str2double(get(handles.IN_duration, 'string'));
    hMain = getappdata(0,'hMain');
else
    dur = round(str2double(get(handles.IN_duration, 'string')));
    hMain = getappdata(0,'hMain');
end
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
        xlim(handles.IN_axes,[0 length(handles.outputdata.audio)/...
            handles.outputdata.fs + handles.addtime])
        xlim(handles.OUT_axes,[0 length(handles.outputdata.audio)/...
            handles.outputdata.fs + handles.addtime])
    else
        xlim(handles.IN_axes,[0 handles.dur])
        xlim(handles.OUT_axes,[0 handles.dur])
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
    %     set(handles.IN_numchsout,'String','1')
    if size(handles.outputdata.audio,2) == 1
        set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','on')
    else
        set([handles.text25,handles.IN_numchsout,handles.sim_chk],'Visible','off')
        set(handles.sim_chk,'Value', 1)
    end
    mainHandles = guidata(handles.main_stage1);
    selectednode = mainHandles.mytree.getSelectedNodes;
    set(handles.IN_name,'String',['rec_' selectednode(1).getName.char])
    set(handles.IN_numchs,'String',num2str(getappdata(hMain,'audio_recorder_numchs')));
    set(handles.IN_numchsout,'String',num2str(getappdata(hMain,'audio_recorder_numchsout')));
    set(handles.text1,'String','Add time');
    set(handles.IN_duration,'String',num2str(getappdata(hMain,'audio_recorder_duration')));
    set(handles.IN_fs,'Enable','off');
    set(handles.IN_fs,'String','-');
    %set(handles.IN_nbits,'Enable','off');
    %set(handles.IN_nbits,'String','-');
    set(handles.recorder_instructions,'String','Instructions');
    set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
    set(handles.audio_recorder,'Position',handles.position);
    set(handles.OUT_axes,'Visible','on');
    children = get(handles.OUT_axes,'Children');
    set(children,'Visible','on');
    handles.numchs = str2num(get(handles.IN_numchs,'String'));
    handles.addtime = str2double(get(handles.IN_duration,'String'));
    handles.fs = handles.outputdata.fs;
    %handles.nbits = handles.outputdata.nbits;
    xlim(handles.IN_axes,[0 handles.dur+handles.addtime])
    xlim(handles.OUT_axes,[0 handles.dur+handles.addtime])
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
    set(handles.recorder_instructions,'String','Instructions');
    set(handles.IN_buffer,'String',num2str(getappdata(hMain,'audio_recorder_buffer')));
    set(handles.audio_recorder,'Position',handles.position-[0 0 0 handles.OUT_axes_position(4)]);
    set(handles.OUT_axes,'Visible','off');
    children = get(handles.OUT_axes,'Children');
    set(children,'Visible','off');
    handles.numchs = str2num(get(handles.IN_numchs,'String'));
    handles.dur = str2double(get(handles.IN_duration,'String'));
    handles.fs = str2double(get(handles.IN_fs,'String'));
    %handles.nbits = str2double(get(handles.IN_nbits,'String'));
    xlim(handles.IN_axes,[0 handles.dur])
    xlim(handles.OUT_axes,[0 handles.dur])
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
            end
        end
    else
        handles.savenewsyscal = 1;
        handles.syscalstats(1).cal = syscalstats.cal;
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
    end
    set(handles.cal_chk,'Enable','on','Value',1)
    set(handles.caltext,'String',[num2str(handles.syscalstats.cal) ' dB'])
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



function recorder_instructions_Callback(hObject, ~, handles) %#ok : queue duration input box
% hObject    handle to recorder_instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recorder_instructions as text
%        str2double(get(hObject,'String')) returns contents of recorder_instructions as a double
message{1} = 'Instructions:';
message{2} = '1. Select Audio Device from drop down menu. (Note: stabilty is not guarenteed when using Windows directsound drivers or built in in/out on Mac OS)';
message{3} = '2. When using ASIO devices on windows machines, you must manually set the buffersize on the device control software as well as the text box in audio recorder';
message{4} = '3. You can test for sample drops by running the routine with error check box ticked. "Retry x N" is the number of times to retry along output channels. If an error is detected during recording and retries >0, the play/record routine is aborted and a second attempt is run. An error dialogue box is shown if there are errors remaining.';
message{5} = '4. Select Audio Device and set device settings in audio recorder BEFORE calibration';
message{6} = '5. If Delay Comp is ticked when recording, latency is removed and the compensated audio is plotted. If Delay Comp is ticked after recording, latency is removed when "Load Recording" is pressed. Latency is only removed once';

helpdlg(message,'AARAE info')

% --- Executes during object creation, after setting all properties.
function recorder_instructions_CreateFcn(hObject, ~, ~) %#ok : queue dueration input box creation
% hObject    handle to recorder_instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
% %     set(hObject,'BackgroundColor','white');




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



function IN_numchsout_Callback(hObject, ~, handles) %#ok : Executed when number of output channels changes
% hObject    handle to IN_numchsout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_numchsout as text
%        str2double(get(hObject,'String')) returns contents of IN_numchsout as a double
hMain = getappdata(0,'hMain');
numchsout = str2num(get(handles.IN_numchsout, 'String'));
hMain = getappdata(0,'hMain');
% Check user's input
numchsout = numchsout(~isnan(numchsout));
numchsout = numchsout(numchsout>0);
numchsout = round(numchsout);
numchsout = unique(numchsout);
if isempty(numchsout)
    set(hObject,'String',num2str(handles.numchsout));
    warndlg('All inputs MUST be real positive numbers!');
else
    handles.numchsout = numchsout;
    setappdata(hMain,'audio_recorder_numchsout',numchsout)
end
guidata(hObject, handles);


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



function Playback_delay_Callback(hObject, ~, handles)
% hObject    handle to Playback_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Playback_delay as text
%        str2double(get(hObject,'String')) returns contents of Playback_delay as a double
hMain = getappdata(0,'hMain');
playbackdelay = str2double(get(hObject,'String'));
if (isnan(playbackdelay) || playbackdelay<0)
    set(hObject,'String',num2str(handles.playbackdelay))
    warndlg('All inputs MUST be real numbers >=0!');
else
    handles.playbackdelay = playbackdelay;
    setappdata(hMain,'playbackdelay',playbackdelay)
end
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function Playback_delay_CreateFcn(hObject, ~, ~)
% hObject    handle to Playback_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LogTextfromRecorder.
function LogTextfromRecorder_Callback(hObject, ~, handles)
% hObject    handle to LogTextfromRecorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = inputdlg('Add comment to the log file and the selection''s history field(s):','Log comment',[15 150]);
% write to log file
if ~isempty(answer)
    answer = char(answer);
    logtext('%% **************************************************\n');
    logtext(['%% User comment ' datestr(now,16) ' \n']);
    logtext('%% \n');
    for n = 1:size(answer,1)
        logtext(['%% ' answer(n,:) '\n']);
    end
    logtext('%% \n');
    logtext('%% **************************************************\n');
end

row = cell(1,4);
row{1,1} = datestr(now);
row{1,2} = 'COMMENT';
row{1,4} = answer;
if isfield(handles,'comment')
    handles.comment = [handles.comment;row];
else
    handles.comment = row;
end

guidata(hObject,handles)


% --- Executes on button press in retry.
function retry_Callback(hObject, ~, handles)
% hObject    handle to retry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.retry = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of retry


% % % % % % % % % % --- Executes on button press in recorder_instructions.
% % % % % % % % % function recorder_instructions_Callback(hObject, eventdata, handles)
% % % % % % % % % % hObject    handle to recorder_instructions (see GCBO)
% % % % % % % % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % % % % % % % handles    structure with handles and user data (see GUIDATA)
% % % % % % % % %
% % % % % % % % %
% % % % % % % % % % --- Executes during object creation, after setting all properties.
% % % % % % % % % function recorder_instructions_CreateFcn(hObject, eventdata, handles)
% % % % % % % % % % hObject    handle to recorder_instructions (see GCBO)
% % % % % % % % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % % % % % % % handles    empty - handles not created until after all CreateFcns called



function nretries_gui_Callback(hObject, ~, handles)
% hObject    handle to nretries_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nretries_gui as text
%        str2double(get(hObject,'String')) returns contents of nretries_gui as a double
hMain = getappdata(0,'hMain');
nretries = round(str2double(get(hObject,'String')));
if (isnan(nretries) || nretries<0)
    set(hObject,'String',num2str(handles.nretries_gui))
    warndlg('All inputs MUST be real numbers >=0!');
else
    set(hObject,'String', nretries)
    setappdata(hMain,'nretries', handles.nretries_gui)
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nretries_gui_CreateFcn(hObject, ~ , handles)
% hObject    handle to nretries_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hMain = getappdata(0,'hMain');
try
    handles.nretries_gui = getappdata(hMain,'nretries',nretries)
    set(hObject, 'String', handles.nretries_gui)
catch
    handles.nretries_gui = '1';
    set(hObject, 'String', handles.nretries_gui)
end


% --- Executes on button press in turn_chk.
function turn_chk_Callback(hObject, eventdata, handles)
% hObject    handle to turn_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of turn_chk
if get(hObject, 'Value') == 1
    set(handles.turn_abs, 'Visible', 'On')
    set(handles.turn_abs, 'String', '0:5:355')
    set(handles.text27, 'String', 'deg')
    set(handles.PlaybackDelayText, 'String', 'Turn increment')
    set(handles.ttAccel, 'Visible', 'On')
    set(handles.ttAccel_txt, 'Visible', 'On')
else
    set(handles.turn_abs, 'Visible', 'off')
    set(handles.text27, 'String', 's')
    set(handles.PlaybackDelayText, 'String', 'Playback Delay')
    set(handles.ttAccel, 'Visible', 'Off')
    set(handles.ttAccel_txt, 'Visible', 'Off')
end


function turn_abs_Callback(hObject, eventdata, handles)
% hObject    handle to turn_abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of turn_abs as text
%        str2double(get(hObject,'String')) returns contents of turn_abs as a double


% --- Executes during object creation, after setting all properties.
function turn_abs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to turn_abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject, 'String', '0:5:355')
end



function ttAccel_Callback(hObject, eventdata, handles)
% hObject    handle to ttAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttAccel as text
%        str2double(get(hObject,'String')) returns contents of ttAccel as a double

% add check that end rotation is multiple of angle increment


% --- Executes during object creation, after setting all properties.
function ttAccel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', '350')

% --- Executes on button press in turn_cont_chk.
function turn_cont_chk_Callback(hObject, eventdata, handles)
% hObject    handle to turn_cont_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of turn_cont_chk
if get(hObject, 'Value') == 1
    set(handles.turn_chk, 'Value', 1)
    set(handles.turn_abs, 'Visible', 'On')
    set(handles.turn_abs, 'String', '30')
    set(handles.text27, 'String', 's/rev')
    set(handles.PlaybackDelayText, 'String', '(22.7-720) s/rev')
    set(handles.ttAccel, 'Visible', 'Off')
    set(handles.ttAccel_txt, 'Visible', 'Off')
else
    if get(handles.turn_chk, 'Value') == 0
        set(handles.turn_abs, 'Visible', 'off')
        set(handles.text27, 'String', 's')
        set(handles.PlaybackDelayText, 'String', 'Playback Delay')
        set(handles.ttAccel, 'Visible', 'Off')
        set(handles.ttAccel_txt, 'Visible', 'Off')
    else
        set(handles.turn_abs, 'Visible', 'On')
        set(handles.turn_abs, 'String', '0:5:355')
        set(handles.text27, 'String', 'deg')
        set(handles.PlaybackDelayText, 'String', 'Turn increment')
        set(handles.ttAccel, 'Visible', 'On')
    end
end


% --- Executes on button press in stack_btn.
function stack_btn_Callback(hObject, eventdata, handles)
% hObject    handle to stack_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[len,chans,~,cycles,outchans,dim6] = size(handles.rec);
% consider locking IO settings or trust user knows what they're doing?
%jono: TODO: check for unintended consequences... 

if cycles >1 
    if ~isfield(handles, 'stack')
    dimwarning = warndlg('dimension 4 > 1, stacking in dimension 6');
    uiwait(dimwarning)
    end
    handles.stackdim = 6;
elseif isfield(handles, 'outputdata') &&...
        isfield(handles.outputdata,'properties') &&...
        isfield(handles.outputdata.properties,'startflag')
    if ~isfield(handles, 'stack')
    dimwarning = warndlg('dimension 4 reserved for silent cycle, stacking in dimension 6');
    uiwait(dimwarning)
    end
    handles.stackdim = 6;
else
    handles.stackdim = 4;
end
if ~isfield(handles, 'stack') && ~isempty(handles.rec)
    handles.stack = handles.rec;
    handles.stackrecordtime = handles.recordtime;
    handles.stackname = {'manualstack001'}; % or ask user...
    count = 1;
elseif ~isempty(handles.rec)
    [d1,d2,d3,d4,d5,d6] = size(handles.stack);
    switch handles.stackdim
        case 4
           if d1==len && d2 == chans && d5 == outchans && d6 == dim6
            handles.stack(:,:,:,d4+1,:,:) = handles.rec;
            count = d4+1;
           else
               error('cannot stack as array sizes differ')
           end
        case 6
            if d1==len && d2 == chans && d4 == cycles && d5 == outchans 
            handles.stack(:,:,:,:,:,d6+1) = handles.rec;
            count = d6+1;
           else
               error('cannot stack as array sizes differ')
            end
    end
end 
set(handles.stack_btn, 'Visible', 'off')
set(handles.stack_btn, 'String', ['Stack ', num2str(count), ' + 1'])
set(handles.stack_edit_btn, 'Visible', 'on', 'Enable', 'on')
handles.stack_size_string = ['stack size: [',...
    num2str(size(handles.stack,1),'%10.3e\n'),...
    '] x [', num2str(size(handles.stack,2)),...
    '] x [', num2str(size(handles.stack,3)),...
    '] x [', num2str(size(handles.stack,4)),...
    '] x [', num2str(size(handles.stack,5)),...
    '] x [', num2str(size(handles.stack,6)),']'];
set(handles.stack_size_text,'string',handles.stack_size_string, 'Enable','on');

guidata(hObject, handles)


% --- Executes on button press in stack_edit_btn.
function stack_edit_btn_Callback(hObject, eventdata, handles)
% hObject    handle to stack_edit_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dimsize = size(handles.stack,handles.stackdim);
param = inputdlg(['select dim ', num2str(handles.stackdim),...
    ' indicies to keep'],'Edit Dims',...
    [1 35],...
    {['1:',num2str(size(handles.stack,handles.stackdim))]});
switch handles.stackdim
    case 4
        handles.stack = handles.stack(:,:,:,str2num(param{:}),:,:);
    case 6
        handles.stack = handles.stack(:,:,:,:,:,str2num(param{:}));
end
set(handles.stack_btn, 'String', ['Stack ', num2str(size(handles.stack,handles.stackdim)), ' + 1'])
 handles.stack_size_string = ['stack size: [',...
    num2str(size(handles.stack,1),'%10.3e\n'),...
    '] x [', num2str(size(handles.stack,2)),...
    '] x [', num2str(size(handles.stack,3)),...
    '] x [', num2str(size(handles.stack,4)),...
    '] x [', num2str(size(handles.stack,5)),...
    '] x [', num2str(size(handles.stack,6)),']'];
set(handles.stack_size_text,'string',handles.stack_size_string, 'Enable','on');
guidata(hObject, handles)
