function varargout = syscal(varargin)
% SYSCAL MATLAB code for syscal.fig
%      SYSCAL, by itself, creates a new SYSCAL or raises the existing
%      singleton*.
%
%      H = SYSCAL returns the handle to a new SYSCAL or the handle to
%      the existing singleton*.
%
%      SYSCAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYSCAL.M with the given input arguments.
%
%      SYSCAL('Property','Value',...) creates a new SYSCAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before syscal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to syscal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help syscal

% Last Modified by GUIDE v2.5 05-Nov-2015 09:06:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @syscal_OpeningFcn, ...
                   'gui_OutputFcn',  @syscal_OutputFcn, ...
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


% --- Executes just before syscal is made visible.
function syscal_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to syscal (see VARARGIN)

% This next couple of lines checks if the GUI is being called from the main
% window, otherwise it doesn't run.
dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'audio_recorder'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.main_stage1 = varargin{mainGuiInput+1};

    fontsize

end

if dontOpen
   disp('-----------------------------------------------------');
   disp('This function is part of the AARAE framework, it is') 
   disp('not a standalone function. To call this function,')
   disp('click on the appropriate calling button on the Audio');
   disp('Recording Window. E.g.:');
   disp('   Recording > SYSCAL');
   disp('-----------------------------------------------------');
else
    handles.mainHandles = guidata(handles.main_stage1);
    if get(handles.mainHandles.pb_enable,'Value') == 0, set(handles.evaldelay_btn,'Enable','off'); end
    outputs = cellstr(get(handles.mainHandles.outputdev_popup,'String'));
    outputdevname = outputs{get(handles.mainHandles.outputdev_popup,'Value')};
    handles.hap = dsp.AudioPlayer('DeviceName',outputdevname,'SampleRate',handles.mainHandles.fs,'QueueDuration',handles.mainHandles.qdur,'BufferSizeSource','Property','BufferSize',handles.mainHandles.buffer);
    inputs = cellstr(get(handles.mainHandles.inputdev_popup,'String'));
    inputdevname = inputs{get(handles.mainHandles.inputdev_popup,'Value')};
    handles.har = dsp.AudioRecorder('DeviceName',inputdevname,'SampleRate',handles.mainHandles.fs,'QueueDuration',handles.mainHandles.qdur,'OutputDataType','double','NumChannels',1,'BufferSizeSource','Property','BufferSize',handles.mainHandles.buffer);
    handles.hsr1 = dsp.SignalSource;
    handles.hsr1.SamplesPerFrame = 1024;
    set(handles.duration_IN,'String','10')
    set(handles.sperframe_IN,'String',num2str(handles.mainHandles.fs*0.1))
    set(handles.percentage_IN,'String','40')
    set(handles.threshold_IN,'String','0.1')
    set(handles.tonelevel_IN,'String','94')
    set(handles.units,'String','Pa')
    set(handles.units_ref,'String','2e-5')
    set(handles.units_type,'String',{'Amplitude units';'Power units'})
    stats{1} = ['Sampling frequency: ' num2str(handles.mainHandles.fs) ' samples/s'];
    set(handles.statstext,'String',stats);
    UserData.state = false;
    set(handles.stop_btn,'UserData',UserData);
    handles.output = struct;
    xlabel(handles.IRaxes,'Time [s]')
    ylim(handles.IRaxes,[-60 10])
    deviceinfo = dspAudioDeviceInfo('inputs');
    set(handles.devnametext,'String',inputdevname)
    inputsel = get(handles.mainHandles.inputdev_popup,'Value');
    channels = cellstr([repmat('IN - Channel',deviceinfo(inputsel).maxInputs,1) num2str((1:deviceinfo(inputsel).maxInputs).')]);
    set(handles.channum_popup,'String',channels)
    handles.output.cal(1,1:deviceinfo(inputsel).maxInputs) = NaN;
    handles.output.units = 'Pa';
    handles.output.units_ref = 2e-5;
    handles.output.units_type = 1;
    guidata(hObject, handles);
    uiwait(hObject);
end

% UIWAIT makes syscal wait for user response (see UIRESUME)
% uiwait(handles.syscal);

% --- Executes when user attempts to close syscal.
function syscal_CloseRequestFcn(hObject, ~, handles) %#ok : Window close request function
% hObject    handle to syscal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = struct;
guidata(hObject,handles)
% Hint: delete(hObject) closes the figure
uiresume(hObject);


% --- Executes on button press in done_btn.
function done_btn_Callback(~, ~, handles) %#ok : Executed when Done bunnos is clicked
% hObject    handle to done_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.syscal);



% --- Outputs from this function are returned to the command line.
function varargout = syscal_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);



% *************************************************************************
% *************************************************************************
% SYSTEM LATENCY MEASUREMENT
% *************************************************************************
% *************************************************************************

% --- Executes on selection change in stimulus_popup.
function stimulus_popup_Callback(~, ~, ~) %#ok : Executed when latency stimulus menu changes
% hObject    handle to stimulus_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stimulus_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stimulus_popup


% --- Executes during object creation, after setting all properties.
function stimulus_popup_CreateFcn(hObject, ~, ~) %#ok : Latency stimulus menu creation
% hObject    handle to stimulus_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in evaldelay_btn.
function evaldelay_btn_Callback(hObject, eventdata, handles) %#ok : Executed when evaluate button is clicked in the latency panel
% hObject    handle to evaldelay_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
duration = str2double(get(handles.slength_IN,'String'));
set(hObject,'Enable','off');
stimulus = get(handles.stimulus_popup,'Value');
switch stimulus
    case 1
        S = impulse1(0,duration,1,handles.mainHandles.fs);
        handles.hsr1.Signal = S.audio(1:duration*handles.mainHandles.fs);
    case 2
        S = linear_sweep(duration,1,handles.mainHandles.fs/2,0,handles.mainHandles.fs);
        handles.hsr1.Signal = S.audio;
    case 3
        S = noise(-1,duration,handles.mainHandles.fs,1,handles.mainHandles.fs/2,1,0,0);
        handles.hsr1.Signal = S.audio;
    case 4
        S = OATSP(duration,0.5,handles.mainHandles.fs,1,0);
        handles.hsr1.Signal = S.audio./max(abs(S.audio));
end
handles.hsr1.Signal = [handles.hsr1.Signal;zeros(handles.hap.QueueDuration*handles.mainHandles.fs,1)];
set(hObject,'BackgroundColor','red');
pause on
pause(0.000001)
pause off
try
    ncycles = ceil(length(handles.hsr1.Signal)/handles.hsr1.SamplesPerFrame);
    audio = zeros(ncycles*handles.hsr1.SamplesPerFrame,1);
    for i = 1:ncycles
        audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:) = step(handles.har);
        step(handles.hap,step(handles.hsr1));
    end
catch sthgwrong
    syswarning = sthgwrong.message;
    set(hObject,'Enable','on');
    warndlg(syswarning,'AARAE info')
    return
end
rec = audio;
if ~isempty(rec)
    qd = handles.mainHandles.fs*handles.har.QueueDuration;
    rec = rec(qd:end);
    rec = rec(1:length(S.audio));
    rec = [rec;zeros(length(handles.hsr1.Signal),1)];
    if stimulus == 1
        ixy = rec;
    else
        ixy = ifft(fft(rec) .* fft([S.audio2;zeros(length(rec)-length(S.audio2),1)]));
        ixy = fftshift(ixy(1:duration*handles.mainHandles.fs*2));
    end
    ixy = ixy(1:ceil(end/2));
    LIxy = (mean((20.*log10(abs(fft(ixy)))).^2)).^0.5;
    ixy = ixy./10^(LIxy/20);
    t = linspace(0,length(ixy)/handles.mainHandles.fs,length(ixy));
    IRlevel = (10.*log10((abs(ixy)./max(abs(ixy))).^2));
    handles.latthresh = str2double(get(handles.latthresh_IN,'String'));
    abovethresh = find(IRlevel > abs(handles.latthresh)*-1);
    [~,I1] = max(IRlevel(abovethresh));
    handles.maxIR = abovethresh(I1);
    I = abovethresh(1);
    plot(handles.IRaxes,t,IRlevel,t,ones(size(t)).*handles.latthresh.*-1,'r',handles.maxIR/handles.mainHandles.fs,IRlevel(handles.maxIR),'or')
    hold(handles.IRaxes,'on')
    plot(handles.IRaxes,I/handles.mainHandles.fs,IRlevel(I),'o','Color',[0 .6 0])
    hold(handles.IRaxes,'off')
    xlabel(handles.IRaxes,'Time [s]')
    ylim(handles.IRaxes,[-60 10])
    set(handles.latency_IN,'String',num2str(I),'Enable','on');
    set(handles.latthresh_IN,'Enable','on')
    handles.output.audio2 = ixy;
    handles.output.fs = handles.mainHandles.fs;
    handles.output.latency = I;
    handles.sysIR = ixy(I:end); % check this - Daniel did not use the index limit
end
release(handles.hap)
release(handles.har)
release(handles.hsr1)
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(handles.latencytext,'String',[num2str(I) ' samples = ~' num2str(I/handles.mainHandles.fs) ' s'])
set([hObject handles.invfdesign_btn handles.invf_popup handles.invfsmooth_popup handles.invfLF_IN handles.invfHF_IN handles.invfIB_IN handles.invfOB_IN handles.invfnfft_IN handles.invflength_IN handles.preproc_popup handles.IRlength_IN handles.postproc_popup],'Enable','on');
set(handles.IRlength_IN,'String',num2str(length(ixy)-I))
set(handles.invfpreview_btn,'Enable','off')
set(handles.invftext,'Visible','off')
set(handles.invfpreview_btn,'Enable','off')
if isfield(handles,'IRwinline'), handles = rmfield(handles,'IRwinline'); end
guidata(hObject,handles)
preproc_popup_Callback(handles.preproc_popup,eventdata,handles)



function slength_IN_Callback(hObject, ~, ~) %#ok : Executes when signal duration changes on the latency panel
% hObject    handle to slength_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slength_IN as text
%        str2double(get(hObject,'String')) returns contents of slength_IN as a double
slength = str2double(get(hObject,'String'));
if isnan(slength) || slength <=0
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String','1')
end


% --- Executes during object creation, after setting all properties.
function slength_IN_CreateFcn(hObject, ~, ~) %#ok : Signal duration input box creation
% hObject    handle to slength_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function latthresh_IN_Callback(hObject, eventdata, handles) %#ok : Executed when latency threshold input box changes
% hObject    handle to latthresh_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latthresh_IN as text
%        str2double(get(hObject,'String')) returns contents of latthresh_IN as a double
latthresh = abs(str2double(get(hObject,'String')));
if ~isnan(latthresh) && 0 < latthresh && latthresh <= 60
    ixy = handles.sysIR;
    handles.latthresh = latthresh;
    t = linspace(0,length(ixy)/handles.mainHandles.fs,length(ixy));
    IRlevel = (10.*log10((abs(ixy)./max(abs(ixy))).^2));
    abovethresh = find(IRlevel > latthresh*-1);
    [~,I1] = max(IRlevel(abovethresh));
    handles.maxIR = abovethresh(I1);
    I = abovethresh(1);
    plot(handles.IRaxes,t,IRlevel,t,ones(size(t)).*latthresh.*-1,'r',handles.maxIR/handles.mainHandles.fs,IRlevel(handles.maxIR),'or')
    hold(handles.IRaxes,'on')
    plot(handles.IRaxes,I/handles.mainHandles.fs,IRlevel(I),'o','Color',[0 .6 0])
    hold(handles.IRaxes,'off')
    ylim(handles.IRaxes,[-60 10])
    handles.output.latency = I;
    set(handles.latency_IN,'String',num2str(I))
    set(handles.IRlength_IN,'String',num2str(length(handles.sysIR)-I))
    set(handles.latencytext,'String',[num2str(I) ' samples = ~' num2str(I/handles.mainHandles.fs) ' s'])
else
    set(hObject,'String',num2str(handles.latthresh))
    warndlg('Threshhold out of boundaries','AARAE info')
end
handles = rmfield(handles,'IRwinline');
guidata(hObject,handles)
preproc_popup_Callback(handles.preproc_popup,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function latthresh_IN_CreateFcn(hObject, ~, ~) %#ok : latency threshold input box creation
% hObject    handle to latthresh_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function latency_IN_Callback(hObject, eventdata, handles) %#ok : Executed when latency input box changes
% hObject    handle to latency_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency_IN as text
%        str2double(get(hObject,'String')) returns contents of latency_IN as a double
latency = str2double(get(hObject,'String'));
if ~isnan(latency) && latency <= handles.maxIR && latency > 0
    handles.output.latency = round(latency);
    set(hObject,'String',round(latency))
    ixy = handles.sysIR;
    t = linspace(0,length(ixy)/handles.mainHandles.fs,length(ixy));
    IRlevel = (10.*log10((abs(ixy)./max(abs(ixy))).^2));
    plot(handles.IRaxes,t,IRlevel,t,ones(size(t)).*handles.latthresh.*-1,'r',handles.maxIR/handles.mainHandles.fs,IRlevel(handles.maxIR),'or')
    hold(handles.IRaxes,'on')
    plot(handles.IRaxes,handles.output.latency/handles.mainHandles.fs,IRlevel(handles.output.latency),'o','Color',[0 .6 0])
    hold(handles.IRaxes,'off')
    ylim(handles.IRaxes,[-60 10])
    set(handles.IRlength_IN,'String',num2str(length(ixy)-handles.output.latency))
    set(handles.latencytext,'String',[num2str(handles.output.latency) ' samples = ~' num2str(latency/handles.mainHandles.fs) ' s'])
else
    set(hObject,'String',num2str(handles.output.latency))
    warndlg('Input cannot be greater than the auto-detected latency or zero','AARAE info')
end
handles = rmfield(handles,'IRwinline');
guidata(hObject,handles)
preproc_popup_Callback(handles.preproc_popup,eventdata,handles)


% --- Executes during object creation, after setting all properties.
function latency_IN_CreateFcn(hObject, ~, ~) %#ok : latency input box creation
% hObject    handle to latency_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in latinstructions_btn.
function latinstructions_btn_Callback(~, ~, ~) %#ok : Executed when instructions button is pressed in the latency panel
% hObject    handle to latinstructions_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message{1} = 'Instructions:';
message{2} = '1. Make sure your audio interface is correctly plugged into your machine.';
message{3} = '2. Make a feedback loop between INPUT 1 and OUTPUT 1';
message{4} = '3. Press the -Evaluate- button to test the system latency.';

helpdlg(message,'AARAE info')










% *************************************************************************
% *************************************************************************
% INVERSE FILTER OF SYSTEM RESPONSE
% *************************************************************************
% *************************************************************************



function IRlength_IN_Callback(hObject, eventdata, handles) %#ok : Executed when the IR length input box changes
% hObject    handle to IRlength_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IRlength_IN as text
%        str2double(get(hObject,'String')) returns contents of IRlength_IN as a double
IRlength = str2double(get(hObject,'String'));
if ~isnan(IRlength) && IRlength <= (length(handles.sysIR)-handles.output.latency) && IRlength > 0
    preproc_popup_Callback(handles.preproc_popup,eventdata,handles)
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
else
    set(hObject,'String',num2str(length(handles.sysIR)-handles.output.latency))
    preproc_popup_Callback(handles.preproc_popup,eventdata,handles)
    warndlg('Invalid window length','AARAE info')
end
%guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function IRlength_IN_CreateFcn(hObject, ~, ~) %#ok : IRlength input box creation
% hObject    handle to IRlength_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in invf_popup.
function invf_popup_Callback(hObject, ~, handles) %#ok : Executed when inverse filter type menu changes
% hObject    handle to invf_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns invf_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from invf_popup
set(handles.invftext,'Visible','off')
set(handles.invfpreview_btn,'Enable','off')
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function invf_popup_CreateFcn(hObject, ~, ~) %#ok : Inverse filter type menu creation
% hObject    handle to invf_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in invfsmooth_popup.
function invfsmooth_popup_Callback(hObject, ~, handles) %#ok : Executed when smoothing type for inverse filter menu changes
% hObject    handle to invfsmooth_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns invfsmooth_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from invfsmooth_popup
set(handles.invftext,'Visible','off')
set(handles.invfpreview_btn,'Enable','off')
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function invfsmooth_popup_CreateFcn(hObject, ~, ~) %#ok : inverse filter smoothing type menu creation
% hObject    handle to invfsmooth_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function invfLF_IN_Callback(hObject, ~, handles) %#ok : inverse filter low cutoff frequency input box creation
% hObject    handle to invfLF_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invfLF_IN as text
%        str2double(get(hObject,'String')) returns contents of invfLF_IN as a double
num = str2double(get(hObject,'String'));
if isnan(num) || num <= 0 || num >= handles.invfHF
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.invfLF))
else
    handles.invfLF = num;
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function invfLF_IN_CreateFcn(hObject, ~, handles) %#ok : inverse filter low cutoff frequency input box creation
% hObject    handle to invfLF_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.invfLF = str2double(get(hObject,'String'));
guidata(hObject,handles)


function invfHF_IN_Callback(hObject, ~, handles) %#ok : Executed when inverse filter high frequency cutoff input box changes
% hObject    handle to invfHF_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invfHF_IN as text
%        str2double(get(hObject,'String')) returns contents of invfHF_IN as a double
num = str2double(get(hObject,'String'));
if isnan(num) || num <= 0 || num <= handles.invfLF
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.invfHF))
else
    handles.invfHF = num;
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function invfHF_IN_CreateFcn(hObject, ~, handles) %#ok : Inverse filter low cutoff frequency input box creation
% hObject    handle to invfHF_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.invfHF = str2double(get(hObject,'String'));
guidata(hObject,handles)


function invfIB_IN_Callback(hObject, ~, handles) %#ok : Executes when inverse filter in-bound level input box changes
% hObject    handle to invfIB_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invfIB_IN as text
%        str2double(get(hObject,'String')) returns contents of invfIB_IN as a double
num = str2double(get(hObject,'String'));
if isnan(num)
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.invfIB))
else
    handles.invfIB = num;
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function invfIB_IN_CreateFcn(hObject, ~, handles) %#ok : Inverse filter in-bound level input box creation
% hObject    handle to invfIB_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.invfIB = str2double(get(hObject,'String'));
guidata(hObject,handles)


function invfOB_IN_Callback(hObject, ~, handles) %#ok : Executes when inverse filter out-of-bound level input box changes
% hObject    handle to invfOB_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invfOB_IN as text
%        str2double(get(hObject,'String')) returns contents of invfOB_IN as a double
num = str2double(get(hObject,'String'));
if isnan(num)
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.invfOB))
else
    handles.invfOB = num;
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function invfOB_IN_CreateFcn(hObject, ~, handles) %#ok : Inverse filter out-of-bound level input box creation
% hObject    handle to invfOB_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.invfOB = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes on button press in invfdesign_btn.
function invfdesign_btn_Callback(hObject, ~, handles) %#ok : Executed when design button is clicken on inferse filter panel
% hObject    handle to invfdesign_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');
pause on
pause(0.000001)
pause off
fs = handles.mainHandles.fs;
f = [handles.invfLF,handles.invfHF];
reg = [handles.invfIB,handles.invfOB];
nfft = handles.invfnfft;
flength = handles.invflength;

% Thresholding
if f(1) > 0 && f(2) < fs/2
    freq=(0:fs/(nfft-1):fs/2)'; 
    f1e=f(1)-f(1)/3;
    f2e=f(2)+f(2)/3;
    if f1e < freq(1)
        f1e=f(1);
        f(1)=f(1)+1;
    end
    if f2e > freq(end)
        f2e=f(2)+1;
    end
    % thersholding B with 1/3 octave interpolated transient edges 
    B=interp1([0 f1e f(1) f(2) f2e freq(end)],[reg(2) reg(2) reg(1) reg(1) reg(2) reg(2)],freq,'pchip');
    B=10.^(-B./20); % from dB to linear
    B=vertcat(B,B(end:-1:1)); 
    b=ifft(B,'symmetric');
    b=circshift(b,nfft/2);
    b=0.5*(1-cos(2*pi*(1:nfft)'/(nfft+1))).*b; % window function
    b=minph(b); % make minimum phase thresholding
    B=fft(b,nfft);
else
    B=0;
end
% Inverse filter design
sysIR = handles.sysIR;%(str2num(get(handles.latency_IN,'String')):end);
if get(handles.preproc_popup,'Value') ~= 1
    IRwindow = [handles.IRwindow;zeros(length(sysIR)-length(handles.IRwindow),1)];
    sysIR = sysIR.*IRwindow;
end
if get(handles.invf_popup,'Value') ~= 4
    H = abs(fft(sysIR,nfft));
else
    H = fft(sysIR,nfft);
end
smoothfactor = get(handles.invfsmooth_popup,'Value');
if smoothfactor == 2, octsmooth = 1; end
if smoothfactor == 3, octsmooth = 3; end
if smoothfactor == 4, octsmooth = 6; end
if smoothfactor == 5, octsmooth = 12; end
if smoothfactor == 6, octsmooth = 24; end
if smoothfactor ~= 1, H = octavesmoothing(H, octsmooth, fs); end
iH=conj(H)./((conj(H).*H)+(conj(B).*B)); % calculating regulated spectral inverse
iH = circshift(ifft(iH,'symmetric'),nfft/2); % could use ifftshift instead
if get(handles.invf_popup,'Value') == 1, handles.invfilter=minph(iH); end
if get(handles.invf_popup,'Value') == 3, handles.invfilter = flipud(minph(iH)); end
if get(handles.invf_popup,'Value') == 4 || get(handles.invf_popup,'Value') == 2, handles.invfilter = iH; end
if get(handles.invf_popup,'Value') == 5
    iHspec = fft(iH);
    phase = angle(iHspec);
    rmsmag = mean(abs(iHspec).^2)^0.5;
    changed_spectrum = ones(length(iHspec),1)*rmsmag .* exp(1i * phase);
    changed_spectrum(1) = 0; % make DC zero
    changed_spectrum(round(length(iHspec)/2)+1) = 0; % make Nyquist zero
    handles.invfilter = ifft(changed_spectrum);
end
if strcmp(get(handles.invflength_IN,'Visible'),'on')
    invfwintype = get(handles.postproc_popup,'Value');
    switch invfwintype
        case 2
            wintype = @rectwin;
        case 3
            wintype = @hann;
    end
    invftype = get(handles.invf_popup,'Value');
    switch invftype
        case 1
            invfwindow = window(wintype,2*flength);
            invfwindow = invfwindow(end/2+1:end);
            invfwindow = [invfwindow;zeros(nfft-flength,1)];
        case 2
            invfwindow = window(wintype,flength);
            invfwindow = [zeros((nfft-flength)/2,1);invfwindow;zeros((nfft-flength)/2,1)];
        case 3
            invfwindow = window(wintype,2*flength);
            invfwindow = invfwindow(end/2+1:end);
            invfwindow = flipud([invfwindow;zeros(nfft-flength,1)]);
        case 4
            invfwindow = window(wintype,flength);
            invfwindow = [zeros((nfft-flength)/2,1);invfwindow;zeros((nfft-flength)/2,1)];
        case 5
            invfwindow = window(wintype,flength);
            invfwindow = [zeros((nfft-flength)/2,1);invfwindow;zeros((nfft-flength)/2,1)];
    end
    handles.invfilter = handles.invfilter.*invfwindow;
end
handles.output.invfilter = handles.invfilter;
set(handles.invftext,'Visible','on')
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set([hObject handles.invfpreview_btn],'Enable','on');
guidata(hObject,handles)


function [h_min] = minph(h)
n = length(h);
h_cep = real(ifft(log(abs(fft(h(:,1))))));
odd = fix(rem(n,2));
wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
h_min = zeros(size(h(:,1)));
h_min(:) = real(ifft(exp(fft(wn.*h_cep(:)))));
% Find the last zero-crossing & apply circular shift.
% Although this compromises the minimum phase response, it reduces the
% artefact from linear convolution of a circularly derived response
s=find((abs(h_min)<abs([h_min(2:end);h_min(1)]))...
    &(abs(h_min(:))<abs([h_min(end);h_min(1:end-1)])),1,'last');
h_min = circshift(h_min,n-s);


% --- Executes on button press in invfpreview_btn.
function invfpreview_btn_Callback(hObject, ~, handles) %#ok : Executed when clicking on the Preview button in the inverse filter panel
% hObject    handle to invfpreview_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');
pause on
pause(0.000001)
pause off

origIR = handles.sysIR;%(handles.output.latency:end);
invfilter = handles.invfilter;

sysIRspec = fft(origIR,length(invfilter));
invfilterspec = fft(invfilter);
convspec = sysIRspec.*invfilterspec;

t = linspace(0,length(invfilter)/handles.mainHandles.fs,length(invfilter));

sysIRflevel = 10.*log10(abs(sysIRspec).^2);
invfilterflevel = 10.*log10(abs(invfilterspec).^2);
convflevel = 10.*log10(abs(convspec).^2);

origIRgd = -diff(unwrap(angle(sysIRspec))).*length(sysIRspec)/(handles.mainHandles.fs*2*pi).*1000;
invfiltergd = -diff(unwrap(angle(invfilterspec))).*length(invfilterspec)/(handles.mainHandles.fs*2*pi).*1000;
timeconvgd = -diff(unwrap(angle(convspec))).*length(convspec)/(handles.mainHandles.fs*2*pi).*1000;

IRtime = abs(ifft(sysIRspec));
convIRtime = abs(ifft(convspec));

freq=0:handles.mainHandles.fs/(length(invfilter)-1):handles.mainHandles.fs/2;
set(0,'Units','characters')
screen = get(0,'Screensize');
set(0,'Units','pixels')
figure;set(gcf,'Units','characters')
set(gcf, 'Position', screen);
subplot(2,2,1);plot(t,IRtime,'k',t,convIRtime,'r')
               hold on
               plot(t,invfilter,'Color',[0 .6 0])
               hold off
               title('Absolute amplitude');xlabel('Time [s]');ylabel('Amplitude')
               leg = legend('System IR','sysIR*filtIR','Filter IR');
               set(leg,'Units','characters')
               posleg = get(leg,'Position');
               set(leg,'Position',[screen(3:4)./2-5 posleg(3:4)])
subplot(2,2,3);plot(t,10.*log10(IRtime.^2),'k',t,10.*log10(convIRtime.^2),'r')
               hold on
               plot(t,10.*log10(invfilter.^2),'Color',[0 .6 0])
               hold off
               title('Squared amplitude');xlabel('Time [s]');ylabel('Amplitude [dB]')
subplot(2,2,2);semilogx(freq,sysIRflevel(1:end/2),'k',freq,convflevel(1:end/2),'r');
               hold on
               semilogx(freq,invfilterflevel(1:end/2),'Color',[0 .6 0])
               hold off
               title('Magnitude spectrum');xlabel('Frequency [Hz]');ylabel('Amplitude [dB]')
subplot(2,2,4);semilogx(freq,origIRgd(1:length(freq)),'k',freq,timeconvgd(1:length(freq)),'r')
               hold on
               semilogx(freq,invfiltergd(1:length(freq)),'Color',[0 .6 0])
               hold off
               title('Group delay');xlabel('Frequency [Hz]');ylabel('Time [ms]')
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');


function invfnfft_IN_Callback(hObject, ~, handles) %#ok : Executed when the nfft input box changes
% hObject    handle to invfnfft_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invfnfft_IN as text
%        str2double(get(hObject,'String')) returns contents of invfnfft_IN as a double
num = str2double(get(hObject,'String'));
if isnan(num) || num <= 0
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.invfnfft))
else
    if mod(num,2) ~= 0, num = num + 1; set(hObject,'String',num2str(num)); end
    handles.invfnfft = num;
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
end
if num < handles.invflength
    set(handles.invflength_IN,'String',num2str(num))
    handles.invflength = num;
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function invfnfft_IN_CreateFcn(hObject, ~, handles) %#ok : nfft input box creation
% hObject    handle to invfnfft_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.invfnfft = str2double(get(hObject,'String'));
guidata(hObject,handles)


function invflength_IN_Callback(hObject, ~, handles) %#ok : Executed when inverse filter lenght input box changes
% hObject    handle to invflength_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invflength_IN as text
%        str2double(get(hObject,'String')) returns contents of invflength_IN as a double
num = str2double(get(hObject,'String'));
if isempty(num) || num <= 0 || num > handles.invfnfft
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.invflength))
else
    if mod(num,2) ~= 0, num = num + 1; set(hObject,'String',num2str(num)); end
    handles.invflength = num;
    set(handles.invftext,'Visible','off')
    set(handles.invfpreview_btn,'Enable','off')
end
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function invflength_IN_CreateFcn(hObject, ~, handles) %#ok : inverse filter length input box creation
% hObject    handle to invflength_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.invflength = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes on selection change in preproc_popup.
function preproc_popup_Callback(hObject, ~, handles)
% hObject    handle to preproc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns preproc_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from preproc_popup
if get(hObject,'Value') == 1
    if isfield(handles,'IRwinline'), delete(handles.IRwinline); handles = rmfield(handles,'IRwinline'); end
    set([handles.prelengthtext handles.IRlength_IN handles.presamplestext],'Visible','off')
else
    set([handles.prelengthtext handles.IRlength_IN handles.presamplestext],'Visible','on')
end
IRlength = str2double(get(handles.IRlength_IN,'String'));
if get(hObject,'Value') == 2
    if isfield(handles,'IRwinline'), delete(handles.IRwinline); handles = rmfield(handles,'IRwinline'); end
    hold(handles.IRaxes,'on')
    handles.IRwindow = window(@rectwin,IRlength);
    t = (linspace(0,IRlength,IRlength)+handles.output.latency)./handles.mainHandles.fs;
    handles.IRwinline = plot(handles.IRaxes,t,10.*log10(handles.IRwindow),'Color',[1 .6 0]);
    hold(handles.IRaxes,'off')
end
if get(hObject,'Value') == 3
    if isfield(handles,'IRwinline'), delete(handles.IRwinline); handles = rmfield(handles,'IRwinline'); end
    hold(handles.IRaxes,'on')
    handles.IRwindow = window(@hann,2*IRlength);
    handles.IRwindow = handles.IRwindow(end/2+1:end);
    t = (linspace(0,IRlength,IRlength)+handles.output.latency)./handles.mainHandles.fs;
    handles.IRwinline = plot(handles.IRaxes,t,10.*log10(handles.IRwindow),'Color',[1 .6 0]);
    hold(handles.IRaxes,'off')
end
set(handles.invftext,'Visible','off')
set(handles.invfpreview_btn,'Enable','off')
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function preproc_popup_CreateFcn(hObject, ~, ~) %#ok : Windowing menu creation for preprocessing subpanel
% hObject    handle to preproc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in postproc_popup.
function postproc_popup_Callback(hObject, ~, handles) %#ok : Executes when windowing menu in the postprocessing subpanel changes
% hObject    handle to postproc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns postproc_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from postproc_popup
if get(hObject,'Value') == 1
    set([handles.postlengthtext handles.invflength_IN handles.postsamplestext],'Visible','off')
else
    set([handles.postlengthtext handles.invflength_IN handles.postsamplestext],'Visible','on')
end
set(handles.invftext,'Visible','off')
set(handles.invfpreview_btn,'Enable','off')
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function postproc_popup_CreateFcn(hObject, ~, ~) %#ok : Windowing menu on the post processing panel creation
% hObject    handle to postproc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in invfiltinstructions_btn.
function invfiltinstructions_btn_Callback(~, ~, ~) %#ok : Executed when instruction button is clicked on the inverse filter panel
% hObject    handle to invfiltinstructions_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message{1} = 'Instructions:';
message{2} = '1. Once you have captured the system impulse response you may set the parameters relevant to generate the inverse filter.';
message{3} = '2. Select from the popup menu the type of filter.';
message{4} = '3. Select the type of smoothing on your filter.';
message{5} = '4. Set the frequency limits for your filter.';
message{6} = '5. Set the in-band and out-of-band gain thresholds.';
message{7} = '6. Set the FFT length and filter length in samples.';
message{7} = 'NOTE: The filter length may not be longer than the FFT length.';
message{8} = ' ';
message{9} = '7. Click on -Design- to design the inverse filter.';
message{10} = '8. Click on -Preview- to view the designed filter compared to the system impulse response and the convolution of your filter with the impulse response.';

helpdlg(message,'AARAE info')












% *************************************************************************
% *************************************************************************
% GAIN CALIBRATION
% *************************************************************************
% *************************************************************************




% --- Executes on button press in load_btn.
function load_btn_Callback(hObject, ~, handles) %#ok : Executed when Load button is clicked
% hObject    handle to load_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'filtaudio'), handles = rmfield(handles,'filtaudio'); end
[filename,pathname] = uigetfile(...
    {'*.wav;*.mat;.*WAV;*.MAT','Calibration tone (*.wav,*.mat)'},...
    'Select audio file',[cd '/Audio']);

if filename ~= 0
    [~,~,ext] = fileparts(filename);
    % Check type of file. First 'if' is for .mat, second is for .wav
    if strcmp(ext,'.mat') || strcmp(ext,'.MAT')
        file = importdata(fullfile(pathname,filename));
        if isstruct(file)
            handles.audio = file.audio(:,1);
%            handles.fs = file.fs;
        else
%            fs = inputdlg('Please specify the sampling frequency','Sampling frequency',1);
%            if (isempty(specs))
%                warndlg('Input field is blank, cannot load data!');
%                handles.audio = [];
%            else
%                fs = str2num(specs{1,1});
%                if (isempty(fs) || fs<=0)
%                    warndlg('Input MUST be a real positive number, cannot load data!');
%                    handles.audio = [];
%                else
                    handles.audio = file(:,1);
%                    handles.fs = fs;
%                end
%            end
        end
    end
    if strcmp(ext,'.wav') || strcmp(ext,'.WAV')
        [handles.audio] = audioread(fullfile(pathname,filename));
        handles.audio = handles.audio(:,1);
    end;
    t = linspace(0,length(handles.audio)/handles.mainHandles.fs,length(handles.audio));
    plot(handles.dispaxes,t,handles.audio)
    xlabel(handles.dispaxes,'Time [s]');
    recrms = (mean(handles.audio.^2))^0.5/2;
    handles.calthresh = line([0 length(t)/handles.mainHandles.fs],[recrms recrms],'Color','k');
    set(handles.threshold_IN,'String',num2str(recrms*0.5))
    set([handles.sperframe_IN,handles.percentage_IN,handles.threshold_IN,...
        handles.tonelevel_IN,handles.units,handles.units_ref,...
        handles.units_type,handles.evalcal_btn,handles.filter_btn],'Enable','on')
else
    return
end
guidata(hObject,handles)


function sperframe_IN_Callback(~, ~, ~) %#ok : Executed when samples per frame input box changes
% hObject    handle to sperframe_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sperframe_IN as text
%        str2double(get(hObject,'String')) returns contents of sperframe_IN as a double


% --- Executes during object creation, after setting all properties.
function sperframe_IN_CreateFcn(hObject, ~, ~) %#ok : Samples per frame input box creation
% hObject    handle to sperframe_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function percentage_IN_Callback(hObject, ~, ~) %#ok : Executed when percentage input box changes
% hObject    handle to percentage_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of percentage_IN as text
%        str2double(get(hObject,'String')) returns contents of percentage_IN as a double
percentage = str2double(get(hObject,'String'));
if isnan(percentage) || percentage <=0 || percentage >= 100
    warndlg('Invalid entry. Percentage must be a scalar between 0 and 100.','AARAE info')
    set(hObject,'String','40')
end

% --- Executes during object creation, after setting all properties.
function percentage_IN_CreateFcn(hObject, ~, ~) %#ok : percentage input box creation
% hObject    handle to percentage_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_IN_Callback(hObject, ~, handles) %#ok : Executed when threshold input box changes 
% hObject    handle to threshold_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_IN as text
%        str2double(get(hObject,'String')) returns contents of threshold_IN as a double
delete(handles.calthresh)
Y = str2double(get(hObject,'String'));
handles.calthresh = line([0 length(handles.audio)/handles.mainHandles.fs],[Y Y],'Color','k');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function threshold_IN_CreateFcn(hObject, ~, ~) %#ok : Threshold input box creation
% hObject    handle to threshold_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in evalcal_btn.
function evalcal_btn_Callback(hObject, ~, handles) %#ok : Executed when Evaluate button is clicked on the gain calibration panel
% hObject    handle to evalcal_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');
pause on
pause(0.000001)
pause off
if isfield(handles,'filtaudio')
    if ~isempty(handles.filtaudio)
        data = handles.filtaudio;
    else
        data = handles.audio;
    end
else
    data = handles.audio;
end
hsr1 = dsp.SignalSource;
hsr1.Signal = data;
hsr1.SamplesPerFrame = str2double(get(handles.sperframe_IN,'String'));
percentage = str2double(get(handles.percentage_IN,'String'));
threshold = str2double(get(handles.threshold_IN,'String'));
tonelevel = str2double(get(handles.tonelevel_IN,'String'));
handles.output.units = get(handles.units,'String');
handles.output.units_ref = str2double(get(handles.units_ref,'String'));
handles.output.units_type = get(handles.units_type,'Value');
ncycles = ceil(length(data)/hsr1.SamplesPerFrame);
ms = zeros(ncycles*hsr1.SamplesPerFrame,1);
for i = 1:ncycles
    chunk = step(hsr1);
    ms((i-1)*hsr1.SamplesPerFrame+1:i*hsr1.SamplesPerFrame,:) = ones(size(chunk)).*mean(chunk.^2);
end
release(hsr1)
ms = ms(1:length(data));
rec = ones(size(ms));
rec(ms<threshold^2) = 0;
m = trimmean(ms(ms>threshold^2),percentage);
rec(ms>m*1.2) = 0;

t = linspace(0,length(data)/handles.mainHandles.fs,length(data));
trimdata = data.*rec(1:length(data));
plot(handles.dispaxes,t,data,'c')
xlabel(handles.dispaxes,'Time [s]');
YLim = get(handles.dispaxes,'YLim');
hold on
plot(handles.dispaxes,t,trimdata,'b')
plot(handles.dispaxes,t,(rec(1:length(data)).*10)-5,'Color','r','LineWidth',1)
set(handles.dispaxes,'YLim',YLim)
hold off
handles.calthresh = line([0 length(t)/handles.mainHandles.fs],[str2double(get(handles.threshold_IN,'String')) str2double(get(handles.threshold_IN,'String'))],'Color','k');
trim = find(trimdata);
dur = length(trim)/handles.mainHandles.fs;
trimlevel = tonelevel - 10 .* log10(mean(trimdata(trim).^2,1));
stats{1} = ['Sampling frequency: ' num2str(handles.mainHandles.fs) ' samples/s'];
stats{2} = ['Trimmed recording length: ' num2str(dur) ' s'];
stats{3} = ['Calibration offset: ' num2str(trimlevel) ' dB'];
if ~isnan(trimlevel), set(handles.caltonetext,'Visible','on'); end
set(handles.statstext,'String',stats);
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');
handles.output.cal(1,get(handles.channum_popup,'Value')) = trimlevel;
guidata(hObject,handles)


function tonelevel_IN_Callback(~, ~, ~) %#ok : Executed when level input box changes
% hObject    handle to tonelevel_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tonelevel_IN as text
%        str2double(get(hObject,'String')) returns contents of tonelevel_IN as a double


% --- Executes during object creation, after setting all properties.
function tonelevel_IN_CreateFcn(hObject, ~, ~) %#ok : level input box creation
% hObject    handle to tonelevel_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in filter_btn.
function filter_btn_Callback(hObject, ~, handles) %#ok : Executed when filter button is clicked
% hObject    handle to filter_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');
pause on
pause(0.000001)
pause off
filter = get(get(handles.filter_panel,'SelectedObject'),'String');
switch filter
    case 'None'
        handles.filtaudio = [];
    case '1 kHz'
        [handles.filtaudio,~] = thirdoctbandfilter(handles.audio,handles.mainHandles.fs,1000);
    case '250 Hz'
        [handles.filtaudio,~] = thirdoctbandfilter(handles.audio,handles.mainHandles.fs,250);
end
t = linspace(0,length(handles.audio)/handles.mainHandles.fs,length(handles.audio));
if ~isempty(handles.filtaudio)
    plot(handles.dispaxes,t,handles.filtaudio)
    xlabel(handles.dispaxes,'Time [s]');
else
    plot(handles.dispaxes,t,handles.audio)
    xlabel(handles.dispaxes,'Time [s]');
end
recrms = (mean(handles.audio.^2))^0.5;
handles.calthresh = line([0 length(t)/handles.mainHandles.fs],[recrms recrms],'Color','k');
set(handles.threshold_IN,'String',num2str(recrms*0.5))
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');
guidata(hObject,handles)



% --- Executes on button press in record_btn.
function record_btn_Callback(hObject, ~, handles) %#ok : Executed when clicking the Record button on the gain calibration panel
% hObject    handle to record_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'filtaudio'), handles = rmfield(handles,'filtaudio'); end
set([handles.sperframe_IN,handles.percentage_IN,handles.threshold_IN,...
    handles.tonelevel_IN,handles.units,handles.units_ref,handles.units_type,...
    handles.evalcal_btn,handles.filter_btn],'Enable','off')
dur = str2double(get(handles.duration_IN,'String'))*handles.mainHandles.fs;
set(hObject,'Enable','off');
set(handles.stop_btn,'Visible','on');
pause on
pause(0.000001)
pause off
% Set record object
set(handles.har,'NumChannels',get(handles.channum_popup,'Value'))
guidata(hObject,handles)
ncycles = ceil(dur/handles.har.SamplesPerFrame);
audio = zeros(ncycles*handles.har.SamplesPerFrame,get(handles.channum_popup,'Value'));
if isfield(handles,'filtaudio'), handles = rmfield(handles,'filtaudio'); end
% Initialize record routine
set(hObject,'BackgroundColor','red');
try
    UserData = get(handles.stop_btn,'UserData');
    for i = 1:ncycles
       UserData = get(handles.stop_btn,'UserData');
       if UserData.state == false
           audio((i-1)*handles.har.SamplesPerFrame+1:i*handles.har.SamplesPerFrame,:) = step(handles.har);
       else
           break
       end
       pause on
       pause(0.000001)
       pause off
    end
catch sthgwrong
    UserData.state = true;
    syswarning = sthgwrong.message;
    warndlg(syswarning,'AARAE info')
    return
end
rec = audio;
% Check recording and adjust for Duration
if ~isempty(rec)
    rec = rec(:,get(handles.channum_popup,'Value'));
    if UserData.state == false
        rec = rec(1:dur);
    else
        UserData.state = false;
        set(handles.stop_btn,'UserData',UserData);
    end
    handles.audio = rec;
    % Plot recording
    time = linspace(0,size(handles.audio,1)/handles.mainHandles.fs,length(handles.audio));
    plot(handles.dispaxes,time,handles.audio);
    xlabel(handles.dispaxes,'Time [s]');
    recrms = (mean(rec.^2))^0.5/2;
    handles.calthresh = line([0 length(time)/handles.mainHandles.fs],[recrms recrms],'Color','k');
    set(handles.threshold_IN,'String',num2str(recrms*0.5))
    set([handles.sperframe_IN,handles.percentage_IN,handles.threshold_IN,...
        handles.tonelevel_IN,handles.units,handles.units_ref,handles.units_type,...
        handles.evalcal_btn,handles.filter_btn],'Enable','on')
end
set(handles.record_btn,'BackgroundColor',[0.94 0.94 0.94]);
set(handles.record_btn,'Enable','on');
set(handles.stop_btn,'Visible','off');
% Release record object
release(handles.har)
guidata(hObject,handles);

% --- Executes on button press in stop_btn.
function stop_btn_Callback(hObject, ~, handles) %#ok : Executed when Stop button is clicked on the gain calibration panel
% hObject    handle to stop_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData.state = true;
set(hObject,'UserData',UserData);
pause on
pause(0.000001)
pause off
guidata(hObject,handles)



function duration_IN_Callback(hObject, ~, handles) %#ok : Executed when duration input box changes
% hObject    handle to duration_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration_IN as text
%        str2double(get(hObject,'String')) returns contents of duration_IN as a double
num = str2double(get(hObject,'String'));
if isnan(num) || num <= 0
    warndlg('Invalid entry','AARAE info')
    set(hObject,'String',num2str(handles.dur))
else
    handles.dur = num;
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function duration_IN_CreateFcn(hObject, ~, handles) %#ok : duration input box creation
% hObject    handle to duration_IN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dur = str2double(get(hObject,'String'));
guidata(hObject,handles)



function units_Callback(hObject, eventdata, handles)
% PHYSICAL UNITS
% hObject    handle to units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of units as text
%        str2double(get(hObject,'String')) returns contents of units as a double


% --- Executes during object creation, after setting all properties.
function units_CreateFcn(hObject, eventdata, handles)
% PHYSICAL UNITS
% hObject    handle to units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function units_ref_Callback(hObject, eventdata, handles)
% REFERENCE VALUE FOR dB
% hObject    handle to units_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of units_ref as text
%        str2double(get(hObject,'String')) returns contents of units_ref as a double


% --- Executes during object creation, after setting all properties.
function units_ref_CreateFcn(hObject, eventdata, handles)
% REFERENCE VALUE FOR dB
% hObject    handle to units_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in units_type.
function units_type_Callback(hObject, eventdata, handles)
% AMPLITUDE OR POWER UNITS
% hObject    handle to units_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns units_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from units_type


% --- Executes during object creation, after setting all properties.
function units_type_CreateFcn(hObject, eventdata, handles)
% AMPLITUDE OR POWER UNITS
% hObject    handle to units_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in channum_popup.
function channum_popup_Callback(~, ~, ~) %#ok : Executes when channel menu changes in the gain calibration panel
% hObject    handle to channum_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channum_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channum_popup


% --- Executes during object creation, after setting all properties.
function channum_popup_CreateFcn(hObject, ~, ~) %#ok : Channel selection menu creation
% hObject    handle to channum_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calfield_btn.
function calfield_btn_Callback(~, ~, handles) %#ok : Executed when cal field values button is clicked
% hObject    handle to calfield_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
deviceinfo = dspAudioDeviceInfo('inputs');
inputsel = get(handles.mainHandles.inputdev_popup,'Value');
msgbox(num2str([(1:deviceinfo(inputsel).maxInputs);handles.output.cal]'))



% --- Executes on button press in calinstructions_btn.
function calinstructions_btn_Callback(~, ~, ~) %#ok : Executed when instruction button is clicked on the gain calibration panel
% hObject    handle to calinstructions_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
message{1} = 'Instructions:';
message{2} = '1. You may record or load the calibration tone before or after evaluating the system latency and generating the corresponding inverse filter (if required).';
message{3} = '2. Click on the -Record- button to capture the calibration tone or alternatively click on -Load from file- to select a prerecorded calibration tone from your computer.';
message{4} = '3. When recording you may choose a set recording time in seconds or stop the recording when you consider you have a long enough recording. Make sure you set enough time to capture the calibration tone';
message{5} = '4. The automated trimming of your calibartion tone requires you to set a frame length in samples to evaluate the calibration tone. This is done by evaluating the recording envelope';
message{6} = '5. The trim mean ratio acts as a filter referenced by the envelope threshold.';
message{7} = '6. Set envelope threshold based on the absolute amplitude value of the recorded or loaded audio.';
message{8} = '7. If you are using a microphone calibrator, provide the calibrator''s sound pressure level in dB set by the manufaturer. Otherwise enter the known level in dB.';
message{9} = '8. Modify the units fields if appropriate (the default reference of 2e-5 Pa should be appropriate for sound pressure level).';
message{10} = ' ';
message{11} = '9. Click on -Filter- to apply an octave band filter around 1 kHz or 250 Hz if required to clean the recording.';
message{12} = '10. Click on -Evaluate- to capture the calibration tone level and usable recording length.';
helpdlg(message,'AARAE info')