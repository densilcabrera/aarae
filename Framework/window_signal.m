% DO NOT EDIT THIS INITIALIZATION FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!
function varargout = window_signal(varargin)
% WINDOW_SIGNAL MATLAB code for window_signal.fig
%      WINDOW_SIGNAL, by itself, creates a new WINDOW_SIGNAL or raises the existing
%      singleton*.
%
%      H = WINDOW_SIGNAL returns the handle to a new WINDOW_SIGNAL or the handle to
%      the existing singleton*.
%
%      WINDOW_SIGNAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINDOW_SIGNAL.M with the given input arguments.
%
%      WINDOW_SIGNAL('Property','Value',...) creates a new WINDOW_SIGNAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before window_signal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to window_signal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help window_signal

% Last Modified by GUIDE v2.5 10-Aug-2015 18:05:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @window_signal_OpeningFcn, ...
    'gui_OutputFcn',  @window_signal_OutputFcn, ...
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


% --- Executes just before window_signal is made visible.
function window_signal_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to window_signal (see VARARGIN)

axis([0 10 -1 1]);
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
end

if dontOpen
    disp('-----------------------------------------------------');
    disp('This function is part of the AARAE framework, it is')
    disp('not a standalone function. To call this function,')
    disp('click on the appropriate calling button on the main');
    disp('Window. E.g.:');
    disp('   Convolve with audio2');
    disp('-----------------------------------------------------');
else
    hMain = getappdata(0,'hMain');
    % Find the IR signal being sent to window
    impulse = find(strcmp(varargin, 'IR'));
    handles.IR = varargin{impulse+1};
    fs = find(strcmp(varargin, 'fs'));
    handles.fs = varargin{fs+1};
    chans = find(strcmp(varargin, 'chans'));
    handles.chans = 1:varargin{chans+1};
    set(handles.ChannelsTextBox,'String',num2str(handles.chans))
    bands = find(strcmp(varargin, 'bands'));
    handles.bands = 1:varargin{bands+1};
    set(handles.BandsTextBox,'String',num2str(handles.bands))
    cycles = find(strcmp(varargin, 'cycles'));
    handles.cycles = 1:varargin{cycles+1};
    set(handles.CyclesTextBox,'String',num2str(handles.cycles))
    outchans = find(strcmp(varargin, 'outchans'));
    handles.outchans = 1:varargin{outchans+1};
    set(handles.OutputChansTextbox,'String',num2str(handles.outchans))
    dim6 = find(strcmp(varargin, 'dim6'));
    handles.dim6 = 1:varargin{dim6+1};
    set(handles.Dim6TextBox,'String',num2str(handles.dim6))
    audio2len = find(strcmp(varargin, 'audio2len'));
    handles.audio2len = varargin{audio2len+1};
    t = linspace(0,size(handles.IR,1),size(handles.IR,1));
    plot(handles.IN_axes,t,10*log10(handles.IR.^2),'Color',[0 0.7 0])
    xlabel(handles.IN_axes,'Samples');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
    [~, id] = max(abs(handles.IR));
    IRlength = max(id);
    %handles.IRlength = IRlength;
    %set(handles.IN_length, 'string',num2str(max(id))); % Get the IR length from input
    set(handles.trimmethod_popup,'Value',getappdata(hMain,'trim_method_after_convolution'));
    
    trimmethod = get(handles.trimmethod_popup,'Value');
    
    
    
    
    switch trimmethod
        case 1
            % 1. Symetrically trim to half length around peak
            trimsamp_low = max(id)-round(IRlength./2);
            trimsamp_high = trimsamp_low + IRlength -1;
        case 2
            % 2. Trim from 10 ms before peak to end
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = size(handles.IR,1);
        case 3
            % 3. Trim from 10 ms before peak, limited to 8 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*8;
        case 4
            % 4. Trim from 10 ms before peak, limited to 4 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*4;
        case 5
            % 5. Trim from 10 ms before peak, limited to 2 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*2;
        case 6
            % 6. Trim from 10 ms before peak, limited to 1 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs;
        case 7
            % 7. Trim from 1 ms before peak to end
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = size(handles.IR,1);
        case 8
            % 8. Trim from 1 ms before peak, limited to 8 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*8;
        case 9
            % 9. Trim from 1 ms before peak, limited to 4 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*4;
        case 10
            % 10. Trim from 1 ms before peak, limited to 2 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*2;
        case 11
            % 11. Trim from 1 ms before peak, limited to 1 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs;
        case 12
            % 12. Causal part to end
            trimsamp_low = handles.audio2len;
            trimsamp_high = size(handles.IR,1);
        case 13
            % 13. Causal part, limited to 8 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs*8-1;
        case 14
            % 14. Causal part, limited to 4 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs*4-1;
        case 15
            % 15. Causal part, limited to 2 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs*2-1;
        case 16
            % 16. Causal part, limited to 1 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs-1;
        case 17
            % 17. Causal part, limited to 0.5 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + round(handles.fs/2)-1;
        case 18
            % 18. Causal part, limited to 0.25 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + round(handles.fs/4)-1;
        case 19
            % 19. Causal part, limited to 0.125 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + round(handles.fs/8)-1;
        otherwise
            % 20. Do not trim
            trimsamp_low = 1;
            trimsamp_high = size(handles.IR,1);
    end
    if trimsamp_low < 1, trimsamp_low = 1; end
    if trimsamp_high > size(handles.IR,1), trimsamp_high = size(handles.IR,1); end
    if trimsamp_high <= trimsamp_low % something went wrong!
        trimsamp_low = 1;
        trimsamp_high = size(handles.IR,1);
    end
    handles.slow = trimsamp_low;
    handles.shigh = trimsamp_high; 
    set(handles.trimlow,'String',num2str(trimsamp_low))
    set(handles.trimhigh,'String',num2str(trimsamp_high))
    L = 10*log10(handles.IR(:).^2);
    minL = min(L);
    maxL = max(L);
    % we get an error here if handles.IR is very short ( e.g. 1 sample
    % long). Probably need to ensure that it is at least 5 samples long.
    B=interp1([0 trimsamp_low trimsamp_low+1 trimsamp_high-2 trimsamp_high-1 t(end)],...
        [minL minL maxL maxL minL minL],t,'linear');
    hold(handles.IN_axes,'on')
    handles.win = plot(handles.IN_axes,B,'r');
    hold(handles.IN_axes,'off')
    trimIR = handles.IR(trimsamp_low:trimsamp_high,:); % Crop IR
    plot(handles.OUT_axes,10*log10(trimIR.^2),'Color',[0 0.7 0]) % Plot cropped IR
    xlabel(handles.OUT_axes,'Samples');
    set(handles.OUT_axes,'XTickLabel',num2str(get(handles.OUT_axes,'XTick').'))
    guidata(hObject, handles);
    uiwait(hObject);
end

% UIWAIT makes window_signal wait for user response (see UIRESUME)
% uiwait(handles.window_signal);


% --- Outputs from this function are returned to the command line.
function varargout = window_signal_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.slow;
varargout{2} = handles.shigh;
varargout{3} = handles.chans;
varargout{4} = handles.bands;
varargout{5} = handles.cycles;
varargout{6} = handles.outchans;
varargout{7} = handles.dim6;
delete(hObject);


% --- Executes on button press in done_btn.
function done_btn_Callback(~, ~, handles) %#ok : Executed when Done button is clicked
% hObject    handle to done_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.window_signal);


% --- Executes when user attempts to close window_signal.
function window_signal_CloseRequestFcn(hObject, ~, ~) %#ok : Executed upon close request function on window
% hObject    handle to window_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(hObject);


function trimlow_Callback(hObject, ~, handles) %#ok : Executed when lower limit input box changes
% hObject    handle to trimlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trimlow as text
%        str2double(get(hObject,'String')) returns contents of trimlow as a double

trimsamp_low = round(str2double(get(hObject,'String')));
% Check user's input
if isnan(trimsamp_low) || trimsamp_low <= 0
    set(hObject,'String',num2str(handles.slow))
    warndlg('Invalid lower limit!','AARAE info');
else
    delete(handles.win)
    handles = rmfield(handles,'win');
    trimsamp_high = round(str2double(get(handles.trimhigh,'String')));
    handles.slow = trimsamp_low;
    t = linspace(0,size(handles.IR,1),size(handles.IR,1));
    L = 10*log10(handles.IR(:).^2);
    minL = min(L);
    maxL = max(L);
    B=interp1([0 trimsamp_low trimsamp_low+1 trimsamp_high-2 trimsamp_high-1 t(end)],...
        [minL minL maxL maxL minL minL],t,'linear');
    hold(handles.IN_axes,'on')
    handles.win = plot(handles.IN_axes,B,'r');
    hold(handles.IN_axes,'off')
    trimIR = handles.IR(trimsamp_low:trimsamp_high,:);
    plot(handles.OUT_axes,10*log10(trimIR.^2))
    xlabel(handles.OUT_axes,'Samples');
    set(handles.OUT_axes,'XTickLabel',num2str(get(handles.OUT_axes,'XTick').'))
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function trimlow_CreateFcn(hObject, ~, ~) %#ok : Lower limit input box creation
% hObject    handle to trimlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trimhigh_Callback(hObject, ~, handles) %#ok : Executed when higher limit input box changes
% hObject    handle to trimhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trimhigh as text
%        str2double(get(hObject,'String')) returns contents of trimhigh as a double

trimsamp_high = round(str2double(get(hObject,'String')));
% Check user's input
if isnan(trimsamp_high) || trimsamp_high >= length(handles.IR)
    set(hObject,'String',num2str(handles.shigh))
    warndlg('Invalid upper limit!','AARAE info');
else
    delete(handles.win)
    handles = rmfield(handles,'win');
    trimsamp_low = round(str2double(get(handles.trimlow,'String')));
    handles.shigh = trimsamp_high;
    t = linspace(0,size(handles.IR,1),size(handles.IR,1));
    L = 10*log10(handles.IR(:).^2);
    minL = min(L);
    maxL = max(L);
    B=interp1([0 trimsamp_low trimsamp_low+1 trimsamp_high-2 trimsamp_high-1 t(end)],...
        [minL minL maxL maxL minL minL],t,'linear');
    hold(handles.IN_axes,'on')
    handles.win = plot(handles.IN_axes,B,'r');
    hold(handles.IN_axes,'off')
    trimIR = handles.IR(trimsamp_low:trimsamp_high,:);
    plot(handles.OUT_axes,10*log10(trimIR.^2))
    xlabel(handles.OUT_axes,'Samples');
    set(handles.OUT_axes,'XTickLabel',num2str(get(handles.OUT_axes,'XTick').'))
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function trimhigh_CreateFcn(hObject, ~, ~) %#ok : Higher limit input box creation
% hObject    handle to trimhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trimmethod_popup.
function trimmethod_popup_Callback(hObject, eventdata, handles)
% hObject    handle to trimmethod_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trimmethod_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trimmethod_popup
hMain = getappdata(0,'hMain');
setappdata(hMain,'trim_method_after_convolution',get(hObject,'Value'));

trimmethod = get(handles.trimmethod_popup,'Value');

[~, id] = max(abs(handles.IR));
IRlength = max(id);




    switch trimmethod
        case 1
            % 1. Symetrically trim to half length around peak
            trimsamp_low = max(id)-round(IRlength./2);
            trimsamp_high = trimsamp_low + IRlength -1;
        case 2
            % 2. Trim from 10 ms before peak to end
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = size(handles.IR,1);
        case 3
            % 3. Trim from 10 ms before peak, limited to 8 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*8;
        case 4
            % 4. Trim from 10 ms before peak, limited to 4 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*4;
        case 5
            % 5. Trim from 10 ms before peak, limited to 2 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*2;
        case 6
            % 6. Trim from 10 ms before peak, limited to 1 s duration
            trimpad = round(handles.fs*0.01);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs;
        case 7
            % 7. Trim from 1 ms before peak to end
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = size(handles.IR,1);
        case 8
            % 8. Trim from 1 ms before peak, limited to 8 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*8;
        case 9
            % 9. Trim from 1 ms before peak, limited to 4 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*4;
        case 10
            % 10. Trim from 1 ms before peak, limited to 2 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs*2;
        case 11
            % 11. Trim from 1 ms before peak, limited to 1 s duration
            trimpad = round(handles.fs*0.001);
            trimsamp_low = max(id)-trimpad;
            trimsamp_high = trimsamp_low + trimpad + handles.fs;
        case 12
            % 12. Causal part to end
            trimsamp_low = handles.audio2len;
            trimsamp_high = size(handles.IR,1);
        case 13
            % 13. Causal part, limited to 8 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs*8-1;
        case 14
            % 14. Causal part, limited to 4 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs*4-1;
        case 15
            % 15. Causal part, limited to 2 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs*2-1;
        case 16
            % 16. Causal part, limited to 1 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + handles.fs-1;
        case 17
            % 17. Causal part, limited to 0.5 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + round(handles.fs/2)-1;
        case 18
            % 18. Causal part, limited to 0.25 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + round(handles.fs/4)-1;
        case 19
            % 19. Causal part, limited to 0.125 s duration
            trimsamp_low = handles.audio2len;
            trimsamp_high = trimsamp_low + round(handles.fs/8)-1;
        otherwise
            % 20. Do not trim
            trimsamp_low = 1;
            trimsamp_high = size(handles.IR,1);
    end
    if trimsamp_low < 1, trimsamp_low = 1; end
    if trimsamp_high > size(handles.IR,1), trimsamp_high = size(handles.IR,1); end
    if trimsamp_high <= trimsamp_low % something went wrong!
        trimsamp_low = 1;
        trimsamp_high = size(handles.IR,1);
    end

handles.slow = trimsamp_low;
handles.shigh = trimsamp_high;
delete(handles.win)
handles = rmfield(handles,'win');
% trimsamp_low = round(str2double(get(handles.trimlow,'String')));
% handles.shigh = trimsamp_high;
t = linspace(0,size(handles.IR,1),size(handles.IR,1));
L = 10*log10(handles.IR(:).^2);
minL = min(L);
maxL = max(L);
B=interp1([0 trimsamp_low trimsamp_low+1 trimsamp_high-2 trimsamp_high-1 t(end)],...
    [minL minL maxL maxL minL minL],t,'linear');
hold(handles.IN_axes,'on')
handles.win = plot(handles.IN_axes,B,'r');
hold(handles.IN_axes,'off')
trimIR = handles.IR(trimsamp_low:trimsamp_high,:);
plot(handles.OUT_axes,10*log10(trimIR.^2),'Color',[0 0.7 0])
xlabel(handles.OUT_axes,'Samples');
set(handles.OUT_axes,'XTickLabel',num2str(get(handles.OUT_axes,'XTick').'))
set(handles.trimlow,'String',num2str(trimsamp_low))
set(handles.trimhigh,'String',num2str(trimsamp_high))
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function trimmethod_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trimmethod_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ChannelsTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to ChannelsTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChannelsTextBox as text
%        str2double(get(hObject,'String')) returns contents of ChannelsTextBox as a double
input = round(str2num(get(hObject,'String')));
input = unique(input(input>=1 & input<=max(handles.chans) & ~isnan(input)));
if isempty(input)
    set(hObject,'String',['1:' num2str(handles.chans)])
    warndlg('Invalid entry!','AARAE info');
else
    handles.chans = input;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ChannelsTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelsTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BandsTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to BandsTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BandsTextBox as text
%        str2double(get(hObject,'String')) returns contents of BandsTextBox as a double
input = round(str2num(get(hObject,'String')));
input = unique(input(input>=1 & input<=max(handles.bands) & ~isnan(input)));
if isempty(input)
    set(hObject,'String',['1:' num2str(handles.bands)])
    warndlg('Invalid entry!','AARAE info');
else
    handles.bands = input;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function BandsTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BandsTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CyclesTextBox_Callback(hObject, eventdata, handles)
% hObject    handle to CyclesTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CyclesTextBox as text
%        str2double(get(hObject,'String')) returns contents of CyclesTextBox as a double
input = round(str2num(get(hObject,'String')));
input = unique(input(input>=1 & input<=max(handles.cycles) & ~isnan(input)));
if isempty(input)
    set(hObject,'String',['1:' num2str(handles.cycles)])
    warndlg('Invalid entry!','AARAE info');
else
    handles.cycles = input;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function CyclesTextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CyclesTextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OutputChansTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to OutputChansTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputChansTextbox as text
%        str2double(get(hObject,'String')) returns contents of OutputChansTextbox as a double
input = round(str2num(get(hObject,'String')));
input = unique(input(input>=1 & input<=max(handles.outchans) & ~isnan(input)));
if isempty(input)
    set(hObject,'String',['1:' num2str(handles.outchans)])
    warndlg('Invalid entry!','AARAE info');
else
    handles.outchans = input;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function OutputChansTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputChansTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dim6TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to Dim6TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dim6TextBox as text
%        str2double(get(hObject,'String')) returns contents of Dim6TextBox as a double
input = round(str2num(get(hObject,'String')));
input = unique(input(input>=1 & input<=max(handles.dim6) & ~isnan(input)));
if isempty(input)
    set(hObject,'String',['1:' num2str(handles.dim6)])
    warndlg('Invalid entry!','AARAE info');
else
    handles.dim6 = input;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function Dim6TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dim6TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Cancel_pushbutton.
function Cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.slow = [];
guidata(hObject, handles);
uiresume(handles.window_signal);
