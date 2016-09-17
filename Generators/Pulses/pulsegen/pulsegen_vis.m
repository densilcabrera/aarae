function varargout = pulsegen_vis(varargin)
% PULSEGEN_VIS MATLAB code for pulsegen_vis.fig
%      PULSEGEN_VIS, by itself, creates a new PULSEGEN_VIS or raises the existing
%      singleton*.
%
%      H = PULSEGEN_VIS returns the handle to a new PULSEGEN_VIS or the handle to
%      the existing singleton*.
%
%      PULSEGEN_VIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PULSEGEN_VIS.M with the given input arguments.
%
%      PULSEGEN_VIS('Property','Value',...) creates a new PULSEGEN_VIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pulsegen_vis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pulsegen_vis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pulsegen_vis

% Last Modified by GUIDE v2.5 29-Nov-2010 19:46:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pulsegen_vis_OpeningFcn, ...
    'gui_OutputFcn',  @pulsegen_vis_OutputFcn, ...
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


% --- Executes just before pulsegen_vis is made visible.
function pulsegen_vis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pulsegen_vis (see VARARGIN)

% Choose default command line output for pulsegen_vis
handles.output = hObject;

handles.edge_handles=[handles.edge_label,handles.edge_slider,handles.edge_edit];
handles.f1_handles=[handles.f1_label,handles.f1_slider,handles.f1_edit,handles.f1_units];
handles.f2_handles=[handles.f2_label,handles.f2_slider,handles.f2_edit,handles.f2_units];
handles.param_handles=[handles.param_label,handles.param_slider,handles.param_edit];
% Update handles structure
guidata(hObject, handles);
update_plot(handles);
% UIWAIT makes pulsegen_vis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pulsegen_vis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pulse_type_menu.
function pulse_type_menu_Callback(hObject, eventdata, handles)
% hObject    handle to pulse_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pulse_type_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pulse_type_menu

if (get(handles.auto_update,'Value'))
    update_plot(handles)
end
val=get(hObject,'Value');
switch(val)
    case {1,2,3,4,5,6,7,8,11}
        set([handles.param_handles,handles.f1_handles,handles.f2_handles],'Visible','off');
        set(handles.edge_handles,'Visible','on');
    case 9
        set([handles.f1_handles,handles.f2_handles],'Visible','off');
        set([handles.param_handles,handles.edge_handles],'Visible','on');
    case {10, 12}
        set([handles.f1_handles,handles.f2_handles],'Visible','on');
        set([handles.param_handles,handles.edge_handles],'Visible','off');
end


% --- Executes during object creation, after setting all properties.
function edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_plot(handles);

% --- Executes on button press in export_button.
function export_button_Callback(hObject, eventdata, handles)
% hObject    handle to export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fs=get(handles.fs_edit,'Value');
T=get(handles.T_edit,'Value');
edge=get(handles.edge_edit,'Value');
f1=get(handles.f1_edit,'Value');
f2=get(handles.f2_edit,'Value');
param=get(handles.param_edit,'Value');
win=get(handles.window_edit,'Value');
modulation=get(handles.modulation_edit,'Value');
type_str=get(handles.pulse_type_menu,'String');
type_val=get(handles.pulse_type_menu,'Value');
type=type_str{type_val};
lp=get(handles.low_pass_edit,'Value');
hp=get(handles.high_pass_edit,'Value');
dispersion=get(handles.dispersion_edit,'Value');

p=pulsegen(Fs,T,edge,type,'start_frequency',f1,'stop_frequency',f2,'args',param,'window',win,'modulation',modulation,'low_pass',lp,...
    'high_pass',hp,'dispersion',dispersion);
t=-T/2:1/Fs:T/2;
switch(type_val)
    case {1,2,3,4,5,6,7,8,11}
str=sprintf('y=pulsegen(%f,%f,%f,''%s'',''window'',%f,''modulation'',%f,''low_pass'',%f,''high_pass'',%f,''dispersion'',%f);',...
    Fs,T,edge,type,win,modulation,lp,hp,dispersion);
    case 9
str=sprintf('y=pulsegen(%f,%f,%f,''%s'',''window'',%f,''modulation'',%f,''low_pass'',%f,''high_pass'',%f,''dispersion'',%f,''args'',%f);',...
    Fs,T,edge,type,win,modulation,lp,hp,dispersion,param);
    case {10,12}
str=sprintf('y=pulsegen(%f,%f,%f,''%s'',''window'',%f,''modulation'',%f,''low_pass'',%f,''high_pass'',%f,''dispersion'',%f,''low_pass'',%f,''high_pass'',%f);',...
    Fs,T,0,type,win,modulation,lp,hp,dispersion,lp,hp);   
end
assignin('base','t',t);
fprintf('%s\n',str);
evalin('base',str);
% --- Executes during object creation, after setting all properties.
function def_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));


function update_plot(handles)

Fs=get(handles.fs_edit,'Value');
T=get(handles.T_edit,'Value');
edge=get(handles.edge_edit,'Value');
f1=get(handles.f1_edit,'Value');
f2=get(handles.f2_edit,'Value');
param=get(handles.param_edit,'Value');
win=get(handles.window_edit,'Value');
modulation=get(handles.modulation_edit,'Value');
type_str=get(handles.pulse_type_menu,'String');
type_val=get(handles.pulse_type_menu,'Value');
type=type_str{type_val};
lp=get(handles.low_pass_edit,'Value');
hp=get(handles.high_pass_edit,'Value');
dispersion=get(handles.dispersion_edit,'Value');

p=pulsegen(Fs,T,edge,type,'start_frequency',f1,'stop_frequency',f2,'args',param,'window',win,'modulation',modulation,'low_pass',lp,...
    'high_pass',hp,'dispersion',dispersion);
t=-T/2:1/Fs:T/2;

pnum=get(handles.plot_menu,'Value');
if (pnum==1)
    plot(handles.axes,t,p);
    xlabel('Time(s)');
else
    lp=length(p);
    Y=fft(p);
    
    sig=abs(Y(1:ceil(lp/2)));
    f=linspace(0,Fs/2,ceil(lp/2));
    semilogy(handles.axes,f(sig>1e-4),sig(sig>1e-4));
    xlabel('Frequency (Hz)');
end

% --- Executes on selection change in plot_menu.
function plot_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_menu
update_plot(handles);


% --- Executes on button press in auto_update.
function auto_update_Callback(hObject, eventdata, handles)
% hObject    handle to auto_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_update
val=get(hObject,'Value');
if (val)
    set(handles.plot_button,'Visible','off');
else
    set(handles.plot_button,'Visible','on');
end




% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to high_pass_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
tag=get(hObject,'tag');
etag=[tag(1:end-6),'edit'];
set(handles.(etag),'String',num2str(get(hObject,'Value')));
set(handles.(etag),'Value',get(hObject,'Value'));
if (get(handles.auto_update,'Value'))
    update_plot(handles)
end



function edit_Callback(hObject, eventdata, handles)
% hObject    handle to f1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1_edit as text
%        str2double(get(hObject,'String')) returns contents of f1_edit as a double
cval=get(hObject,'String');
val=str2double(cval);
tag=get(hObject,'tag');
stag=[tag(1:end-4),'slider'];
if (~isnan(val)&&(~isempty(val)))
   
    set(hObject,'Value',val);
    set(handles.(stag),'Value',val);
    if (get(handles.auto_update,'Value'))
        update_plot(handles)
    end

else
    set(hObject,'String',num2str(get(hObject,'Value')));
end



function fs_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_edit as text
%        str2double(get(hObject,'String')) returns contents of fs_edit as a double
cval=get(hObject,'String');
val=str2double(cval);
if (~isnan(val)&&(~isempty(val)))
   
    set(hObject,'Value',val);
    set([handles.f1_slider,handles.f2_slider,handles.modulation],'Max',val/2);
    if (get(handles.auto_update,'Value'))
        update_plot(handles)
    end

else
    set(hObject,'String',num2str(get(hObject,'Value')));
end


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edge_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
au=get(handles.auto_update,'Value');
set(handles.auto_update,'Value',0);
set(handles.edge_slider,'Value',1);
slider_Callback(handles.edge_slider,eventdata,handles);
set(handles.window_slider,'Value',0);
slider_Callback(handles.window_slider,eventdata,handles);
set(handles.modulation_slider,'Value',0);
slider_Callback(handles.modulation_slider,eventdata,handles);
set(handles.f1_slider,'Value',10);
slider_Callback(handles.f1_slider,eventdata,handles);
set(handles.f2_slider,'Value',20);
slider_Callback(handles.f2_slider,eventdata,handles);
set(handles.param_slider,'Value',0.5);
slider_Callback(handles.param_slider,eventdata,handles);
set(handles.low_pass_slider,'Value',1);
slider_Callback(handles.low_pass_slider,eventdata,handles);
set(handles.high_pass_slider,'Value',0);
slider_Callback(handles.high_pass_slider,eventdata,handles);
set(handles.dispersion_slider,'Value',0);
slider_Callback(handles.dispersion_slider,eventdata,handles);
set(handles.fs_edit,'String','48000');
set(handles.fs_edit,'Value',48000);
set(handles.auto_update,'Value',au);
update_plot(handles);
