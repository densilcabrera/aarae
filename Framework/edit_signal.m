% DO NOT EDIT THIS INITIALIZATION FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!
function varargout = edit_signal(varargin)
% EDIT_SIGNAL MATLAB code for edit_signal.fig
%      EDIT_SIGNAL, by itself, creates a new EDIT_SIGNAL or raises the existing
%      singleton*.
%
%      H = EDIT_SIGNAL returns the handle to a new EDIT_SIGNAL or the handle to
%      the existing singleton*.
%
%      EDIT_SIGNAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDIT_SIGNAL.M with the given input arguments.
%
%      EDIT_SIGNAL('Property','Value',...) creates a new EDIT_SIGNAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before edit_signal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to edit_signal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help edit_signal

% Last Modified by GUIDE v2.5 10-Sep-2014 14:29:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @edit_signal_OpeningFcn, ...
                   'gui_OutputFcn',  @edit_signal_OutputFcn, ...
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


% --- Executes just before edit_signal is made visible.
function edit_signal_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to edit_signal (see VARARGIN)

% Choose default command line output for edit_signal
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
    mainHandles = guidata(handles.main_stage1);
    handles.fid = mainHandles.fid;

    fontsize

end



if dontOpen
   disp('-----------------------------------------------------');
   disp('This function is part of the AARAE framework, it is') 
   disp('not a standalone function. To call this function,')
   disp('click on the appropriate calling button on the main');
   disp('Window. E.g.:');
   disp('   Edit');
   disp('-----------------------------------------------------');
else
    % Call the 'desktop'
    hMain = getappdata(0,'hMain');
    handles.version = 1;
    handles.testsignal(handles.version) = getappdata(hMain,'testsignal');
    
    % *********************************************************************    
    % Note that it is not possible for a basic processor to add a field
    % when used from the Edit GUI, so it is necessary to add the following
    % fields if missing (with null values) so that they can be edited:
    % * cal
    % * properties.units, properties.units_ref, properties.units_type
    % * another possible candidate could be chanID (not implemented) or any
    % other field that a basic processor works with.
    % It would be possible to automatically remove these fields from the
    % output of the function if they still have null values (not
    % implemented).
    if ~isfield(handles.testsignal(handles.version),'cal')
        handles.testsignal(handles.version).cal = ...
            nan(1,size(handles.testsignal(handles.version).audio,2));
    end
    if ~isfield(handles.testsignal(handles.version),'properties')
        handles.testsignal(handles.version).properties.units = '';
        handles.testsignal(handles.version).properties.units_ref = 1;
        handles.testsignal(handles.version).properties.units_type = 1;
    else
        if ~isfield(handles.testsignal(handles.version).properties,'units')
            handles.testsignal(handles.version).properties.units = '';
        end
        if ~isfield(handles.testsignal(handles.version).properties,'units_ref')
            handles.testsignal(handles.version).properties.units_ref = 1;
        end
        if ~isfield(handles.testsignal(handles.version).properties,'units_type')
            handles.testsignal(handles.version).properties.units_type = 1;
        end
    end
    % *********************************************************************
    
    % Bring up the data from the selected leaf to be edited
    audiodata = handles.testsignal(handles.version);
    mainHandles = guidata(handles.main_stage1);
    selectedNodes = mainHandles.mytree.getSelectedNodes;
    handles.selNodeName = selectedNodes(1).getName.char;
    set(handles.IN_name,'String',handles.selNodeName)
    audiodatatext = evalc('audiodata');
    set(handles.audiodatatext,'String',audiodatatext);
    handles.fs = audiodata.fs;
    dur = length(handles.testsignal(handles.version).audio)/handles.fs;
    % Allocate memory space for the edited signal
    handles.rel_time = linspace(0,dur,length(handles.testsignal(handles.version).audio));
    handles.xi(handles.version) = min(handles.rel_time);
    handles.xf(handles.version) = max(handles.rel_time);
    set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
    set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    % Plot signal to be cropped
    if ~ismatrix(handles.testsignal(handles.version).audio)
        set(handles.channel_panel,'Visible','on');
        set(handles.IN_nchannel,'String','1');
        set(handles.tchannels,'String',['/ ' num2str(size(handles.testsignal(handles.version).audio,2))]);
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
        if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
        if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
        if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    else
        set(handles.channel_panel,'Visible','off');
        linea = handles.testsignal(handles.version).audio;
        cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea)
    xlabel(handles.IN_axes,'Time');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
    % Update handles structure
    guidata(hObject, handles);
    % UIWAIT makes edit_signal wait for user response (see UIRESUME)
%    if ndims(audiodata.audio) <= 4
        uiwait(hObject)
%    else
%        warndlg('Edition of 4-Dimensional audio or greater not yet enabled, sorry!','AARAE info','modal')
%    end
end



% --- Outputs from this function are returned to the command line.
function varargout = edit_signal_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.xi(handles.version);
varargout{2} = handles.xf(handles.version);
delete(hObject);


% --- Executes on button press in oo_btn.
function oo_btn_Callback(hObject, ~, handles) %#ok : Executed when Overwrite Original button is clicked
% hObject    handle to oo_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the main window handles (Consider removing)
mainHandles = guidata(handles.main_stage1);

% Check if there's a chunk or not
if ~isempty(handles.testsignal(handles.version))
    aarae_fig = findobj('type','figure','tag','aarae');
    selectedNodes = mainHandles.mytree.getSelectedNodes;
    removefield = matlab.lang.makeValidName(selectedNodes(1).getName.char);
    
    % Save as you go
    if exist([cd '/Utilities/Backup/' selectedNodes(1).getName.char '.mat'],'file')
        delete([cd '/Utilities/Backup/' selectedNodes(1).getName.char '.mat'])
    end
    signaldata = handles.testsignal(handles.version);
    signaldata.name = handles.selNodeName;
    save([cd '/Utilities/Backup/' handles.selNodeName '.mat'], 'signaldata');
    
    set(mainHandles.(matlab.lang.makeValidName(removefield)),'Name',handles.selNodeName);
    set(mainHandles.(matlab.lang.makeValidName(removefield)),'Value',handles.selNodeName);
    mainHandles.(matlab.lang.makeValidName(handles.selNodeName)) = mainHandles.(matlab.lang.makeValidName(removefield));
    if ~strcmp(selectedNodes(1).getName.char,handles.selNodeName), mainHandles = rmfield(mainHandles,removefield); end
    mainHandles.(matlab.lang.makeValidName(handles.selNodeName)).UserData = handles.testsignal(handles.version);
    mainHandles.mytree.reloadNode(mainHandles.(matlab.lang.makeValidName(handles.selNodeName)).getParent);
    mainHandles.mytree.setSelectedNode(mainHandles.(matlab.lang.makeValidName(handles.selNodeName)));
    guidata(aarae_fig, mainHandles);
end
guidata(hObject,handles);
uiresume(handles.edit_signal);


% --- Executes on button press in cancel_btn.
function cancel_btn_Callback(~, ~, handles) %#ok : Executed when Cancel button is clicked
% hObject    handle to cancel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.OUT_start,'String','-');
set(handles.OUT_end,'String','-');
fprintf(handles.fid, '%% CANCEL EDIT\n');
uiresume(handles.edit_signal);


% --- Executes when user attempts to close edit_signal.
function edit_signal_CloseRequestFcn(hObject, ~, ~) %#ok : Executed when window is closed
% hObject    handle to edit_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(hObject);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function edit_signal_WindowButtonDownFcn(hObject, ~, handles) %#ok : Executed when user clicks on the window, used for click and drag cropping
% hObject    handle to edit_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

click = get(hObject,'CurrentObject');
obj = get(click,'Type');
if ((click == handles.IN_axes) || strcmp(obj,'line')) % If user clicks anywhere on the axes where current signal is displayed...
    point1 = get(handles.IN_axes,'CurrentPoint');    % button down detected
    rbbox; % return figure units
    point2 = get(handles.IN_axes,'CurrentPoint'); % button up detected
    xi = min(point1(1,1),point2(1,1));
    xf = max(point1(1,1),point2(1,1));
    if (xi >= min(handles.rel_time) && xf <= max(handles.rel_time) && xi ~= xf) % Check if selection is valid
        handles.version = handles.version + 1;
        set([handles.undo_btn handles.reset_btn],'Enable','on');
        set(handles.redo_btn,'Enable','off');
        set(handles.OUT_start,'String',num2str(xi));
        set(handles.OUT_end,'String',num2str(xf));
        % Save the selected chunk in a diferent variable
        handles.testsignal(handles.version) = handles.testsignal(handles.version - 1);
        if handles.timescale == 1
            indices = cat(2,{round((xi-min(handles.rel_time))*handles.fs)+1:round((xf-min(handles.rel_time))*handles.fs)},repmat({':'},1,ndims(handles.testsignal(handles.version - 1).audio)-1));
            handles.testsignal(handles.version).audio = handles.testsignal(handles.version - 1).audio(indices{:});
        elseif handles.timescale == 2
            indices = cat(2,{round((xi-min(handles.rel_time)))+1:round((xf-min(handles.rel_time)))},repmat({':'},1,ndims(handles.testsignal(handles.version - 1).audio)-1));
            handles.testsignal(handles.version).audio = handles.testsignal(handles.version - 1).audio(indices{:});
        end
        handles.rel_time = linspace(xi,xf,length(handles.testsignal(handles.version).audio));
        handles.xi(handles.version) = xi;
        handles.xf(handles.version) = xf;
        handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
        handles.testsignal = handles.testsignal(1:handles.version);
        handles.xi = handles.xi(1:handles.version);
        handles.xf = handles.xf(1:handles.version);
        handles.timescale = handles.timescale(1:handles.version);
        if ~ismatrix(handles.testsignal(handles.version).audio)
            set(handles.channel_panel,'Visible','on');
            linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
            if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
            if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
            if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
            set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
        else
            set(handles.channel_panel,'Visible','off');
            linea = handles.testsignal(handles.version).audio;
            cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
            set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
        end
        pixels = get_axes_width(handles.IN_axes);
        [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
        plot(handles.IN_axes,t,linea);
        xlabel(handles.IN_axes,'Time');
        set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
        audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
        %mainHandles = guidata(handles.main_stage1);
        %selectedNodes = mainHandles.mytree.getSelectedNodes;
        audiodatatext = evalc('audiodata');
        set(handles.audiodatatext,'String',audiodatatext);
        guidata(hObject,handles);
    else % Display out of boundaries warnings
        warndlg('Data selection out of boundaries','WARNING');
        set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
        set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    end
end



function OUT_start_Callback(hObject, ~, handles) %#ok : Executed when start input box changes
% hObject    handle to OUT_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OUT_start as text
%        str2double(get(hObject,'String')) returns contents of OUT_start as a double
xi = str2double(get(hObject,'String'));
xf = str2double(get(handles.OUT_end,'String'));
if ~isnan(xi) && ~isnan(xf) && xi >= min(handles.rel_time) && xi ~= xf
    handles.version = handles.version + 1;
    set([handles.undo_btn handles.reset_btn],'Enable','on');
    set(handles.redo_btn,'Enable','off');
    handles.testsignal(handles.version) = handles.testsignal(handles.version - 1);
    if get(handles.timescale_popup,'Value') == 1
        indices = cat(2,{ceil((xi-min(handles.rel_time))*handles.fs)+1:length(handles.testsignal(handles.version - 1).audio)},repmat({':'},1,ndims(handles.testsignal(handles.version - 1).audio)-1));
        handles.testsignal(handles.version).audio = handles.testsignal(handles.version - 1).audio(indices{:});
    elseif get(handles.timescale_popup,'Value') == 2
        indices = cat(2,{ceil((xi-min(handles.rel_time)))+1:length(handles.testsignal(handles.version - 1).audio)},repmat({':'},1,ndims(handles.testsignal(handles.version - 1).audio)-1));
        handles.testsignal(handles.version).audio = handles.testsignal(handles.version - 1).audio(indices{:});
    end
    handles.rel_time = linspace(xi,xf,length(handles.testsignal(handles.version).audio));
    handles.xi(handles.version) = xi;
    handles.xf(handles.version) = handles.xf(handles.version - 1);
    handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
    handles.testsignal = handles.testsignal(1:handles.version);
    handles.xi = handles.xi(1:handles.version);
    handles.xf = handles.xf(1:handles.version);
    handles.timescale = handles.timescale(1:handles.version);
    if ~ismatrix(handles.testsignal(handles.version).audio)
        set(handles.channel_panel,'Visible','on');
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
        if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
        if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
        if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    else
        set(handles.channel_panel,'Visible','off');
        linea = handles.testsignal(handles.version).audio;
        cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea);
    xlabel(handles.IN_axes,'Time');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
    audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
    %mainHandles = guidata(handles.main_stage1);
    %selectedNodes = mainHandles.mytree.getSelectedNodes;
    audiodatatext = evalc('audiodata');
    set(handles.audiodatatext,'String',audiodatatext);
else % Display out of boundaries warnings
    addsilence = questdlg('Would you like to add silence before the audio data displayed?',...
                          'Data selection out of boundaries',...
                          'Yes', 'No', 'No');
    switch addsilence
        case 'Yes'
            handles.version = handles.version + 1;
            set([handles.undo_btn handles.reset_btn],'Enable','on');
            set(handles.redo_btn,'Enable','off');
            handles.testsignal(handles.version) = handles.testsignal(handles.version - 1);
            sizeprever = size(handles.testsignal(handles.version-1).audio);
            if get(handles.timescale_popup,'Value') == 1
                handles.testsignal(handles.version).audio = cat(1,zeros([round(abs(xi-min(handles.rel_time))*handles.fs),sizeprever(2:end)]),handles.testsignal(handles.version - 1).audio);
            elseif get(handles.timescale_popup,'Value') == 2
                handles.testsignal(handles.version).audio = cat(1,zeros([round(abs(xi-min(handles.rel_time))),sizeprever(2:end)]),handles.testsignal(handles.version - 1).audio);
            end
            handles.rel_time = linspace(xi,xf,length(handles.testsignal(handles.version).audio));
            handles.xi(handles.version) = xi;
            handles.xf(handles.version) = handles.xf(handles.version - 1);
            handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
            handles.testsignal = handles.testsignal(1:handles.version);
            handles.xi = handles.xi(1:handles.version);
            handles.xf = handles.xf(1:handles.version);
            handles.timescale = handles.timescale(1:handles.version);
            if ~ismatrix(handles.testsignal(handles.version).audio)
                set(handles.channel_panel,'Visible','on');
                linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
                if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
                if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
                if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
                set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
            else
                set(handles.channel_panel,'Visible','off');
                linea = handles.testsignal(handles.version).audio;
                cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
                set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
            end
            pixels = get_axes_width(handles.IN_axes);
            [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
            plot(handles.IN_axes,t,linea);
            xlabel(handles.IN_axes,'Time');
            set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
            audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
            %mainHandles = guidata(handles.main_stage1);
            %selectedNodes = mainHandles.mytree.getSelectedNodes;
            audiodatatext = evalc('audiodata');
            set(handles.audiodatatext,'String',audiodatatext);
        case 'No'
            set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
            set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    end
end
guidata(hObject,handles);



function OUT_end_Callback(hObject, ~, handles) %#ok : Executed when end point input box changes
% hObject    handle to OUT_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OUT_end as text
%        str2double(get(hObject,'String')) returns contents of OUT_end as a double
xi = str2double(get(handles.OUT_start,'String'));
xf = str2double(get(hObject,'String'));
if ~isnan(xi) && ~isnan(xf) && xf <= max(handles.rel_time) && xi ~= xf
    handles.version = handles.version + 1;
    set([handles.undo_btn handles.reset_btn],'Enable','on');
    set(handles.redo_btn,'Enable','off');
    handles.testsignal(handles.version) = handles.testsignal(handles.version - 1);
    if handles.timescale == 1
        indices = cat(2,{1:round((xf-min(handles.rel_time))*handles.fs)},repmat({':'},1,ndims(handles.testsignal(handles.version - 1).audio)-1));
        handles.testsignal(handles.version).audio = handles.testsignal(handles.version - 1).audio(indices{:});
    elseif handles.timescale == 2
        indices = cat(2,{1:round((xf-min(handles.rel_time)))},repmat({':'},1,ndims(handles.testsignal(handles.version - 1).audio)-1));
        handles.testsignal(handles.version).audio = handles.testsignal(handles.version - 1).audio(indices{:});
    end
    handles.rel_time = linspace(xi,xf,length(handles.testsignal(handles.version).audio));
    handles.xi(handles.version) = handles.xi(handles.version - 1);
    handles.xf(handles.version) = xf;
    handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
    handles.testsignal = handles.testsignal(1:handles.version);
    handles.xi = handles.xi(1:handles.version);
    handles.xf = handles.xf(1:handles.version);
    handles.timescale = handles.timescale(1:handles.version);
    if ~ismatrix(handles.testsignal(handles.version).audio)
        set(handles.channel_panel,'Visible','on');
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
        if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
        if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
        if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    else
        set(handles.channel_panel,'Visible','off');
        linea = handles.testsignal(handles.version).audio;
        cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea);
    xlabel(handles.IN_axes,'Time');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
    audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
    %mainHandles = guidata(handles.main_stage1);
    %selectedNodes = mainHandles.mytree.getSelectedNodes;
    audiodatatext = evalc('audiodata');
    set(handles.audiodatatext,'String',audiodatatext);
else % Display out of boundaries warnings
    addsilence = questdlg('Would you like to add silence after the audio data displayed?',...
                          'Data selection out of boundaries',...
                          'Yes', 'No', 'No');
    switch addsilence
        case 'Yes'
            handles.version = handles.version + 1;
            set([handles.undo_btn handles.reset_btn],'Enable','on');
            set(handles.redo_btn,'Enable','off');
            handles.testsignal(handles.version) = handles.testsignal(handles.version - 1);
            sizeprever = size(handles.testsignal(handles.version-1).audio);
            if get(handles.timescale_popup,'Value') == 1
                handles.testsignal(handles.version).audio = cat(1,handles.testsignal(handles.version - 1).audio,zeros([round(abs(xf-max(handles.rel_time))*handles.fs),sizeprever(2:end)]));
            elseif get(handles.timescale_popup,'Value') == 2
                handles.testsignal(handles.version).audio = cat(1,handles.testsignal(handles.version - 1).audio,zeros([round(abs(xf-max(handles.rel_time))),sizeprever(2:end)]));
            end
            handles.rel_time = linspace(xi,xf,length(handles.testsignal(handles.version).audio));
            handles.xi(handles.version) = handles.xi(handles.version - 1);
            handles.xf(handles.version) = xf;
            handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
            handles.testsignal = handles.testsignal(1:handles.version);
            handles.xi = handles.xi(1:handles.version);
            handles.xf = handles.xf(1:handles.version);
            handles.timescale = handles.timescale(1:handles.version);
            if ~ismatrix(handles.testsignal(handles.version).audio)
                set(handles.channel_panel,'Visible','on');
                linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
                if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
                if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
                if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
                set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
            else
                set(handles.channel_panel,'Visible','off');
                linea = handles.testsignal(handles.version).audio;
                cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
                set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
            end
            pixels = get_axes_width(handles.IN_axes);
            [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
            plot(handles.IN_axes,t,linea);
            xlabel(handles.IN_axes,'Time');
            set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
            audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
            %mainHandles = guidata(handles.main_stage1);
            %selectedNodes = mainHandles.mytree.getSelectedNodes;
            audiodatatext = evalc('audiodata');
            set(handles.audiodatatext,'String',audiodatatext);
        case 'No'
            set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
            set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    end
end
guidata(hObject,handles);


% --- Executes on button press in reset_btn.
function reset_btn_Callback(hObject, ~, handles) %#ok : Executed when reset button is clicked
% hObject    handle to reset_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
% Call the 'desktop'
hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');
% Bring up the data from the selected leaf to be edited
handles.version = handles.version + 1;
set(handles.undo_btn,'Enable','on');
set([handles.redo_btn hObject],'Enable','off');
handles.testsignal(handles.version) = audiodata;
handles.fs = audiodata.fs;
dur = length(handles.testsignal(handles.version).audio)/handles.fs;
% Allocate memory space for the edited signal
handles.rel_time = linspace(0,dur,length(handles.testsignal(handles.version).audio));
handles.xi(handles.version) = min(handles.rel_time);
handles.xf(handles.version) = max(handles.rel_time);
set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
set(handles.timescale_popup,'Value',1);
handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
handles.testsignal = handles.testsignal(1:handles.version);
handles.xi = handles.xi(1:handles.version);
handles.xf = handles.xf(1:handles.version);
handles.timescale = handles.timescale(1:handles.version);
% Plot signal to be cropped
if ~ismatrix(handles.testsignal(handles.version).audio)
    set(handles.channel_panel,'Visible','on');
    set(handles.IN_nchannel,'String','1');
    set(handles.tchannels,'String',['/ ' num2str(size(handles.testsignal(handles.version).audio,2))]);
    linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
    if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
    if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
    if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
    set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
else
    set(handles.channel_panel,'Visible','off');
    linea = handles.testsignal(handles.version).audio;
    cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
    set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
end
pixels = get_axes_width(handles.IN_axes);
[t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
plot(handles.IN_axes,t,linea)
xlabel(handles.IN_axes,'Time');
set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
%mainHandles = guidata(handles.main_stage1);
%selectedNodes = mainHandles.mytree.getSelectedNodes;
audiodatatext = evalc('audiodata');
set(handles.audiodatatext,'String',audiodatatext);
%    handles.line = findobj(gcf,'type','line');
% Update handles structure
guidata(hObject, handles);
fprintf(handles.fid, '%% RESET EDIT\n');


% --- Executes on selection change in edit_box.
function edit_box_Callback(hObject, ~, handles) %#ok : Executed when selection changes in the edit function box
% hObject    handle to edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns edit_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_box
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
[~,funname] = fileparts(selection);
if ~strcmp(selection,'crop.m')
    handles.funname = funname;
    helptext = evalc(['help ' funname]);
    set(hObject,'Tooltip',helptext);
    set(handles.apply_btn,'Enable','on');
    set(handles.apply_btn,'BackgroundColor',[0.94 0.94 0.94]);
else
    handles.funname = [];
    set(hObject,'Tooltip','Click where you want to begin the cropped selection on the axes, hold down the left mouse button while dragging the pointer over the region that you want to crop, let go of the mouse button when you finish the selection.');
    set(handles.apply_btn,'Enable','off');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_box_CreateFcn(hObject, ~, handles) %#ok : Funtion listing box creation
% hObject    handle to edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

curdir = cd;
tools = what([curdir '/Processors/Basic']);
if ~isempty(tools.m)
    set(hObject,'String',['crop.m';cellstr(tools.m)],'Value',1);
else
    set(hObject,'String','crop.m','Value',1);
end
guidata(hObject,handles)


% --- Executes on button press in apply_btn.
function apply_btn_Callback(hObject, ~, handles) %#ok : Executed when apply button is clicked
% hObject    handle to apply_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');
handles.version = handles.version + 1;
processed = feval(handles.funname,handles.testsignal(handles.version - 1));
if ~isempty(processed)
    set([handles.undo_btn handles.reset_btn],'Enable','on');
    set(handles.redo_btn,'Enable','off');
    if isstruct(processed)
        % The following prevents new fields from being added - maybe this
        % is too limiting
        newdata = handles.testsignal(handles.version - 1); 
        dif = intersect(fieldnames(handles.testsignal(handles.version - 1)),fieldnames(processed));
        for i = 1:size(dif,1)
            newdata.(dif{i,1}) = processed.(dif{i,1});
        end
        newdata = addhistory(newdata,'Edited');
    else
        newdata = handles.testsignal(handles.version - 1);
        newdata.audio = processed;
    end
    handles.testsignal(handles.version) = newdata;
    handles.xi(handles.version) = handles.xi(handles.version - 1);
    handles.xf(handles.version) = handles.xf(handles.version - 1);
    handles.timescale(handles.version) = get(handles.timescale_popup,'Value');
    handles.testsignal = handles.testsignal(1:handles.version);
    handles.xi = handles.xi(1:handles.version);
    handles.xf = handles.xf(1:handles.version);
    handles.timescale = handles.timescale(1:handles.version);
    handles.rel_time = linspace(handles.xi(handles.version),handles.xf(handles.version),length(handles.testsignal(handles.version).audio));
    if ~ismatrix(handles.testsignal(handles.version).audio)
        set(handles.channel_panel,'Visible','on');
        set(handles.tchannels,'String',['/ ' num2str(size(handles.testsignal(handles.version).audio,2))]);
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
        if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
        if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
        if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    else
        set(handles.channel_panel,'Visible','off');
        linea = handles.testsignal(handles.version).audio;
        cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea)
    xlabel(handles.IN_axes,'Time');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
    audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
    %mainHandles = guidata(handles.main_stage1);
    %selectedNodes = mainHandles.mytree.getSelectedNodes;
    audiodatatext = evalc('audiodata');
    set(handles.audiodatatext,'String',audiodatatext);
    fprintf(handles.fid, ['%% ' datestr(now,16) ' - Processed in Edit window using ' handles.funname '\n']);
    %fprintf(handles.fid,[audiodatatext,'\n']);
    % Log verbose metadata
    logaudioleaffields(processed);
    if isfield(handles,'choosefromhigherdims')
         handles.choosefromhigherdims = [];
    end
else
    handles.version = handles.version - 1;
end
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');
guidata(hObject, handles);



function IN_nchannel_Callback(hObject, ~, handles) %#ok :  Executed when channel number input box changes
% hObject    handle to IN_nchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_nchannel as text
%        str2double(get(hObject,'String')) returns contents of IN_nchannel as a double
channel = str2double(get(handles.IN_nchannel,'String'));

if (channel <= size(handles.testsignal(handles.version).audio,2)) && (channel > 0) && ~isnan(channel)
    handles.channel = channel;
    linea(:,:) = handles.testsignal(handles.version).audio(:,channel,:);
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea)
    xlabel(handles.IN_axes,'Time');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
else
    warndlg('Invalid channel');
    set(handles.IN_nchannel,'String',num2str(handles.channel));
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function IN_nchannel_CreateFcn(hObject, ~, handles) %#ok : Channel number input box creation
% hObject    handle to IN_nchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.channel = 1;
guidata(hObject,handles)


% --- Executes on selection change in timescale_popup.
function timescale_popup_Callback(hObject, ~, handles) %#ok : Executed when timescale menu changes (samples or seconds)
% hObject    handle to timescale_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns timescale_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from timescale_popup
contents = cellstr(get(hObject,'String'));
timescale = contents{get(hObject,'Value')};

if handles.timescale(handles.version) == 1 && strcmp(timescale,'Samples')
    handles.xi(handles.version) = round(handles.xi(handles.version)*handles.testsignal(handles.version).fs);
    handles.xf(handles.version) = round(handles.xf(handles.version)*handles.testsignal(handles.version).fs);
    set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
    set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    handles.rel_time = linspace(handles.xi(handles.version),handles.xf(handles.version),length(handles.testsignal(handles.version).audio));
    if ~ismatrix(handles.testsignal(handles.version).audio)
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
    else
        linea = handles.testsignal(handles.version).audio;
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea)
    xlabel(handles.IN_axes,'Samples');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
elseif handles.timescale(handles.version) == 2 && strcmp(timescale,'Seconds')
    handles.xi(handles.version) = handles.xi(handles.version)/handles.testsignal(handles.version).fs;
    handles.xf(handles.version) = handles.xf(handles.version)/handles.testsignal(handles.version).fs;
    set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
    set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    handles.rel_time = linspace(handles.xi(handles.version),handles.xf(handles.version),length(handles.testsignal(handles.version).audio));
    if ~ismatrix(handles.testsignal(handles.version).audio)
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
    else
        linea = handles.testsignal(handles.version).audio;
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea)
    xlabel(handles.IN_axes,'Time');
    set(handles.IN_axes,'XTickLabel',num2str(get(handles.IN_axes,'XTick').'))
end
handles.timescale(handles.version) = get(handles.timescale_popup,'Value');

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function timescale_popup_CreateFcn(hObject, ~, handles) %#ok : Timescale menu popup creation
% hObject    handle to timescale_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.timescale(1) = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in undo_btn.
function undo_btn_Callback(hObject, ~, handles) %#ok : Executed when undo button is clicked
% hObject    handle to undo_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.version = handles.version - 1;
set([handles.redo_btn handles.reset_btn],'Enable','on');
set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
set(handles.timescale_popup,'Value',handles.timescale(handles.version));
handles.rel_time = linspace(handles.xi(handles.version),handles.xf(handles.version),length(handles.testsignal(handles.version).audio));
if ~ismatrix(handles.testsignal(handles.version).audio)
    set(handles.channel_panel,'Visible','on');
    linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
    if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
    if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
    if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
    set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
else
    set(handles.channel_panel,'Visible','off');
    linea = handles.testsignal(handles.version).audio;
    cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
    set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
end
pixels = get_axes_width(handles.IN_axes);
[t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
plot(handles.IN_axes,t,linea);
audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
handles.testsignal(handles.version) = addhistory(handles.testsignal(handles.version),'UNDO PREVIOUS EDIT');
%mainHandles = guidata(handles.main_stage1);
%selectedNodes = mainHandles.mytree.getSelectedNodes;
audiodatatext = evalc('audiodata');
set(handles.audiodatatext,'String',audiodatatext);
fprintf(handles.fid,'%% UNDO PREVIOUS EDIT\n');
if handles.version == 1, set(hObject,'Enable','off'); end
guidata(hObject,handles);

% --- Executes on button press in redo_btn.
function redo_btn_Callback(hObject, ~, handles) %#ok : Executed when redo button is clicked
% hObject    handle to redo_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.version = handles.version + 1;
set([handles.undo_btn handles.reset_btn],'Enable','on');
if handles.version <= length(handles.testsignal)
    set(handles.OUT_start,'String',num2str(handles.xi(handles.version)));
    set(handles.OUT_end,'String',num2str(handles.xf(handles.version)));
    set(handles.timescale_popup,'Value',handles.timescale(handles.version));
    handles.rel_time = linspace(handles.xi(handles.version),handles.xf(handles.version),length(handles.testsignal(handles.version).audio));
    if ~ismatrix(handles.testsignal(handles.version).audio)
        set(handles.channel_panel,'Visible','on');
        linea(:,:) = handles.testsignal(handles.version).audio(:,str2double(get(handles.IN_nchannel,'String')),:);
        if ndims(handles.testsignal(handles.version).audio) == 3, cmap = colormap(hsv(size(handles.testsignal(handles.version).audio,3))); end
        if ndims(handles.testsignal(handles.version).audio) == 4, cmap = colormap(copper(size(handles.testsignal(handles.version).audio,4))); end
        if ndims(handles.testsignal(handles.version).audio) >= 5, cmap = colormap(cool(size(handles.testsignal(handles.version).audio,5))); end
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    else
        set(handles.channel_panel,'Visible','off');
        linea = handles.testsignal(handles.version).audio;
        cmap = colormap(lines(size(handles.testsignal(handles.version).audio,2)));
        set(handles.edit_signal,'DefaultAxesColorOrder',cmap)
    end
    pixels = get_axes_width(handles.IN_axes);
    [t, linea] = reduce_to_width(handles.rel_time, real(linea), pixels, [-inf inf]);
    plot(handles.IN_axes,t,linea);
    audiodata = handles.testsignal(handles.version); %#ok : Used in evalc below
    %mainHandles = guidata(handles.main_stage1);
    %selectedNodes = mainHandles.mytree.getSelectedNodes;
    audiodatatext = evalc('audiodata');
    set(handles.audiodatatext,'String',audiodatatext);
    if handles.version == 1, set(hObject,'Enable','off'); end
    if handles.version == length(handles.testsignal), set(hObject,'Enable','off'); end
    guidata(hObject,handles);
    fprintf(handles.fid,'%% REDO PREVIOUS EDIT\n');
end


% --- Executes on selection change in device_popup.
function device_popup_Callback(hObject, ~, handles) %#ok : Eexcuted when output device selection menu changes
% hObject    handle to device_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns device_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from device_popup

selection = get(hObject,'Value');
handles.odeviceid = handles.odeviceidlist(selection);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function device_popup_CreateFcn(hObject, ~, handles) %#ok : Output device selection creation
% hObject    handle to device_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

devinfo = audiodevinfo; % Get available device information
odevicelist = {devinfo.output.Name}; % Populate list
handles.odeviceidlist =  cell2mat({devinfo.output.ID});
handles.odeviceid = handles.odeviceidlist(1,1);
set(hObject,'String',odevicelist);
guidata(hObject,handles);


% --- Executes on button press in play_btn.
function play_btn_Callback(hObject, ~, handles) %#ok : Executed when play button is clicked
% hObject    handle to play_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the 'desktop'
hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');

if isempty(audiodata)
    warndlg('No signal loaded!');
else
    % Retrieve information from the selected leaf
    testsignal = handles.testsignal(handles.version).audio./max(max(max(abs(handles.testsignal(handles.version).audio))));
    fs = handles.testsignal(handles.version).fs;
    nbits = 16;
    doesSupport = audiodevinfo(0, handles.odeviceid, fs, nbits, size(testsignal,2));
    if doesSupport && ismatrix(testsignal)
        % Play signal
        handles.player = audioplayer(testsignal,fs,nbits,handles.odeviceid);
        play(handles.player);
        set(handles.stop_btn,'Enable','on');
    else
        warndlg('Device not supported for playback!');
    end
end
guidata(hObject, handles);


% --- Executes on button press in stop_btn.
function stop_btn_Callback(hObject, ~, handles) %#ok : Executed when stop button is clicked
% hObject    handle to stop_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isplaying(handles.player)
    stop(handles.player);
end
guidata(hObject,handles);


% write to new
% --- Executes on button press in wn_btn.
function wn_btn_Callback(hObject, ~, handles) %#ok : Executed when Write to new button is clicked
% hObject    handle to wn_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mainHandles = guidata(handles.main_stage1);
if ~isempty(handles.testsignal(handles.version))
    aarae_fig = findobj('type','figure','tag','aarae');
    if strcmp(handles.testsignal(handles.version).datatype,'syscal')
        iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/boardicon.gif');
    else
        iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
    end
    handles.selNodeName = [handles.selNodeName '_edit'];
    handles.testsignal(handles.version).name = handles.selNodeName;
    mainHandles.(matlab.lang.makeValidName(handles.selNodeName)) = uitreenode('v0', handles.selNodeName, handles.selNodeName,  iconPath, true);
    switch handles.testsignal(handles.version).datatype
        case 'testsignals', ivalue = 1;
        case 'measurements', ivalue = 2;
        case 'syscal', ivalue = 2;
        case 'processed', ivalue = 3;
        case 'results', ivalue = 4;
    end
    [branch,ok] = listdlg('ListString',{'Test signals','Measurements','Processed','Results'},'SelectionMode','single','PromptString','Save edited audio in:','InitialValue',ivalue);
    datatype = handles.testsignal(handles.version).datatype;
    if ok == 0
        mainHandles.(matlab.lang.makeValidName(handles.testsignal(handles.version).datatype)).add(mainHandles.(matlab.lang.makeValidName(handles.selNodeName)));
    else
        if branch == 1, mainHandles.testsignals.add(mainHandles.(matlab.lang.makeValidName(handles.selNodeName))); handles.testsignal(handles.version).datatype = 'testsignals'; end
        if branch == 2, mainHandles.measurements.add(mainHandles.(matlab.lang.makeValidName(handles.selNodeName))); handles.testsignal(handles.version).datatype = 'measurements'; end
        if branch == 3, mainHandles.processed.add(mainHandles.(matlab.lang.makeValidName(handles.selNodeName))); handles.testsignal(handles.version).datatype = 'processed'; end
        if branch == 4, mainHandles.results.add(mainHandles.(matlab.lang.makeValidName(handles.selNodeName))); handles.testsignal(handles.version).datatype = 'results'; end
    end
    if strcmp(datatype,'syscal'), handles.testsignal(handles.version).datatype = 'syscal'; end
    
    % Save as you go
    signaldata = handles.testsignal(handles.version); 
    signaldata.name = [handles.selNodeName '_edit'];
    save([cd '/Utilities/Backup/' handles.selNodeName '.mat'], 'signaldata');
    
    mainHandles.(matlab.lang.makeValidName(handles.selNodeName)).UserData = handles.testsignal(handles.version);
    mainHandles.mytree.reloadNode(mainHandles.(matlab.lang.makeValidName(handles.selNodeName)).getParent);
    mainHandles.mytree.setSelectedNode(mainHandles.(matlab.lang.makeValidName(handles.selNodeName)));
    guidata(aarae_fig, mainHandles);
    fprintf(handles.fid,'%% WRITE TO NEW\n');
end
guidata(hObject,handles);
uiresume(handles.edit_signal);



function IN_name_Callback(hObject, ~, handles) %#ok : Executed when the content of the name input changes
% hObject    handle to IN_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_name as text
%        str2double(get(hObject,'String')) returns contents of IN_name as a double
if ~isempty(get(hObject,'String'))
    handles.selNodeName = matlab.lang.makeValidName(get(hObject,'String'));
else
    set(hObject,'String',matlab.lang.makeValidName(handles.selNodeName))
end
% handles.testsignal.name = [handles.selNodeName '_edit'];
% fprintf(handles.fid,['%% New name: ', [handles.selNodeName '_edit'],'\n']);
handles.testsignal(2).name = matlab.lang.makeValidName(handles.selNodeName);
fprintf(handles.fid,['%% New name: ', matlab.lang.makeValidName(handles.selNodeName),'\n']);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function IN_name_CreateFcn(hObject, ~, ~) %#ok : Creation of name input box
% hObject    handle to IN_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
