% DO NOT EDIT THIS INITIALIZATION FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!
function varargout = calculator(varargin)
% CALCULATOR MATLAB code for calculator.fig
%      CALCULATOR, by itself, creates a new CALCULATOR or raises the existing
%      singleton*.
%
%      H = CALCULATOR returns the handle to a new CALCULATOR or the handle to
%      the existing singleton*.
%
%      CALCULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCULATOR.M with the given input arguments.
%
%      CALCULATOR('Property','Value',...) creates a new CALCULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calculator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calculator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calculator

% Last Modified by GUIDE v2.5 01-Nov-2013 14:38:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calculator_OpeningFcn, ...
                   'gui_OutputFcn',  @calculator_OutputFcn, ...
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


% --- Executes just before calculator is made visible.
function calculator_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calculator (see VARARGIN)


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


% Update handles structure
guidata(hObject, handles);

if dontOpen
   disp('-----------------------------------------------------');
   disp('This function is part of the AARAE framework, it is') 
   disp('not a standalone function. To call this function,')
   disp('click on the appropriate calling button on the main');
   disp('Window. E.g.:');
   disp('   Calculators');
   disp('-----------------------------------------------------');
else
   uiwait(hObject);
end


% --- Outputs from this function are returned to the command line.
function varargout = calculator_OutputFcn(hObject, ~, ~) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [];
delete(hObject);


% --- Executes on button press in calc_btn.
function calc_btn_Callback(hObject, ~, handles) %#ok : Executes on Calculate button click
% hObject    handle to calc_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainHandles = guidata(handles.main_stage1);
set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');
already_logged = false;
if nargout(handles.funname) == 1
    signaldata = feval(handles.funname);
    if ~isempty(signaldata)
        mainHandles.mytree.setSelectedNode(mainHandles.root);
        if isfield(signaldata,'audio')
            signaldata = checkcal(signaldata);
            %signaldata.nbits = 16;
            iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
            if isfield(signaldata,'tables')
                % if tables are ALSO present, then write to the log file
                % before removing the tables field
                % In a possible future revision, we could split the output
                % into two leaves (one with tables, the other with audio)
                fprintf(mainHandles.fid, ['%% ' datestr(now,16) ' - Used calculator ' handles.funname '\n']);  
                logaudioleaffields(signaldata,0); % Log verbose metadata including tables
                signaldata = rmfield(signaldata,'tables');
                already_logged = true;
            end
        else
            iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/notesicon.gif');
        end

        if length(fieldnames(signaldata)) ~= 1
            signaldata.datatype = 'results';
            leafname = isfield(mainHandles,matlab.lang.makeValidName(handles.funname));
            
            % THE FOLLOWING HAS BEEN REPLACED BY THE CODE THAT FOLLOWS IT
            % (SEE COMMENT ABOUT BUG FROM CHANGING handles.funname)
%             if leafname == 1
%                 index = 1;
%                 % This while cycle is just to make sure no signals are
%                 % overwriten
%                 if length(matlab.lang.makeValidName(handles.funname)) >= namelengthmax, handles.funname = handles.funname(1:round(end/2)); end
%                 while isfield(handles,matlab.lang.makeValidName([handles.funname,'_',num2str(index)])) == 1
%                     index = index + 1;
%                 end
%                 % the following causes a bug if you run the function 3
%                 % times without closing Calculators
%                 handles.funname = matlab.lang.makeValidName([handles.funname,'_',num2str(index)]);
%             end
%             signaldata.name = handles.funname;
%             
%             % Save as you go
%             save([cd '/Utilities/Backup/' handles.funname '.mat'], 'signaldata');
%             
%             mainHandles.(matlab.lang.makeValidName(handles.funname)) = uitreenode('v0', handles.funname,  handles.funname,  iconPath, true);
%             mainHandles.(matlab.lang.makeValidName(handles.funname)).UserData = signaldata;
%             mainHandles.results.add(mainHandles.(matlab.lang.makeValidName(handles.funname)));

            signaldata.name = handles.funname;
            if leafname == 1
                index = 1;
                % This while cycle is just to make sure no signals are
                % overwriten
                
                if length(matlab.lang.makeValidName(signaldata.name)) >= namelengthmax, signaldata.name = signaldata.name(1:round(end/2)); end
                while isfield(handles,matlab.lang.makeValidName([signaldata.name,'_',num2str(index)])) == 1
                    index = index + 1;
                end
                signaldata.name = matlab.lang.makeValidName([signaldata.name,'_',num2str(index)]);
            end
            
            % Save as you go
            save([cd '/Utilities/Backup/' signaldata.name '.mat'], 'signaldata');
            
            mainHandles.(matlab.lang.makeValidName(handles.funname)) = uitreenode('v0', handles.funname,  handles.funname,  iconPath, true);
            mainHandles.(matlab.lang.makeValidName(handles.funname)).UserData = signaldata;
            mainHandles.results.add(mainHandles.(matlab.lang.makeValidName(handles.funname)));

            mainHandles.mytree.reloadNode(mainHandles.results);
            mainHandles.mytree.expand(mainHandles.results);
            % The following line causes a bug when figures and audio are
            % both output by the calculator. This is particularly an issue
            % with colormap figures (such as image). Probably selection of
            % the new audio leaf should be done when we return to the aarae
            % main gui.
           % mainHandles.mytree.setSelectedNode(mainHandles.(matlab.lang.makeValidName(handles.funname)));
            set([mainHandles.clrall_btn,mainHandles.export_btn],'Enable','on')
        end
        if ~already_logged
        fprintf(mainHandles.fid, ['%% ' datestr(now,16) ' - Used calculator ' handles.funname '\n']);
        
        % Log verbose metadata
        logaudioleaffields(signaldata,0);
        end
    end
else
    feval(handles.funname);
end
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');
h = findobj('type','figure','-not','tag','aarae','-not','tag','calculators');
aarae_fig = findobj('type','figure','tag','aarae');
index = 1;
filename = dir([cd '/Utilities/Temp/' handles.funname num2str(index) '.fig']);
if ~isempty(filename)
    while isempty(dir([cd '/Utilities/Temp/' handles.funname num2str(index) '.fig'])) == 0
        index = index + 1;
    end
end
if ~isempty(h)
    for i = 1:length(h)
        % There may be a timing bug here if big-data figures are closed too early
        saveas(h(i),[cd '/Utilities/Temp/' handles.funname num2str(index) '.fig']);
        fprintf(mainHandles.fid,['%% Result figure name: ', handles.funname num2str(index), '.fig, temporarily stored in /Utilities/Temp/\n']);
        index = index + 1;
    end
end
fprintf(mainHandles.fid,'\n');
results = dir([cd '/Utilities/Temp']);
guidata(aarae_fig, mainHandles);
guidata(hObject, handles);


% --- Executes on button press in close_btn.
function close_btn_Callback(~, ~, handles) %#ok : Executes on Close button click
% hObject    handle to close_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.calculators);


% --- Executes when user attempts to close calculator.
function calculators_CloseRequestFcn(hObject, ~, ~) %#ok : Window close request function
% hObject    handle to calculator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(hObject);


% --- Executes on selection change in signalcat_box.
function signalcat_box_Callback(hObject, ~, handles) %#ok : Calculator categories selection box
% hObject    handle to signalcat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns signalcat_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from signalcat_box

% Displays the available calculators for the selected processing category
contents = cellstr(get(hObject,'String'));
calccat = contents{get(hObject,'Value')};
calculators = what([cd '/Calculators/' calccat]);
try
if ~isempty(cellstr(calculators.m))
    set(handles.signal_box,'Visible','on','String',[' ';cellstr(calculators.m)],'Value',1,'Tooltip','');
    set(handles.calc_btn,'Visible','off');
else
    set(handles.signal_box,'Visible','off');
    set(handles.calc_btn,'Visible','off');
end
guidata(hObject,handles);
catch
end

% --- Executes during object creation, after setting all properties.
function signalcat_box_CreateFcn(hObject, ~, handles) %#ok : Calculator categories box creation
% hObject    handle to signalcat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Populate function box with the function available in the folder 'Calculators'
curdir = cd;
signals = removedotfiles(dir([curdir '/Calculators']));
%set(hObject,'String',[' ';cellstr({signals(3:length(signals)).name}')]);
set(hObject,'String',[' ';cellstr({signals.name}')]);
handles.funname = [];
guidata(hObject,handles)


% --- Executes on selection change in signal_box.
function signal_box_Callback(hObject, ~, handles) %#ok : Calculator list selection box
% hObject    handle to signal_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns signal_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from signal_box
% Allows 'help' display when mouse is hoovered over fun_box
try
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
[~,funname] = fileparts(selection);
if ~strcmp(selection,' ')
    handles.funname = funname;
    helptext = evalc(['help ' funname]);
    set(hObject,'Tooltip',helptext);
    set(handles.calc_btn,'Visible','on');
else
    handles.funname = [];
    set(hObject,'Tooltip','');
    set(handles.calc_btn,'Visible','off');
end
guidata(hObject,handles);
catch
end

% --- Executes during object creation, after setting all properties.
function signal_box_CreateFcn(hObject, ~, ~) %#ok : Calculator selection box creation
% hObject    handle to signal_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% *************************************************************************
% Local utility functions
% *************************************************************************

function OUT = removedotfiles(IN)
% This function takes a structure that was created by Matlab's dir function
% and removes all items that have a period at the start of the name field.
% The purpose of this is to remove system files from the list
try
    include = true(length(IN),1);
    for n = 1:length(IN)
        ind = strfind(IN(n).name,'.');
        if ~isempty(ind)
            if ind(1) == 1
                include(n)= false;
            end
        end
    end
    OUT = IN(include);
catch
    OUT = IN;
end