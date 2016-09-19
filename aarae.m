% DO NOT EDIT THIS INITIALIZATION FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!
function varargout = aarae(varargin)
% AARAE MATLAB code for aarae.fig
%      AARAE, by itself, creates a new AARAE or raises the existing
%      singleton*.
%
%      H = AARAE returns the handle to a new AARAE or the handle to
%      the existing singleton*.
%
%      AARAE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AARAE.M with the given input arguments.
%
%      AARAE('Property','Value',...) creates a new AARAE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aarae_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aarae_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aarae

% Last Modified by GUIDE v2.5 15-Sep-2016 21:02:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @aarae_OpeningFcn, ...
    'gui_OutputFcn',  @aarae_OutputFcn, ...
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



% *************************************************************************
% *************************************************************************
%                         STARTING AARAE
% *************************************************************************
% *************************************************************************




% *************************************************************************
% INITIALIZE AARAE THINGS
% *************************************************************************

% --- Executes just before aarae is made visible.
function aarae_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to aarae
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aarae (see VARARGIN)

AARAEversion = 'AARAE Release 9'; % *** UPDATE THIS WITH EACH RELEASE ***

% Choose default command line output for aarae
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Setup the 'desktop'
setappdata(0, 'hMain', gcf);
hMain = getappdata(0,'hMain');
setappdata(hMain,'testsignal',[]);
setappdata(hMain,'audio_recorder_input',1)
setappdata(hMain,'audio_recorder_output',1)
setappdata(hMain,'audio_recorder_numchs',1)
setappdata(hMain,'audio_recorder_duration',1)
setappdata(hMain,'audio_recorder_fs',48000)
%setappdata(hMain,'audio_recorder_nbits',16)
setappdata(hMain,'audio_recorder_qdur',0.1) % set by testing on Matlab 2015b
setappdata(hMain,'audio_recorder_buffer',8192) % set by testing on Matlab 2015b
setappdata(hMain,'audio_recorder_silencerequest',0)
setappdata(hMain,'trim_method_after_convolution',1)
setappdata(hMain,'AARAEversion',AARAEversion)
set(hObject,'Name','AARAE')
% Read settings file
Settings = [];
if ~isempty(dir([cd '/Settings.mat']))
    load([cd '/Settings.mat']);
    handles.Settings = Settings;
else
    Settings.maxtimetodisplay = 10;
    Settings.frequencylimits = 'Default';
    Settings.calibrationtoggle = 1;
    Settings.maxlines = 100;
    Settings.colormap = 'Jet';
    Settings.specmagscale = 'Raw';
    Settings.username = 'AARAE User';
    handles.Settings = Settings;
    save([cd '/Settings.mat'],'Settings')
end

handles.compareaudio = -1; % prevents a timing bug when large files are loaded and the compare button pressed immediately

if ~isdir([cd '/Log']), mkdir([cd '/Log']); end
if ~isdir([cd '/Utilities/Temp'])
    mkdir([cd '/Utilities/Temp']);
else
    results = dir([cd '/Utilities/Temp']);
    set(handles.result_box,'String',[' ';cellstr({results(3:length(results)).name}')]);
end
if ~isdir([cd '/Utilities/Backup'])
    mkdir([cd '/Utilities/Backup']);
    handles.defaultaudiopath = [cd '/Audio'];
else
    handles.defaultaudiopath = [cd '/Utilities/Backup'];
    set(handles.recovery_txt,'Visible','on')
end

% Add folder paths for filter functions and signal analyzers
addpath(genpath(cd));
handles.player = audioplayer(0,48000);

if ~isdir([cd '/Audio/REQUIRED_AUDIO'])
    mkdir([cd '/Audio/REQUIRED_AUDIO']);
end

try
    [handles.reference_audio.audio, handles.reference_audio.fs] = audioread('/Audio/REQUIRED_AUDIO/REFERENCE_AUDIO.wav');
catch
    h = warndlg('REFERENCE_AUDIO.wav is missing from the REQUIRED_AUDIO directory (within AARAE''s Audio directory). Please select an alternative reference audio recording.');
    uiwait(h)
    audiochoice = importaudio;
    if ~isempty(audiochoice)
        handles.reference_audio = mean(audiochoice(:,:,1,1,1,1),2); %mixdown if multichan
    else
        handles.reference_audio.audio = sin(2*pi*1000*(0:48000)'./48000); % 1 kHz sine
        handles.reference_audio.fs = 48000;
    end
end

try
    [handles.silenceplease.audio, handles.silenceplease.fs] = audioread('/Audio/REQUIRED_AUDIO/SILENCE_PLEASE.wav');
catch
    modulator = 10*sin(2*pi*8*(0:48000)'./48000); % 8 Hz freq modulation
    handles.silenceplease.audio = sin(modulator+2*pi*1000*(0:48000)'./48000);
    handles.silenceplease.fs = 48000;
end

try
    [handles.thankyou.audio, handles.thankyou.fs] = audioread('/Audio/REQUIRED_AUDIO/THANKYOU.wav');
catch
    handles.thankyou.audio = sin(2*pi*1000*(0:48000)'./48000);
    handles.thankyou.fs = 48000;
end

% Setup for Densil's tree
iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/matlabicon.gif');
handles.root = uitreenode('v0', 'project', 'My project', iconPath, false);

iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/foldericon.gif');
handles.testsignals = uitreenode('v0', 'testsignals', 'Test signals', iconPath, false);
handles.measurements = uitreenode('v0', 'measurements', 'Measurements', iconPath, false);
handles.processed = uitreenode('v0', 'processed', 'Processed', iconPath, false);
handles.results = uitreenode('v0', 'results', 'Results', iconPath, false);

handles.root.add(handles.testsignals);
handles.root.add(handles.measurements);
handles.root.add(handles.processed);
handles.root.add(handles.results);

[handles.mytree,container] = uitree('v0','Root', handles.root,'SelectionChangeFcn',@mySelectFcn);
set(container, 'Parent', hObject);
treeheight_char = get(handles.process_panel,'Position')+get(handles.analysis_panel,'Position')+get(handles.uipanel1,'Position');
treewidth_char = get(handles.analysis_panel,'Position');
set(handles.analysis_panel,'Units','pixels');
treewidth_pix = get(handles.analysis_panel,'Position');
factor = treewidth_pix./treewidth_char;
set(handles.mytree,'Position',[0,treewidth_char(1,2)*factor(1,2),treewidth_char(1,1)*factor(1,1),treeheight_char(1,4)*factor(1,4)]);
handles.mytree.expand(handles.root);
handles.mytree.setSelectedNode(handles.root);
handles.mytree.setMultipleSelectionEnabled(true);

% Generate activity log
activity = dir([cd '/Log' '/activity log 1.txt']);
if isempty(activity)
    activitylog = '/activity log 1.txt';
    handles.fid = fopen([cd '/Log' activitylog], 'w');
    handles.activitylog = activitylog; % used for export all
else
    index = 2;
    % This while cycle is just to make sure no files are overwriten
    while isempty(dir([cd '/Log' '/activity log ',num2str(index),'.txt'])) == 0
        index = index + 1;
    end
    activitylog = ['/activity log ',num2str(index),'.txt'];
    handles.fid = fopen([cd '/Log' activitylog], 'w');
    handles.activitylog = activitylog; % used for export all
end

handles.alternate = 0;
fprintf(handles.fid, ['%% AARAE session started ' datestr(now) '\n']);
fprintf(handles.fid, ['%% User: ' handles.Settings.username '\n\n']);
fprintf(handles.fid, '%% Audio and Acoustical Response Environment (AARAE) for Matlab\n');
fprintf(handles.fid, ['%% ' AARAEversion  '\n']);
fprintf(handles.fid, '%% http://aarae.org\n\n');
fprintf(handles.fid, '%% Reference:\n');
fprintf(handles.fid, '%% Cabrera, D., Jimenez, D., & Martens, W. L. (2014, November).\n');
fprintf(handles.fid, '%% Audio and Acoustical Response Analysis Environment (AARAE): \n');
fprintf(handles.fid, '%% a tool to support education and research in acoustics. \n');
fprintf(handles.fid, '%% In Proceedings of Internoise. Melbourne, Australia.\n\n');
fprintf(handles.fid, '%% For other publications about AARAE refer to aarae.org\n\n');
fprintf(handles.fid, '%% This AARAE log file contains a record of activity, including:\n');
fprintf(handles.fid, '%%   * descriptions of audio and other data;\n');
fprintf(handles.fid, '%%   * descriptions of activity;\n');
fprintf(handles.fid, '%%   * results tables from analysers;\n');
fprintf(handles.fid, '%%   * names of exported files; and\n');
fprintf(handles.fid, '%%   * function call equivalents to AARAE activity.\n\n');
fprintf(handles.fid, '%% When opening the log file in a spreadsheet program, use comma delimiting to distribute table data into the appropriate spreadsheet cells.\n');
fprintf(handles.fid, '%% (However, bear in mind that the commas required in the log file''s function calls will be removed when you do that.)\n\n');
fprintf(handles.fid, '%% The code written to this log file may be useful for adapting in writing an AARAE workflow function.\n');
fprintf(handles.fid, '%% Examples of workflows are in AARAE''s Workflows folder.\n\n');
fprintf(handles.fid, '%% **************************************************\n\n\n\n');

guidata(hObject, handles);

% AARAE function for dealing with differences between Windows and Mac font
% size
fontsize

% Set waiting flag in appdata
setappdata(handles.aarae,'waiting',1)
% UIWAIT makes aarae wait for user response (see UIRESUME)
uiwait(handles.aarae);






% *************************************************************************
% SETTINGS BUTTON
% *************************************************************************

% --- Executes on button press in settings_btn.
function settings_btn_Callback(hObject, ~, handles) %#ok : Executed when Settings button is clicked
% hObject    handle to settings_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Settings = settings('main_stage1', handles.aarae);%inputdlg('Maximum time period to display','AARAE settings',[1 50],{num2str(handles.Settings.maxtimetodisplay)});
if ~isempty(Settings)
    %newpref = cell2struct(newpref,{'maxtimetodisplay'});
    %newpref.maxtimetodisplay = str2double(newpref.maxtimetodisplay);
    handles.Settings = Settings;
    save([cd '/Settings.mat'],'Settings')
    guidata(hObject,handles)
    selectedNodes = handles.mytree.getSelectedNodes;
    handles.mytree.setSelectedNode(handles.root);
    handles.mytree.setSelectedNode(selectedNodes(1));
end



% *************************************************************************
% LOGTEXT BUTTON
% *************************************************************************

% --- Executes on button press in logtextbutton.
function logtextbutton_Callback(hObject, ~, handles)
% hObject    handle to logtextbutton (see GCBO)
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
% write to selected AARAE structures
selectedNodes = handles.mytree.getSelectedNodes;
for nleafs = 1:length(selectedNodes)
    handles.nleafs = nleafs;
    guidata(hObject,handles)
    signaldata = selectedNodes(nleafs).handle.UserData;
    if ~isempty(signaldata)
        row = cell(1,4);
        row{1,1} = datestr(now);
        row{1,2} = 'COMMENT';
        row{1,4} = answer;
        if isfield(signaldata,'history')
            signaldata.history = [signaldata.history;row];
        else
            signaldata.history = row;
        end
        % Save as you go
        delete([cd '/Utilities/Backup/' selectedNodes(nleafs).getName.char '.mat'])
        save([cd '/Utilities/Backup/' selectedNodes(nleafs).getName.char '.mat'], 'signaldata','-v7.3');
        try
            handles.(matlab.lang.makeValidName(char(selectedNodes(nleafs).getName))).UserData = signaldata;
        catch
            warndlg('Sorry - an error occured in writing comment to the tree. Please let Densil know about this.','Bug!');
        end
        selectedParent = selectedNodes(nleafs).getParent;
        handles.mytree.reloadNode(selectedParent);
    end
    handles.mytree.setSelectedNodes(selectedNodes)
end
guidata(hObject,handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over logtextbutton.
function logtextbutton_ButtonDownFcn(hObject, ~, ~)
% hObject    handle to logtextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% *************************************************************************
% KEYBOARD SHORTCUTS
% *************************************************************************

% --- Executes on key press with focus on aarae or any of its controls.
function aarae_WindowKeyPressFcn(hObject, eventdata, handles) %#ok
% hObject    handle to aarae (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if strcmp(eventdata.Modifier,'alt')
    handles.alternate = 1;
else
    handles.alternate = 0;
end
%guidata(hObject,handles)
selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;
if ~isempty(signaldata)
    if strcmp(eventdata.Key,'l') && ~isfield(handles,'legend')
        if ismatrix(signaldata.audio)
            if isfield(signaldata,'chanID')
                handles.legend = legend(handles.axestime,signaldata.chanID);
            end
        end
        if ~ismatrix(signaldata.audio)
            if isfield(signaldata,'bandID')
                handles.legend = legend(handles.axestime,cellstr(num2str(signaldata.bandID')));
            end
        end
    elseif strcmp(eventdata.Key,'l') && isfield(handles,'legend')
        legend(handles.axestime,'off');
        handles = rmfield(handles,'legend');
    end
end
if ~isempty(eventdata.Modifier)
    if strcmp(eventdata.Modifier,'control') == 1
        switch eventdata.Key
            case 'r'
                rec_btn_Callback(hObject, eventdata, handles)
                handles = guidata(hObject);
            case 'g'
                genaudio_btn_Callback(hObject, eventdata, handles)
                handles = guidata(hObject);
            case 'l'
                load_btn_Callback(hObject, eventdata, handles)
                handles = guidata(hObject);
            case 'c'
                calc_btn_Callback(hObject, eventdata, handles)
                handles = guidata(hObject);
            case 'e'
                if strcmp(get(handles.tools_panel,'Visible'),'on')
                    edit_btn_Callback(hObject, eventdata, handles)
                    handles = guidata(hObject);
                end
            case 's'
                if strcmp(get(handles.tools_panel,'Visible'),'on')
                    save_btn_Callback(hObject, eventdata, handles)
                    handles = guidata(hObject);
                end
            case 'delete'
                if strcmp(get(handles.tools_panel,'Visible'),'on')
                    delete_btn_Callback(hObject, eventdata, handles)
                    handles = guidata(hObject);
                end
            case 't'
                logtextbutton_Callback(hObject, eventdata, handles)
                %handles = guidata(hObject);
        end
    end
end
guidata(hObject,handles)













% *************************************************************************
% *************************************************************************
%                  ENDING, CLEARING & EXPORTING FROM AARAE
% *************************************************************************
% *************************************************************************


% *************************************************************************
% FINISH AARAE SESSION
% *************************************************************************

% --- Executes on button press in finish_btn.
function finish_btn_Callback(~, eventdata, handles) %#ok
% hObject    handle to finish_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.export_btn,'Enable'),'on')
    choice = questdlg('Are you sure to want to finish this AARAE session? Unexported data will be lost.',...
        'Exit AARAE',...
        'Yes','No','Export all & exit','Yes');
    switch choice
        case 'Yes'
            uiresume(handles.aarae);
        case 'Export all & exit'
            export_btn_Callback(handles.export_btn,eventdata,handles)
            uiresume(handles.aarae);
    end
else
    uiresume(handles.aarae);
end





% *************************************************************************
% EXPORT ALL DATA FROM AARAE (CREATE A 'PROJECT')
% *************************************************************************

% --- Executes on button press in export_btn.
function export_btn_Callback(hObject, ~, handles)
% hObject    handle to export_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uitreedescription

root = handles.root; % Get selected leaf
root = root(1);
first = root.getFirstChild;
nbranches = root.getChildCount;
branches = cell(nbranches,1);
branches{1,1} = char(first.getValue);
nleaves = 0;
nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{1,1}))(1).getChildCount;
next = first.getNextSibling;
for n = 2:nbranches
    branches{n,1} = char(next.getValue);
    nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{n,1}))(1).getChildCount;
    next = next.getNextSibling;
end
if nleaves == 0
    warndlg('Nothing to export!','AARAE info');
else
    %    leaves = cell(nleaves,1);
    %    i = 0;
    %    for n = 1:size(branches,1)
    %        currentbranch = handles.(matlab.lang.makeValidName(branches{n,1}));
    %        if currentbranch.getChildCount ~= 0
    %            i = i + 1;
    %            first = currentbranch.getFirstChild;
    %            %leafnames(i,:) = first.getName;
    %            leaves{i,1} = char(first.getValue);
    %            next = first.getNextSibling;
    %            if ~isempty(next)
    %                for m = 1:currentbranch.getChildCount-1
    %                    i = i + 1;
    %                    %leafnames(i,:) = next.getName;
    %                    leaves{i,1} = char(next.getValue);
    %                    next = next.getNextSibling;
    %                end
    %            end
    %        end
    %    end
    if ~isdir([cd '/Projects']), mkdir([cd '/Projects']); end
    folder = uigetdir([cd '/Projects'],'Export all');
    if ischar(folder)
        set(hObject,'BackgroundColor','red');
        set(hObject,'Enable','off');
        pause on
        pause(0.001)
        pause off
        %        h = waitbar(0,['1 of ' num2str(size(leaves,1))],'Name','Saving files...');
        %        steps = size(leaves,1);
        %        for i = 1:size(leaves,1)
        %            current = handles.(matlab.lang.makeValidName(leaves{i,:}));
        %            current = current(1);
        %            data = current.handle.UserData; %#ok : used in save
        %            if ~exist([folder '/' leaves{i,:} '.mat'],'file')
        %                try
        %                    save([folder '/' leaves{i,:} '.mat'], 'data');
        %                catch error
        %                    warndlg(error.message,'AARAE info')
        %                end
        %            else
        %                button = questdlg(['A file called ' leaves{i,:} '.mat is already in the destination folder, would you like to replace it?'],...
        %                                  'AARAE info','Yes','No','Append','Append');
        %                switch button
        %                    case 'Yes'
        %                        save([folder '/' leaves{i,:} '.mat'], 'data');
        %                    case 'Append'
        %                        index = 1;
        %                        % This while cycle is just to make sure no signals are
        %                        % overwriten
        %                        while exist([folder '/' leaves{i,:} '_' num2str(index) '.mat'],'file')
        %                            index = index + 1;
        %                        end
        %                        try
        %                            save([folder '/' leaves{i,:} '_' num2str(index) '.mat'], 'data');
        %                        catch error
        %                            warndlg(error.message,'AARAE info')
        %                        end
        %                end
        %            end
        %            waitbar(i / steps,h,sprintf('%d of %d',i,size(leaves,1)))
        %        end
        %        close(h)
        if isdir([cd '/Utilities/Backup'])
            leaves = dir([cd '/Utilities/Backup/*.mat']);
            copyfile([cd '/Utilities/Backup'],folder);
        end
        if isdir([cd '/Utilities/Temp'])
            nfigs = dir([cd '/Utilities/Temp/*.fig']);
            copyfile([cd '/Utilities/Temp'],[folder '/figures']);
        end
        if isfield(handles,'activitylog')
            if isdir([cd '/Log'])
                copyfile([cd '/Log' handles.activitylog],folder);
            end
        end
        addpath(genpath([cd '/Projects']))
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Exported ' num2str(size(leaves,1)) ' data files and ' num2str(size(nfigs,1)) ' figures to "%s" \n\n'],folder);
        set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
        set(hObject,'Enable','off');
    else
        addpath(genpath([cd '/Projects']))
    end
end
guidata(hObject,handles)






% *************************************************************************
% CLEAR ALL FROM THE AARAE WORKSPACE
% *************************************************************************

% --- Executes on button press in clrall_btn.
function clrall_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to clrall_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
root = handles.root; % Get selected leaf
root = root(1);
first = root.getFirstChild;
nbranches = root.getChildCount;
branches = cell(nbranches,1);
branches{1,1} = char(first.getValue);
nleaves = 0;
nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{1,1}))(1).getChildCount;
next = first.getNextSibling;
for n = 2:nbranches
    branches{n,1} = char(next.getValue);
    nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{n,1}))(1).getChildCount;
    next = next.getNextSibling;
end
leaves = cell(nleaves,1);
i = 0;
for n = 1:size(branches,1)
    currentbranch = handles.(matlab.lang.makeValidName(branches{n,1}));
    if currentbranch.getChildCount ~= 0
        i = i + 1;
        first = currentbranch.getFirstChild;
        %leafnames(i,:) = first.getName;
        leaves{i,:} = char(first.getValue);
        next = first.getNextSibling;
        if ~isempty(next)
            for m = 1:currentbranch.getChildCount-1
                i = i + 1;
                %leafnames(i,:) = next.getName;
                leaves{i,:} = char(next.getValue);
                next = next.getNextSibling;
            end
        end
    end
end
if nleaves == 0
    warndlg('Nothing to delete!','AARAE info');
else
    %leafnames = char(leafnames);
    %leaves = char(leaves);
    delete = questdlg('Current workspace will be cleared, would you like to proceed?',...
        'Warning',...
        'Yes','No','Yes');
    switch delete
        case 'Yes'
            set(hObject,'BackgroundColor','red');
            set(hObject,'Enable','off');
            for i = 1:size(leaves,1)
                current = handles.(matlab.lang.makeValidName(leaves{i,1}));
                handles.mytree.remove(current);
                handles = rmfield(handles,matlab.lang.makeValidName(leaves{i,1}));
            end
            handles.mytree.reloadNode(handles.root);
            handles.mytree.setSelectedNode(handles.root);
            rmpath([cd '/Utilities/Temp']);
            rmdir([cd '/Utilities/Temp'],'s');
            rmpath([cd '/Utilities/Backup']);
            rmdir([cd '/Utilities/Backup'],'s');
            mkdir([cd '/Utilities/Temp']);
            mkdir([cd '/Utilities/Backup']);
            addpath([cd '/Utilities/Temp']);
            addpath([cd '/Utilities/Backup']);
            set(handles.result_box,'Value',1);
            set(handles.result_box,'String',cell(1,1));
            fprintf(handles.fid, ['%% ' datestr(now,16) ' - Cleared workspace \n\n']);
            set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
            set(hObject,'Enable','off');
            set(handles.export_btn,'Enable','off');
    end
end
guidata(hObject,handles)


% --- Executes on button press in CloseFiguresButton.
function CloseFiguresButton_Callback(hObject, eventdata, handles)
% hObject    handle to CloseFiguresButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('type','figure','-not','tag','aarae');
delete(h)



% *************************************************************************
% CLOSE AARAE
% *************************************************************************

% --- Executes when user attempts to close aarae.
function aarae_CloseRequestFcn(hObject,eventdata,handles) %#ok
% hObject    handle to aarae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.export_btn,'Enable'),'on')
    choice = questdlg('Are you sure to want to finish this AARAE session? Unexported data will be lost.',...
        'Exit AARAE',...
        'Yes','No','Export all & exit','Yes');
    switch choice
        case 'Yes'
            if getappdata(handles.aarae,'waiting')
                % The GUI is still in UIWAIT, so call UIRESUME and return
                uiresume(hObject);
                setappdata(handles.aarae,'waiting',0);
            else
                % The GUI is no longer waiting, so destroy it now.
                delete(hObject);
            end
            %   uiresume(handles.aarae);
        case 'Export all & exit'
            export_btn_Callback(handles.export_btn,eventdata,handles)
            if getappdata(handles.aarae,'waiting')
                % The GUI is still in UIWAIT, so call UIRESUME and return
                uiresume(hObject);
                setappdata(handles.aarae,'waiting',0);
            else
                % The GUI is no longer waiting, so destroy it now.
                delete(hObject);
            end
            %uiresume(handles.aarae);
    end
else
    uiresume(handles.aarae);
end
% Check appdata flag to see if the main GUI is in a wait state
%if getappdata(handles.aarae,'waiting')
% The GUI is still in UIWAIT, so call UIRESUME and return
%    uiresume(hObject);
%    setappdata(handles.aarae,'waiting',0);
%else
%    % The GUI is no longer waiting, so destroy it now.
%    delete(hObject);
%end





% *************************************************************************
% CLEAN UP AFTER CLOSING GUI
% *************************************************************************

% --- Outputs from this function are returned to the command line.
function varargout = aarae_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to aarae
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clean up the mess after closing the GUI
rmappdata(0,'hMain');
rmpath(genpath(cd));
rmdir([cd '/Utilities/Temp'],'s');
rmdir([cd '/Utilities/Backup'],'s');
fprintf(handles.fid, ['\n%% - End of AARAE session - ' datestr(now)]);
fclose('all');
if ~isempty(handles.aarae)
    delete(handles.aarae);
end
java.lang.Runtime.getRuntime.gc % Java 'garbage collection'
% Get default command line output from handles structure
varargout{1} = [];









% *************************************************************************
% *************************************************************************
%                THE 'START' BUTTONS: ACQUIRING AUDIO DATA
% *************************************************************************
% *************************************************************************



% *************************************************************************
% GENERATE AUDIO
% *************************************************************************

% --- Executes on button press in genaudio_btn.
function genaudio_btn_Callback(hObject, ~, handles)
% hObject    handle to genaudio_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.recovery_txt,'Visible'),'on'), set(handles.recovery_txt,'Visible','off'); end
% Call the 'desktop'
hMain = getappdata(0,'hMain');
setappdata(hMain,'testsignal',[]);

% Call the window that allows signal generation
newleaf = genaudio('main_stage1', handles.aarae);
%set(handles.aarae,'CurrentObject',[])
% Update the tree with the generated signal
handles.mytree.setSelectedNode(handles.root);
if ~isempty(getappdata(hMain,'testsignal'))
    signaldata = getappdata(hMain,'testsignal');
    if isfield(signaldata,'tag'), signaldata = rmfield(signaldata,'tag'); end
    signaldata.datatype = 'testsignals';
    if isfield(signaldata,'audio')
        %signaldata.nbits = 16;
        iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
    else
        iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/notesicon.gif');
    end
    signaldata = checkcal(signaldata);
    signaldata = addhistory(signaldata,['Generated by ' handles.Settings.username ' using ' getappdata(hMain,'AARAEversion')]);
    leafname = isfield(handles,matlab.lang.makeValidName(newleaf));
    if leafname == 1
        index = 1;
        % This while cycle is just to make sure no signals are
        % overwriten
        if length(matlab.lang.makeValidName([newleaf,'_',num2str(index)])) >= namelengthmax-2, newleaf = newleaf(1:round(end/2)); end
        while isfield(handles,matlab.lang.makeValidName([newleaf,'_',num2str(index)])) == 1
            index = index + 1;
        end
        newleaf = [newleaf,'_',num2str(index)];
    end
    
    % Save as you go
    save([cd '/Utilities/Backup/' newleaf '.mat'], 'signaldata','-v7.3');
    
    % Generate new leaf
    signaldata.name = matlab.lang.makeValidName(newleaf);
    handles.(signaldata.name) = uitreenode('v0', newleaf,  newleaf,  iconPath, true);
    handles.(signaldata.name).UserData = signaldata;
    handles.testsignals.add(handles.(signaldata.name));
    handles.mytree.reloadNode(handles.testsignals);
    handles.mytree.expand(handles.testsignals);
    handles.mytree.setSelectedNode(handles.(signaldata.name));
    set([handles.clrall_btn,handles.export_btn],'Enable','on')
    % Log event
    fprintf(handles.fid, ['%% ' datestr(now,16) ' - Generated ' newleaf ': duration = ' num2str(length(signaldata.audio)/signaldata.fs) ' s ; fs = ' num2str(signaldata.fs) ' Hz; size = ' num2str(size(signaldata.audio)) '\n']);
    % Log verbose metadata
    logaudioleaffields(signaldata,0);
end
java.lang.Runtime.getRuntime.gc % Java garbage collection
guidata(hObject, handles);






% *************************************************************************
% LOAD DATA INTO AARAE
% *************************************************************************

% --- Executes on button press in load_btn.
function load_btn_Callback(hObject, ~, handles)
% hObject    handle to load_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call the 'desktop'
hMain = getappdata(0,'hMain');
if strcmp(get(handles.recovery_txt,'Visible'),'on'), set(handles.recovery_txt,'Visible','off'); end
% Get path of the file to load
[filename,handles.defaultaudiopath,filterindex] = uigetfile(...
    {'*.wav;*.mat;.WAV;.MAT','Test Signals (*.wav,*.mat)';...
    '*.wav;*.mat;.WAV;.MAT','Measurement files (*.wav,*.mat)';...
    '*.wav;*.mat;.WAV;.MAT','Processed files (*.wav,*.mat)';...
    '*.wav;*.mat;.WAV;.MAT','Result files (*.wav,*.mat)'},...
    'Select audio file',handles.defaultaudiopath,...
    'MultiSelect','on');

if ~iscell(filename)
    if ischar(filename), filename = cellstr(filename);
    else return
    end
end
h = waitbar(0,['Loading 1 of ' num2str(length(filename)) ' files'],'Name','AARAE info','WindowStyle','modal','Tag','progressbar');
steps = length(filename);
for i = 1:length(filename)
    if filename{i} ~= 0
        newleaf = cell(1,1);
        [~,newleaf{1,1},ext] = fileparts(filename{i});
        % Check type of file. First 'if' is for .mat, second is for .wav
        if strcmp(ext,'.mat') || strcmp(ext,'.MAT')
            file = importdata(fullfile(handles.defaultaudiopath,filename{i}));
            if isstruct(file)
                signaldata = file;
            else
                specs = inputdlg('Please specify the sampling frequency','Sampling frequency',1);
                if (isempty(specs))
                    warndlg('Input field is blank, cannot load data!','AARAE info');
                    %signaldata = [];
                    return;
                else
                    fs = str2double(specs{1,1});
                    %nbits = str2double(specs{2,1});
                    if isnan(fs) || fs<=0 % || isnan(nbits) || nbits<=0)
                        warndlg('Input MUST be a real positive number, cannot load data!','AARAE info');
                        %signaldata = [];
                        return;
                    else
                        signaldata = [];
                        signaldata.audio = file;
                        signaldata.fs = fs;
                        %signaldata.nbits = nbits;
                    end
                end
            end
        end
        if strcmp(ext,'.wav') || strcmp(ext,'.WAV')
            signaldata = [];
            [signaldata.audio,signaldata.fs] = audioread(fullfile(handles.defaultaudiopath,filename{i}));
            %signaldata.nbits = 16;
        end;
        
        % Generate new leaf and update the tree
        if ~isempty(signaldata)
            if ~isfield(signaldata,'chanID') && isfield(signaldata,'audio')
                signaldata.chanID = cellstr([repmat('Chan',size(signaldata.audio,2),1) num2str((1:size(signaldata.audio,2))')]);
            end
            if ~isfield(signaldata,'datatype') || (isfield(signaldata,'datatype') && strcmp(signaldata.datatype,'IR'))
                if filterindex == 1, signaldata.datatype = 'testsignals'; end;
                if filterindex == 2, signaldata.datatype = 'measurements'; end;
                if filterindex == 3, signaldata.datatype = 'processed'; end;
                if filterindex == 4, signaldata.datatype = 'results'; end;
            end
            signaldata = checkcal(signaldata);
            signaldata = addhistory(signaldata,['Loaded using ' getappdata(hMain,'AARAEversion')]);
            if isfield(signaldata,'audio') && ~strcmp(signaldata.datatype,'syscal')
                iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
            elseif strcmp(signaldata.datatype,'syscal')
                iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/boardicon.gif');
            elseif isfield(signaldata,'data')
                iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/help_fx.png');
            elseif isfield(signaldata,'tables')
                iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/HDF_grid.gif');
            else
                iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/notesicon.gif');
            end
            leafname = isfield(handles,matlab.lang.makeValidName(newleaf{1,1}));
            if leafname == 1
                index = 1;
                % This while cycle is just to make sure no duplicate names
                if length(matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) >= namelengthmax-2, newleaf{1,1} = newleaf{1,1}(1:round(end/2)); end
                while isfield(handles,matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) == 1
                    index = index + 1;
                end
                newleaf{1,1} = [newleaf{1,1},'_',num2str(index)];
            end
            signaldata.name = matlab.lang.makeValidName(newleaf{1,1});
            
            % Save as you go
            save([cd '/Utilities/Backup/' newleaf{1,1} '.mat'], 'signaldata','-v7.3');
            
            handles.(signaldata.name) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
            handles.(signaldata.name).UserData = signaldata;
            if strcmp(signaldata.datatype,'syscal'), signaldata.datatype = 'measurements'; end
            handles.(matlab.lang.makeValidName(signaldata.datatype)).add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
            handles.mytree.reloadNode(handles.(matlab.lang.makeValidName(signaldata.datatype)));
            handles.mytree.expand(handles.(matlab.lang.makeValidName(signaldata.datatype)));
            %handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(newleaf{1,1})));
            set([handles.clrall_btn,handles.export_btn],'Enable','on')
            fprintf(handles.fid, ['%% ' datestr(now,16) ' - Loaded "' filename{i} '" to branch "' char(handles.(matlab.lang.makeValidName(signaldata.datatype)).getName) '"\n']);
            % Log verbose metadata
            logaudioleaffields(signaldata);
        end
        guidata(hObject, handles);
    end
    try
        waitbar(i / steps,h,sprintf('Loading %d of %d files',i,length(filename)));
    catch
        warndlg('Loading interrupted by user!','AARAE info')
        break
    end
end
if i == steps, close(h); end
handles.mytree.setSelectedNode(handles.(signaldata.name));




% *************************************************************************
% RECORD AUDIO
% *************************************************************************
% Launch the audio_recorder & deal with the output

% --- Executes on button press in rec_btn.
function rec_btn_Callback(hObject, ~, handles)
% hObject    handle to rec_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.recovery_txt,'Visible'),'on'), set(handles.recovery_txt,'Visible','off'); end
% Call the 'desktop'
hMain = getappdata(0,'hMain');

%set(handles.datatypetext,'String','No signal loaded');

% Call the audio recorder window
audiodata = audio_recorder('main_stage1', handles.aarae);

% Generate new leaf and update tree with the recording
handles.mytree.setSelectedNode(handles.root);
newleaf = getappdata(hMain,'signalname');
savenewsyscalstats = getappdata(hMain,'savenewsyscalstats');
if savenewsyscalstats == 1
    handles.syscalstats = getappdata(hMain,'syscalstats');
    handles.mytree.setSelectedNode(handles.root);
    iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/boardicon.gif');
    handles.syscalstats.datatype = 'syscal';
    funname = ['System_calibration ' strrep(datestr(rem(now,1)),':','_')];
    leafname = isfield(handles,matlab.lang.makeValidName(funname));
    if leafname == 1
        index = 1;
        % This while cycle is just to make sure no signals are
        % overwriten
        while isfield(handles,matlab.lang.makeValidName([funname,'_',num2str(index)])) == 1
            index = index + 1;
        end
        funname = [funname,'_',num2str(index)];
    end
    
    % Save as you go
    tempsyscalstats = handles.syscalstats; %#ok : Used in followring line
    save([cd '/Utilities/Backup/' funname '.mat'], 'tempsyscalstats','-v7.3');
    
    handles.(matlab.lang.makeValidName(funname)) = uitreenode('v0', funname,  funname,  iconPath, true);
    handles.(matlab.lang.makeValidName(funname)).UserData = handles.syscalstats;
    handles.measurements.add(handles.(matlab.lang.makeValidName(funname)));
    handles.mytree.reloadNode(handles.measurements);
    handles.mytree.expand(handles.measurements);
    handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(funname)));
    set([handles.clrall_btn,handles.export_btn],'Enable','on')
    fprintf(handles.fid, ['%% ' datestr(now,16) ' - Saved system calibration data ' handles.funname '\n\n']);
end
if ~isempty(audiodata)
    audiodata.datatype = 'measurements';
    iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
    leafname = isfield(handles,matlab.lang.makeValidName(newleaf));
    if leafname == 1
        index = 1;
        % This while cycle is just to make sure no signals are
        % overwriten
        while isfield(handles,matlab.lang.makeValidName([newleaf,'_',num2str(index)])) == 1
            index = index + 1;
        end
        newleaf = [newleaf,'_',num2str(index)];
    end
    
    audiodata_fields = fieldnames(audiodata);
    audiodata_emptyfields = structfun(@isempty,audiodata);
    audiodata = rmfield(audiodata,audiodata_fields(audiodata_emptyfields));
    audiodata = checkcal(audiodata);
    audiodata = addhistory(audiodata,['Loaded into ' getappdata(hMain,'AARAEversion')]);
    audiodata.name = matlab.lang.makeValidName(newleaf);
    % Save as you go
    save([cd '/Utilities/Backup/' newleaf '.mat'], 'audiodata','-v7.3');
    
    handles.(audiodata.name) = uitreenode('v0', newleaf,  newleaf,  iconPath, true);
    handles.(audiodata.name).UserData = audiodata;
    handles.measurements.add(handles.(audiodata.name));
    handles.mytree.reloadNode(handles.measurements);
    handles.mytree.expand(handles.measurements);
    handles.mytree.setSelectedNode(handles.(audiodata.name));
    set([handles.clrall_btn,handles.export_btn],'Enable','on')
    fprintf(handles.fid, ['%% ' datestr(now,16) ' - Recorded "' newleaf '": duration = ' num2str(length(audiodata.audio)/audiodata.fs) 's\n']);
    % Log verbose metadata
    logaudioleaffields(audiodata);
end
handles.alternate=0; %
java.lang.Runtime.getRuntime.gc % Java garbage collection
guidata(hObject, handles);





% *************************************************************************
% 123 CALCULATORS BUTTON
% *************************************************************************
% Launch the Calculators GUI.

% --- Executes on button press in calc_btn.
function calc_btn_Callback(~, ~, handles)
% hObject    handle to calc_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.recovery_txt,'Visible'),'on'), set(handles.recovery_txt,'Visible','off'); end
% Call the window that displays calculators
calculator('main_stage1', handles.aarae);

% Handles update is done inside calculator.m












% *************************************************************************
% *************************************************************************
%                           THE 'TOOLS' BUTTONS
% *************************************************************************
% *************************************************************************
% edit, save, delete, convolve, calibrate and compare



% *************************************************************************
% EDIT AUDIO
% *************************************************************************

% --- Executes on button press in edit_btn.
function edit_btn_Callback(~, ~, handles)
% hObject    handle to edit_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the 'desktop'
hMain = getappdata(0,'hMain');
signaldata = getappdata(hMain,'testsignal');

if isempty(signaldata), warndlg('No signal loaded!','Whoops...!');
    %elseif ndims(signaldata.audio) > 2, warndlg('Cannot edit file!','Whoops...!');
else
    % Call editing window
    selectedNodes = handles.mytree.getSelectedNodes;
    selectedNodes = selectedNodes(1);
    fprintf(handles.fid, '%% ***********************************************\n');
    fprintf(handles.fid, ['%% ' datestr(now,16) ' - Opened Edit window with "' char(selectedNodes.getName) '"\n\n']);
    [xi,xf] = edit_signal('main_stage1', handles.aarae);
    % Update tree with edited signal
    if ~isempty(xi) && ~isempty(xf)
        signaldata = getappdata(hMain,'testsignal');
        signaldata = checkcal(signaldata);
        signaldata = addhistory(signaldata,'Edited');
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Edited "' char(selectedNodes.getName) '": cropped from %fs to %fs (at input fs); new duration = ' num2str(length(signaldata.audio)/signaldata.fs) ' s\n\n'],xi,xf);
        fprintf(handles.fid, '%% ***********************************************\n');
        signaldata.name = char(selectedNodes.getName); % name is set
        %within edit_signal, so we should not need to set it here.
    else
        handles.mytree.setSelectedNode(handles.root);
    end
end








% *************************************************************************
% SAVE DATA FROM AARAE
% *************************************************************************

% --- Executes on button press in save_btn.
function save_btn_Callback(hObject, ~, handles)
% hObject    handle to save_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the 'desktop'
%hMain = getappdata(0,'hMain');
%audiodata = getappdata(hMain,'testsignal');
selectedNodes = handles.mytree.getSelectedNodes;
for nleafs = 1:length(selectedNodes)
    audiodata = selectedNodes(nleafs).handle.UserData;
    if ~isempty(audiodata) %Check if there's data to save
        name = inputdlg('File name: (Please specify .wav for wave files)','Save as MATLAB File',1,{[char(selectedNodes(nleafs).getName) '.mat']}); %Request file name
        if ~isempty(name)
            %name = name{1,1};
            [~,name{1,1},ext]=fileparts(name{1,1});
            if strcmp(ext,'.mat'), ensave = 1;
            elseif strcmp(ext,'.wav') && ismatrix(audiodata.audio), ensave = 1;
            elseif strcmp(ext,'.wav') && ~ismatrix(audiodata.audio), ensave = 1; ext = '.mat';
            elseif isempty(ext), ensave = 1;
            else ensave = 0;
            end
        else
            return
        end
        if isempty(name{1,1}) || ensave == 0
            warndlg('No data saved','AARAE info');
            return
        else
            if isempty(ext), ext = '.mat'; end
            if strcmp(ext,'.wav') && (~isfield(audiodata,'audio') || ~isfield(audiodata,'fs')), ext = '.mat'; end
            folder_name = uigetdir(handles.defaultaudiopath,'Save AARAE file');
            if ischar(folder_name)
                handles.defaultaudiopath = folder_name;
                listing = dir([folder_name '/' name{1,1} ext]);
                if isempty(listing)
                    if strcmp(ext,'.mat'), save([folder_name '/' name{1,1} ext],'audiodata','-v7.3'); end
                    if strcmp(ext,'.wav'), audiowrite([folder_name '/' name{1,1} ext],audiodata.audio,audiodata.fs); end
                else
                    index = 1;
                    % This while cycle is just to make sure no signals are
                    % overwriten
                    while isempty(dir([name{1,1},'_',num2str(index),ext])) == 0
                        index = index + 1;
                    end
                    name{1,1} = [name{1,1},'_',num2str(index),ext];
                    if strcmp(ext,'.mat'), save([folder_name '/' name{1,1} ext],'audiodata','-v7.3'); end
                    if strcmp(ext,'.wav'), audiowrite(audiodata.audio,audiodata.fs,[folder_name '/' name{1,1}]); end
                end
                %current = cd;
                fprintf(handles.fid, ['%% ' datestr(now,16) ' - Saved "' char(selectedNodes(nleafs).getName) '" to file "' name{1,1} ext '" in folder "%s"' '\n'],folder_name);
            end
        end
    end
end
guidata(hObject, handles);




% *************************************************************************
% DELETE DATA FROM AARAE
% *************************************************************************

% --- Executes on button press in delete_btn.
function delete_btn_Callback(hObject, ~, handles)
% hObject    handle to delete_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the 'desktop'
hMain = getappdata(0,'hMain');
deleteans = questdlg('Current data will be lost, would you like to proceed?',...
    'Warning',...
    'Yes','No','Yes');
switch deleteans
    case 'Yes'
        % Deletes selected leaf from the tree
        setappdata(hMain,'testsignal',[]);
        selectedNodes = handles.mytree.getSelectedNodes;
        for nleafs = 1:length(selectedNodes)
            audiodata = selectedNodes(nleafs).handle.UserData;
            if strcmp(audiodata.datatype,'syscal')
                handles = rmfield(handles,'syscalstats');
                set(handles.signaltypetext,'String',[])
            end
            if ~isempty(audiodata)
                
                % Delete from backup files
                filename = [cd '/Utilities/Backup/' selectedNodes(nleafs).getName.char '.mat'];
                delete(filename)
                
                selectedParent = selectedNodes(nleafs).getParent;
                handles.mytree.remove(selectedNodes(nleafs));
                handles.mytree.reloadNode(selectedParent);
                handles.mytree.setSelectedNode(handles.root);
                handles = rmfield(handles,matlab.lang.makeValidName(char(selectedNodes(nleafs).getName)));
                fprintf(handles.fid, ['%% ' datestr(now,16) ' - Deleted "' char(selectedNodes(nleafs).getName) '" from branch "' char(selectedParent.getName) '"\n\n']);
            end
        end
        guidata(hObject, handles);
    case 'No'
        guidata(hObject, handles);
end
java.lang.Runtime.getRuntime.gc % Java garbage collection







% *************************************************************************
% CONVOLVE AUDIO WITH AUDIO2
% *************************************************************************

% --- Executes on button press in IR_btn (convolve audio with audio2).
function IR_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to IR_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'BackgroundColor','red');
set(hObject,'Enable','off');

hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');
selectedNodes = handles.mytree.getSelectedNodes;
if ~isfield(audiodata,'audio2'), audiodata.audio2 = []; end
if isempty(audiodata.audio2)
    % there is no secondary audio field, so let's get something as a
    % substitute - this is envisaged to be expecially useful for transfer
    % function based calculations (where the recorded audio is compared
    % with the original audio). Currently this is restricted to only one
    % column of audio.
    selection = choose_audio;
    if isempty(selection)
        handles.alternate = 0;
        set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
        set(hObject,'Enable','on');
        java.lang.Runtime.getRuntime.gc % Java garbage collection
        guidata(hObject, handles);
        return
    else
        [~,chans,bands,dim4,dim5,dim6]= size(selection.audio);
        if chans + bands + dim4 + dim5 + dim6 > 5
            param = inputdlg({['Select channel (up to ' num2str(chans) ')'];...
                ['Select band (up to ' num2str(bands) ')'];...
                ['Select dim 4 cycle (up to ' num2str(dim4) ')'];...
                ['Select dim 5 (up to ' num2str(dim5) ')'];...
                ['Select dim 6 (up to ' num2str(dim6) ')']},...
                'Select just 1 column of audio',... % window title.
                [1 60],...
                {'1';'1';'1';'1';'1'}); % Preset answers
            
            param = round(abs(str2num(char(param))));
            if length(param) ~= 5
                handles.alternate = 0;
                set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
                set(hObject,'Enable','on');
                java.lang.Runtime.getRuntime.gc % Java garbage collection
                guidata(hObject, handles);
                return
            else
                if param(1) < 1 || param(1) > chans, param(1) = 1; end
                if param(2) < 1 || param(2) > bands, param(2) = 1; end
                if param(3) < 1 || param(3) > dim4, param(3) = 1; end
                if param(4) < 1 || param(4) > dim5, param(4) = 1; end
                if param(5) < 1 || param(5) > dim6, param(5) = 1; end
                selection.audio = selection.audio(:,param(1),param(2),param(3),param(4),param(5)); % additional audio data
            end
        end
        if selection.fs ~= audiodata.fs
            gcd_fs = gcd(audiodata.fs,selection.fs); % greatest common denominator
            audiodata.audio2 = resample(selection.audio,audiodata.fs/gcd_fs,selection.fs/gcd_fs);
        else
            audiodata.audio2 = selection.audio;
        end
    end
    handles.alternate = 1; % we must choose the processing method.
end

if handles.alternate ~= 1
    % convolveaudiowithaudio2 is an AARAE utility
    h = msgbox('Computing impulse response, please wait...','AARAE info','modal');
    [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,0,1);
    calcmethod = 1;
    nameprefix = 'IR_';
else
    % select an alternative operation
    str = {'1. Convolve audio with audio2';... % same as default
        '2. Cross-correlate audio with audio2';...
        '3. Circular convolution of audio with audio2';...
        '4. Circular cross-correlation of audio with audio2';...
        '5. Transfer function from audio2 to audio (-200 dB threshold)';...
        '6. Transfer function from audio2 to audio (-90 dB threshold)';...
        '7. Transfer function from audio2 to audio (-80 dB threshold)';...
        '8. Transfer function from audio2 to audio (-70 dB threshold)';...
        '9. Transfer function from reversed audio2 to audio (-200 dB threshold)';...
        '10. Transfer function from reversed audio2 to audio (-90 dB threshold)';...
        '11. Transfer function from reversed audio2 to audio (-80 dB threshold)';...
        '12. Transfer function from reversed audio2 to audio (-70 dB threshold)';...
        '13. Transfer function from audio to audio2 (-200 dB threshold)';...
        '14. Transfer function from audio to audio2 (-90 dB threshold)';...
        '15. Transfer function from audio to audio2 (-80 dB threshold)';...
        '16. Transfer function from audio to audio2 (-70 dB threshold)';...
        '17. Transfer function from reversed audio to audio2 (-200 dB threshold)';...
        '18. Transfer function from reversed audio to audio2 (-90 dB threshold)';...
        '19. Transfer function from reversed audio to audio2 (-80 dB threshold)';...
        '20. Transfer function from reversed audio to audio2 (-70 dB threshold)';...
        '21. Time domain deconvolution of audio2 from audio';...
        '22. Time domain deconvolution of audio from audio2';...
        '23. Time domain deconvolution of time-reversed audio2 from audio';...
        '24. Time domain deconvolution of time-reversed audio from audio2';...
        '25. Time domain convolution of audio with audio2';...
        '26. Time domain convolution of audio with time-reversed audio2'};
    
    
    [calcmethod,ok] = listdlg('PromptString','Select the processing method',...
        'SelectionMode','single',...
        'ListString',str,...
        'ListSize', [400,400],...
        'CancelString','Alternative audio2');
    if ~ok
        selection = choose_audio; % select alternative to audio2
        selection = choose_from_higher_dimensions(selection,1,1);
        audiodata.audio2 = selection.audio(:,1,1,1,1,1);
        if selection.fs ~= audiodata.fs
            audiodata.audio2 = resample(audiodata.audio2,audiodata.fs,selection.fs);
        end
        [calcmethod,ok] = listdlg('PromptString','Select the processing method',...
            'SelectionMode','single',...
            'ListString',str,...
            'ListSize', [400,400]);
    end
    if ok
        h = msgbox('Computing impulse response (or other cross function), please wait...','AARAE info','modal');
        % The following 'switch' is unnecessary at present, but may be
        % useful if each method is customised further (e.g. the name
        % prefix)
        switch calcmethod
            case 1
                % same as normal method
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,0,1);
                nameprefix = 'IR_';
            case 2
                % cross-correlate audio with audio2
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,2,1);
                nameprefix = 'X_';
            case 3
                % circular convolution
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,3,1);
                nameprefix = 'IR_';
            case 4
                % circularly cross-correlate audio with audio2
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,4,1);
                nameprefix = 'X_';
            case 5
                % Transfer function from audio2 to audio (-200 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,5,1);
                nameprefix = 'IR_';
            case 6
                % Transfer function from audio2 to audio (-90 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,6,1);
                nameprefix = 'IR_';
            case 7
                % Transfer function from audio2 to audio (-80 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,7,1);
                nameprefix = 'IR_';
            case 8
                % Transfer function from audio2 to audio (-70 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,8,1);
                nameprefix = 'IR_';
            case 9
                % Transfer function from reversed audio2 to audio (-200 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,9,1);
                nameprefix = 'IR_';
            case 10
                % Transfer function from reversed audio2 to audio (-90 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,10,1);
                nameprefix = 'IR_';
            case 11
                % Transfer function from reversed audio2 to audio (-80 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,11,1);
                nameprefix = 'IR_';
            case 12
                % Transfer function from reversed audio2 to audio (-70 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,12,1);
                nameprefix = 'IR_';
            case 13
                % 13. Transfer function from audio to audio2 (-200 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,13,1);
                nameprefix = 'IR_';
            case 14
                % 14. Transfer function from audio to audio2 (-90 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,14,1);
                nameprefix = 'IR_';
            case 15
                % 15. Transfer function from audio to audio2 (-80 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,15,1);
                nameprefix = 'IR_';
            case 16
                % 16. Transfer function from audio to audio2 (-70 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,16,1);
                nameprefix = 'IR_';
            case 17
                % 17. Transfer function from reversed audio to audio2 (-200 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,17,1);
                nameprefix = 'IR_';
            case 18
                % 18. Transfer function from reversed audio to audio2 (-90 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,18,1);
                nameprefix = 'IR_';
            case 19
                % 19. Transfer function from reversed audio to audio2 (-80 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,19,1);
                nameprefix = 'IR_';
            case 20
                % 20. Transfer function from reversed audio to audio2 (-70 dB threshold)
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,20,1);
                nameprefix = 'IR_';
            case 21
                % 21. Time domain deconvolution of audio2 from audio
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,21,1);
                nameprefix = 'X21_';
            case 22
                % 22. Time domain deconvolution of audio from audio2
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,22,1);
                nameprefix = 'X22_';
            case 23
                % 23. Time domain deconvolution of time-reversed audio2 from audio
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,23,1);
                nameprefix = 'X23_';
            case 24
                % 24. Time domain deconvolution of time-reversed audio from audio2
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,24,1);
                nameprefix = 'X24_';
            case 25
                % 25. Time domain convolution of audio with audio2
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,25,1);
                nameprefix = 'X25_';
            case 26
                % 26. Time domain convolution of audio with time-reversed audio2
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,26,1);
                nameprefix = 'X26_';
                
                % Auto-functions are not currently implemented here
                %             case 27
                %                 [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,27,1);
                %                 nameprefix = 'X27_';
                %             case 28
                %                 [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,28,1);
                %                 nameprefix = 'X28_';
                %             case 29
                %                 [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,29,1);
                %                 nameprefix = 'X29_';
                %             case 30
                %                 [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,30,1);
                %                 nameprefix = 'X30_';
                %             case 31
                %                 [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,31,1);
                %                 nameprefix = 'X31_';
            otherwise
                [IR,method,scalingmethod] = convolveaudiowithaudio2(audiodata,[],0,0,1);
                nameprefix = 'IR_';
        end
    else
        handles.alternate = 0;
        set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
        set(hObject,'Enable','on');
        java.lang.Runtime.getRuntime.gc % Java garbage collection
        guidata(hObject, handles);
        return
    end
end
if ishandle(h), close(h); end
if ~isempty(IR)
    % Trim the IR. tempIR is used for visual display in window_signal
    if ~ismatrix(IR), tempIR(:,:) = IR(:,1,:,end,1,1); else tempIR = IR; end
    [trimsamp_low,trimsamp_high,chanind,bandind,cycind,outchanind,dim6ind] = ...
        window_signal('main_stage1', handles.aarae,...
        'IR',tempIR,...
        'fs',audiodata.fs,...
        'audio2len',size(audiodata.audio2,1),...
        'chans',size(IR,2),...
        'bands',size(IR,3),...
        'cycles',size(IR,4),...
        'outchans',size(IR,5),...
        'dim6',size(IR,6));  % Calls the trimming GUI window to trim the IR
    if isempty(trimsamp_low)
        % the 'cancel' button makes trimsamp_low empty
        handles.alternate = 0;
        set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
        set(hObject,'Enable','on');
        java.lang.Runtime.getRuntime.gc % Java garbage collection
        guidata(hObject, handles);
        return
    end
    IR = IR(trimsamp_low:trimsamp_high,chanind,bandind,cycind,outchanind,dim6ind);
    IRlength = length(IR);
    
    % Create new leaf and update the tree
    handles.mytree.setSelectedNode(handles.root);
    newleaf = matlab.lang.makeValidName([nameprefix selectedNodes(1).getName.char]);
    leafname = isfield(handles,matlab.lang.makeValidName(newleaf));
    if leafname == 1
        index = 1;
        % This while cycle is just to make sure no signals are
        % overwriten
        if length(matlab.lang.makeValidName([newleaf,'_',num2str(index)])) >= namelengthmax-2, newleaf = newleaf(1:round(end/2)); end
        while isfield(handles,matlab.lang.makeValidName([newleaf,'_',num2str(index)])) == 1
            index = index + 1;
        end
        newleaf = matlab.lang.makeValidName([newleaf,'_',num2str(index)]);
    end
    if ~isempty(getappdata(hMain,'testsignal'))
        signaldata = audiodata;
        signaldata.name = newleaf;
        signaldata = rmfield(signaldata,'audio2');
        signaldata.audio = IR;
        %signaldata.fs = fs; % This should be unnecessary
        %signaldata.nbits = 16;
        
        if isfield(signaldata,'chanID')
            try
                signaldata.chanID = signaldata.chanID{chanind,1};
                if length(signaldata.chanID{:,1}) ~= size(signaldata.audio,2)
                    signaldata.chanID = cellstr([repmat('Chan',size(signaldata.audio,2),1) num2str((1:size(signaldata.audio,2))')]);
                end
            catch
                signaldata.chanID = cellstr([repmat('Chan',size(signaldata.audio,2),1) num2str((1:size(signaldata.audio,2))')]);
            end
        else
            signaldata.chanID = cellstr([repmat('Chan',size(signaldata.audio,2),1) num2str((1:size(signaldata.audio,2))')]);
        end
        
        if isfield(signaldata,'bandID')
            try
                signaldata.bandID = signaldata.bandID(bandind);
                if length(signaldata.bandID) ~= size(signaldata.audio,3)
                    signaldata.bandID = 1:size(signaldata.audio,3);
                end
            catch
                signaldata.bandID = 1:size(signaldata.audio,3);
            end
        end
        
        if isfield(signaldata,'properties')
            if isfield(signaldata.properties,'relgain')
                if method == 1 || method >= 5
                    signaldata.properties = rmfield(signaldata.properties,'relgain');
                end
            end
            if isfield(signaldata.properties,'startflag')
                if method == 1 || method >= 5
                    signaldata.properties = rmfield(signaldata.properties,'startflag');
                end
            end
        end
        
        signaldata.datatype = 'measurements';
        iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
        
        % Save as you go
        save([cd '/Utilities/Backup/' newleaf '.mat'], 'signaldata','-v7.3');
        
        handles.(matlab.lang.makeValidName(newleaf)) = uitreenode('v0', newleaf,  newleaf,  iconPath, true);
        handles.(matlab.lang.makeValidName(newleaf)).UserData = signaldata;
        handles.measurements.add(handles.(matlab.lang.makeValidName(newleaf)));
        handles.mytree.reloadNode(handles.measurements);
        handles.mytree.expand(handles.measurements);
        handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(newleaf)));
        set([handles.clrall_btn,handles.export_btn],'Enable','on')
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Processed "' char(selectedNodes(1).getName) '" to generate an impulse response of ' num2str(IRlength) ' points\n']);
        switch calcmethod
            case 1
                fprintf(handles.fid,'calcmethod = 1; %% Convolve audio with audio2\n');
            case 2
                fprintf(handles.fid,'calcmethod = 2; %% Cross-correlate audio with audio2\n');
            case 3
                fprintf(handles.fid,'calcmethod = 3; %% Circular convolution of audio with audio2 (based on the length of audio2)\n');
            case 4
                fprintf(handles.fid,'calcmethod = 4; %% Circular cross-correlation of audio with audio2 (based on the length of audio2)\n');
            case 5
                fprintf(handles.fid,'calcmethod = 5; %% Transfer function from audio2 to audio (-200 dB threshold)\n');
            case 6
                fprintf(handles.fid,'calcmethod = 6; %% Transfer function from audio2 to audio (-90 dB threshold)\n');
            case 7
                fprintf(handles.fid,'calcmethod = 7; %% Transfer function from audio2 to audio (-80 dB threshold)\n');
            case 8
                fprintf(handles.fid,'calcmethod = 8; %% Transfer function from audio2 to audio (-70 dB threshold)\n');
            case 9
                fprintf(handles.fid,'calcmethod = 9; %% Transfer function from reversed audio2 to audio (-200 dB threshold)\n');
            case 10
                fprintf(handles.fid,'calcmethod = 10; %% Transfer function from reversed audio2 to audio (-90 dB threshold)\n');
            case 11
                fprintf(handles.fid,'calcmethod = 11; %% Transfer function from reversed audio2 to audio (-80 dB threshold)\n');
            case 12
                fprintf(handles.fid,'calcmethod = 12; %% Transfer function from reversed audio2 to audio (-70 dB threshold)\n');
            case 13
                fprintf(handles.fid,'calcmethod = 13; %% Transfer function from audio to audio2 (-200 dB threshold)\n');
            case 14
                fprintf(handles.fid,'calcmethod = 14; %% Transfer function from audio to audio2 (-90 dB threshold)\n');
            case 15
                fprintf(handles.fid,'calcmethod = 15; %% Transfer function from audio to audio2 (-80 dB threshold)\n');
            case 16
                fprintf(handles.fid,'calcmethod = 16; %% Transfer function from audio to audio2 (-70 dB threshold)\n');
            case 17
                fprintf(handles.fid,'calcmethod = 17; %% Transfer function from reversed audio to audio2 (-200 dB threshold)\n');
            case 18
                fprintf(handles.fid,'calcmethod = 18; %% Transfer function from reversed audio to audio2 (-90 dB threshold)\n');
            case 19
                fprintf(handles.fid,'calcmethod = 19; %% Transfer function from reversed audio to audio2 (-80 dB threshold)\n');
            case 20
                fprintf(handles.fid,'calcmethod = 20; %% Transfer function from reversed audio to audio2 (-70 dB threshold)\n');
            case 21
                fprintf(handles.fid,'calcmethod = 21; %% Time domain deconvolution of audio2 from audio\n');
            case 22
                fprintf(handles.fid,'calcmethod = 22; %% Time domain deconvolution of audio from audio2\n');
            case 23
                fprintf(handles.fid,'calcmethod = 23; %% Time domain deconvolution of time-reversed audio2 from audio\n');
            case 24
                fprintf(handles.fid,'calcmethod = 24; %% Time domain deconvolution of time-reversed audio from audio2\n');
            case 25
                fprintf(handles.fid,'calcmethod = 25; %% Time domain convolution of audio with audio2\n');
            case 26
                fprintf(handles.fid,'calcmethod = 26; %% Time domain convolution of audio with time-reversed audio2\n');
        end
        switch method
            case 1
                fprintf(handles.fid,'method = 1; %% Synchronous average of cycles (excluding silent cycle)\n');
            case 2
                fprintf(handles.fid,'method = 2; %% Stack multicycle IR measurements in dimension 4\n');
            case 3
                fprintf(handles.fid,'method = 3; %% Reshape higher dimensions (>3) to channels\n');
            case 4
                fprintf(handles.fid,'method = 4; %% Simply convolve (without averaging, stacking or selecting)\n');
            case 5
                fprintf(handles.fid,'method = 5; %% Select the cleanest cycle\n');
            case 6
                fprintf(handles.fid,'method = 6; %% Select the cleanest IR (multichannel)\n');
            case 7
                fprintf(handles.fid,'method = 7; %% Select the cleanest single IR (best channel)\n');
            case 8
                fprintf(handles.fid,'method = 8; %% Select the silent cycle or the IR with the lowest SNR (multichannel)\n');
            case 9
                fprintf(handles.fid,'method = 9; %% Exclude the IR with the lowest SNR (multichannel)\n');
            case 10
                fprintf(handles.fid,'method = 10; %% Stack of IRs cumulatively averaged from best to worst SNR, with silent cycle (if available)\n');
            case 11
                fprintf(handles.fid,'method = 11; %% Stack of IRs cumulatively averaged from best to worst SNR, with silent cycle (if available), and visualisation the stack\n');
        end
        
        fprintf(handles.fid,['X = convolveaudiowithaudio2(X,','method',',',num2str(scalingmethod),',calcmethod);\n']);
        
        fprintf(handles.fid,['X.audio = X.audio(',num2str(trimsamp_low),':',num2str(trimsamp_high),...
            ',[',num2str(chanind),...
            '],[',num2str(bandind),...
            '],[',num2str(cycind),...
            '],[',num2str(outchanind),...
            '],[',num2str(dim6ind),...
            ']);\n']);
        
        fprintf(handles.fid,'\n');
        % Log verbose metadata (not necessary here)
        % logaudioleaffields(signaldata);
    end
end
handles.alternate = 0;
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');
java.lang.Runtime.getRuntime.gc % Java garbage collection
guidata(hObject, handles);




% *************************************************************************
% CAL (CALIBRATION) BUTTON
% *************************************************************************


% --- Executes on button press in cal_btn.
function cal_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to cal_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMain = getappdata(0,'hMain');
signaldata = getappdata(hMain,'testsignal');
selectedNodes = handles.mytree.getSelectedNodes;
%selectedNodes = selectedNodes(1);

method = menu('Calibration',...
    'Choose from AARAE',...
    'Locate file on disc',...
    'Input value',...
    'Specify Leq',...
    'Specify weighted Leq',...
    'Cancel');
cal_level = [];
switch method
    case 1
        root = handles.root; % Get selected leaf
        root = root(1);
        first = root.getFirstChild;
        nbranches = root.getChildCount;
        branches = cell(nbranches,1);
        branches{1,1} = char(first.getValue);
        nleaves = 0;
        nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{1,1}))(1).getChildCount;
        next = first.getNextSibling;
        for n = 2:nbranches
            branches{n,1} = char(next.getValue);
            nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{n,1}))(1).getChildCount;
            next = next.getNextSibling;
        end
        leaves = cell(nleaves,1);
        i = 0;
        for n = 1:size(branches,1)
            currentbranch = handles.(matlab.lang.makeValidName(branches{n,1}));
            if currentbranch.getChildCount ~= 0
                i = i + 1;
                first = currentbranch.getFirstChild;
                %leafnames(i,:) = first.getName;
                leaves{i,:} = char(first.getValue);
                next = first.getNextSibling;
                if ~isempty(next)
                    for m = 1:currentbranch.getChildCount-1
                        i = i + 1;
                        %leafnames(i,:) = next.getName;
                        leaves{i,:} = char(next.getValue);
                        next = next.getNextSibling;
                    end
                end
            end
        end
        [s,ok] = listdlg('PromptString','Select a file:',...
            'SelectionMode','single',...
            'ListString',leaves);
        if ok == 1
            caldata = handles.(matlab.lang.makeValidName(leaves{s,1})).handle.UserData;
            if ~isfield(caldata,'audio')
                warndlg('Incompatible calibration file','Warning!');
            else
                cal_level = 10 .* log10(mean(caldata.audio.^2,1));
                %calibration = 1./(10.^((cal_level)./20));
                %signaldata.audio = signaldata.audio.*calibration;
                if (size(signaldata.audio,2) == size(cal_level,2) || size(cal_level,2) == 1) && ismatrix(caldata.audio)
                    cal_offset = inputdlg({'Calibration tone RMS level',...
                        'Underlying units',...
                        'Reference value for decibels',...
                        'Units type (1 for amplitude, 2 for power)'},...
                        'Reference audio values',...
                        [1 50],{'94','Pa','2e-5','1'});
                    if isnan(str2double(char(cal_offset(1))))
                        return
                    else
                        cal_level = str2double(char(cal_offset(1))) - cal_level;
                        units = char(cal_offset(2));
                        units_ref = str2double(char(cal_offset(3)));
                        units_type = str2double(char(cal_offset(4)));
                    end
                    if size(cal_level,2) == 1, cal_level = repmat(cal_level,1,size(signaldata.audio,2)); end
                else
                    warndlg('Incompatible calibration file','Warning!');
                end
            end
        else
            warndlg('No signal loaded!','Whoops...!');
        end
    case 2
        [filename,handles.defaultaudiopath] = uigetfile(...
            {'*.wav;*.mat;.WAV;.MAT','Calibration file (*.wav,*.mat)'},...
            'Select audio file',handles.defaultaudiopath);
        if ~ischar(filename)
            return
        else
            [~,~,ext] = fileparts(filename);
        end
        if filename ~= 0
            % Check type of file. First 'if' is for .mat, second is for .wav
            if strcmp(ext,'.mat') || strcmp(ext,'.MAT')
                file = importdata(fullfile(handles.defaultaudiopath,filename));
                if isstruct(file)
                    caltone = file.audio;
                else
                    caltone = file;
                end
            elseif strcmp(ext,'.wav') || strcmp(ext,'.WAV')
                caltone = audioread(fullfile(handles.defaultaudiopath,filename));
            else
                caltone = [];
            end
            if size(caltone,1) < size(caltone,2), caltone = caltone'; end
            cal_level = 10 * log10(mean(caltone.^2,1));
            if (size(signaldata.audio,2) == size(cal_level,2) || size(cal_level,2) == 1) && ismatrix(caltone)
                cal_offset = inputdlg({'Calibration tone RMS level',...
                    'Underlying units',...
                    'Reference value for decibels',...
                    'Units type (1 for amplitude, 2 for power)'},...
                    'Reference audio values',...
                    [1 50],{'94','Pa','2e-5','1'});
                if isnan(str2double(char(cal_offset(1))))
                    return
                else
                    cal_level = str2double(char(cal_offset(1))) - cal_level;
                    units = char(cal_offset(2));
                    units_ref = str2double(char(cal_offset(3)));
                    units_type = str2double(char(cal_offset(4)));
                end
                if size(cal_level,2) == 1, cal_level = repmat(cal_level,1,size(signaldata.audio,2)); end
            else
                warndlg('Incompatible calibration file!','AARAE info');
            end
        end
    case 3
        chans = size(signaldata.audio,2);
        if isfield(signaldata,'cal')
            def = cellstr(num2str(signaldata.cal'));
        else
            def = cellstr(num2str(zeros(chans,1)));
        end
        if chans > 1
            cal_level = inputdlg(cellstr([repmat('channel ',chans,1) num2str((1:chans)')]),...
                'Calibration value (dB)',[1 60],def);
            cal_level = str2num(char(cal_level))'; %#ok to prevent from spaces introduced in the input boxes
            if size(cal_level,1) > size(cal_level,2), cal_level = cal_level'; end
            if isempty(cal_level) || chans ~= size(cal_level,2)
                warndlg('Calibration values mismatch!','AARAE info');
                return
            end
            cal_units = inputdlg({'Underlying units',...
                'Reference value for decibels',...
                'Units type (1 for amplitude, 2 for power)'},...
                'Reference audio values',...
                [1 50],{'Pa','2e-5','1'});
            if isnan(str2double(char(cal_units(2))))
                return
            else
                units = char(cal_units(1));
                units_ref = str2double(char(cal_units(2)));
                units_type = str2double(char(cal_units(3)));
            end
        else
            cal_offset = inputdlg({'Calibration value (dB)',...
                'Underlying units',...
                'Reference value for decibels',...
                'Units type (1 for amplitude, 2 for power)'},...
                'Reference audio values',...
                [1 50],{'0','Pa','2e-5','1'});
            if isnan(str2double(char(cal_offset(1))))
                return
            else
                cal_level = str2double(char(cal_offset(1)));
                units = char(cal_offset(2));
                units_ref = str2double(char(cal_offset(3)));
                units_type = str2double(char(cal_offset(4)));
            end
        end
    case 4
        caldata = selectedNodes(1).handle.UserData;
        cal_level = 10 .* log10(mean(caldata.audio.^2,1));
        cal_level = repmat(20*log10(mean(10.^(cal_level./20),2)),1,size(caldata.audio,2));
        if size(caldata.audio,2) > 1
            cal_offset = inputdlg('Signal RMS level',...
                'Calibration value',[1 50],cellstr(num2str(70+zeros(size(cal_level)))));
            if isempty(cal_offset)
                return;
            else
                cal_offset = str2num(char(cal_offset)); %#ok : to allow spaces between calibration values
            end
            if (isequal(size(cal_offset),size(cal_level)) || size(cal_offset,2) == 1) && ismatrix(caldata.audio)
                cal_level = cal_offset - cal_level;
            else
                warndlg('Calibration values mismatch!','AARAE info');
                return
            end
            cal_units = inputdlg({'Underlying units',...
                'Reference value for decibels',...
                'Units type (1 for amplitude, 2 for power)'},...
                'Reference audio values',...
                [1 50],{'Pa','2e-5','1'});
            if isnan(str2double(char(cal_units(2))))
                return
            else
                units = char(cal_units(1));
                units_ref = str2double(char(cal_units(2)));
                units_type = str2double(char(cal_units(3)));
            end
        else
            cal_offset = inputdlg({'Signal RMS level (dB)',...
                'Underlying units',...
                'Reference value for decibels',...
                'Units type (1 for amplitude, 2 for power)'},...
                'Reference audio values',...
                [1 50],{'70','Pa','2e-5','1'});
            if isnan(str2double(char(cal_offset(1))))
                return
            else
                cal_level = str2double(char(cal_offset(1))) - cal_level;
                units = char(cal_offset(2));
                units_ref = str2double(char(cal_offset(3)));
                units_type = str2double(char(cal_offset(4)));
            end
        end
    case 5
        caldata = selectedNodes(1).handle.UserData;
        weights = what([cd '/Processors/Filters']);
        if ~isempty(weights.m)
            [selection,ok] = listdlg('ListString',cellstr(weights.m),'SelectionMode','single');
        else
            warndlg('No weighting filters found!','AARAE info')
            return
        end
        if ok ~= 1, return; end
        [~,funname] = fileparts(weights.m{selection,1});
        caldata = feval(funname,caldata);
        cal_level = 10 .* log10(mean(caldata.audio.^2,1));
        cal_level = repmat(20*log10(mean(10.^(cal_level./20),2)),1,size(caldata.audio,2));
        if size(caldata.audio,2) > 1
            cal_offset = inputdlg('Signal RMS level',...
                'Calibration value',[1 50],cellstr(num2str(70+zeros(size(cal_level)))));
            if isempty(cal_offset)
                return;
            else
                cal_offset = str2num(char(cal_offset)); %#ok : to allow spaces between calibration values
            end
            if (isequal(size(cal_offset),size(cal_level)) || size(cal_offset,2) == 1) && ismatrix(caldata.audio)
                cal_level = cal_offset - cal_level;
            else
                warndlg('Calibration values mismatch!','AARAE info');
                return
            end
            cal_units = inputdlg({'Underlying units',...
                'Reference value for decibels',...
                'Units type (1 for amplitude, 2 for power)'},...
                'Reference audio values',...
                [1 50],{'Pa','2e-5','1'});
            if isnan(str2double(char(cal_units(2))))
                return
            else
                units = char(cal_units(1));
                units_ref = str2double(char(cal_units(2)));
                units_type = str2double(char(cal_units(3)));
            end
        else
            cal_offset = inputdlg({'Signal RMS level (dB)',...
                'Underlying units',...
                'Reference value for decibels',...
                'Units type (1 for amplitude, 2 for power)'},...
                'Reference audio values',...
                [1 50],{'70','Pa','2e-5','1'});
            if isnan(str2double(char(cal_offset(1))))
                return
            else
                cal_level = str2double(char(cal_offset(1))) - cal_level;
                units = char(cal_offset(2));
                units_ref = str2double(char(cal_offset(3)));
                units_type = str2double(char(cal_offset(4)));
            end
        end
    case 6
        return
end
if ~isempty(cal_level)
    for i = 1:length(selectedNodes)
        signaldata = selectedNodes(i).handle.UserData;
        callevel = cal_level;
        if size(signaldata.audio,2) < length(cal_level), callevel = cal_level(1:size(signaldata.audio,2)); end
        if size(signaldata.audio,2) > length(cal_level), callevel = [cal_level NaN(1,size(signaldata.audio,2)-length(cal_level))]; end
        signaldata.cal = callevel;
        if exist('units','var')
            signaldata.properties.units = units;
        end
        if exist('units_ref','var')
            signaldata.properties.units_ref = units_ref;
        end
        if exist('units_type','var')
            signaldata.properties.units_type = units_type;
        end
        signaldata = checkcal(signaldata);
        signaldata = addhistory(signaldata,'Calibrated');
        % Save as you go
        delete([cd '/Utilities/Backup/' selectedNodes(i).getName.char '.mat'])
        save([cd '/Utilities/Backup/' selectedNodes(i).getName.char '.mat'], 'signaldata','-v7.3');
        try
            %             if ~isfield(signaldata,'name')
            %                 disp('name field does not exist - see warning dialog')
            %                 fprintf(handles.fid,'name field does not exist - see below');
            %             end
            % handles.(signaldata.name).UserData = signaldata; % apply cal field directly to the audio loaded to the tree
            % selectedNodes(i).handle.UserData = signaldata; % produces error in MATLAB 2015b
            handles.(matlab.lang.makeValidName(char(selectedNodes(i).getName))).UserData = signaldata;
        catch
            warndlg('Sorry - an error occured in writing cal to the tree. Please let Densil know about this. As an alternative you can calibrate using cal_aarae in Processors-Basic (or in the edit window).','Bug!');
            fprintf(handles.fid,'Sorry - an error occured in writing cal to the tree. Please let Densil know about this. As an alternative you can calibrate using cal_aarae in Processors-Basic (or in the edit window).');
        end
        selectedParent = selectedNodes(i).getParent;
        handles.mytree.reloadNode(selectedParent);
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Calibrated "' char(selectedNodes(i).getName) '": adjusted to ' num2str(cal_level) 'dB \n\n']);
    end
    handles.mytree.setSelectedNodes(selectedNodes)
end
guidata(hObject,handles)





% *************************************************************************
% COMPARE BUTTON
% *************************************************************************


% --- Executes on button press in compare_btn.
function compare_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to compare_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedNodes = handles.mytree.getSelectedNodes;
if handles.compareaudio == 1
    % ******* Determine the size of the audio to be compared
    numberofnodes = length(selectedNodes);
    if numberofnodes > handles.Settings.maxlines
        numberofnodes = handles.Settings.maxlines;
    end
    
    [len,chans,bands,cycles,outchans,dim6] = deal(ones(numberofnodes,1));
    %         bandIDs = cell(numberofnodes,1);
    for i = 1:numberofnodes
        if handles.alternate==1 && isfield(selectedNodes(i).handle.UserData,'audio2')
            [len(i), chans(i), bands(i),cycles(i),outchans(i),dim6(i)] = ...
                size(selectedNodes(i).handle.UserData.audio2);
        else
            [len(i), chans(i), bands(i),cycles(i),outchans(i),dim6(i)] = ...
                size(selectedNodes(i).handle.UserData.audio);
            %                 if isfield(selectedNodes(i).handle.UserData,'bandID'
            %                     % need to try to match bandIDs
            %                     bandIDs{i} = selectedNodes(i).handle.UserData.bandID(:)';
            %                 end
        end
    end
    
    % get plottype from the 'time' chart (the upper chart) in the AARAE GUI
    axes = 'time';
    plottype = get(handles.(matlab.lang.makeValidName([axes '_popup'])),'Value');
    
    parameterstring10 = '';
    switch plottype
        case 1
            figurename = 'Real amplitude';
            defaultsmoothing = '0';
            plotcategory = 'Timeaxis';
        case 2
            figurename = 'Squared amplitude';
            parameterstring10 = 'Smoothing filter length [samples]';
            defaultsmoothing = '0';
            plotcategory = 'Timeaxis';
        case 3
            figurename = 'Level [dB]';
            parameterstring10 = 'Smoothing filter length [samples]';
            defaultsmoothing = '0';
            plotcategory = 'Timeaxis';
        case 4
            figurename = 'Envelope';
            parameterstring10 = 'Smoothing filter length [samples]';
            defaultsmoothing = '0';
            plotcategory = 'Timeaxis';
        case 5
            figurename = 'Instantaneous frequency [Hz]';
            parameterstring10 = 'Smoothing filter length [samples]';
            defaultsmoothing = '5';
            plotcategory = 'Timeaxis';
        case 6
            figurename = 'Absolute amplitude';
            parameterstring10 = 'Smoothing filter length [samples]';
            defaultsmoothing = '0';
            plotcategory = 'Timeaxis';
        case 7
            figurename = 'Imaginary amplitude';
            defaultsmoothing = '0';
            plotcategory = 'Timeaxis';
        case 8
            figurename = 'Level spectrum [dB]';
            parameterstring10 = 'Fractional octave band smoothing (e.g. 1 for octave, 3 for 1/3 octave, etc)';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 9
            figurename = 'Squared spectrum';
            parameterstring10 = 'Fractional octave band smoothing (e.g. 1 for octave, 3 for 1/3 octave, etc)';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 10
            figurename = 'Absolute spectrum';
            parameterstring10 = 'Fractional octave band smoothing (e.g. 1 for octave, 3 for 1/3 octave, etc)';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 11
            figurename = 'Real spectrum';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 12
            figurename = 'Imaginary spectrum';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 13
            figurename = 'Phase spectrum [radians]';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 14
            figurename = 'Unwrapped phase spectrum [radians]';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 15
            figurename = 'Phase spectrum [degrees]';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 16
            figurename = 'Unwrapped phase spectrum [rad/2pi]';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        case 17
            figurename = 'Group delay [ms]';
            defaultsmoothing = '0';
            plotcategory = 'Freqaxis';
        otherwise
            figurename = '';
            plotcategory = '';
    end
    
    
    dimsize = [numberofnodes;max(chans);max(bands);max(cycles);max(outchans);max(dim6)];
    
    if handles.Settings.calibrationtoggle == 1 && handles.alternate~=1
        default11 = '1';
    else
        default11 = '0';
    end
    
    if max(dimsize) > 1
        % ******* Make a dialog box
        parameterstring1 = 'Generate subplots for Nothing [0]';
        parameterstring2 = 'Distinct line HSV hues for Nothing [0]';
        parameterstring3 = 'Distinct line HSV saturations for Nothing [0]';
        parameterstring4 = 'Distinct line HSV values for Nothing [0]';
        if dimsize(1) > 1
            parameterstring1 = [parameterstring1 ', Audio selection [1]'];
            parameterstring2 = [parameterstring2 ', Audio selection [1]'];
            parameterstring3 = [parameterstring3 ', Audio selection [1]'];
            parameterstring4 = [parameterstring4 ', Audio selection [1]'];
        end
        if dimsize(2) > 1
            parameterstring1 = [parameterstring1 ', Channels [2]'];
            parameterstring2 = [parameterstring2 ', Channels [2]'];
            parameterstring3 = [parameterstring3 ', Channels [2]'];
            parameterstring4 = [parameterstring4 ', Channels [2]'];
        end
        if dimsize(3) > 1
            parameterstring1 = [parameterstring1 ', Bands [3]'];
            parameterstring2 = [parameterstring2 ', Bands [3]'];
            parameterstring3 = [parameterstring3 ', Bands [3]'];
            parameterstring4 = [parameterstring4 ', Bands [3]'];
        end
        if dimsize(4) > 1
            parameterstring1 = [parameterstring1 ', Cycles [4]'];
            parameterstring2 = [parameterstring2 ', Cycles [4]'];
            parameterstring3 = [parameterstring3 ', Cycles [4]'];
            parameterstring4 = [parameterstring4 ', Cycles [4]'];
        end
        if dimsize(5) > 1
            parameterstring1 = [parameterstring1 ', Asynchronous output channels [5]'];
            parameterstring2 = [parameterstring2 ', Asynchronous output channels [5]'];
            parameterstring3 = [parameterstring3 ', Asynchronous output channels [5]'];
            parameterstring4 = [parameterstring4 ', Asynchronous output channels [5]'];
        end
        if dimsize(6) > 1
            parameterstring1 = [parameterstring1 ', Dimension 6 [6]'];
            parameterstring2 = [parameterstring2 ', Dimension 6 [6]'];
            parameterstring3 = [parameterstring3 ', Dimension 6 [6]'];
            parameterstring4 = [parameterstring4 ', Dimension 6 [6]'];
        end
        
        parameterstring5 = ['Channel indices (up to ', num2str(max(chans)),')'];
        parameterstring6 = ['Band indices (up to ', num2str(max(bands)),')'];
        parameterstring7 = ['Cycle indices (up to ', num2str(max(cycles)),')'];
        parameterstring8 = ['Asynchonous output indices (up to ', num2str(max(outchans)),')'];
        parameterstring9 = ['Dimension 6 indices (up to ', num2str(max(dim6)),')'];
        
        
        nonsignletondims = find(dimsize>1,4,'first');
        if isempty(nonsignletondims)
            [default1, default2, default3, default4] = deal('0');
        elseif length(nonsignletondims) == 1
            [default1, default2] = deal(num2str(nonsignletondims));
            [default3, default4] = deal('0');
        elseif length(nonsignletondims) == 2
            default1 = num2str(nonsignletondims(1));
            default2 = num2str(nonsignletondims(2));
            [default3, default4] = deal('0');
        elseif length(nonsignletondims) == 3
            default1 = num2str(nonsignletondims(1));
            default2 = num2str(nonsignletondims(2));
            [default3, default4] = deal(num2str(nonsignletondims(3)));
        else
            default1 = num2str(nonsignletondims(1));
            default2 = num2str(nonsignletondims(2));
            default3 = num2str(nonsignletondims(3));
            default4 = num2str(nonsignletondims(4));
        end
        
        
        defaultchans = ['1:' num2str(dimsize(2))];
        defaultbands = ['1:' num2str(dimsize(3))];
        defaultcycles = ['1:' num2str(dimsize(4))];
        defaultoutchans = ['1:' num2str(dimsize(5))];
        defaultdim6 = ['1:' num2str(dimsize(6))];
        
        
        
        % Dialog box for selection
        param = inputdlg({parameterstring1;...
            parameterstring2;...
            parameterstring3;...
            parameterstring4;...
            parameterstring5;...
            parameterstring6;...
            parameterstring7;...
            parameterstring8;...
            parameterstring9;...
            parameterstring10;...
            'Raw values [0]; Use calibration if available [1]; or Normalise [2]; or Normalize to rms [3] (if relevant to plot type)';...
            'Change the plot type [0 | 1] (some extra plot types are also available)'},...% inputdlg window.
            'Data Mapping & Selection',...
            [1 90],...
            {default1;default2;default3;default4;...
            defaultchans;...
            defaultbands;...
            defaultcycles;...
            defaultoutchans;...
            defaultdim6;...
            defaultsmoothing;...
            default11;'0'});
        
        
        
        if length(param) < 12, param = []; end
        if ~isempty(param)
            subplotdim = str2num(char(param(1)));
            Hdim = str2num(char(param(2)));
            Sdim = str2num(char(param(3)));
            Vdim = str2num(char(param(4)));
            chanplot = str2num(char(param(5)));
            bandplot = str2num(char(param(6)));
            cycleplot = str2num(char(param(7)));
            outchanplot = str2num(char(param(8)));
            dim6plot = str2num(char(param(9)));
            smoothingmethod = str2num(char(param(10)));
            cal_or_norm = str2num(char(param(11)));
            changeplottype = str2num(char(param(12)));
        else
            % get out of here if the user presses 'cancel'
            return
        end
    else
        [subplotdim, Hdim, Sdim, Vdim, chanplot, bandplot, cycleplot,...
            outchanplot, dim6plot] = deal(1);
        
        % Dialog box for selection
        param = inputdlg({parameterstring10;...
            'Raw values [0]; Use calibration if available [1]; or Normalise [2]; or Normalize to rms [3] (if relevant to plot type)';...
            'Change the plot type [0 | 1]'},...% inputdlg window.
            'Smoothing & plot type',...
            [1 90],...
            {defaultsmoothing;default11;'0'});
        if ~isempty(param)
            smoothingmethod = str2num(char(param(1)));
            cal_or_norm = str2num(char(param(2)));
            changeplottype = str2num(char(param(3)));
        else
            % cancel
            return
        end
    end
    
    if changeplottype == 1
        str = {'1. Time - Real amplitude';...
            '2. Time  - Squared amplitude';...
            '3. Time - Level';...
            '4. Time - Hilbert envelope';...
            '5. Time - Instantaneous frequency';...
            '6. Time - Absolute amplitude';...
            '7. Time - Imaginary amplitude';...
            '8. Frequency - Level';...
            '9. Frequency - Squared magnitude';...
            '10. Frequency - Absolute magnitude';...
            '11. Frequency - Real';...
            '12. Frequency - Imaginary';...
            '13. Frequency - Phase [radians]';...
            '14. Frequency - Unwrapped phase [radians]';...
            '15. Frequency - Phase [degrees]';...
            '16. Frequency - Unwrapped phase [rad/2pi]';...
            '17. Frequency - Group delay';...
            '18. Cumulative time-distribution of real amplitude';...
            '19. Cumulative time-distribution of level';...
            '20. Cumulative time-distribution of Hilbert envelope';...
            '21. Cumulative sum of energy over time [dB]';...
            '22. Scatter plots of time and frequency power centroid and Leq';...
            '23. Scatter plots (A-weighted) of time and frequency power centroid and Leq'};
        
        
        
        % To Do: Put other special plots on this list. possibilities
        % include:
        % * Frequency - Octave band level
        % * Frequency - 1/3-octave band level
        % * Frequency - phase delay
        % * Quefrency - complex / real cepstrum (might not be useful)
        % * Spectrogram (requires quite different code, so might be
        %      inappropriate)
        % * Complex spectrum (using 3d plot?) - might be too slow/big
        % * analytic signal (same issues as above)
        % * Step response (maybe trivial in LTI)
        % * Simple numerical comparisons (e.g. SPL - like levelstats), eg
        %       cumulative level dist, or barplots
        %
        % * spectrum stats expressed graphically (centroid, stdev, etc)
        % * Spatial analysis of multichannel signals (e.g. 2chan level diff
        %       using a temporal integrator, HOA analysis).
        % * Autocorrelation analysis
        % * Special comparison analyses using a reference: level diff from
        %       ref in time & freq, xcorr, tf, etc. However we need a way
        %       of identifying the reference audio
        %
        
        
        % % PHASE DELAY
        %[signaldata.audio,f] = ...
        %phasedelay(signaldata.audio,1,length(signaldata.audio),signaldata.fs);
        % (just an idea)
        %end
        
        
        
        [plottypeselection,ok] = listdlg('PromptString','Select the plot type',...
            'SelectionMode','single',...
            'ListString',str,...
            'ListSize', [400,400]);
        if ok
            plottype = plottypeselection;
            switch plottype
                case 1
                    figurename = 'Real amplitude';
                    plotcategory = 'Timeaxis';
                case 2
                    figurename = 'Squared amplitude';
                    plotcategory = 'Timeaxis';
                case 3
                    figurename = 'Level [dB]';
                    plotcategory = 'Timeaxis';
                case 4
                    figurename = 'Envelope';
                    plotcategory = 'Timeaxis';
                case 5
                    figurename = 'Instantaneous frequency [Hz]';
                    plotcategory = 'Timeaxis';
                case 6
                    figurename = 'Absolute amplitude';
                    plotcategory = 'Timeaxis';
                case 7
                    figurename = 'Imaginary amplitude';
                    plotcategory = 'Timeaxis';
                case 8
                    figurename = 'Level spectrum [dB]';
                    plotcategory = 'Freqaxis';
                case 9
                    figurename = 'Squared spectrum';
                    plotcategory = 'Freqaxis';
                case 10
                    figurename = 'Absolute spectrum';
                    plotcategory = 'Freqaxis';
                case 11
                    figurename = 'Real spectrum';
                    plotcategory = 'Freqaxis';
                case 12
                    figurename = 'Imaginary spectrum';
                    plotcategory = 'Freqaxis';
                case 13
                    figurename = 'Phase spectrum [radians]';
                    plotcategory = 'Freqaxis';
                case 14
                    figurename = 'Unwrapped phase spectrum [radians]';
                    plotcategory = 'Freqaxis';
                case 15
                    figurename = 'Phase spectrum [degrees]';
                    plotcategory = 'Freqaxis';
                case 16
                    figurename = 'Unwrapped phase spectrum [rad/2pi]';
                    plotcategory = 'Freqaxis';
                case 17
                    figurename = 'Group delay [ms]';
                    plotcategory = 'Freqaxis';
                case 18
                    figurename = 'Cumulative time distribution (real amplitude)';
                    plotcategory = 'Timeaxis';
                case 19
                    figurename = 'Cumulative time distribution [dB]';
                    plotcategory = 'Timeaxis';
                case 20
                    figurename = 'Cumulative time distribution (Hilbert envelope)';
                    plotcategory = 'Timeaxis';
                case 21
                    figurename = 'Cumulative sum of energy over time [dB]';
                    plotcategory = 'Timeaxis';
                case 22
                    figurename = 'Time and freq power centroids and levels';
                    plotcategory = 'Scatter3axis';
                case 23
                    figurename = 'A-weighted time and freq power centroids and levels';
                    plotcategory = 'Scatter3axis';
                otherwise
                    figurename = '';
                    plotcategory = '';
            end
        end
    end
    
    % size of selected audio
    [chansselect,bandsselect,cyclesselect,outchansselect,dim6select] =...
        deal(ones(numberofnodes,1));
    for i = 1:numberofnodes
        if handles.alternate==1 && isfield(selectedNodes(i).handle.UserData,'audio2')
            [~, chansselect(i), bandsselect(i),...
                cyclesselect(i),outchansselect(i),dim6select(i)] = ...
                size(selectedNodes(i).handle.UserData.audio2(:,...
                chanplot(chanplot<=chans(i) & chanplot<=handles.Settings.maxlines),...
                bandplot(bandplot<=bands(i) & bandplot<=handles.Settings.maxlines),...
                cycleplot(cycleplot<=cycles(i) & cycleplot<=handles.Settings.maxlines),...
                outchanplot(outchanplot<=outchans(i) & outchanplot<=handles.Settings.maxlines),...
                dim6plot(dim6plot<=dim6(i) & dim6plot<=handles.Settings.maxlines)));
        else
            [~, chansselect(i), bandsselect(i),...
                cyclesselect(i),outchansselect(i),dim6select(i)] = ...
                size(selectedNodes(i).handle.UserData.audio(:,...
                chanplot(chanplot<=chans(i) & chanplot<=handles.Settings.maxlines),...
                bandplot(bandplot<=bands(i) & bandplot<=handles.Settings.maxlines),...
                cycleplot(cycleplot<=cycles(i) & cycleplot<=handles.Settings.maxlines),...
                outchanplot(outchanplot<=outchans(i) & outchanplot<=handles.Settings.maxlines),...
                dim6plot(dim6plot<=dim6(i) & dim6plot<=handles.Settings.maxlines)));
        end
    end
    
    dimsizeselect = [numberofnodes;max(chansselect);...
        max(bandsselect);max(cyclesselect);max(outchansselect);...
        max(dim6select)];
    
    % number of subplots, hues, saturations and values (in HSV color-space)
    if subplotdim >= 1 && subplotdim <= 6
        numberofsubplots = dimsizeselect(subplotdim);
    else
        numberofsubplots = 1;
    end
    if Hdim >= 1 && Hdim <= 6
        numberofH = dimsizeselect(Hdim);
    else
        numberofH = 1;
    end
    if Sdim >= 1 && Sdim <= 6
        numberofS = dimsizeselect(Sdim);
    else
        numberofS = 1;
    end
    if Vdim >= 1 && Vdim <= 6
        numberofV = dimsizeselect(Vdim);
    else
        numberofV = 1;
    end
    
    [r, c] = subplotpositions(numberofsubplots, 0.5);
    linecolor = HSVplotcolours2(numberofH, numberofS, numberofV);
    
    
    
    
    % make the figure
    compplot = figure('Name',figurename);
    pointmark = {'o';'x';'^';'s';'p';'hexagram';'d';'v';'>';'<';'+';'*'}; % if markers are needed (note that there are only 12 here)
    for i = 1:numberofnodes
        signaldata = selectedNodes(i).handle.UserData;
        if handles.alternate==1 && isfield(selectedNodes(i).handle.UserData,'audio2')
            signaldata.audio = signaldata.audio2;
        end
        % default units settings, changed below if relevant
        units = '';
        units_ref = 1;
        units_type = 1;
        if ~isempty(signaldata) && isfield(signaldata,'audio')
            % apply calibration (if requested and if possible)
            if isfield(signaldata,'cal') && cal_or_norm == 1  && handles.alternate~=1
                if size(signaldata.audio,2) == length(signaldata.cal)
                    signaldata.cal(isnan(signaldata.cal)) = 0;
                    if isfield(signaldata,'properties')
                        if isfield(signaldata.properties,'units')
                            units = signaldata.properties.units;
                        else
                            units = '';
                        end
                        if isfield(signaldata.properties,'units_ref')
                            units_ref = signaldata.properties.units_ref;
                        else
                            units_ref = 1;
                        end
                        if isfield(signaldata.properties,'units_type')
                            units_type = signaldata.properties.units_type;
                        else
                            units_type = 1;
                        end
                        %                     else
                        %                         units = '';
                        %                         units_ref = 1;
                        %                         units_type = 1;
                    end
                    if units_type == 1
                        signaldata.audio = signaldata.audio * units_ref;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref/20);
                    else
                        signaldata.audio = signaldata.audio * units_ref.^0.5;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
                    end
                    signaldata.audio = signaldata.audio .* ...
                        repmat(10.^(signaldata.cal(:)'./20),[size(signaldata.audio,1),1,...
                        size(signaldata.audio,3),size(signaldata.audio,4),...
                        size(signaldata.audio,5),size(signaldata.audio,6)]);
                end
            end
            To = floor(str2double(get(handles.To_time,'String'))*signaldata.fs)+1;
            Tf = floor(str2double(get(handles.Tf_time,'String'))*signaldata.fs);
            if Tf > length(signaldata.audio), Tf = length(signaldata.audio); end
            signaldata.audio = signaldata.audio(To:Tf,...
                chanplot(chanplot<=chans(i)),...
                bandplot(bandplot<=bands(i)),...
                cycleplot(cycleplot<=cycles(i)),...
                outchanplot(outchanplot<=outchans(i)),...
                dim6plot(dim6plot<=dim6(i)));
            t = linspace(0,length(signaldata.audio),length(signaldata.audio))./signaldata.fs;
            f = signaldata.fs .* ((1:length(signaldata.audio))-1) ./ length(signaldata.audio);
            switch handles.Settings.specmagscale;
                case {'Divided by length'}
                    spectscale = 1./len(i);
                case {'x sqrt2/length'}
                    spectscale = 2.^0.5./len(i);
                case {'x 2/length'}
                    spectscale = 2./len(i);
                otherwise
                    spectscale = 1;
            end
            
            
            if plottype == 1
                % REAL AMPLITUDE AS A FUNCTION OF TIME
                signaldata.audio = real(signaldata.audio);
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 2
                % SQUARED AMPLITUDE AS A FUNCTION OF TIME
                if smoothingmethod > 1
                    signaldata.audio = filter(ones(1,smoothingmethod)/smoothingmethod,...
                        1,signaldata.audio.^2);
                else
                    signaldata.audio = signaldata.audio.^2;
                end
                if cal_or_norm == 1
                    if units_type == 2
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^2]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(mean(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 3
                % LEVEL AS A FUNCTION OF TIME
                if smoothingmethod > 1
                    signaldata.audio = 10.*log10(...
                        filter(ones(1,smoothingmethod)/smoothingmethod,...
                        1,signaldata.audio.^2));
                else
                    signaldata.audio = 10.*log10(signaldata.audio.^2);
                end
                if cal_or_norm == 1
                    if units_type == 1
                        signaldata.audio = signaldata.audio - 20*log10(units_ref);
                    else
                        signaldata.audio = signaldata.audio - 10*log10(units_ref);
                    end
                    units_string = ' [dB]';
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        - repmat(max(signaldata.audio),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        - repmat(pow2db(mean(db2pow(signaldata.audio))),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 4
                % TIME ENVELOPE (VIA HILBERT)
                for b = 1:bandsselect(i)
                    for d4 = 1:cyclesselect(i)
                        for d5 = 1:outchansselect(i)
                            for d6 = 1:dim6select(i)
                                signaldata.audio(:,:,b,d4,d5,d6) = ...
                                    abs(hilbert(real(signaldata.audio(:,:,b,d4,d5,d6))));
                            end
                        end
                    end
                end
                if smoothingmethod > 1
                    signaldata.audio = filter(ones(1,smoothingmethod)/smoothingmethod,...
                        1,signaldata.audio);
                end
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 5
                % INSTANTANEOUS FREQUENCY
                for b = 1:bandsselect(i)
                    for d4 = 1:cyclesselect(i)
                        for d5 = 1:outchansselect(i)
                            for d6 = 1:dim6select(i)
                                signaldata.audio(:,:,b,d4,d5,d6) = ...
                                    diff([angle(hilbert(real(...
                                    signaldata.audio(:,:,b,d4,d5,d6)))); ...
                                    zeros(1,chansselect(i))])*signaldata.fs/2/pi;
                            end
                        end
                    end
                end
                if smoothingmethod > 1
                    signaldata.audio = ...
                        medfilt1(signaldata.audio,...
                        smoothingmethod);
                end
            end
            
            if plottype == 6
                % ABSOLUTE AMPLITUDE AS A FUNCTION OF TIME
                if smoothingmethod > 1
                    signaldata.audio = filter(ones(1,smoothingmethod)/smoothingmethod,...
                        1,abs(signaldata.audio));
                else
                    signaldata.audio = abs(signaldata.audio);
                end
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 7
                % IMAGINARY AMPLITUDE AS A FUNCTION OF TIME
                signaldata.audio = imag(signaldata.audio);
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2 && mean(mean(mean(mean(mean(max(abs(signaldata.audio))))))) > 0
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            % DO THE FFT HERE FOR THE SPECTRUM PLOTS
            if strcmp(plotcategory,'Freqaxis')
                % try to avoid out-of-memory error by limiting the maximum
                % size of the fft).
                if numel(signaldata.audio) < 1e6
                    signaldata.audio = fft(signaldata.audio);
                else
                    for ch = 1:chansselect(i)
                        for b = 1:bandsselect(i)
                            for d4 = 1:cyclesselect(i)
                                for d5 = 1:outchansselect(i)
                                    for d6 = 1:dim6select(i)
                                        signaldata.audio(:,ch,b,d4,d5,d6) = ...
                                            fft(signaldata.audio(:,ch,b,d4,d5,d6));
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            if plottype == 8
                % LEVEL AS A FUNCTION OF FREQUENCY
                if smoothingmethod < 1
                    signaldata.audio =...
                        10*log10(abs(signaldata.audio.*spectscale).^2);
                else
                    for b = 1:bandsselect(i)
                        for d4 = 1:cyclesselect(i)
                            for d5 = 1:outchansselect(i)
                                for d6 = 1:dim6select(i)
                                    % need to vectorize octavesmoothing!
                                    signaldata.audio(:,:,b,d4,d5,d6) = octavesmoothing(abs(signaldata.audio(:,:,b,d4,d5,d6)).^2,...
                                        smoothingmethod, signaldata.fs);
                                end
                            end
                        end
                    end
                    lowlimit = 128/(len(i)/signaldata.fs); % avoid very low freq hump error
                    signaldata.audio = signaldata.audio(f>lowlimit,:,:,:,:,:);
                    f = f(f>lowlimit);
                    signaldata.audio = 10*log10(signaldata.audio);
                end
                if cal_or_norm == 1
                    units_string = ' [dB]';
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        - repmat(max(signaldata.audio),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        - repmat(pow2db(mean(db2pow(signaldata.audio))),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 9
                % SQUARED MAGNITUDE AS A FUNCTION OF FREQUENCY (POWER SPECTRUM)
                signaldata.audio = (abs(signaldata.audio).*spectscale).^2;
                if smoothingmethod >= 1
                    for b = 1:bandsselect(i)
                        for d4 = 1:cyclesselect(i)
                            for d5 = 1:outchansselect(i)
                                for d6 = 1:dim6select(i)
                                    % need to vectorize octavesmoothing!
                                    signaldata.audio(:,:,b,d4,d5,d6) = octavesmoothing(signaldata.audio(:,:,b,d4,d5,d6),...
                                        smoothingmethod, signaldata.fs);
                                end
                            end
                        end
                    end
                    lowlimit = 128/(len(i)/signaldata.fs); % avoid very low freq hump error
                    signaldata.audio = signaldata.audio(f>lowlimit,:,:,:,:,:);
                    f = f(f>lowlimit);
                end
                if cal_or_norm == 1
                    if units_type == 2
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^2]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(mean(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 10
                % MAGNITUDE AS A FUNCTION OF FREQUENCY
                signaldata.audio = abs(signaldata.audio).*spectscale;
                if smoothingmethod >= 1
                    for b = 1:bandsselect(i)
                        for d4 = 1:cyclesselect(i)
                            for d5 = 1:outchansselect(i)
                                for d6 = 1:dim6select(i)
                                    % need to vectorize octavesmoothing!
                                    signaldata.audio(:,:,b,d4,d5,d6) = octavesmoothing(signaldata.audio(:,:,b,d4,d5,d6),...
                                        smoothingmethod, signaldata.fs);
                                end
                            end
                        end
                    end
                    lowlimit = 128/(len(i)/signaldata.fs); % avoid very low freq hump error
                    signaldata.audio = signaldata.audio(f>lowlimit,:,:,:,:,:);
                    f = f(f>lowlimit);
                end
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 11
                % REAL SPECTRUM
                signaldata.audio = real(signaldata.audio).*spectscale;
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 12
                % IMAGINARY SPECTRUM
                signaldata.audio = imag(signaldata.audio).*spectscale;
                if cal_or_norm == 1
                    if units_type == 1
                        units_string = [' [' units ']'];
                    else
                        if ~isempty(units)
                            units_string = [' [(' units ')^0.5]'];
                        else
                            units_string = '';
                        end
                    end
                elseif cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            % VARIOUS PHASE SPECTRA
            if plottype == 13, signaldata.audio = angle(signaldata.audio); end
            if plottype == 14, signaldata.audio = unwrap(angle(signaldata.audio)); end
            if plottype == 15, signaldata.audio = angle(signaldata.audio) .* 180/pi; end
            if plottype == 16
                signaldata.audio = unwrap(angle(signaldata.audio)) ./(2*pi); end
            
            if plottype == 17
                % GROUP DELAY
                signaldata.audio = -diff(unwrap(angle(signaldata.audio))).*length(signaldata.audio)/(signaldata.fs*2*pi).*1000;
                f = f(1:end-1);
            end
            
            if plottype == 18
                % CUMULATIVE TIME DISTRIBUTION OF AMPLITUDE
                signaldata.audio = sort(real(signaldata.audio));
                if cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 19
                % CUMULATIVE TIME DISTRIBUTION OF LEVEL
                signaldata.audio = sort(10*log10(abs(signaldata.audio).^2));
                if cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        - repmat(max(signaldata.audio),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        - repmat(pow2db(mean(db2pow(signaldata.audio))),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 20
                % CUMULATIVE TIME DISTRIBUTION OF ENVELOPE
                for b = 1:bandsselect(i)
                    for d4 = 1:cyclesselect(i)
                        for d5 = 1:outchansselect(i)
                            for d6 = 1:dim6select(i)
                                signaldata.audio(:,:,b,d4,d5,d6) = ...
                                    hilbert(real(signaldata.audio(:,:,b,d4,d5,d6)));
                            end
                        end
                    end
                end
                signaldata.audio = sort(abs(signaldata.audio));
                if cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        ./ repmat(max(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        ./ repmat(rms(abs(signaldata.audio)),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            if plottype == 21
                % CUMULATIVE SUM OF ENERGY OVER TIME (in dB)
                signaldata.audio = 10*log10(cumsum(abs(signaldata.audio).^2)./signaldata.fs);
                if cal_or_norm == 2
                    signaldata.audio = signaldata.audio...
                        - repmat(max(signaldata.audio),[size(signaldata.audio,1),1,1,1,1,1]);
                elseif cal_or_norm == 3
                    signaldata.audio = signaldata.audio...
                        - repmat(pow2db(mean(db2pow(signaldata.audio))),[size(signaldata.audio,1),1,1,1,1,1]);
                end
            end
            
            if plottype == 23, signaldata.audio = Aweight(signaldata.audio,signaldata.fs); end
            if plottype == 22 || plottype == 23
                % Leq
                Leq = 10*log10(mean(signaldata.audio.^2));
                Lmax = 10*log10(max(signaldata.audio.^2));
                L50 = 10*log10(median(signaldata.audio.^2));
                
                % power temporal centroid (centre time)
                t=repmat(t',[1,size(signaldata.audio,2),...
                    size(signaldata.audio,3),size(signaldata.audio,4),...
                    size(signaldata.audio,5),size(signaldata.audio,6)]);
                tcentroid = sum(t.*signaldata.audio.^2) ./ ...
                    sum(signaldata.audio.^2);
                
                
                % power spectral centroid
                if numel(signaldata.audio) < 1e6
                    signaldata.audio = fft(signaldata.audio);
                else
                    for ch = 1:chansselect(i)
                        for b = 1:bandsselect(i)
                            for d4 = 1:cyclesselect(i)
                                for d5 = 1:outchansselect(i)
                                    for d6 = 1:dim6select(i)
                                        signaldata.audio(:,ch,b,d4,d5,d6) = ...
                                            fft(signaldata.audio(:,ch,b,d4,d5,d6));
                                    end
                                end
                            end
                        end
                    end
                end
                f=repmat(f(1:round(end/2))',[1,size(signaldata.audio,2),...
                    size(signaldata.audio,3),size(signaldata.audio,4),...
                    size(signaldata.audio,5),size(signaldata.audio,6)]);
                spectrum = abs(signaldata.audio(1:length(f),:,:,:,:,:)).^2;
                fcentroid = sum(f.*spectrum(1:size(f,1),:,:,:,:,:)) ./ ...
                    sum(spectrum(1:size(f,1),:,:,:,:,:));
                
                
            end
            
            
            
            
            %             if strcmp(get(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Visible'),'on')
            %                 smoothfactor = get(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Value');
            %                 if smoothfactor == 2, octsmooth = 1; end
            %                 if smoothfactor == 3, octsmooth = 3; end
            %                 if smoothfactor == 4, octsmooth = 6; end
            %                 if smoothfactor == 5, octsmooth = 12; end
            %                 if smoothfactor == 6, octsmooth = 24; end
            %                 if smoothfactor ~= 1, signaldata.audio = octavesmoothing(signaldata.audio, octsmooth, signaldata.fs); end
            %             end
            
            for ch = 1:chansselect(i)
                for b = 1:bandsselect(i)
                    for d4 = 1:cyclesselect(i)
                        for d5 = 1:outchansselect(i)
                            for d6 = 1:dim6select(i)
                                labelstring = '';
                                titlestring = '';
                                switch subplotdim
                                    case 1
                                        plotnum = i;
                                        if isfield(signaldata,'name')
                                            titlestring = strrep(signaldata.name,'_',' ');
                                        end
                                        if handles.alternate==1 && isfield(signaldata,'audio2')
                                            titlestring = [titlestring ' audio2'];
                                            if max(chansselect(i)) > 1
                                                labelstring = [labelstring ' chan ' num2str(ch)];
                                            end
                                        else
                                            if max(chansselect(i)) > 1
                                                if isfield(signaldata,'chanID')
                                                    try
                                                        labelstring = [labelstring signaldata.chanID{ch,1}];
                                                    catch
                                                    end
                                                end
                                            end
                                            if max(bandsselect(i)) > 1
                                                if isfield(signaldata,'bandID')
                                                    try
                                                        labelstring = [labelstring num2str(signaldata.bandID(b))  ' Hz'];
                                                    catch
                                                    end
                                                end
                                            end
                                        end
                                    case 2
                                        plotnum = ch;
                                        if handles.alternate==1 && isfield(signaldata,'audio2')
                                            titlestring = ['audio2 chan ' num2str(ch)];
                                        else
                                            if isfield(signaldata,'chanID')
                                                try
                                                    titlestring = signaldata.chanID{ch,1};
                                                catch
                                                end
                                            end
                                            if max(bandsselect(i)) > 1
                                                if isfield(signaldata,'bandID')
                                                    try
                                                        labelstring = [labelstring num2str(signaldata.bandID(b)) ' Hz'];
                                                    catch
                                                    end
                                                end
                                            end
                                        end
                                    case 3
                                        plotnum = b;
                                        if isfield(signaldata,'bandID')
                                            try
                                                titlestring = [num2str(signaldata.bandID(b)) ' Hz'];
                                            catch
                                            end
                                        end
                                        if max(chansselect(i)) > 1
                                            if isfield(signaldata,'chanID')
                                                try
                                                    labelstring = [labelstring signaldata.chanID{ch,1}];
                                                catch
                                                end
                                            end
                                        end
                                    case 4
                                        plotnum = d4;
                                        if isfield(signaldata,'properties')
                                            if isfield(signaldata.properties,'relgain')
                                                try
                                                    titlestring = [num2str(signaldata.properties.relgain(d4)) ' dB'];
                                                catch
                                                end
                                            end
                                        end
                                        if max(chansselect(i)) > 1
                                            if isfield(signaldata,'chanID')
                                                try
                                                    labelstring = [labelstring signaldata.chanID{ch,1}];
                                                catch
                                                end
                                            end
                                        end
                                        if max(bandsselect(i)) > 1
                                            if isfield(signaldata,'bandID')
                                                try
                                                    labelstring = [labelstring num2str(signaldata.bandID(b)) ' Hz'];
                                                catch
                                                end
                                            end
                                        end
                                    case 5
                                        plotnum = d5;
                                        if max(chansselect(i)) > 1
                                            if isfield(signaldata,'chanID')
                                                try
                                                    labelstring = [labelstring signaldata.chanID{ch,1}];
                                                catch
                                                end
                                            end
                                        end
                                        if max(bandsselect(i)) > 1
                                            if isfield(signaldata,'bandID')
                                                try
                                                    labelstring = [labelstring num2str(signaldata.bandID(b)) ' Hz'];
                                                catch
                                                end
                                            end
                                        end
                                    case 6
                                        plotnum = d6;
                                        if max(chansselect(i)) > 1
                                            if isfield(signaldata,'chanID')
                                                try
                                                    labelstring = [labelstring signaldata.chanID{ch,1}];
                                                catch
                                                end
                                            end
                                        end
                                        if max(bandsselect(i)) > 1
                                            if isfield(signaldata,'bandID')
                                                try
                                                    labelstring = [labelstring num2str(signaldata.bandID(b)) ' Hz'];
                                                catch
                                                end
                                            end
                                        end
                                    otherwise
                                        plotnum = 1;
                                        if max(chansselect(i)) > 1
                                            if isfield(signaldata,'chanID')
                                                try
                                                    labelstring = [labelstring signaldata.chanID{ch,1}];
                                                catch
                                                end
                                            end
                                        end
                                        if max(bandsselect(i)) > 1
                                            if isfield(signaldata,'bandID')
                                                try
                                                    labelstring = [labelstring num2str(signaldata.bandID(b)) ' Hz'];
                                                catch
                                                end
                                            end
                                        end
                                end
                                switch Hdim
                                    case 1
                                        Hind = i;
                                    case 2
                                        Hind = ch;
                                    case 3
                                        Hind = b;
                                    case 4
                                        Hind = d4;
                                    case 5
                                        Hind = d5;
                                    case 6
                                        Hind = d6;
                                    otherwise
                                        Hind = 1;
                                end
                                switch Sdim
                                    case 1
                                        Sind = i;
                                    case 2
                                        Sind = ch;
                                    case 3
                                        Sind = b;
                                    case 4
                                        Sind = d4;
                                    case 5
                                        Sind = d5;
                                    case 6
                                        Sind = d6;
                                    otherwise
                                        Sind = 1;
                                end
                                switch Vdim
                                    case 1
                                        Vind = i;
                                    case 2
                                        Vind = ch;
                                    case 3
                                        Vind = b;
                                    case 4
                                        Vind = d4;
                                    case 5
                                        Vind = d5;
                                    case 6
                                        Vind = d6;
                                    otherwise
                                        Vind = 1;
                                end
                                if strcmp(plotcategory,'Timeaxis')
                                    % TIME X-AXIS
                                    subplot(r,c,plotnum)
                                    plot(t,real(signaldata.audio(:,ch,b,d4,d5,d6)), ...
                                        'color',permute(linecolor(Hind,Sind,Vind,:),[1,4,2,3]),...
                                        'DisplayName',labelstring);
                                    if numberofsubplots-c < plotnum
                                        xlabel('Time [s]');
                                    end
                                    if exist('units_string','var')
                                        if rem(plotnum-1,c) == 0
                                            ylabel(units_string);
                                        end
                                    end
                                elseif strcmp(plotcategory,'Freqaxis')
                                    % FREQ X-AXIS
                                    h=subplot(r,c,plotnum);
                                    plot(f,real(signaldata.audio(:,ch,b,d4,d5,d6)), ...
                                        'color',permute(linecolor(Hind,Sind,Vind,:),[1,4,2,3]),...
                                        'DisplayName',labelstring);
                                    if numberofsubplots-c < plotnum
                                        xlabel('Frequency [Hz]');
                                    end
                                    if ischar(handles.Settings.frequencylimits)
                                        xlim([f(2) signaldata.fs/2])
                                    else
                                        xlim(handles.Settings.frequencylimits)
                                    end
                                    log_check = get(handles.(matlab.lang.makeValidName(['log' axes '_chk'])),'Value');
                                    if log_check == 1
                                        set(h,'XScale','log')
                                    else
                                        set(h,'XScale','linear','XTickLabelMode','auto')
                                    end
                                    if exist('units_string','var')
                                        if rem(plotnum-1,c) == 0
                                            ylabel(units_string);
                                        end
                                    end
                                elseif plottype == 21 || plottype == 22 % strcmp(plotcategory,'Scatter3axis')
                                    mark = pointmark{mod(plotnum-1,12)+1};
                                    colr = permute(linecolor(Hind,Sind,Vind,:),[1,4,2,3]);
                                    subplot(2,2,1);
                                    plot(tcentroid(1,ch,b,d4,d5,d6),fcentroid(1,ch,b,d4,d5,d6), ...
                                        'Marker',mark,...
                                        'MarkerSize',8,...
                                        'color',colr,...
                                        'MarkerFaceColor',colr,...
                                        'DisplayName',labelstring);
                                    xlabel('Time [s]');
                                    ylabel('Frequency [Hz]');
                                    hold on
                                    
                                    subplot(2,2,2);
                                    errorbar(tcentroid(1,ch,b,d4,d5,d6),Leq(1,ch,b,d4,d5,d6), ...
                                        Leq(1,ch,b,d4,d5,d6)-L50(1,ch,b,d4,d5,d6),...
                                        Lmax(1,ch,b,d4,d5,d6)-Leq(1,ch,b,d4,d5,d6),...
                                        'Marker',mark,...
                                        'MarkerSize',8,...
                                        'color',colr,...
                                        'MarkerFaceColor',colr,...
                                        'DisplayName',labelstring);
                                    xlabel('Time [s]');
                                    ylabel('Leq (L50 & Lmax) [dB]');
                                    
                                    hold on
                                    
                                    hsubplot3 = subplot(2,2,3);
                                    errorbar(fcentroid(1,ch,b,d4,d5,d6),Leq(1,ch,b,d4,d5,d6), ...
                                        Leq(1,ch,b,d4,d5,d6)-L50(1,ch,b,d4,d5,d6),...
                                        Lmax(1,ch,b,d4,d5,d6)-Leq(1,ch,b,d4,d5,d6),...
                                        'Marker',mark,...
                                        'MarkerSize',8,...
                                        'color',colr,...
                                        'MarkerFaceColor',colr,...
                                        'DisplayName',labelstring);
                                    xlabel('Frequency [Hz]');
                                    ylabel('Leq (L50 & Lmax) [dB]');
                                    hold on
                                    
                                    hsubplot4 = subplot(2,2,4); % legend plot
                                    ylim([0 1]);
                                    xlim([0 1]);
                                    if ~isempty(labelstring)
                                        colr1 = [0 0 0];
                                    else
                                        colr1 = colr;
                                    end
                                    plot(0.03, plotnum/(numberofsubplots+1),...
                                        'Marker',mark,...
                                        'MarkerSize',8,...
                                        'MarkerFaceColor',colr1,...
                                        'color','k')
                                    hold on
                                    text(0.06,plotnum/(numberofsubplots+1),titlestring)
                                    text(0.6,Hind/(numberofH+1),labelstring,'color',colr)
                                    set(hsubplot4,'YTickLabel','',...
                                        'YTick',zeros(1,0),...
                                        'XTickLabel','',...
                                        'XTick',zeros(1,0))
                                    title('Legend');
                                end
                                if plottype ~= 22 && plottype ~= 23
                                    title(titlestring)
                                    hold on
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if strcmp(plotcategory,'Timeaxis') || strcmp(plotcategory,'Freqaxis')
        iplots = get(compplot,'Children');
        if length(iplots) > 1
            xlims = cell2mat(get(iplots,'Xlim'));
            set(iplots,'Xlim',[min(xlims(:,1)) max(xlims(:,2))])
            ylims = cell2mat(get(iplots,'Ylim'));
            set(iplots,'Ylim',[min(ylims(:,1)) max(ylims(:,2))])
            uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
                'Position', [0 0 65 30],...
                'Callback', 'setaxeslimits');
        end
    end
    if plottype == 22 || plottype == 23 % strcmp(plotcategory,'Scatter3axis')
        iplots = get(compplot,'Children');
        xlims = cell2mat(get(iplots,'Xlim'));
        ylims = cell2mat(get(iplots,'Ylim'));
        %Levellims = ylims(2,:);
        Timelims = xlims(4,:);
        Freqlims = xlims(2,:);
        if Timelims(1) < 0
            set(iplots(3),'Xlim',[0,Timelims(2)])
            set(iplots(4),'Xlim',[0,Timelims(2)])
        end
        if Freqlims(1) < 0
            set(iplots(2),'Xlim',[0,Freqlims(2)])
            set(iplots(4),'Ylim',[0,Freqlims(2)])
        end
        uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
            'Position', [0 0 65 30],...
            'Callback', 'set3axeslimits');
    end
    
    
    % **************************************************
    % THE FOLLOWING IS THE PREVIOUS AUDIO COMPARISON PLOT
    
    %     compplot = figure;
    %     for i = 1:length(selectedNodes)
    %         linea = [];
    %         axes = 'time';
    %         signaldata = selectedNodes(i).handle.UserData;
    %         if ~isempty(signaldata) && isfield(signaldata,'audio')
    %             plottype = get(handles.(matlab.lang.makeValidName([axes '_popup'])),'Value');
    %             t = linspace(0,length(signaldata.audio),length(signaldata.audio))./signaldata.fs;
    %             f = signaldata.fs .* ((1:length(signaldata.audio))-1) ./ length(signaldata.audio);
    %             if ~ismatrix(signaldata.audio)
    %                 if ndims(signaldata.audio) == 3, cmap = colormap(hsv(size(signaldata.audio,3))); end
    %                 if ndims(signaldata.audio) >= 4, cmap = colormap(copper(size(signaldata.audio,4))); end
    %                 try
    %                     linea(:,:) = signaldata.audio(:,str2double(get(handles.IN_nchannel,'String')),:);
    %                 catch
    %                     linea = zeros(size(t));
    %                 end
    %             else
    %                 cmap = colormap(lines(size(signaldata.audio,2)));
    %                 linea = signaldata.audio;
    %             end
    %             if isfield(signaldata,'cal') && handles.Settings.calibrationtoggle == 1
    %                 if size(linea,2) == length(signaldata.cal)
    %                     signaldata.cal(isnan(signaldata.cal)) = 0;
    %                     linea = linea.*repmat(10.^(signaldata.cal(:)'./20),length(linea),1);
    %                 elseif ~ismatrix(signaldata.audio) && size(signaldata.audio,2) == length(signaldata.cal)
    %                     signaldata.cal(isnan(signaldata.cal)) = 0;
    %                     cal = repmat(signaldata.cal(str2double(get(handles.IN_nchannel,'String'))),1,size(linea,2));
    %                     linea = linea.*repmat(10.^(cal(:)'./20),length(linea),1);
    %                 end
    %             end
    %             switch handles.Settings.specmagscale;
    %                 case {'Divided by length'}
    %                     spectscale = 1./length(linea);
    %                 case {'Times sqrt2/length'}
    %                     spectscale = 2.^0.5./length(linea);
    %                 otherwise
    %                     spectscale = 1;
    %             end
    %             if plottype == 1, linea = real(linea); end
    %             if plottype == 2, linea = linea.^2; end
    %             if plottype == 3, linea = 10.*log10(linea.^2); end
    %             if plottype == 4, linea = abs(hilbert(real(linea))); end
    %             if plottype == 5, linea = medfilt1(diff([angle(hilbert(real(linea))); zeros(1,size(linea,2))])*signaldata.fs/2/pi, 5); end
    %             if plottype == 6, linea = abs(linea); end
    %             if plottype == 7, linea = imag(linea); end
    %             if plottype == 8, linea = 10*log10(abs(fft(linea).*spectscale).^2); end %freq
    %             if plottype == 9, linea = (abs(fft(linea)).*spectscale).^2; end
    %             if plottype == 10, linea = abs(fft(linea)).*spectscale; end
    %             if plottype == 11, linea = real(fft(linea)).*spectscale; end
    %             if plottype == 12, linea = imag(fft(linea)).*spectscale; end
    %             if plottype == 13, linea = angle(fft(linea)); end
    %             if plottype == 14, linea = unwrap(angle(fft(linea))); end
    %             if plottype == 15, linea = angle(fft(linea)) .* 180/pi; end
    %             if plottype == 16, linea = unwrap(angle(fft(linea))) ./(2*pi); end
    %             if plottype == 17, linea = -diff(unwrap(angle(fft(linea)))).*length(linea)/(signaldata.fs*2*pi).*1000; end
    %             if strcmp(get(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Visible'),'on')
    %                 smoothfactor = get(handles.(matlab.lang.makeValidName(['smooth' axes '_popup'])),'Value');
    %                 if smoothfactor == 2, octsmooth = 1; end
    %                 if smoothfactor == 3, octsmooth = 3; end
    %                 if smoothfactor == 4, octsmooth = 6; end
    %                 if smoothfactor == 5, octsmooth = 12; end
    %                 if smoothfactor == 6, octsmooth = 24; end
    %                 if smoothfactor ~= 1, linea = octavesmoothing(linea, octsmooth, signaldata.fs); end
    %             end
    %             if length(selectedNodes) == 1
    %                 [r, c] = subplotpositions(size(linea,2), 0.5);
    %                 for j = 1:size(linea,2)
    %                     if plottype <= 7
    %                         subplot(r,c,j);
    %                         set(gca,'NextPlot','replacechildren','ColorOrder',cmap(j,:))
    %                         plot(t,real(linea(:,j))) % Plot signal in time domain
    %                         if ismatrix(signaldata.audio) && isfield(signaldata,'chanID'), title(signaldata.chanID{j,1}); end
    %                         if ~ismatrix(signaldata.audio) && isfield(signaldata,'bandID'), title(num2str(signaldata.bandID(1,j))); end
    %                         xlabel('Time [s]');
    %                     end
    %                     if plottype >= 8
    %                         h = subplot(r,c,j);
    %                         set(gca,'NextPlot','replacechildren','ColorOrder',cmap(j,:))
    %                         plot(f(1:length(linea(:,j))),linea(:,j));% Plot signal in frequency domain
    %                         if ismatrix(signaldata.audio) && isfield(signaldata,'chanID'), title(signaldata.chanID{j,1}); end
    %                         if ~ismatrix(signaldata.audio) && isfield(signaldata,'bandID'), title(num2str(signaldata.bandID(1,j))); end
    %                         xlabel('Frequency [Hz]');
    %                         if ischar(handles.Settings.frequencylimits)
    %                             xlim([f(2) signaldata.fs/2])
    %                         else
    %                             xlim(handles.Settings.frequencylimits)
    %                         end
    %                         log_check = get(handles.(matlab.lang.makeValidName(['log' axes '_chk'])),'Value');
    %                         if log_check == 1
    %                             set(h,'XScale','log')
    %                         else
    %                             set(h,'XScale','linear','XTickLabelMode','auto')
    %                         end
    %                     end
    %                 end
    %             else
    %                 if plottype <= 7
    %                     subplot(length(selectedNodes),1,i);
    %                     set(gca,'NextPlot','replacechildren','ColorOrder',cmap)
    %                     plot(t,real(linea)) % Plot signal in time domain
    %                     title(strrep(selectedNodes(i).getName.char,'_',' '))
    %                     xlabel('Time [s]');
    %                 end
    %                 if plottype >= 8
    %                     h = subplot(length(selectedNodes),1,i);
    %                     set(gca,'NextPlot','replacechildren','ColorOrder',cmap)
    %                     plot(f(1:length(linea)),linea);% Plot signal in frequency domain
    %                     title(strrep(selectedNodes(i).getName.char,'_',' '))
    %                     xlabel('Frequency [Hz]');
    %                     if ischar(handles.Settings.frequencylimits)
    %                         xlim([f(2) signaldata.fs/2])
    %                     else
    %                         xlim(handles.Settings.frequencylimits)
    %                     end
    %                     log_check = get(handles.(matlab.lang.makeValidName(['log' axes '_chk'])),'Value');
    %                     if log_check == 1
    %                         set(h,'XScale','log')
    %                     else
    %                         set(h,'XScale','linear','XTickLabelMode','auto')
    %                     end
    %                 end
    %             end
    %         end
    %     end
    %     iplots = get(compplot,'Children');
    %     if length(iplots) > 1
    %         xlims = cell2mat(get(iplots,'Xlim'));
    %         set(iplots,'Xlim',[min(xlims(:,1)) max(xlims(:,2))])
    %         ylims = cell2mat(get(iplots,'Ylim'));
    %         set(iplots,'Ylim',[min(ylims(:,1)) max(ylims(:,2))])
    %         uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
    %             'Position', [0 0 65 30],...
    %             'Callback', 'setaxeslimits');
    %     end
    
    % **************************************************
    
    
elseif handles.compareaudio == 0
    comparedata('main_stage1', handles.aarae);
end
handles.alternate = 0;
guidata(hObject, handles);
















% *************************************************************************
% *************************************************************************
%                               PLAY AUDIO
% *************************************************************************
% *************************************************************************



% *************************************************************************
% PLAY AUDIO
% *************************************************************************

% --- Executes on button press in play_btn.
function play_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to play_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the 'desktop'
hMain = getappdata(0,'hMain');
% Retrieve information from the selected leaf
audiodata = getappdata(hMain,'testsignal');
if isempty(audiodata)
    warndlg('No signal loaded!');
else
    if handles.alternate==1 && isfield(audiodata,'audio2')
        testsignal = real(audiodata.audio2);
    else
        testsignal = real(audiodata.audio);
    end
    testsignal = sum(sum(sum(sum(testsignal,3),4),5),6); % sum higher dimensions
    if size(testsignal,2) > 2, testsignal = sum(testsignal,2); end % mix channels if > 2
    testsignal = testsignal./max(max(abs(testsignal)));
    fs = audiodata.fs;
    nbits = 16;
    doesSupport = audiodevinfo(0, handles.odeviceid, fs, nbits, size(testsignal,2));
    if doesSupport && ismatrix(testsignal)
        % Play signal
        handles.player = audioplayer(testsignal,fs,nbits,handles.odeviceid);
        play(handles.player);
        selectedNodes = handles.mytree.getSelectedNodes;
        contents = cellstr(get(handles.device_popup,'String'));
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Played "' char(selectedNodes(1).getName) '" using ' contents{get(hObject,'Value')} '\n\n']);
    else
        warndlg('Device not supported for playback!');
    end
end
handles.alternate = 0;
guidata(hObject, handles);





% *************************************************************************
% STOP AUDIO
% *************************************************************************


% --- Executes on button press in stop_btn.
function stop_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to stop_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isplaying(handles.player)
    stop(handles.player);
end
guidata(hObject,handles);





% *************************************************************************
% PLAY AUDIO TIME-REVERSED
% *************************************************************************


% --- Executes on button press in playreverse_btn.
function playreverse_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to playreverse_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');

if isempty(audiodata)
    warndlg('No signal loaded!');
else
    if handles.alternate==1 && isfield(audiodata,'audio2')
        testsignal = real(audiodata.audio2);
    else
        testsignal = real(audiodata.audio);
    end
    testsignal = sum(sum(sum(sum(testsignal,3),4),5),6); % sum higher dimensions
    if size(testsignal,2) > 2, testsignal = sum(testsignal,2); end % mix channels if > 2
    testsignal = flipud(testsignal)./max(max(abs(testsignal)));
    fs = audiodata.fs;
    nbits = 16;
    doesSupport = audiodevinfo(0, handles.odeviceid, fs, nbits, size(testsignal,2));
    if doesSupport && ismatrix(testsignal)
        % Play signal
        handles.player = audioplayer(testsignal,fs,nbits,handles.odeviceid);
        play(handles.player);
        selectedNodes = handles.mytree.getSelectedNodes;
        contents = cellstr(get(handles.device_popup,'String'));
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Played "' char(selectedNodes(1).getName) '" using ' contents{get(hObject,'Value')} '\n\n']);
    else
        warndlg('Device not supported for playback!');
    end
end
handles.alternate = 0;
guidata(hObject, handles);






% *************************************************************************
% PLAY AUDIO IN RANDOM PHASE (STEADY-STATE RENDER)
% *************************************************************************


% --- Executes on button press in randphaseplay_btn.
function randphaseplay_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to randphaseplay_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');

if isempty(audiodata)
    warndlg('No signal loaded!');
else
    if handles.alternate==1 && isfield(audiodata,'audio2')
        testsignal = real(audiodata.audio2);
    else
        testsignal = real(audiodata.audio);
    end
    testsignal = sum(sum(sum(sum(testsignal,3),4),5),6); % sum higher dimensions
    if size(testsignal,2) > 2, testsignal = sum(testsignal,2); end % mix channels if > 2
    testsignal = testsignal./max(max(abs(testsignal)));
    len = length(testsignal);
    len = 2 * ceil(len/2);
    spectrum = fft(testsignal,len);
    magnitude = abs(spectrum);
    randphase = rand(len/2-1,1) .* 2 * pi;
    randphase = repmat([0; randphase; 0; flipud(-randphase)],1,size(testsignal,2));
    changed_spectrum = magnitude .* exp(1i * randphase);
    testsignal = ifft(changed_spectrum);
    fs = audiodata.fs;
    nbits = 16;
    doesSupport = audiodevinfo(0, handles.odeviceid, fs, nbits, size(testsignal,2));
    if doesSupport && ismatrix(testsignal)
        % Play signal
        handles.player = audioplayer(testsignal,fs,nbits,handles.odeviceid);
        play(handles.player);
        selectedNodes = handles.mytree.getSelectedNodes;
        contents = cellstr(get(handles.device_popup,'String'));
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Played "' char(selectedNodes(1).getName) '" using ' contents{get(hObject,'Value')} '\n\n']);
    else
        warndlg('Device not supported for playback!');
    end
end
handles.alternate = 0;
guidata(hObject, handles);







% *************************************************************************
% PLAY AUDIO WITH FLATTENED MAGNITUDE SPECTRUM
% *************************************************************************


% --- Executes on button press in flatmagplay_btn.
function flatmagplay_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to flatmagplay_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');

if isempty(audiodata)
    warndlg('No signal loaded!');
else
    if handles.alternate==1 && isfield(audiodata,'audio2')
        testsignal = real(audiodata.audio2);
    else
        testsignal = real(audiodata.audio);
    end
    testsignal = sum(sum(sum(sum(testsignal,3),4),5),6); % sum higher dimensions
    if size(testsignal,2) > 2, testsignal = sum(testsignal,2); end % mix channels if > 2
    testsignal = testsignal./max(max(abs(testsignal)));
    len = length(testsignal);
    %len = 2 .* ceil(len./2);
    spectrum = fft(testsignal,len);
    phase = angle(spectrum);
    rmsmag = mean(abs(spectrum).^2).^0.5; % root mean square magnitude
    
    % combine magnitude with phase
    changed_spectrum = ones(len,size(testsignal,2)).*repmat(rmsmag,size(testsignal,1),1) .* exp(1i .* phase);
    changed_spectrum(1) = 0; % make DC zero
    changed_spectrum(ceil(len/2)) = 0; % make Nyquist zero
    testsignal = ifft(changed_spectrum);
    fs = audiodata.fs;
    nbits = 16;
    doesSupport = audiodevinfo(0, handles.odeviceid, fs, nbits, size(testsignal,2));
    if doesSupport && ismatrix(testsignal)
        % Play signal
        handles.player = audioplayer(testsignal,fs,nbits,handles.odeviceid);
        play(handles.player);
        selectedNodes = handles.mytree.getSelectedNodes;
        contents = cellstr(get(handles.device_popup,'String'));
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Played "' char(selectedNodes(1).getName) '" using ' contents{get(hObject,'Value')} '\n\n']);
    else
        warndlg('Device not supported for playback!');
    end
end
handles.alternate = 0;
guidata(hObject, handles);







% *************************************************************************
% PLAY AUDIO CONVOLVED WITH REFERENCE AUDIO
% *************************************************************************


% --- Executes on button press in convplay_btn.
function convplay_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to convplay_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMain = getappdata(0,'hMain');
audiodata = getappdata(hMain,'testsignal');

if isempty(audiodata)
    warndlg('No signal loaded!');
else
    if handles.alternate==1 && isfield(audiodata,'audio2')
        testsignal = real(audiodata.audio2);
    else
        testsignal = real(audiodata.audio);
    end
    testsignal = sum(sum(sum(sum(testsignal,3),4),5),6); % sum higher dimensions
    if size(testsignal,2) > 2, testsignal = sum(testsignal,2); end % mix channels if > 2
    testsignal = testsignal./max(max(abs(testsignal)));
    fs = audiodata.fs;
    if handles.reference_audio.fs ~= fs
        gcd_fs = gcd(handles.reference_audio.fs,fs); % greatest common denominator
        reference_audio = resample(handles.reference_audio.audio,fs/gcd_fs,handles.reference_audio.fs/gcd_fs);
    else
        reference_audio = handles.reference_audio.audio;
    end
    reference_audio = repmat(reference_audio,1,size(testsignal,2));
    len1 = length(testsignal);
    len2 = length(reference_audio);
    outputlength = len1 + len2 - 1;
    testsignal = ifft(fft(testsignal, outputlength) .* fft(reference_audio, outputlength));
    testsignal = testsignal./max(max(max(abs(testsignal))));
    nbits = 16;
    doesSupport = audiodevinfo(0, handles.odeviceid, fs, nbits, size(testsignal,2));
    if doesSupport && ismatrix(testsignal)
        % Play signal
        handles.player = audioplayer(testsignal,fs,nbits,handles.odeviceid);
        play(handles.player);
        selectedNodes = handles.mytree.getSelectedNodes;
        contents = cellstr(get(handles.device_popup,'String'));
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Played "' char(selectedNodes(1).getName) '" using ' contents{get(hObject,'Value')} '\n\n']);
    else
        warndlg('Device not supported for playback!');
    end
end
handles.alternate = 0;
guidata(hObject, handles);



% *************************************************************************
% AUDIO PLAYBACK DEVICE POPUP
% *************************************************************************

% --- Executes on selection change in device_popup.
function device_popup_Callback(hObject, ~, handles) %#ok
% hObject    handle to device_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns device_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from device_popup

selection = get(hObject,'Value');
handles.odeviceid = handles.odeviceidlist(selection);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function device_popup_CreateFcn(hObject, ~, handles) %#ok
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









% *************************************************************************
% *************************************************************************
%                               WORKFLOWS
% *************************************************************************
% *************************************************************************
% This should be almost the same as the main Processors function


% --- Executes on button press in RunWorkflowButton.
function RunWorkflowButton_Callback(hObject, ~, handles)
% hObject    handle to RunWorkflowButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedNodes = handles.mytree.getSelectedNodes;
funcallback = [];
for nleafs = 1:length(selectedNodes)
    handles.nleafs = nleafs;
    guidata(hObject,handles)
    signaldata = selectedNodes(nleafs).handle.UserData;
    if ~isempty(signaldata)
        set(hObject,'BackgroundColor','red');
        set(hObject,'Enable','off');
        set(handles.CloseFiguresButton,'Enable','off');
        name = selectedNodes(nleafs).getName.char;
        errorflag = false;
        funname = 'AARAE_workflow_processor';
        try
            if ~isempty(funcallback) && strcmp(funname,funcallback.name)
                processed = feval(funname,signaldata,funcallback.inarg{:});
            else
                processed = feval(funname,signaldata);
            end
            if isfield(processed,'funcallback')
                funcallback = processed.funcallback;
                [~,funcallback.name] = fileparts(funcallback.name);
            end
        catch err
            processed = [];
            msgString = getReport(err);
            disp(['AARAE processor error running ' funname '.']);
            disp(msgString); % displays the error message without creating an error.
            fprintf(handles.fid,['AARAE processor error running ' funname 'with ' name(nleafs)'.\n']);
            fprintf(handles.fid,msgString);
            fprintf(handles.fid,'.\n');
            errorflag = true;
        end
    else
        % no audio input
        set(hObject,'BackgroundColor','red');
        set(hObject,'Enable','off');
        set(handles.CloseFiguresButton,'Enable','off');
        name = [];
        errorflag = false;
        funname = 'AARAE_workflow_processor';
        guidata(hObject,handles)
        try
            processed = feval(funname,[]);
            %handles = guidata(hObject,handles);
            if isfield(processed,'funcallback')
                funcallback = processed.funcallback;
                [~,funcallback.name] = fileparts(funcallback.name);
            end
        catch err
            processed = [];
            msgString = getReport(err);
            disp(['AARAE processor error running ' funname '.']);
            disp(msgString); % displays the error message without creating an error.
            fprintf(handles.fid,['AARAE processor error running ' funname '.\n']);
            fprintf(handles.fid,msgString);
            fprintf(handles.fid,'.\n');
            errorflag = true;
        end
    end
    aarae_fig = findobj('Tag','aarae');
    handles = guidata(aarae_fig);
    if ~isempty(processed)
        % Generate new leaf and update tree
        newleaf = cell(1,1);
        newleaf{1,1} = matlab.lang.makeValidName([name '_' char(processed.funcallback.inarg{1,1})]);
        if ~isempty(signaldata)
            if isstruct(processed)
                dif = intersect(fieldnames(signaldata),fieldnames(processed));
                newdata = signaldata;
                for i = 1:size(dif,1)
                    newdata.(dif{i,1}) = processed.(dif{i,1});
                end
                newfields = setxor(fieldnames(signaldata),fieldnames(processed));
                for i = 1:size(newfields,1)
                    if ~isfield(newdata,newfields{i,1})
                        newdata.(newfields{i,1}) = processed.(newfields{i,1});
                    end
                end
                
                % remove audio and audio2 if there are tables (from
                % AARAE_workflow_processor). Usually tables should
                % only be created by analysers (not processors).
                if isfield(newdata,'tables') && isfield(newdata,'audio')
                    newdata = rmfield(newdata,'audio');
                end
                if isfield(newdata,'tables') && isfield(newdata,'audio2')
                    newdata = rmfield(newdata,'audio2');
                end
                
            else
                % should not be necessary
                newdata = signaldata;
                newdata.audio = processed;
            end
        else
            newdata = processed;
        end
        iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
        leafname = isfield(handles,matlab.lang.makeValidName(newleaf{1,1}));
        if leafname == 1
            index = 1;
            % This while cycle is just to make sure no signals are
            % overwriten
            if length(matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) >= namelengthmax-2, newleaf{1,1} = newleaf{1,1}(1:round(end/2)); end
            while isfield(handles,matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) == 1
                index = index + 1;
            end
            newleaf{1,1} = matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)]);
        end
        
        newdata_fields = fieldnames(newdata);
        newdata_emptyfields = structfun(@isempty,newdata);
        newdata = rmfield(newdata,newdata_fields(newdata_emptyfields));
        newdata.name = matlab.lang.makeValidName(newleaf{1,1});
        newdata = checkcal(newdata);
        % Save as you go
        save([cd '/Utilities/Backup/' newleaf{1,1} '.mat'], 'newdata','-v7.3');
        
        if strcmp(newdata.datatype,'testsignals')
            handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
            handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
            handles.testsignals.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
            handles.mytree.reloadNode(handles.testsignals);
            handles.mytree.expand(handles.testsignals);
        elseif strcmp(newdata.datatype,'measurements')
            handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
            handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
            handles.measurements.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
            handles.mytree.reloadNode(handles.measurements);
            handles.mytree.expand(handles.measurements);
        elseif strcmp(newdata.datatype,'results')
            if isfield(newdata,'audio')
                iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
            elseif isfield(newdata,'tables')
                iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/HDF_grid.gif');
            else
                iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/notesicon.gif');
            end
            handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
            handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
            handles.results.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
            handles.mytree.reloadNode(handles.results);
            handles.mytree.expand(handles.results);
        else
            handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
            handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
            handles.processed.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
            handles.mytree.reloadNode(handles.processed);
            handles.mytree.expand(handles.processed);
        end
        fprintf(handles.fid, ['%% ' datestr(now,16) ' - Ran ' newdata.funcallback.inarg{1,1} ' (AARAE workflow)\n']);
        % Log verbose metadata
        %callbackstring = logaudioleaffields(newdata);
        logaudioleaffields(newdata);
        if isfield(handles,'choosefromhigherdims')
            handles.choosefromhigherdims = [];
        end
    else
        newleaf{1,1} = [];
    end
    if errorflag
        h=warndlg('One or more attempts to run a workflow processor failed due to an error. Please refer to the error report in the Matlab Command Window.','AARAE info','modal');
        uiwait(h)
    end
    h = findobj('type','figure','-not','tag','aarae');
    if ~isempty(h)
        %         index = 1;
        %         filename = dir([cd '/Utilities/Temp/' newdata.funcallback.inarg{1,1} num2str(index) '.fig']);
        %         if ~isempty(filename)
        %             while isempty(dir([cd '/Utilities/Temp/' newdata.funcallback.inarg{1,1} num2str(index) '.fig'])) == 0
        %                 index = index + 1;
        %             end
        %         end
        
        
        for i = 1:length(h)
            % Write figure's UserData property
            UserData = {'Environment','AARAE';...
                'Function name',funname;...
                'Input',char(selectedNodes(nleafs).getName);...
                'Callback string', 'Workflow';...
                'Time created',datestr(now,31)};
            set(h(i),'UserData',UserData)
            try
                Figtag = get(h(i),'Tag');
                if ~isempty(Figtag)
                    if strcmp(Figtag,'AARAE table')
                        Figtag = 'T'; % indicate table in filename
                    else
                        Figtag = 'C'; % indicate chart in filename
                    end
                else
                    Figtag = 'C'; % indicate chart in filename
                end
                index = 1;
                filename = dir([cd '/Utilities/Temp/' newdata.funcallback.inarg{1,1} ' ' funname Figtag num2str(index) '.fig']);
                if ~isempty(filename)
                    while isempty(dir([cd '/Utilities/Temp/' newdata.funcallback.inarg{1,1} ' ' funname Figtag num2str(index) '.fig'])) == 0
                        index = index + 1;
                    end
                end
                filename = [newdata.funcallback.inarg{1,1} ' ' funname Figtag num2str(index) '.fig'];
                saveas(h(i),[cd '/Utilities/Temp/' filename]);
                fprintf(handles.fid,['%% Result figure name: ', filename, ', temporarily stored in /Utilities/Temp/\n']);
            catch
            end
            
            
            %             saveas(h(i),[cd '/Utilities/Temp/' newdata.funcallback.inarg{1,1} num2str(index) '.fig']);
            %             fprintf(handles.fid,['%% Result figure name: ', newdata.funcallback.inarg{1,1} num2str(index), '.fig, temporarily stored in /Utilities/Temp/\n']);
            %             index = index + 1;
        end
        results = dir([cd '/Utilities/Temp']);
        set(handles.result_box,'String',[' ';cellstr({results(3:length(results)).name}')]);
    end
    if length(selectedNodes) > 1, delete(h); end
end
set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
set(hObject,'Enable','on');
% The following causes problems, so we won't change the selected node
% if ~isempty(newleaf{1,1})
%     handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(newleaf{1,1})));
% end
set([handles.clrall_btn,handles.export_btn],'Enable','on')
set(handles.CloseFiguresButton,'Enable','on');
java.lang.Runtime.getRuntime.gc % Java garbage collection
guidata(hObject,handles);










% *************************************************************************
% *************************************************************************
%                              PROCESSORS
% *************************************************************************
% *************************************************************************





% *************************************************************************
% PROCESSORS BOX GUI FUNCTIONS
% *************************************************************************

% --- Executes on selection change in procat_box.
function procat_box_Callback(hObject, ~, handles) %#ok
% hObject    handle to procat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns procat_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from procat_box

% Displays the processes available for the selected process category
%hMain = getappdata(0,'hMain');
%signaldata = getappdata(hMain,'testsignal');

contents = cellstr(get(hObject,'String'));
handles.procat = contents{get(hObject,'Value')};
processes = what([cd '/Processors/' handles.procat]);

try
    if ~isempty(processes.m)
        set(handles.proc_box,'Visible','on','String',[' ';cellstr(processes.m)],'Value',1);
        set([handles.proc_btn,handles.proc_help_btn],'Visible','off');
    else
        set(handles.proc_box,'Visible','off');
        set([handles.proc_btn,handles.proc_help_btn],'Visible','off');
    end
catch
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function procat_box_CreateFcn(hObject, ~, handles) %#ok
% hObject    handle to procat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Displays the process categories available in the Processors folder

curdir = cd;
processes = removedotfiles(dir([curdir '/Processors']));
%set(hObject,'String',[' ';cellstr({processes(3:length(processes)).name}')]);
set(hObject,'String',[' ';cellstr({processes.name}')]);
handles.funname = [];
guidata(hObject,handles)

% --- Executes on selection change in proc_box.
function proc_box_Callback(hObject, ~, handles) %#ok
% hObject    handle to proc_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns proc_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from proc_box
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
[~,funname] = fileparts(selection);
if ~strcmp(selection,' ')
    handles.funname = funname;
    helptext = evalc(['help ' funname]);
    set(hObject,'Tooltip',helptext);
    set([handles.proc_btn,handles.proc_help_btn],'Visible','on');
    set(handles.proc_btn,'BackgroundColor',[0.94 0.94 0.94]);
    set([handles.proc_btn,handles.proc_help_btn],'Enable','on');
else
    handles.funname = [];
    set(hObject,'Tooltip','');
    set([handles.proc_btn,handles.proc_help_btn],'Visible','off');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function proc_box_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to proc_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% *************************************************************************
% PROCESS AUDIO - MAIN FUNCTION FOR PROCESSING
% *************************************************************************

% --- Executes on button press in proc_btn.
function proc_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to proc_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedNodes = handles.mytree.getSelectedNodes;
funcallback = [];
for nleafs = 1:length(selectedNodes)
    handles.nleafs = nleafs;
    guidata(hObject,handles)
    signaldata = selectedNodes(nleafs).handle.UserData;
    if ~isempty(signaldata)
        set(hObject,'BackgroundColor','red');
        set(hObject,'Enable','off');
        set(handles.CloseFiguresButton,'Enable','off');
        % Processes the selected leaf using the selected process from proc_box
        contents = cellstr(get(handles.procat_box,'String'));
        category = contents{get(handles.procat_box,'Value')};
        contents = cellstr(get(handles.proc_box,'String'));
        file = contents(get(handles.proc_box,'Value'));
        name = selectedNodes(nleafs).getName.char;
        errorflag = false;
        for multi = 1:size(file,1)
            processed = [];
            [~,~,ext] = fileparts(file{multi,1});
            if strcmp(ext,'.mat')
                warndlg('Option not yet available','AARAE info','modal')
                %                content = load([cd '/Processors/' category '/' num2str(signaldata.fs) 'Hz/' char(file(multi,:))]);
                %                filterbank = content.filterbank;
                %                w = whos('filterbank');
                %                if strcmp(w.class,'dfilt.df2sos')
                %                    for i = 1:length(filterbank)
                %                        for j = 1:min(size(signaldata.audio))
                %                            processed(:,j,i) = filter(filterbank(1,i),signaldata.audio(:,j));
                %                        end
                %                    end
                %                    bandID = [];
                %                elseif strcmp(w.class,'double')
                %                    processed = filter(filterbank,1,signaldata.audio);
                %                end
            elseif strcmp(ext,'.m')
                [~,funname] = fileparts(char(file(multi,:)));
                try
                    if ~isempty(funcallback) && strcmp(funname,funcallback.name)
                        processed = feval(funname,signaldata,funcallback.inarg{:});
                    else
                        processed = feval(funname,signaldata);
                    end
                    if isfield(processed,'funcallback')
                        funcallback = processed.funcallback;
                        [~,funcallback.name] = fileparts(funcallback.name);
                    end
                catch err
                    processed = [];
                    msgString = getReport(err);
                    disp(['AARAE processor error running ' funname '.']);
                    disp(msgString); % displays the error message without creating an error.
                    fprintf(handles.fid,['AARAE processor error running ' funname '.\n']);
                    fprintf(handles.fid,msgString);
                    fprintf(handles.fid,'.\n');
                    errorflag = true;
                end
            else
                processed = [];
            end
            aarae_fig = findobj('Tag','aarae');
            handles = guidata(aarae_fig);
            if ~isempty(processed)
                % Generate new leaf and update tree
                newleaf = cell(1,1);
                % remove suffix if present
                C = strsplit(char(file(multi,:)),'.');
                file1 = C{1};
                newleaf{1,1} = matlab.lang.makeValidName([name '_' file1]);
                if ~isempty(signaldata)
                    if isstruct(processed)
                        dif = intersect(fieldnames(signaldata),fieldnames(processed));
                        newdata = signaldata;
                        for i = 1:size(dif,1)
                            newdata.(dif{i,1}) = processed.(dif{i,1});
                        end
                        newfields = setxor(fieldnames(signaldata),fieldnames(processed));
                        for i = 1:size(newfields,1)
                            if ~isfield(newdata,newfields{i,1})
                                newdata.(newfields{i,1}) = processed.(newfields{i,1});
                            end
                        end
                        
                        % remove audio and audio2 if there are tables (from
                        % AARAE_workflow_processor). Usually tables should
                        % only be created by analysers (not processors).
                        if isfield(newdata,'tables') && isfield(newdata,'audio')
                            newdata = rmfield(newdata,'audio');
                        end
                        if isfield(newdata,'tables') && isfield(newdata,'audio2')
                            newdata = rmfield(newdata,'audio2');
                        end
                        
                    else
                        newdata = signaldata;
                        newdata.audio = processed;
                    end
                    if ~strcmp(funname,'AARAE_workflow_processor') || ~isfield(newdata,'datatype')
                        newdata.datatype = 'processed';
                    end
                    iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
                    leafname = isfield(handles,matlab.lang.makeValidName(newleaf{1,1}));
                    if leafname == 1
                        index = 1;
                        % This while cycle is just to make sure no signals are
                        % overwriten
                        if length(matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) >= namelengthmax-2, newleaf{1,1} = newleaf{1,1}(1:round(end/2)); end
                        while isfield(handles,matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) == 1
                            index = index + 1;
                        end
                        newleaf{1,1} = matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)]);
                    end
                    
                    newdata_fields = fieldnames(newdata);
                    newdata_emptyfields = structfun(@isempty,newdata);
                    newdata = rmfield(newdata,newdata_fields(newdata_emptyfields));
                    newdata.name = matlab.lang.makeValidName(newleaf{1,1});
                    newdata = checkcal(newdata);
                    newdata = addhistory(newdata,'Processed');
                    % Save as you go
                    save([cd '/Utilities/Backup/' newleaf{1,1} '.mat'], 'newdata','-v7.3');
                    
                    
                    if strcmp(newdata.datatype,'testsignals')
                        handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
                        handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
                        handles.testsignals.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
                        handles.mytree.reloadNode(handles.testsignals);
                        handles.mytree.expand(handles.testsignals);
                        set([handles.clrall_btn,handles.export_btn],'Enable','on')
                    elseif strcmp(newdata.datatype,'measurements')
                        handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
                        handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
                        handles.measurements.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
                        handles.mytree.reloadNode(handles.measurements);
                        handles.mytree.expand(handles.measurements);
                        set([handles.clrall_btn,handles.export_btn],'Enable','on')
                    elseif strcmp(newdata.datatype,'results')
                        if isfield(newdata,'audio')
                            iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
                        elseif isfield(newdata,'tables')
                            iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/HDF_grid.gif');
                        else
                            iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/notesicon.gif');
                        end
                        handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
                        handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
                        handles.results.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
                        handles.mytree.reloadNode(handles.results);
                        handles.mytree.expand(handles.results);
                        set([handles.clrall_btn,handles.export_btn],'Enable','on')
                    else
                        handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
                        handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = newdata;
                        handles.processed.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
                        handles.mytree.reloadNode(handles.processed);
                        handles.mytree.expand(handles.processed);
                        set([handles.clrall_btn,handles.export_btn],'Enable','on')
                    end
                end
                fprintf(handles.fid, ['%% ' datestr(now,16) ' - Processed "' name '" using ' funname ' in ' category '\n']);
                % Log verbose metadata
                callbackstring=logaudioleaffields(newdata);
                if isfield(handles,'choosefromhigherdims')
                    handles.choosefromhigherdims = [];
                end
            else
                newleaf{1,1} = [];
            end
        end
        if errorflag
            h=warndlg('One or more attempts to run a processor failed due to an error. Please refer to the error report in the Matlab Command Window.','AARAE info','modal');
            uiwait(h)
        end
    end
    h = findobj('type','figure','-not','tag','aarae');
    if ~isempty(h)
        %         index = 1;
        %         filename = dir([cd '/Utilities/Temp/' handles.funname num2str(index) '.fig']);
        %         if ~isempty(filename)
        %             while isempty(dir([cd '/Utilities/Temp/' handles.funname num2str(index) '.fig'])) == 0
        %                 index = index + 1;
        %             end
        %         end
        %         for i = 1:length(h)
        %             saveas(h(i),[cd '/Utilities/Temp/' handles.funname num2str(index) '.fig']);
        %             fprintf(handles.fid,['%% Result figure name: ', handles.funname num2str(index), '.fig, temporarily stored in /Utilities/Temp/\n']);
        %             index = index + 1;
        %         end
        for i = 1:length(h)
            % Write figure's UserData property
            if ~exist('callbackstring','var')
                callbackstring = '';
            end
            UserData = {'Environment','AARAE';...
                'Function name',funname;...
                'Input',char(selectedNodes(nleafs).getName);...
                'Callback string', callbackstring;...
                'Time created',datestr(now,31)};
            set(h(i),'UserData',UserData)
            try
                Figtag = get(h(i),'Tag');
                if ~isempty(Figtag)
                    if strcmp(Figtag,'AARAE table')
                        Figtag = 'T'; % indicate table in filename
                    else
                        Figtag = 'C'; % indicate chart in filename
                    end
                else
                    Figtag = 'C'; % indicate chart in filename
                end
                index = 1;
                filename = dir([cd '/Utilities/Temp/' char(selectedNodes(nleafs).getName) ' ' funname Figtag num2str(index) '.fig']);
                if ~isempty(filename)
                    while isempty(dir([cd '/Utilities/Temp/' char(selectedNodes(nleafs).getName) ' ' funname Figtag num2str(index) '.fig'])) == 0
                        index = index + 1;
                    end
                end
                filename = [char(selectedNodes(nleafs).getName) ' ' funname Figtag num2str(index) '.fig'];
                saveas(h(i),[cd '/Utilities/Temp/' filename]);
                fprintf(handles.fid,['%% Result figure name: ', filename, ', temporarily stored in /Utilities/Temp/\n']);
            catch
            end
        end
        results = dir([cd '/Utilities/Temp']);
        set(handles.result_box,'String',[' ';cellstr({results(3:length(results)).name}')]);
    end
    if length(selectedNodes) > 1 || size(file,1) > 1, delete(h); end
    set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
    set(hObject,'Enable','on');
end
if ~isempty(newleaf{1,1})
    handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(newleaf{1,1})));
end
set(handles.CloseFiguresButton,'Enable','on');
java.lang.Runtime.getRuntime.gc % Java garbage collection
guidata(hObject,handles);





% *************************************************************************
% PROCESSORS HELP BUTTON
% *************************************************************************

% --- Executes on button press in proc_help_btn.
function proc_help_btn_Callback(~, ~, handles) %#ok : Executed when help button for processors is clicked
% hObject    handle to proc_help_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.proc_box,'String'));
selection = contents{get(handles.proc_box,'Value')};
eval(['doc ' selection])











% *************************************************************************
% *************************************************************************
%                              ANALYSERS
% *************************************************************************
% *************************************************************************



% *************************************************************************
% SELECTION CHANGE IN ANALYSERS FUNCAT BOX
% *************************************************************************

% --- Executes on selection change in funcat_box.
function funcat_box_Callback(hObject, ~, handles) %#ok
% hObject    handle to funcat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns funcat_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from funcat_box

% Displays the available analysers for the selected processing category
contents = cellstr(get(hObject,'String'));
handles.funcat = contents{get(hObject,'Value')};
analyzers = what([cd '/Analysers/' handles.funcat]);
try
    if ~isempty(cellstr(analyzers.m))
        set(handles.fun_box,'Visible','on','String',[' ';cellstr(analyzers.m)],'Value',1,'Tooltip','');
        set([handles.analyze_btn,handles.analyser_help_btn],'Visible','off');
    else
        set(handles.fun_box,'Visible','off');
        set([handles.analyze_btn,handles.analyser_help_btn],'Visible','off');
    end
catch
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function funcat_box_CreateFcn(hObject, ~, handles) %#ok
% hObject    handle to funcat_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Populate function box with the function available in the folder 'Analysers'
curdir = cd;
analyzers = removedotfiles(dir([curdir '/Analysers']));
%set(hObject,'String',[' ';cellstr({analyzers(3:length(analyzers)).name}')]);
set(hObject,'String',[' ';cellstr({analyzers.name}')]);
handles.funname = [];
guidata(hObject,handles)

% To do: review whether these spaces at the top of the list are of any use
% - perhaps they should not be added.


% *************************************************************************
% ANALYSE AUDIO
% *************************************************************************

% --- Executes on button press in analyze_btn.
function analyze_btn_Callback(hObject, ~, handles) %#ok
% hObject    handle to analyze_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Analyses the selected leaf using the function selected from fun_box
% Call the 'desktop'
% hMain = getappdata(0,'hMain');
% audiodata = getappdata(hMain,'testsignal');
selectedNodes = handles.mytree.getSelectedNodes;
funcallback = [];
h = findobj('type','figure','-not','tag','aarae');
delete(h)
set(handles.CloseFiguresButton,'Enable','off');
for nleafs = 1:length(selectedNodes)
    handles.nleafs = nleafs;
    guidata(hObject,handles)
    audiodata = selectedNodes(nleafs).handle.UserData;
    % Evaluate selected function for the leaf selected from the tree
    if ~isempty(handles.funname) && ~isempty(audiodata)
        set(hObject,'BackgroundColor','red');
        set(hObject,'Enable','off');
        contents = cellstr(get(handles.fun_box,'String'));
        file = contents(get(handles.fun_box,'Value'));
        errorflag = false;
        for multi = 1:size(file,1)
            [~,funname] = fileparts(char(file(multi,:)));
            if nargout(funname) == 1 || nargout(funname) == -2
                try
                    if ~isempty(funcallback) && strcmp(funname,funcallback.name)
                        out = feval(funname,audiodata,funcallback.inarg{:});
                    else
                        out = feval(funname,audiodata);
                    end
                    if isfield(out,'funcallback')
                        funcallback = out.funcallback;
                        [~,funcallback.name] = fileparts(funcallback.name);
                    end
                catch err
                    out = [];
                    msgString = getReport(err);
                    disp(['AARAE analyser error running ' funname '.']);
                    disp(msgString); % displays the error message without creating an error.
                    fprintf(handles.fid,['AARAE analyser error running ' funname '.\n']);
                    fprintf(handles.fid,msgString);
                    fprintf(handles.fid,'.\n');
                    errorflag = true;
                end
            else
                out = [];
                try
                    feval(funname,audiodata);
                catch err
                    out = [];
                    msgString = getReport(err);
                    disp(['AARAE analyser error running ' funname '.']);
                    disp(msgString); % displays the error message without creating an error.
                    fprintf(handles.fid,['AARAE analyser error running ' funname '.\n']);
                    fprintf(handles.fid,msgString);
                    fprintf(handles.fid,'.\n');
                    errorflag = true;
                end
            end
            aarae_fig = findobj('Tag','aarae');
            handles = guidata(aarae_fig);
            newleaf = cell(1,1);
            newleaf{1,1} = matlab.lang.makeValidName([char(selectedNodes(nleafs).getName) '_' funname]);
            if ~isempty(out)
                signaldata = out;
                signaldata.datatype = 'results';
                if isfield(signaldata,'audio')
                    iconPath = fullfile(matlabroot,'/toolbox/fixedpoint/fixedpointtool/resources/plot.png');
                elseif isfield(signaldata,'tables')
                    iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/HDF_grid.gif');
                else
                    iconPath = fullfile(matlabroot,'/toolbox/matlab/icons/notesicon.gif');
                end
                if length(fieldnames(out)) ~= 1
                    leafname = isfield(handles,matlab.lang.makeValidName(newleaf{1,1}));
                    if leafname == 1
                        index = 1;
                        % This while cycle is just to make sure no signals are
                        % overwriten
                        if length(matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) >= namelengthmax-2, newleaf{1,1} = newleaf{1,1}(1:round(end/2)); end
                        while isfield(handles,matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)])) == 1
                            index = index + 1;
                        end
                        newleaf{1,1} = matlab.lang.makeValidName([newleaf{1,1},'_',num2str(index)]);
                    end
                    
                    signaldata_fields = fieldnames(signaldata);
                    signaldata_emptyfields = structfun(@isempty,signaldata);
                    signaldata = rmfield(signaldata,signaldata_fields(signaldata_emptyfields));
                    signaldata.name = newleaf{1,1};
                    signaldata = checkcal(signaldata);
                    signaldata = addhistory(signaldata,'Analysed');
                    % Save as you go
                    save([cd '/Utilities/Backup/' newleaf{1,1} '.mat'], 'signaldata','-v7.3');
                    
                    handles.(matlab.lang.makeValidName(newleaf{1,1})) = uitreenode('v0', newleaf{1,1},  newleaf{1,1},  iconPath, true);
                    handles.(matlab.lang.makeValidName(newleaf{1,1})).UserData = signaldata;
                    handles.results.add(handles.(matlab.lang.makeValidName(newleaf{1,1})));
                    handles.mytree.reloadNode(handles.results);
                    handles.mytree.expand(handles.results);
                    %handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(newleaf)));
                    set([handles.clrall_btn,handles.export_btn],'Enable','on')
                end
                fprintf(handles.fid, ['%% ' datestr(now,16) ' - Analysed "' char(selectedNodes(nleafs).getName) '" using ' funname ' in ' handles.funcat '\n']);% In what category???
                % Log verbose metadata
                callbackstring = logaudioleaffields(signaldata);
                if isfield(handles,'choosefromhigherdims')
                    handles.choosefromhigherdims = [];
                end
                %                 % Log contents of results tables - THIS IS NOW IN logaudioleaffields(signaldata);
                %                 if isfield(signaldata,'tables')
                %                     for tt = 1:size(signaldata.tables,2)
                %                         try
                %                         fprintf(handles.fid, ['%%  RESULTS TABLE: ', signaldata.tables(tt).Name,'\n']);
                %                         rownamelen = zeros(length(signaldata.tables(tt).RowName),1);
                %                         for rw = 1:length(signaldata.tables(tt).RowName)
                %                             rownamelen(rw) = length(signaldata.tables(tt).RowName{rw});
                %                         end
                %                         maxrownamelen = max(rownamelen);
                %                         if maxrownamelen < 15, maxrownamelen =15; end
                %                         fprintf(handles.fid, ['%%,', repmat(' ',[1,maxrownamelen+1])]);
                %                         for cl = 1:length(signaldata.tables(tt).ColumnName)
                %                             if cl < length(signaldata.tables(tt).ColumnName)
                %                                 colnamestring = char(signaldata.tables(tt).ColumnName(cl));
                %                                 if length(colnamestring)>14,colnamestring = colnamestring(1:14); end
                %                                 fprintf(handles.fid, [colnamestring,',',repmat(' ',[1,15-length(colnamestring)])]);
                %                             else
                %                                 colnamestring = char(signaldata.tables(tt).ColumnName(cl));
                %                                 if length(colnamestring)>14,colnamestring = colnamestring(1:14); end
                %                                 fprintf(handles.fid, [colnamestring,'\n']);
                %                             end
                %                         end
                %
                %                         for rw = 1:length(signaldata.tables(tt).RowName)
                %                             fprintf(handles.fid, ['%% ', char(signaldata.tables(tt).RowName(rw)),',',repmat(' ',[1,maxrownamelen-rownamelen(rw)+1])]);
                %                             for cl = 1:size(signaldata.tables(tt).Data,2)
                %                                 if cl < size(signaldata.tables(tt).Data,2)
                %                                     fprintf(handles.fid, [num2str(signaldata.tables(tt).Data(rw,cl)),',',repmat(' ',[1,15-length(num2str(signaldata.tables(tt).Data(rw,cl)))])]);
                %                                 else
                %                                     fprintf(handles.fid, [num2str(signaldata.tables(tt).Data(rw,cl)),'\n']);
                %                                 end
                %                             end
                %                         end
                %                         catch
                %                             fprintf(handles.fid,'Table format not interpreted for logging');
                %                         end
                %                         fprintf(handles.fid,'\n');
                %                     end
                %                 end
                
                h = findobj('type','figure','-not','tag','aarae');
                if ~isempty(h)
                    %                 index = 1;
                    %                 filename = dir([cd '/Utilities/Temp/' char(selectedNodes(nleafs).getName) funname num2str(index) '.fig']);
                    %                 if ~isempty(filename)
                    %                     while isempty(dir([cd '/Utilities/Temp/' char(selectedNodes(nleafs).getName) funname num2str(index) '.fig'])) == 0
                    %                         index = index + 1;
                    %                     end
                    %                 end
                    for i = 1:length(h)
                        % Write figure's UserData property
                        if ~exist('callbackstring','var')
                            callbackstring = '';
                        end
                        UserData = {'Environment','AARAE';...
                            'Function name',funname;...
                            'Input',char(selectedNodes(nleafs).getName);...
                            'Callback string', callbackstring;...
                            'Time created',datestr(now,31)};
                        set(h(i),'UserData',UserData)
                        try
                            Figtag = get(h(i),'Tag');
                            if ~isempty(Figtag)
                                if strcmp(Figtag,'AARAE table')
                                    Figtag = 'T'; % indicate table in filename
                                else
                                    Figtag = 'C'; % indicate chart in filename
                                end
                            else
                                Figtag = 'C'; % indicate chart in filename
                            end
                            index = 1;
                            filename = dir([cd '/Utilities/Temp/' char(selectedNodes(nleafs).getName) ' ' funname Figtag num2str(index) '.fig']);
                            if ~isempty(filename)
                                while isempty(dir([cd '/Utilities/Temp/' char(selectedNodes(nleafs).getName) ' ' funname Figtag num2str(index) '.fig'])) == 0
                                    index = index + 1;
                                end
                            end
                            filename = [char(selectedNodes(nleafs).getName) ' ' funname Figtag num2str(index) '.fig'];
                            saveas(h(i),[cd '/Utilities/Temp/' filename]);
                            fprintf(handles.fid,['%% Result figure name: ', filename, ', temporarily stored in /Utilities/Temp/\n']);
                        catch
                        end
                    end
                end
                fprintf(handles.fid,'\n');
                results = dir([cd '/Utilities/Temp']);
                % could sort the results here in terms of the parts of the
                % filenames (leafname, funname, T|C, index)
                %                 D = dir('directory_name') returns the results in an M-by-1
                %                   structure with the fields:
                %                       name    -- Filename
                %                       date    -- Modification date
                %                       bytes   -- Number of bytes allocated to the file
                %                       isdir   -- 1 if name is a directory and 0 if not
                %                       datenum -- Modification date as a MATLAB serial date number.
                %                       This value is locale-dependent.
                
                
                set(handles.result_box,'String',[' ';cellstr({results(3:length(results)).name}')]);
                if length(selectedNodes) > 1 || size(file,1) > 1, delete(h); end % close figures
            end
        end
        if errorflag
            h=warndlg('One or more attempts to run an analyser failed due to an error. Please refer to the error report in the Matlab Command Window.','AARAE info','modal');
            uiwait(h)
        end
        set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
        set(hObject,'Enable','on');
    end
end
set(handles.CloseFiguresButton,'Enable','on');
java.lang.Runtime.getRuntime.gc % Java garbage collection
%handles.mytree.setSelectedNode(handles.(matlab.lang.makeValidName(newleaf)));
guidata(hObject,handles)



% *************************************************************************
% SELECTION CHANGE IN FUN ANALYSERS BOX
% *************************************************************************

% --- Executes on selection change in fun_box.
function fun_box_Callback(hObject, ~, handles) %#ok
% hObject    handle to fun_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fun_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fun_box

% Allows 'help' display when mouse is hoovered over fun_box
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
[~,funname] = fileparts(selection);
if ~strcmp(selection,' ')
    handles.funname = funname;
    helptext = evalc(['help ' funname]);
    set(hObject,'Tooltip',helptext);
    set([handles.analyze_btn,handles.analyser_help_btn],'Visible','on');
    set(handles.analyze_btn,'BackgroundColor',[0.94 0.94 0.94]);
    set([handles.analyze_btn,handles.analyser_help_btn],'Enable','on');
else
    handles.funname = [];
    set(hObject,'Tooltip','');
    set([handles.analyze_btn,handles.analyser_help_btn],'Visible','off');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function fun_box_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to fun_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% *************************************************************************
% ANALYSERS HELP BUTTON
% *************************************************************************
% --- Executes on button press in analyser_help_btn.
function analyser_help_btn_Callback(~, ~, handles) %#ok : Executed when help button for analysers is clicked
% hObject    handle to analyser_help_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.fun_box,'String'));
selection = contents{get(handles.fun_box,'Value')};
eval(['doc ' selection])















% *************************************************************************
% *************************************************************************
%                           MAIN PLOTS IN GUI
% *************************************************************************
% *************************************************************************




% *************************************************************************
% SELECT CHANNEL TO DISPLAY
% *************************************************************************


function IN_nchannel_Callback(hObject, ~, handles) %#ok
% hObject    handle to IN_nchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IN_nchannel as text
%        str2double(get(hObject,'String')) returns contents of IN_nchannel as a double
hMain = getappdata(0,'hMain');
signaldata = getappdata(hMain,'testsignal'); % Get leaf contents from the 'desktop'
channel = str2double(get(handles.IN_nchannel,'String'));

if (channel <= size(signaldata.audio,2)) && (channel > 0) && ~isnan(channel)
    refreshplots(handles,'time')
    refreshplots(handles,'freq')
else
    warndlg('Invalid channel');
    set(handles.IN_nchannel,'String',num2str(handles.channel));
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function IN_nchannel_CreateFcn(hObject, ~, handles) %#ok
% hObject    handle to IN_nchannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.channel = 1;
guidata(hObject,handles);







% *************************************************************************
% CLICK ON A GUI PLOT TO CREATE A NEW FIGURE
% *************************************************************************

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function aarae_WindowButtonDownFcn(hObject, ~, handles) %#ok
% hObject    handle to aarae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
click = get(hObject,'CurrentObject');
if ~isempty(click) && ((click == handles.axestime) || (get(click,'Parent') == handles.axestime))
    hMain = getappdata(0,'hMain');
    signaldata = getappdata(hMain,'testsignal');
    if ~isempty(signaldata)
        if handles.alternate == 1 && isfield(signaldata,'audio2')
            data = signaldata.audio2;
        else
            data = signaldata.audio;
        end
        To = floor(str2double(get(handles.To_time,'String'))*signaldata.fs)+1;
        Tf = floor(str2double(get(handles.Tf_time,'String'))*signaldata.fs);
        if Tf > length(signaldata.audio), Tf = length(data); end
        if ~ismatrix(data)
            channel = str2double(get(handles.IN_nchannel,'String'));
            if handles.alternate == 0
                linea(:,:) = data(To:Tf,channel,:);
            else
                linea(:,:) = data(:,channel,:);
            end
            if ndims(data) == 3, cmap = colormap(hsv(size(data,3))); end
            if ndims(data) == 4, cmap = colormap(copper(size(data,4))); end
            if ndims(data) >= 5, cmap = colormap(cool(size(data,5))); end
        else
            if handles.alternate == 0
                linea = data(To:Tf,:);
            else
                linea = data;
            end
            cmap = colormap(lines(size(data,2)));
        end
        units = '';
        units_ref = 1;
        units_type = 1;
        if isfield(signaldata,'cal') && handles.Settings.calibrationtoggle == 1
            if size(linea,2) == length(signaldata.cal)
                if isfield(signaldata,'properties')
                    if isfield(signaldata.properties,'units')
                        units = signaldata.properties.units;
                    else
                        units = '';
                    end
                    if isfield(signaldata.properties,'units_ref')
                        units_ref = signaldata.properties.units_ref;
                    else
                        units_ref = 1;
                    end
                    if isfield(signaldata.properties,'units_type')
                        units_type = signaldata.properties.units_type;
                    else
                        units_type = 1;
                    end
                else
                    units = '';
                    units_ref = 1;
                    units_type = 1;
                end
                if units_type == 1
                    if numel(units_ref) == 1
                        linea = linea * units_ref;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    end
                else
                    if numel(units_ref) == 1
                        linea = linea * units_ref.^0.5;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref.^0.5,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./10);
                    end
                end
                signaldata.cal(isnan(signaldata.cal)) = 0;
                linea = linea.*repmat(10.^(signaldata.cal./20),length(linea),1);
            elseif ~ismatrix(signaldata.audio) && size(signaldata.audio,2) == length(signaldata.cal)
                if isfield(signaldata,'properties')
                    if isfield(signaldata.properties,'units')
                        units = signaldata.properties.units;
                    else
                        units = '';
                    end
                    if isfield(signaldata.properties,'units_ref')
                        units_ref = signaldata.properties.units_ref;
                    else
                        units_ref = 1;
                    end
                    if isfield(signaldata.properties,'units_type')
                        units_type = signaldata.properties.units_type;
                    else
                        units_type = 1;
                    end
                else
                    units = '';
                    units_ref = 1;
                    units_type = 1;
                end
                if units_type == 1
                    linea = linea .* units_ref;
                    signaldata.cal = signaldata.cal ./ 10.^(units_ref/20);if numel(units_ref) == 1
                        linea = linea * units_ref;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    end
                else
                    if numel(units_ref) == 1
                        linea = linea * units_ref.^0.5;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref.^0.5,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./10);
                    end
                end
                signaldata.cal(isnan(signaldata.cal)) = 0;
                cal = repmat(signaldata.cal(str2double(get(handles.IN_nchannel,'String'))),1,size(linea,2));
                linea = linea.*repmat(10.^(cal./20),length(linea),1);
            end
        else
            units = '';
            units_ref = 1;
            units_type = 1;
        end
        t = linspace(To,Tf,length(linea))./signaldata.fs;
        f = signaldata.fs .* ((1:length(linea))-1) ./ length(linea);
        switch handles.Settings.specmagscale;
            case {'Divided by length'}
                spectscale = 1./length(linea);
            case {'x sqrt2/length'}
                spectscale = 2.^0.5./length(linea);
            case {'x 2/length'}
                spectscale = 2./length(linea);
            otherwise
                spectscale = 1;
        end
        h = figure;
        set(h,'DefaultAxesColorOrder',cmap);
        plottype = get(handles.time_popup,'Value');
        if plottype == 1
            linea = real(linea);
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 2
            linea = linea.^2;
            if ~isempty(units)
                if units_type == 2
                    units_string = [' [' units ']'];
                else
                    
                    units_string = [' [(' units ')^2]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 3
            if units_type == 1
                if numel(units_ref) == 1
                    linea = 10.*log10(linea.^2 ./  units_ref.^2);
                else
                    linea = 10.*log10(linea.^2 ./...
                        repmat(units_ref.^2,size(linea,1),1));
                end
            else
                if numel(units_ref) == 1
                    linea = 10.*log10(linea.^2 ./  units_ref);
                else
                    linea = 10.*log10(linea.^2 ./...
                        repmat(units_ref,size(linea,1),1));
                end
            end
            units_string = '';
        end
        if plottype == 4
            linea = abs(hilbert(real(linea)));
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 5
            linea = medfilt1(diff([angle(hilbert(real(linea))); zeros(1,size(linea,2))])*signaldata.fs/2/pi, 5);
            units_string = '';
        end
        if plottype == 6
            linea = abs(linea);
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 7
            linea = imag(linea);
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 8
            if units_type == 1
                if numel(units_ref) == 1
                    linea = 10*log10(abs(fft(linea).*spectscale ./...
                        units_ref).^2);
                else
                    linea = 10*log10(abs(fft(linea).*spectscale ./...
                        repmat(units_ref,size(linea,1),1)).^2);
                end
            else
                if numel(units_ref) == 1
                    linea = 10*log10(abs(fft(linea).*spectscale ./...
                        units_ref.^0.5).^2);
                else
                    linea = 10*log10(abs(fft(linea).*spectscale ./...
                        repmat(units_ref.^0.5,size(linea,1),1)).^2);
                end
            end
            units_string = '';
        end %freq
        if plottype == 9
            linea = (abs(fft(linea)).*spectscale).^2;
            if ~isempty(units)
                if units_type == 2
                    units_string = [' [' units ']'];
                else
                    
                    units_string = [' [(' units ')^2]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 10
            linea = abs(fft(linea)).*spectscale;
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 11
            linea = real(fft(linea)).*spectscale;
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 12
            linea = imag(fft(linea)).*spectscale;
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 13
            linea = angle(fft(linea));
            units_string = '';
        end
        if plottype == 14
            linea = unwrap(angle(fft(linea)));
            units_string = '';
        end
        if plottype == 15
            linea = angle(fft(linea)) .* 180/pi;
            units_string = '';
        end
        if plottype == 16
            linea = unwrap(angle(fft(linea))) ./(2*pi);
            units_string = '';
        end
        if plottype == 17
            linea = -diff(unwrap(angle(fft(linea)))).*length(linea)/(signaldata.fs*2*pi).*1000;
            units_string = '';
        end
        if plottype <= 7
            plot(t,real(linea)) % Plot signal in time domain
            xlabel('Time [s]');
            yl = cellstr(get(handles.time_popup,'String'));
            yl = yl{get(handles.time_popup,'Value')};
            ylabel([yl(8:end) units_string]);
        end
        if plottype >= 8
            if plottype == 8
                smoothfactor = get(handles.smoothtime_popup,'Value');
                if smoothfactor == 2, octsmooth = 1; end
                if smoothfactor == 3, octsmooth = 3; end
                if smoothfactor == 4, octsmooth = 6; end
                if smoothfactor == 5, octsmooth = 12; end
                if smoothfactor == 6, octsmooth = 24; end
                if smoothfactor ~= 1, linea = octavesmoothing(linea, octsmooth, signaldata.fs); end
            end
            %if plottype == 17,
            plot(f(1:length(linea)),linea)%,'Marker','None'); end
            %if plottype ~= 17, semilogx(f,linea); end % Plot signal in frequency domain
            xlabel('Frequency [Hz]');
            yl = cellstr(get(handles.time_popup,'String'));
            yl = yl{get(handles.time_popup,'Value')};
            ylabel([yl(9:end) units_string]);
            if ischar(handles.Settings.frequencylimits)
                xlim([f(2) signaldata.fs/2])
            else
                xlim(handles.Settings.frequencylimits)
            end
            log_check = get(handles.logtime_chk,'Value');
            if log_check == 1
                set(gca,'XScale','log')
                %set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
            else
                set(gca,'XScale','linear')%,'XTickLabelMode','auto')
                %set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
            end
        end
        if handles.alternate == 0
            if ~ismatrix(data)
                if isfield(signaldata,'bandID')
                    legend(num2str(signaldata.bandID'));
                end
            else
                if isfield(signaldata,'chanID')
                    legend(signaldata.chanID);
                end
            end
        end
        handles.alternate = 0;
    end
    selectedNodes = handles.mytree.getSelectedNodes;
    selectedNodes = selectedNodes(1);
    title(strrep(selectedNodes.getName.char,'_',' '))
end
units = '';
units_ref = 1;
units_type = 1;
if ~isempty(click) && ((click == handles.axesfreq) || (get(click,'Parent') == handles.axesfreq))
    hMain = getappdata(0,'hMain');
    signaldata = getappdata(hMain,'testsignal');
    if ~isempty(signaldata)
        if handles.alternate == 1 && isfield(signaldata,'audio2')
            data = signaldata.audio2;
        else
            data = signaldata.audio;
        end
        To = floor(str2double(get(handles.To_freq,'String'))*signaldata.fs)+1;
        Tf = floor(str2double(get(handles.Tf_freq,'String'))*signaldata.fs);
        if Tf > length(signaldata.audio), Tf = length(data); end
        if ~ismatrix(data)
            channel = str2double(get(handles.IN_nchannel,'String'));
            if handles.alternate == 0
                linea(:,:) = data(To:Tf,channel,:);
            else
                linea(:,:) = data(:,channel,:);
            end
            if ndims(data) == 3, cmap = colormap(hsv(size(data,3))); end
            if ndims(data) == 4, cmap = colormap(copper(size(data,4))); end
            if ndims(data) >= 5, cmap = colormap(cool(size(data,5))); end
        else
            if handles.alternate == 0
                linea = data(To:Tf,:);
            else
                linea = data;
            end
            cmap = colormap(lines(size(data,2)));
        end
        if isfield(signaldata,'cal') && handles.Settings.calibrationtoggle == 1
            if size(linea,2) == length(signaldata.cal)
                if isfield(signaldata,'properties')
                    if isfield(signaldata.properties,'units')
                        units = signaldata.properties.units;
                    else
                        units = '';
                    end
                    if isfield(signaldata.properties,'units_ref')
                        units_ref = signaldata.properties.units_ref;
                    else
                        units_ref = 1;
                    end
                    if isfield(signaldata.properties,'units_type')
                        units_type = signaldata.properties.units_type;
                    else
                        units_type = 1;
                    end
                else
                    units = '';
                    units_ref = 1;
                    units_type = 1;
                end
                signaldata.cal(isnan(signaldata.cal)) = 0;
                if units_type == 1
                    if numel(units_ref) == 1
                        linea = linea * units_ref;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    end
                else
                    if numel(units_ref) == 1
                        linea = linea * units_ref.^0.5;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref.^0.5,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./10);
                    end
                end
                
                linea = linea.*repmat(10.^(signaldata.cal./20),length(linea),1);
            elseif ~ismatrix(signaldata.audio) && size(signaldata.audio,2) == length(signaldata.cal)
                if isfield(signaldata,'properties')
                    if isfield(signaldata.properties,'units')
                        units = signaldata.properties.units;
                    else
                        units = '';
                    end
                    if isfield(signaldata.properties,'units_ref')
                        units_ref = signaldata.properties.units_ref;
                    else
                        units_ref = 1;
                    end
                    if isfield(signaldata.properties,'units_type')
                        units_type = signaldata.properties.units_type;
                    else
                        units_type = 1;
                    end
                else
                    units = '';
                    units_ref = 1;
                    units_type = 1;
                end
                if units_type == 1
                    if numel(units_ref) == 1
                        linea = linea * units_ref;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./20);
                    end
                else
                    if numel(units_ref) == 1
                        linea = linea * units_ref.^0.5;
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref/10);
                    elseif length(units_ref) == length(signaldata.cal)
                        linea = linea .* repmat(units_ref.^0.5,size(linea,1),1);
                        signaldata.cal = signaldata.cal ./ 10.^(units_ref./10);
                    end
                end
                signaldata.cal(isnan(signaldata.cal)) = 0;
                cal = repmat(signaldata.cal(str2double(get(handles.IN_nchannel,'String'))),1,size(linea,2));
                linea = linea.*repmat(10.^(cal./20),length(linea),1);
            end
        else
            units = '';
            units_ref = 1;
            units_type = 1;
        end
        t = linspace(To,Tf,length(linea))./signaldata.fs;
        f = signaldata.fs .* ((1:length(linea))-1) ./ length(linea);
        switch handles.Settings.specmagscale;
            case {'Divided by length'}
                spectscale = 1./length(linea);
            case {'x sqrt2/length'}
                spectscale = 2.^0.5./length(linea);
            case {'x 2/length'}
                spectscale = 2./length(linea);
            otherwise
                spectscale = 1;
        end
        h = figure;
        set(h,'DefaultAxesColorOrder',cmap);
        plottype = get(handles.freq_popup,'Value');
        if plottype == 1
            linea = real(linea);
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 2
            linea = linea.^2;
            if ~isempty(units)
                if units_type == 2
                    units_string = [' [' units ']'];
                else
                    
                    units_string = [' [(' units ')^2]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 3
            if units_type == 1
                if numel(units_ref) == 1
                    linea = 10.*log10(linea.^2 ./  units_ref.^2);
                else
                    linea = 10.*log10(linea.^2 ./...
                        repmat(units_ref.^2,size(linea,1),1));
                end
            else
                if numel(units_ref) == 1
                    linea = 10.*log10(linea.^2 ./  units_ref);
                else
                    linea = 10.*log10(linea.^2 ./...
                        repmat(units_ref,size(linea,1),1));
                end
            end
            units_string = '';
        end
        if plottype == 4
            linea = abs(hilbert(real(linea)));
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 5
            linea = medfilt1(diff([angle(hilbert(real(linea))); zeros(1,size(linea,2))])*signaldata.fs/2/pi, 5);
            units_string = '';
        end
        if plottype == 6
            linea = abs(linea);
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 7
            linea = imag(linea);
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 8
            if units_type == 1
                if numel(units_ref) == 1
                    linea = 10*log10(abs(fft(linea).*spectscale ./...
                        units_ref).^2);
                else
                    linea = 10*log10(abs(fft(linea).*spectscale ./ ...
                        repmat(units_ref,size(linea,1),1)).^2);
                end
            else
                if numel(units_ref) == 1
                    linea = 10*log10(abs(fft(linea).*spectscale ./...
                        units_ref.^0.5).^2);
                else
                    linea = 10*log10(abs(fft(linea).*spectscale ./ ...
                        repmat(units_ref.^0.5,size(linea,1),1)).^2);
                end
            end
            units_string = '';
        end %freq
        if plottype == 9
            linea = (abs(fft(linea)).*spectscale).^2;
            if ~isempty(units)
                if units_type == 2
                    units_string = [' [' units ']'];
                else
                    
                    units_string = [' [(' units ')^2]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 10
            linea = abs(fft(linea)).*spectscale;
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 11
            linea = real(fft(linea)).*spectscale;
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 12
            linea = imag(fft(linea)).*spectscale;
            if ~isempty(units)
                if units_type == 1
                    units_string = [' [' units ']'];
                else
                    units_string = [' [(' units ')^0.5]'];
                end
            else
                units_string = '';
            end
        end
        if plottype == 13
            linea = angle(fft(linea));
            units_string = '';
        end
        if plottype == 14
            linea = unwrap(angle(fft(linea)));
            units_string = '';
        end
        if plottype == 15
            linea = angle(fft(linea)) .* 180/pi;
            units_string = '';
        end
        if plottype == 16
            linea = unwrap(angle(fft(linea))) ./(2*pi);
            units_string = '';
        end
        if plottype == 17
            linea = -diff(unwrap(angle(fft(linea)))).*length(linea)/(signaldata.fs*2*pi).*1000;
            units_string = '';
        end
        if plottype <= 7
            plot(t,real(linea)) % Plot signal in time domain
            xlabel('Time [s]');
            yl = cellstr(get(handles.freq_popup,'String'));
            yl = yl{get(handles.freq_popup,'Value')};
            ylabel([yl(8:end) units_string]);
        end
        if plottype >= 8
            if plottype == 8
                smoothfactor = get(handles.smoothfreq_popup,'Value');
                if smoothfactor == 2, octsmooth = 1; end
                if smoothfactor == 3, octsmooth = 3; end
                if smoothfactor == 4, octsmooth = 6; end
                if smoothfactor == 5, octsmooth = 12; end
                if smoothfactor == 6, octsmooth = 24; end
                if smoothfactor ~= 1, linea = octavesmoothing(linea, octsmooth, signaldata.fs); end
            end
            %if plottype == 17,
            plot(f(1:length(linea)),linea)%,'Marker','None'); end
            %if plottype ~= 17, semilogx(f,linea); end % Plot signal in frequency domain
            xlabel('Frequency [Hz]');
            yl = cellstr(get(handles.freq_popup,'String'));
            yl = yl{get(handles.freq_popup,'Value')};
            ylabel([yl(9:end) units_string]);
            if ischar(handles.Settings.frequencylimits)
                xlim([f(2) signaldata.fs/2])
            else
                xlim(handles.Settings.frequencylimits)
            end
            log_check = get(handles.logfreq_chk,'Value');
            if log_check == 1
                set(gca,'XScale','log')
                %set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
            else
                set(gca,'XScale','linear')%,'XTickLabelMode','auto')
                %set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
            end
        end
        if handles.alternate == 0
            if ~ismatrix(data)
                if isfield(signaldata,'bandID')
                    legend(num2str(signaldata.bandID'));
                end
            else
                if isfield(signaldata,'chanID')
                    legend(signaldata.chanID);
                end
            end
        end
        handles.alternate = 0;
    end
    selectedNodes = handles.mytree.getSelectedNodes;
    selectedNodes = selectedNodes(1);
    title(strrep(selectedNodes.getName.char,'_',' '))
end
if ~isempty(click) && ((click == handles.axesdata) || isequal(click,findobj('Tag','tempaxes','Parent',hObject)) || (get(click,'Parent') == handles.axesdata) || isequal(get(click,'Parent'),findobj('Tag','tempaxes','Parent',hObject))) && ~strcmp(get(click,'Type'),'text')
    selectedNodes = handles.mytree.getSelectedNodes;
    data = selectedNodes(1).handle.UserData;
    if ~isfield(data,'tables')
        figure;
        haxes = axes;
        doresultplot(handles,haxes)
    else
        figure;
        ntable = get(handles.ntable_popup,'Value');
        Xvalues = get(handles.Xvalues_sel,'SelectedObject');
        Xvalues = get(Xvalues,'tag');
        switch Xvalues
            case 'radiobutton1'
                bar(data.tables(ntable).Data(:,get(handles.Yvalues_box,'Value')),'FaceColor',[0 0.5 0.5])
                set(gca,'Xtick',1:length(data.tables(ntable).RowName),'XTickLabel',data.tables(ntable).RowName)
            case 'radiobutton2'
                bar(data.tables(ntable).Data(get(handles.Yvalues_box,'Value'),:),'FaceColor',[0 0.5 0.5])
                set(gca,'Xtick',1:length(data.tables(ntable).ColumnName),'XTickLabel',data.tables(ntable).ColumnName)
        end
        ycontents = cellstr(get(handles.Yvalues_box,'String'));
        ylabel(gca,ycontents{get(handles.Yvalues_box,'Value')})
    end
end
if ~isempty(click) && strcmp(get(click,'Type'),'text')
    if click == get(get(click,'Parent'),'Xlabel')
        if strcmp(get(get(click,'Parent'),'XScale'),'linear')
            set(get(click,'Parent'),'XScale','log')
        else
            set(get(click,'Parent'),'XScale','linear')
        end
    end
    if click == get(get(click,'Parent'),'Ylabel')
        if strcmp(get(get(click,'Parent'),'YScale'),'linear')
            set(get(click,'Parent'),'YScale','log')
        else
            set(get(click,'Parent'),'YScale','linear')
        end
    end
    if click == get(get(click,'Parent'),'Zlabel')
        if strcmp(get(get(click,'Parent'),'ZScale'),'linear')
            set(get(click,'Parent'),'ZScale','log')
        else
            set(get(click,'Parent'),'ZScale','linear')
        end
    end
end
handles.alternate = 0;
guidata(hObject,handles)



% *************************************************************************
% PLOT CONTROLS
% *************************************************************************

% --- Executes on selection change in smoothfreq_popup.
function smoothfreq_popup_Callback(hObject, ~, handles) %#ok
% hObject    handle to smoothfreq_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns smoothfreq_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from smoothfreq_popup
refreshplots(handles,'freq')
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function smoothfreq_popup_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to smoothfreq_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in smoothtime_popup.
function smoothtime_popup_Callback(hObject, ~, handles) %#ok
% hObject    handle to smoothtime_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns smoothtime_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from smoothtime_popup
refreshplots(handles,'time')
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function smoothtime_popup_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to smoothtime_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in logfreq_chk.
function logfreq_chk_Callback(hObject, ~, handles) %#ok
% hObject    handle to logfreq_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logfreq_chk
refreshplots(handles,'freq')
guidata(hObject,handles)

% --- Executes on button press in logtime_chk.
function logtime_chk_Callback(hObject, ~, handles) %#ok
% hObject    handle to logtime_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logtime_chk
refreshplots(handles,'time')
guidata(hObject,handles)




function To_freq_Callback(hObject, ~, handles) %#ok : Executed when initial time input box changes above the lower axes
% hObject    handle to To_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of To_freq as text
%        str2double(get(hObject,'String')) returns contents of To_freq as a double
To_freq = str2double(get(hObject,'String'));
selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;
if isnan(To_freq) || To_freq < 0 || To_freq == length(signaldata.audio)/signaldata.fs
    %warndlg('Invalid entry','AARAE info','modal')
    set(hObject,'String',num2str(handles.To_freq_IN))
elseif To_freq >= str2double(get(handles.Tf_freq,'String'))
    Tf_freq = To_freq + (handles.Tf_freq_IN - handles.To_freq_IN);
    if Tf_freq > length(signaldata.audio)/signaldata.fs
        Tf_freq = length(signaldata.audio)/signaldata.fs;
    end
    set(handles.Tf_freq,'String',num2str(Tf_freq))
    handles.To_freq_IN = To_freq;
    handles.Tf_freq_IN = Tf_freq;
    refreshplots(handles,'freq')
    guidata(hObject,handles)
else
    handles.To_freq_IN = To_freq;
    refreshplots(handles,'freq')
    guidata(hObject,handles)
end



% --- Executes during object creation, after setting all properties.
function To_freq_CreateFcn(hObject, ~, ~) %#ok : creation of initial time input box for the lower axes
% hObject    handle to To_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tf_freq_Callback(hObject, ~, handles) %#ok : Executed when final time input box for lower axes changes
% hObject    handle to Tf_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tf_freq as text
%        str2double(get(hObject,'String')) returns contents of Tf_freq as a double
Tf_freq = str2double(get(hObject,'String'));
selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;
if isnan(Tf_freq) || Tf_freq <= str2double(get(handles.To_freq,'String')) || Tf_freq > length(signaldata.audio)/signaldata.fs
    %warndlg('Invalid entry','AARAE info','modal')
    set(hObject,'String',handles.Tf_freq_IN)
    refreshplots(handles,'freq')
    guidata(hObject,handles)
else
    handles.Tf_freq_IN = Tf_freq;
    refreshplots(handles,'freq')
    guidata(hObject,handles)
end


% --- Executes during object creation, after setting all properties.
function Tf_freq_CreateFcn(hObject, ~, ~) %#ok : creation of final time input box for the lower axes
% hObject    handle to Tf_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function To_time_Callback(hObject, ~, handles) %#ok : Executed when initial time input box changes above upper axes
% hObject    handle to To_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of To_time as text
%        str2double(get(hObject,'String')) returns contents of To_time as a double
To_time = str2double(get(hObject,'String'));
selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;
if isnan(To_time) || To_time < 0 || To_time == length(signaldata.audio)/signaldata.fs
    set(hObject,'String',num2str(handles.To_time_IN))
elseif To_time >= str2double(get(handles.Tf_time,'String'))
    Tf_time = To_time + (handles.Tf_time_IN - handles.To_time_IN);
    if Tf_time > length(signaldata.audio)/signaldata.fs
        Tf_time = length(signaldata.audio)/signaldata.fs;
    end
    set(handles.Tf_time,'String',num2str(Tf_time))
    handles.To_time_IN = To_time;
    handles.Tf_time_IN = Tf_time;
    refreshplots(handles,'time')
    guidata(hObject,handles)
else
    handles.To_time_IN = To_time;
    refreshplots(handles,'time')
    guidata(hObject,handles)
end


% --- Executes during object creation, after setting all properties.
function To_time_CreateFcn(hObject, ~, ~) %#ok : creation of initial time input box for the upper axes
% hObject    handle to To_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tf_time_Callback(hObject, ~, handles) %#ok : Executed when final time input box for upper axes changes
% hObject    handle to Tf_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tf_time as text
%        str2double(get(hObject,'String')) returns contents of Tf_time as a double
Tf_time = str2double(get(hObject,'String'));
selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;
if isnan(Tf_time) || Tf_time <= str2double(get(handles.To_time,'String')) || Tf_time > length(signaldata.audio)/signaldata.fs
    %warndlg('Invalid entry','AARAE info','modal')
    set(hObject,'String',handles.Tf_time_IN)
    refreshplots(handles,'time')
    guidata(hObject,handles)
else
    handles.Tf_time_IN = Tf_time;
    refreshplots(handles,'time')
    guidata(hObject,handles)
end


% --- Executes during object creation, after setting all properties.
function Tf_time_CreateFcn(hObject, ~, ~) %#ok : creation of final time input box for the upper axes
% hObject    handle to Tf_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in freq_popup.
function freq_popup_Callback(hObject, ~, handles) %#ok
% hObject    handle to freq_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns freq_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from freq_popup
refreshplots(handles,'freq')
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function freq_popup_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to freq_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in time_popup.
function time_popup_Callback(hObject, ~, handles) %#ok
% hObject    handle to time_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns time_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_popup
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
set(handles.compare_btn,'TooltipString',['Compare selected signals in ' selection])
refreshplots(handles,'time')
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function time_popup_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to time_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function time_units_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_units_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function freq_units_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_units_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
guidata(hObject,handles)






% *************************************************************************
% *************************************************************************
% RESULT FIGURES BOX FUNCTIONS
% *************************************************************************
% *************************************************************************


% --- Executes on selection change in result_box.
function result_box_Callback(hObject, ~, handles) %#ok
% hObject    handle to result_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns result_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from result_box
get(handles.aarae,'SelectionType');
contents = cellstr(get(hObject,'String'));
if strcmp(get(handles.aarae,'SelectionType'),'open') && ~isempty(contents{get(hObject,'Value')})
    contents = cellstr(get(hObject,'String'));
    file = contents{get(hObject,'Value')};
    [~,~,ext] = fileparts(file);
    if ~strcmp(file,' ')
        switch ext
            case '.fig'
                openfig(file);
            otherwise
                open(file)
        end
    end
end

% --- Executes during object creation, after setting all properties.
function result_box_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to result_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





















% *************************************************************************
% *************************************************************************
%                           PROPERTIES BUTTON
% *************************************************************************
% *************************************************************************

% --- Executes on button press in properties_btn.
function properties_btn_Callback(~, ~, handles) %#ok
% hObject    handle to properties_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedNodes = handles.mytree.getSelectedNodes;
signaldata = selectedNodes(1).handle.UserData;
if isfield(signaldata,'properties')
    properties = signaldata.properties; %#ok : used for evalc
    msgbox([selectedNodes(1).getName.char evalc('properties')],'AARAE info')
end





% *************************************************************************
% *************************************************************************
%                           AARAE TABLES IN GUI
% *************************************************************************
% *************************************************************************



% --- Executes on selection change in ntable_popup.
function ntable_popup_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to ntable_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ntable_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ntable_popup
eventdata.NewValue = get(handles.Xvalues_sel,'SelectedObject');
Xvalues_sel_SelectionChangeFcn(hObject, eventdata, handles)
Yvalues_box_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ntable_popup_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to ntable_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Xvalues_box.
function Xvalues_box_Callback(~, ~, ~) %#ok
% hObject    handle to Xvalues_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Xvalues_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Xvalues_box


% --- Executes during object creation, after setting all properties.
function Xvalues_box_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to Xvalues_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Yvalues_box.
function Yvalues_box_Callback(~, ~, handles)
% hObject    handle to Yvalues_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Yvalues_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Yvalues_box
selectedNodes = handles.mytree.getSelectedNodes;
data = selectedNodes(1).handle.UserData;
ntable = get(handles.ntable_popup,'Value');
Xvalues = get(handles.Xvalues_sel,'SelectedObject');
Xvalues = get(Xvalues,'tag');
switch Xvalues
    case 'radiobutton1'
        bar(handles.axesdata,data.tables(ntable).Data(:,get(handles.Yvalues_box,'Value')),'FaceColor',[0 0.5 0.5])
        set(handles.axesdata,'Xtick',1:length(data.tables(ntable).RowName),'XTickLabel',data.tables(ntable).RowName)
    case 'radiobutton2'
        bar(handles.axesdata,data.tables(ntable).Data(get(handles.Yvalues_box,'Value'),:),'FaceColor',[0 0.5 0.5])
        set(handles.axesdata,'Xtick',1:length(data.tables(ntable).ColumnName),'XTickLabel',data.tables(ntable).ColumnName)
end
ycontents = cellstr(get(handles.Yvalues_box,'String'));
ylabel(handles.axesdata,ycontents{get(handles.Yvalues_box,'Value')})

% --- Executes during object creation, after setting all properties.
function Yvalues_box_CreateFcn(hObject, ~, ~) %#ok
% hObject    handle to Yvalues_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in Xvalues_sel.
function Xvalues_sel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Xvalues_sel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
selectedNodes = handles.mytree.getSelectedNodes;
data = selectedNodes(1).handle.UserData;
ntable = get(handles.ntable_popup,'Value');
switch get(eventdata.NewValue,'Tag')
    case 'radiobutton1'
        set(handles.Xvalues_box,'String',data.tables(ntable).RowName,'Value',1)
        set(handles.Yvalues_box,'String',data.tables(ntable).ColumnName,'Value',1)
    case 'radiobutton2'
        set(handles.Yvalues_box,'String',data.tables(ntable).RowName,'Value',1)
        set(handles.Xvalues_box,'String',data.tables(ntable).ColumnName,'Value',1)
end
Yvalues_box_Callback(hObject, eventdata, handles)





% --- Executes on button press in wild_btn.
%function wild_btn_Callback(hObject, eventdata, handles)
% hObject    handle to wild_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%!matlab -nodesktop








% *************************************************************************
% *************************************************************************
%                           RESULT PLOTS
% *************************************************************************
% *************************************************************************


% --- Executes on selection change in chartfunc_popup.
function chartfunc_popup_Callback(~, eventdata, handles) %#ok : Executed when selection changes in chart selection popup menu
% hObject    handle to chartfunc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chartfunc_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chartfunc_popup
doresultplot(handles,handles.axesdata)

% --- Executes during object creation, after setting all properties.
function chartfunc_popup_CreateFcn(hObject, ~, ~) %#ok : creation of chart selection type popup menu
% hObject    handle to chartfunc_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes when selected cell(s) is changed in cattable.
function cattable_CellSelectionCallback(hObject, eventdata, handles) %#ok : opens listdlg for changing selection of categorical dimensions
% hObject    handle to cattable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
tabledata = get(hObject,'Data');
if size(eventdata.Indices,1) ~= 0 && eventdata.Indices(1,2) == 2
    chkbox = tabledata{eventdata.Indices(1,1),4};
    if isempty(chkbox), chkbox = false; end
    if chkbox == false
        selectedNodes = handles.mytree.getSelectedNodes;
        audiodata = selectedNodes(1).handle.UserData;
        handles.tabledata = tabledata;
        catname = tabledata{eventdata.Indices(1,1),1};
        liststr = audiodata.(matlab.lang.makeValidName(catname));
        if size(liststr,1) < size(liststr,2), liststr = liststr'; end
        if ~iscellstr(liststr) && ~isnumeric(liststr), liststr = cellstr(num2str(cell2mat(liststr)));
        elseif isnumeric(liststr), liststr = cellstr(num2str(liststr)); end
        [sel,ok] = listdlg('ListString',liststr,'InitialValue',str2num(tabledata{eventdata.Indices(1),eventdata.Indices(2)})); %#ok : necessary for getting selection vector
        if ok == 1
            logsel = ['[' num2str(sel) ']'];
            tabledata{eventdata.Indices(1),eventdata.Indices(2)} = logsel;
        end
        set(hObject,'Data',{''})
        set(hObject,'Data',tabledata)
        guidata(handles.aarae,handles)
        doresultplot(handles,handles.axesdata)
    else
        % Possible code to truncate 'continuous' selection
    end
end





% *************************************************************************
% RESULT PLOT SETTINGS
% *************************************************************************

% --- Executes when entered data in editable cell(s) in cattable.
function cattable_CellEditCallback(hObject, eventdata, handles) %#ok : Executed when cell information on the uitable changes
% hObject    handle to cattable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if size(eventdata.Indices,1) ~= 0 && eventdata.Indices(1,2) == 4
    tabledata = get(handles.cattable,'Data');
    catorcont = tabledata(:,4);
    naxis = length(find([catorcont{:}] == true));
    if naxis < 3
        if islogical(catorcont{eventdata.Indices(1,1),1}) && catorcont{eventdata.Indices(1,1),1} == true
            tabledata{eventdata.Indices(1,1),2} = ':';
        else
            tabledata{eventdata.Indices(1,1),2} = '[1]';
        end
        set(handles.cattable,'Data',tabledata);
        switch naxis
            case 0
                set(handles.chartfunc_popup,'String',{'distributionPlot','boxplot'},'Value',1)
            case 1
                set(handles.chartfunc_popup,'String',{'plot','semilogx','semilogy','loglog','distributionPlot','boxplot'},'Value',1)
            case 2
                set(handles.chartfunc_popup,'String',{'mesh','surf','imagesc'},'Value',1)
        end
        doresultplot(handles,handles.axesdata)
    else
        doresultplot(handles,handles.axesdata)
    end
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
