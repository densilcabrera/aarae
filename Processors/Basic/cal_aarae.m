function OUT = cal_aarae(IN,method,values,fileref)
% This function can be used to calibrate audio, and should be functionally
% equivalent to AARAE's CAL button


if ~exist('method','var')
    method = [];
    if isempty(method)
        method = menu('Calibration',...
            'Choose from AARAE',...
            'Locate file on disc',...
            'Input value',...
            'Specify Leq',...
            'Specify weighted Leq',...
            'Cancel');
    end
end
if ~exist('values','var')
    values = [];
end
if ~exist('fileref','var')
    fileref = [];
end

handles = guidata(findobj('Tag','aarae'));
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
            fileref = s;
            calfile = handles.(matlab.lang.makeValidName(leaves{s,1})).handle.UserData;
            if ~isfield(calfile,'audio')
                warndlg('Incompatible calibration file','Warning!');
                OUT = [];
                return
            else
                cal_level = 10 .* log10(mean(calfile.audio.^2,1));
                %calibration = 1./(10.^((cal_level)./20));
                %IN.audio = IN.audio.*calibration;
                if ~isempty(values) && length(values) > 1 && length(values) ~= size(IN.audio,2)
                        values = [];
                end
                if (size(IN.audio,2) == size(cal_level,2) || size(cal_level,2) == 1) && ismatrix(IN.audio)
                    if isempty(values)
                        cal_offset = inputdlg('Calibration tone RMS level',...
                        'Calibration value',[1 50],{'0'});
                    else
                        cal_offset = values;
                    end
                    if isnan(str2double(char(cal_offset)))
                        return
                    else
                        cal_level = str2double(char(cal_offset)) - cal_level;
                    end
                    if size(cal_level,2) == 1, cal_level = repmat(cal_level,1,size(IN.audio,2)); end
                else
                    warndlg('Incompatible calibration file','Warning!');
                    OUT = [];
                    return
                end
            end
        else
            warndlg('No signal loaded!','Whoops...!');
            OUT = [];
            return
        end
        values = cal_offset;
    case 2
        if isempty(fileref)
        [filename,handles.defaultaudiopath] = uigetfile(...
            {'*.wav;*.mat;.WAV;.MAT','Calibration file (*.wav,*.mat)'},...
            'Select audio file',handles.defaultaudiopath);
        else
            filename = fileref;
            if ~ischar(filename)
                [filename,handles.defaultaudiopath] = uigetfile(...
            {'*.wav;*.mat;.WAV;.MAT','Calibration file (*.wav,*.mat)'},...
            'Select audio file',handles.defaultaudiopath);
            end
        end
        if ~ischar(filename)
            OUT = [];
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
            if ~isempty(values) && length(values) > 1 && length(values) ~= size(IN.audio,2)
                 values = [];
            end
            if (size(IN.audio,2) == size(cal_level,2) || size(cal_level,2) == 1) && ismatrix(caltone)
                if isempty(values)
                cal_offset = inputdlg('Calibration tone RMS level',...
                    'Calibration value',[1 50],{'0'});
                else
                    cal_offset = values;
                end
                if isnan(str2double(char(cal_offset)))
                    return
                else
                    cal_level = str2double(char(cal_offset)) - cal_level;
                end
                if size(cal_level,2) == 1, cal_level = repmat(cal_level,1,size(IN.audio,2)); end
            else
                warndlg('Incompatible calibration file!','AARAE info');
                OUT = [];
                return
            end
        end
        values = cal_offset;
    case 3
        chans = size(IN.audio,2);
        if ~isempty(values) && (length(values) == 1 || length(values) == chans)
            cal_level = values;
        else
            if isfield(IN,'cal')
                def = cellstr(num2str(IN.cal'));
            else
                def = cellstr(num2str(zeros(chans,1)));
            end
            cal_level = inputdlg(cellstr([repmat('channel ',chans,1) num2str((1:chans)')]),...
                'Calibration value',[1 60],def);
            cal_level = str2num(char(cal_level))'; %#ok to prevent from spaces introduced in the input boxes
        end
        if size(cal_level,1) > size(cal_level,2), cal_level = cal_level'; end
        if isempty(cal_level) || chans ~= size(cal_level,2)
            warndlg('Calibration values mismatch!','AARAE info');
            OUT = [];
            return
        end
        values = cal_level;
    case 4
        cal_level = 10 .* log10(mean(IN.audio.^2,1));
        cal_level = repmat(20*log10(mean(10.^(cal_level./20),2)),1,size(IN.audio,2));
        if ~isempty(values) && (length(values) == 1 || length(values) == size(IN.audio,2))
            cal_offset = values;
        else
            cal_offset = inputdlg('Signal RMS level',...
                'Calibration value',[1 50],cellstr(num2str(zeros(size(cal_level)))));
            if isempty(cal_offset)
                OUT = [];
                return;
            else
                cal_offset = str2num(char(cal_offset)); %#ok : to allow spaces between calibration values
            end
        end
        if (isequal(size(cal_offset),size(cal_level)) || size(cal_offset,2) == 1) && ismatrix(IN.audio)
            cal_level = cal_offset - cal_level;
        else
            warndlg('Calibration values mismatch!','AARAE info');
            OUT = [];
            return
        end
        values = cal_offset;
    case 5
        weights = what([cd '/Processors/Filters']);
        if ~isempty(weights.m)
            [selection,ok] = listdlg('ListString',cellstr(weights.m),'SelectionMode','single');
        else
            warndlg('No weighting filters found!','AARAE info')
            OUT = [];
            return
        end
        if ok == 1
            [~,funname] = fileparts(weights.m{selection,1});
            INweight = feval(funname,IN);
            cal_level = 10 .* log10(mean(INweight.audio.^2,1));
            cal_level = repmat(20*log10(mean(10.^(cal_level./20),2)),1,size(IN.audio,2));
            cal_offset = inputdlg('Signal RMS level',...
                'Calibration value',[1 50],cellstr(num2str(zeros(size(cal_level)))));
            if isempty(cal_offset)
                OUT = [];
                return;
            else
                cal_offset = str2num(char(cal_offset)); %#ok : to allow spaces between calibration values
            end
            if (isequal(size(cal_offset),size(cal_level)) || size(cal_offset,2) == 1) && ismatrix(IN.audio)
                cal_level = cal_offset - cal_level;
            else
                warndlg('Calibration values mismatch!','AARAE info');
                OUT = [];
                return
            end
        else
            OUT = [];
            return
        end
    otherwise
        OUT = [];
        return
end
if ~isempty(cal_level)
        callevel = cal_level;
        if size(IN.audio,2) < length(cal_level), callevel = cal_level(1:size(IN.audio,2)); end
        if size(IN.audio,2) > length(cal_level), callevel = [cal_level NaN(1,size(IN.audio,2)-length(cal_level))]; end
        OUT = IN;
        OUT.cal = callevel;
        OUT.funcallback.name = 'cal_aarae.m';
        OUT.funcallback.inarg = {method,values,fileref};
else
    OUT = [];
end

