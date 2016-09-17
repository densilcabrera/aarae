function OUT = editfields_aarae(IN)
% This function allows you to edit some of the fields of an audio leaf.
%
% WARNING: it is best to know what you are doing if you are editing fields
% - some edits could cause errors or make the audio leaf unuseable for a
% particular function.
%
% Currently the explicitly supported fields are:
%
% * fs - audio sampling rate. Editing this allows the audio to be
%   interpreted differently on playback, processing and analysis. Note that
%   some sampling rates may not be possible to use for playback.
%
% * bandID - a vector, specifying the centre frequency of bands in
%   Hz. This is used by some analysers, mainly for display, but sometimes
%   also it may affect the analysis. If bandID is present, and fs is
%   edited, then it may make sense to edit bandIDs accordingly.
%
% * chanID - a cell array of strings describing each channel. Note that
%   aarae has some common formats for chanID (see the function makechanID
%   in the Utilities directory). Other formats can be used if they do not
%   need to be interpreted by readchanID (also in Utilities).
%
% * cal - calibration offset in dB. Note that you can edit this in the main
%   AARAE GUI using the cal button (although it might be easier to edit
%   multichannel cal data here).
%
% Other fields can usually also be edited, using the generic interface
% (which might or might not work for a given field). This includes
% subfields of properties. Note that the functioncallback and datatype
% fields cannot be edited. Also audio fields cannot be edited using this
% function.



% Make a list of editable fields
fnames = fieldnames(IN);
if isfield(IN,'properties')
    fnamesprop = fieldnames(IN.properties);
else
    fnamesprop = [];
end

isproperties = [zeros(length(fnames),1);ones(length(fnamesprop),1)];
iseditable = false(size(isproperties));
numfields = length(iseditable);

fields = [fnames;fnamesprop];


for n = 1:numfields
    if isempty(regexp(char(fields{n}),'audio','once')) ...
            && isempty(regexp(char(fields{n}),'funcallback','once')) ...
        && isempty(regexp(char(fields{n}),'datatype','once')) ...
        && isempty(regexp(char(fields{n}),'properties','once'))
        iseditable(n) = true; 
    end
end

 ind = find(iseditable);
 isproperties = isproperties(ind);
 fields = fields(ind);

[S,ok] = listdlg('Name','Edit Fields',...
    'PromptString','Select field(s) to edit',...
    'ListString',fields);
if ok == 0
    OUT = [];
    return
end
chosenfields = fields(S);
isproperties = isproperties(S);

for n = 1: length(chosenfields)
    
    
    % fs
    if strcmp(chosenfields{n},'fs')
        
        prompt = {['Current sampling rate is ',...
            num2str(IN.fs),...
            ' Hz. (Note that editing this will not result in resampling of the audio.)']};
        dlg_title = 'Edit Sampling Rate';
        num_lines = 1;
        def = {num2str(IN.fs)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if ~isempty(answer)
            fs = abs(round(str2double(answer)));
            if ~isnan(fs)
                OUT.fs = fs;
            end
        end
        
        
        
        % bandID
    elseif  strcmp(chosenfields{n},'bandID')
        f = figure;
        
        dat = IN.bandID(:);
        columnname = {'bandID'};
        columnformat = {'numeric'};
        columneditable = true;
        t = uitable(f,...
            'Units','normalized',...
            'Position',[0 0 0.5 1],...
            'Data',dat,...
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columneditable,...
            'Rowname',[]);
        try
            h = uicontrol(f,...
                'Units','normalized',...
                'Position',[0.6 0.1 0.3 0.2],...
                'String','Continue',...
                'Callback','uiresume(gcbf)');
            
            uiwait(gcf);
            % retrieve handle to uitable
            tH = findobj(gcf,'Type','uitable');
            % retrieve data
            answer = get(tH,'Data');
            if length(answer) == size(IN.audio,3)
                OUT.bandID = answer;
            else
                OUT.bandID=IN.bandID;
            end
            close(f);
        catch
            OUT.bandID=IN.bandID;
        end
        
        
        % chanID
    elseif strcmp(chosenfields{n},'chanID')
        f = figure;
        
        dat = IN.chanID(:);
        columnname = {'chanID'};
        columnformat = {'char'};
        columneditable = true;
        t = uitable(f,...
            'Units','normalized',...
            'Position',[0 0 0.5 1],...
            'Data',dat,...
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columneditable,...
            'Rowname',[]);
        try
            h = uicontrol(f,...
                'Units','normalized',...
                'Position',[0.6 0.1 0.3 0.2],...
                'String','Continue',...
                'Callback','uiresume(gcbf)');
            
            uiwait(gcf);
            % retrieve handle to uitable
            tH = findobj(gcf,'Type','uitable');
            % retrieve data
            answer = get(tH,'Data');
            if length(answer) == size(IN.audio,2)
                OUT.chanID = answer;
            else
                OUT.chanID=IN.chanID;
            end
            
            close(f);
        catch
            OUT.chanID=IN.chanID;
        end
        
        
        % cal
    elseif  strcmp(chosenfields{n},'cal')
        f = figure;
        
        dat = IN.cal(:);
        columnname = {'cal'};
        columnformat = {'numeric'};
        columneditable = true;
        t = uitable(f,...
            'Units','normalized',...
            'Position',[0 0 0.5 1],...
            'Data',dat,...
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columneditable,...
            'Rowname',[]);
        try
            h = uicontrol(f,...
                'Units','normalized',...
                'Position',[0.6 0.1 0.3 0.2],...
                'String','Continue',...
                'Callback','uiresume(gcbf)');
            
            uiwait(gcf);
            % retrieve handle to uitable
            tH = findobj(gcf,'Type','uitable');
            % retrieve data
            answer = get(tH,'Data');
            if length(answer) == size(IN.audio,2)
                OUT.cal = answer;
            else
                OUT.cal=IN.cal;
            end
            close(f);
        catch
            OUT.cal=IN.cal;
        end
        
    else
        % Generic editor (which might not always work)
        
        if isproperties(n)
            dat = IN.properties.(chosenfields{n});
        else
            dat = IN.(chosenfields{n});
        end
        
        
        if length(dat) == 1
            % single value editor
            try
            prompt = ['Edit ',chosenfields{n}];
            dlg_title = 'Single Value Edit';
            num_lines = 1;
            if ischar(dat)
                def = {dat};
            else
                def = {num2str(dat)};
            end
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if ~isempty(answer)
                if ischar(dat)
                    if isproperties(n)
                        OUT.properties(chosenfields{n}) = answer;
                    else
                        OUT.(chosenfields{n}) = answer;
                    end
                else
                    if isproperties(n)
                        OUT.properties(chosenfields{n}) = str2double(answer);
                    else
                        OUT.(chosenfields{n}) = str2double(answer);
                    end
                end
                
            end
            catch
            end
            
        else
            % uitable editor
            
            
            f = figure;
            try
            [~,c] = size(dat);
            
            columnname = {chosenfields(n)};
            if ischar(dat(1,1))
                columnformat = {'char'};
            else
                columnformat = {'numeric'};
            end
            columneditable = true(1,c);
            t = uitable(f,...
                'Units','normalized',...
                'Position',[0 0 0.5 1],...
                'Data',dat,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnEditable', columneditable,...
                'Rowname',[]);
            
                h = uicontrol(f,...
                    'Units','normalized',...
                    'Position',[0.6 0.1 0.3 0.2],...
                    'String','Continue',...
                    'Callback','uiresume(gcbf)');
                
                uiwait(gcf);
                % retrieve handle to uitable
                tH = findobj(gcf,'Type','uitable');
                % retrieve data
                answer = get(tH,'Data');
                
                if ischar(dat(1,1))
                    if isproperties(n)
                        OUT.properties.(chosenfields{n}) = answer;
                    else
                        OUT.(chosenfields{n}) = answer;
                    end
                else
                    if isproperties(n)
                        OUT.properties.(chosenfields{n}) = str2double(answer);
                    else
                        OUT.(chosenfields{n}) = str2double(answer);
                    end
                end
                
                close(f);
            catch
                close(f);
%                 if isproperties(n)
%                     OUT.properties.(chosenfields{n}) = IN.properties.(chosenfields{n});
%                 else
%                     OUT.(chosenfields{n}) = IN.(chosenfields{n});
%                 end
            end
            
            
            
        end
    end
    

if ~exist('OUT','var')
    OUT = [];
end
    
end
end % eof
