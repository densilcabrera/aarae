function indices = partial_selection(in)

handles = guidata(findobj('Tag','aarae'));

audiosize = size(in.audio);
prompt = cell(1,ndims(in.audio));
def = cell(1,ndims(in.audio));
dlgtitle = 'Data selection';

for i = 1:ndims(in.audio)
    if i == 1, prompt{1,i} = ['Time in samples: (max. ' num2str(audiosize(1,i)) ')'];
    elseif i == 2, prompt{1,i} = ['Channels: (max. ' num2str(audiosize(1,i)) ')'];
    elseif i == 3, prompt{1,i} = ['Bands: (max. ' num2str(audiosize(1,i)) ')'];
    else prompt{1,i} = ['Dim' num2str(i) ': (max. ' num2str(audiosize(1,i)) ')'];
    end
    def{1,i} = ['1:' num2str(audiosize(1,i))];
end

answer = inputdlg(prompt,dlgtitle,[1 50],def);
originalindices = cellfun(@eval,def,'UniformOutput',false);
if any(cellfun(@isempty,answer)) || isempty(answer)
    warndlg('Invalid selection, process will be applied to the whole audio file.','AARAE info','modal')
    fprintf(handles.fid,'%% Invalid selection - partial selection was not executed. \n');
    indices = originalindices;
else
    indices = cell(size(answer));
    for i = 1:length(answer), indices{i,1} = ['[' answer{i,1} ']']; end
    indices = cellfun(@eval,indices,'UniformOutput',false);
    partialsel = false;
    for i = 1:length(answer)
        if numel(indices{i}) ~= numel(originalindices{i})
            partialsel = true;
        end
    end
    
    try
        in.audio(indices{:});
        if partialsel
            selectionstring = 'X.audio = X.audio(';
            for i = 1:length(answer)
                dimstring = char(answer{i});
                if ~isempty(regexp(dimstring,',','once')) || ~isempty(regexp(dimstring,' ','once'))
                    dimstring = ['[', dimstring, ']'];
                end
                if i<length(answer)
                    selectionstring = [selectionstring, dimstring,','];
                else
                    selectionstring = [selectionstring, dimstring,');'];
                end
            end
            
            fprintf(handles.fid,'%% Partial selection of audio using AARAE''s partial_selection.m utility function - the following code is equivalent:\n');
            %fprintf(handles.fid,'COMPLETE.audio = X.audio;\n');
            fprintf(handles.fid,[selectionstring,'\n']);
%             handles.partialselindices=answer; % why do we need/want to do this?
%             guidata(findobj('Tag','aarae'),handles); % maybe this is dangerous!
        end
    catch
        warndlg('Invalid selection, process will be applied to the whole audio file.','AARAE info','modal')
        fprintf(handles.fid,'%% Invalid selection - partial selection was not executed. \n');
        indices = originalindices;
    end
end