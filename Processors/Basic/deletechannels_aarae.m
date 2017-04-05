function out = deletechannels_aarae(in,S,S2,ok)
% This function allows you to select the channels that you wish to retain,
% and delete the remaining channels.
% 29/03/2017: update to select output channels (dim5)

%if isfield(in,'chanID')
%    param = in.chanID;
%else
%    param = [];
%end
if isstruct(in)
    data = in.audio;
    if isfield(in,'chanID')
        if length(in.chanID) == size(data,2)
            param = in.chanID;
        else
            param = cellstr([repmat('Chan',size(data,2),1) num2str((1:size(data,2))')]);
        end
    else
        param = cellstr([repmat('Chan',size(data,2),1) num2str((1:size(data,2))')]);
    end
    
    if isfield(in,'OutchanID')
        if length(in.OutchanID) == size(data,5)
            param2 = in.OutchanID;
        else
            param2 = cellstr([repmat('Chan',size(data,5),1) num2str((1:size(data,5))')]);
        end
    else
        param2 = cellstr([repmat('Chan',size(data,5),1) num2str((1:size(data,5))')]);
    end
else
    data = in;
    param = cellstr([repmat('Chan',size(data,2),1) num2str((1:size(data,2))')]);
    param2 = cellstr([repmat('Chan',size(data,5),1) num2str((1:size(data,5))')]);
end




if ~isempty(param) && ~isempty(param2)
    if nargin < 2
        [S,ok] = listdlg('Name','Channel selection',...
            'PromptString','Delete unselected channels',...
            'ListString',param,...
            'ListSize', [160 320]);
        if size(param2,1) >1 
        [S2,ok2] = listdlg('Name','Channel selection',...
            'PromptString','Delete unselected Output channels',...
            'ListString',param2,...
            'ListSize', [180 320]);
        if ok ==1 && ok2 == 1
            ok = 1;
        else
            ok = 0;
        end
        else
            S2 = 1:size(data,5);
        end
           
    end

    if ok == 1 && ~isempty(S) && ~isempty(S2)
        try
            out = in;
            out.audio = in.audio(:,S,:,:,S2,:);
            if isfield(in,'chanID') && length(in.chanID) == size(in.audio,2)
                out.chanID = in.chanID(S);
                if isfield(in,'cal')
                    out.cal = in.cal(S);
                end
            else
                % write new chan IDs (preserving original chan numbers)
                out.chanID = cellstr([repmat('Chan',size(out.audio,2),1) num2str(S')]);
                if isfield(in,'cal')
                    out = rmfield(out,'cal'); % remove cal field if it cannot be interpreted
                end
            end
            if isfield(in,'OutchanID') && length(in.OutchanID) == size(in.audio,5)
                out.OutchanID = in.OutchanID(S2);
            else
                % write new Outchan IDs (preserving original chan numbers)
                out.OutchanID = cellstr([repmat('OutChan',size(out.audio,5),1) num2str(S2')]);
            end
            out.funcallback.name = 'deletechannels_aarae.m';
            out.funcallback.inarg = {S,S2,ok};
        catch
            out = [];
        end
    else
        out = [];
    end
else
    out = [];
    h = warndlg('Audio has a single channel','AARAE info','modal');
    uiwait(h)
end
end