function out = audiofield_extract(in)
% If more than one audio field exists (usually audio and audio2) then this
% function replaces audio with audio2, and most of the fields associated
% with the audio are removed (fs is retained).
%
% If more than two audio fields exist (not usual), the the above is done,
% and audio3 replaces audio2 etc.


% list of possible audio field names
audiofieldnames = {'audio','audio2','audio3','audio4','audio5','audio6',...
    'audio7','audio8','audio9','audio10','audio11','audio12','audio13',...
    'audio14','audio15','audio16'};

% check which audio fields exist
audiofields = [isfield(in,audiofieldnames{1}), ...
    isfield(in,audiofieldnames{2}), ...
    isfield(in,audiofieldnames{3}), ...
    isfield(in,audiofieldnames{4}), ...
    isfield(in,audiofieldnames{5}), ...
    isfield(in,audiofieldnames{6}), ...
    isfield(in,audiofieldnames{7}), ...
    isfield(in,audiofieldnames{8}), ...
    isfield(in,audiofieldnames{9}), ...
    isfield(in,audiofieldnames{10}), ...
    isfield(in,audiofieldnames{11}), ...
    isfield(in,audiofieldnames{12}), ...
    isfield(in,audiofieldnames{13}), ...
    isfield(in,audiofieldnames{14}), ...
    isfield(in,audiofieldnames{15}), ...
    isfield(in,audiofieldnames{16})];

% remove fields that don't exist from the list of field names
audiofieldnames = audiofieldnames(audiofields);

 

if sum(audiofields) == 1
    % if only one audio field exists, then there is nothing to do
    warndlg('Unable to extract another audio field because only one audio field exists','AARAE info','modal');
    out = [];
    return
% elseif sum(audiofields) == 2
%     % if only two audio fields exist, then extract the second field
%     out.(audiofieldnames{1}) = in.(audiofieldnames{2});
%     out.fs = in.fs;
%     out.chanID = cellstr([repmat('Chan',size(out.audio,2),1) num2str((1:size(out.audio,2))')]);
else
    for n = 2:sum(audiofields)
    % shift audio fields down by 1
        out.(audiofieldnames{n-1}) = in.(audiofieldnames{n});
    end
    
    % The following line would circularly shift the main audio field to the
    % highest
    %out.(audiofieldnames{sum(audiofields)}) = in.(audiofieldnames{1});
    % or return an empty field instead:
    out.(audiofieldnames{sum(audiofields)}) = [];
    
    out.fs = in.fs;
    
      
    
    % zap the cal field
    if isfield(in,'cal')
        out.cal = [];
%       out = rmfield(out,'cal') % fields cannot be deleted within a
%       processor because aarae preserves the fields in its framework
    end

    
    % overwrite chanIDs
    %out.chanID = cellstr([repmat('Chan',size(out.audio,2),1) num2str((1:size(out.audio,2))')]);
    out.chanID = makechanID(size(out.audio,2),0);

    
    % zap bandID field
    if isfield(in,'bandID')
        out.bandID = [];
%       out = rmfield(out,'bandID') % fields cannot be deleted within a
%       processor because aarae preserves the fields in its framework
    end
    
    % potentially zap other fields (in future revisions)

end
out.funcallback.name = 'audiofield_extract.m';
out.funcallback.inarg = {};