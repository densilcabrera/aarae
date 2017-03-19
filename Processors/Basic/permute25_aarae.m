function out = permute25_aarae(in)
% The 2nd and 5th dimensions of the input audio are swapped.
%
% Note, cal may need to be applied channels before running this. Consider
% updating all scripts where cal is applied to check for chan||OutChan???
%
% Code by Jonothan Holmes
% version 2.0 (Feb 2017)

    out = in;
if ndims(in.audio) > 1
    out.audio = permute(in.audio,[1,5,3,4,2,6]);
    
    if isfield(in,'OutchanID')
        if iscell(in.OutchanID)
            out.chanID = in.OutchanID;
        else
            out.chanID = {in.OutchanID};
        end
    else
        out.chanID = cellstr([repmat('chan ',size(out.audio,2),1) num2str((1:size(out.audio,2))')]);
    end
    
    if isfield(in,'chanID')
        out.OutchanID = in.chanID;
    end
    %%% cal field
    
    out.funcallback.name = 'permute25_aarae.m';
    out.funcallback.inarg = {};
else
    out = [];
    h = warndlg('Data is one-dimensional','AARAE info','modal');
    uiwait(h)
end