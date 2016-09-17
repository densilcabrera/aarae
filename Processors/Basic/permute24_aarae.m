function out = permute24_aarae(in)
% The 2nd and 4th dimensions of the input audio are swapped.
%
% Code by Densil Cabrera
% version 1.0 (2 August 2014)
if ndims(in.audio) > 1
    out.audio = permute(in.audio,[1,4,3,2,5,6]);
    
    if isfield(in,'dim4ID')
        if iscell(in.dim4ID)
            out.chanID = in.dim4ID;
        else
            out.chanID = {in.dim4ID};
        end
    else
        out.chanID = cellstr([repmat('chan ',size(out.audio,2),1) num2str((1:size(out.audio,2))')]);
    end
    
    if isfield(in,'chanID')
        out.dim4ID = in.chanID;
    end
    
    out.funcallback.name = 'permute24_aarae.m';
    out.funcallback.inarg = {};
else
    out = [];
    h = warndlg('Data is one-dimensional','AARAE info','modal');
    uiwait(h)
end