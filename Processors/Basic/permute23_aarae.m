function out = permute23_aarae(in)
% The 2nd and third dimensions of the input audio are swapped.
% Hence channels become bands and bands become channels.
% ChanID and bandID fields are also swapped (if they exist).
%
% Code by Densil Cabrera & Daniel Jimenez
% version 1.0 (5 November 2013)
if ndims(in.audio) > 2
    out.audio = permute(in.audio,[1,3,2,4,5,6]);

    if isfield(in,'bandID')
        out.chanID = cellstr(num2str(in.bandID'));
    end

    if isfield(in,'chanID')
        out.bandID = str2num(char(in.chanID));
        if ~isrow(out.bandID), out.bandID = out.bandID'; end
    end
    
    if isempty(out.chanID{1,1})
        out.chanID = cellstr([repmat('Chan',size(out.audio,2),1) num2str((1:size(out.audio,2))')]);
    end
    out.funcallback.name = 'permute23_aarae.m';
    out.funcallback.inarg = {};
else
    out = [];
    h = warndlg('Data is two-dimensional','AARAE info','modal');
    uiwait(h)
end