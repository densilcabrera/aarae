function out = audio2data(in)

t = linspace(0,length(in.audio)/in.fs,length(in.audio));
if ismatrix(in.audio)
    doresultleaf(in.audio,'Amplitude',{'Time','channels','Frequency'},...
                 'Time',      t,                   's',           true,...
                 'channels',  in.chanID,           'categorical', [],...
                 'name','audiodata');
elseif ndims(in.audio) == 3
    doresultleaf(in.audio,'Amplitude',{'Time','channels','Frequency'},...
                 'Time',      t,                   's',           true,...
                 'channels',  in.chanID,           'categorical', [],...
                 'Frequency', num2cell(in.bandID), 'Hz',          false,...
                 'name','audiodata');
end
out = [];