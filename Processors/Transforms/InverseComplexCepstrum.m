function out = InverseComplexCepstrum(in)
% This function returns the inverse complex cepstrum of the input audio
% using Matlab's icceps function.

[~,chans,bands] = size(in.audio);
out = in;
for ch = 1:chans
    for b = 1:bands
        out.audio(:,ch,b) = icceps(in.audio(:,ch,b));
    end
end
out.funcallback.name = 'InverseComplexCepstrum.m';
out.funcallback.inarg = {};