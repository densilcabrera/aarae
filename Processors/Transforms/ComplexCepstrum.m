function out = ComplexCepstrum(in)
% This function returns the complex cepstrum of the input audio
% using Matlab's cceps function. If complex valued data are input, then
% only the real part is used to derive the cepstrum.
out = in;
[~,chans,bands] = size(in);
for ch = 1:chans
    for b = 1:bands
        out.audio(:,ch,b) = cceps(real(in.audio(:,ch,b)));
    end
end
out.funcallback.name = 'ComplexCepstrum.m';
out.funcallback.inarg = {};