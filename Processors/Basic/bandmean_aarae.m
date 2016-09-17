function out = bandmean_aarae(in)
% The input audio waveform is averaged across bands (dimension 3), i.e. a
% mixdown of the bands, by averaging them.
    out.audio = mean(in.audio,3);
    out.bandID = 1;
    out.funcallback.name = 'bandmean_aarae.m';
    out.funcallback.inarg = {};
end