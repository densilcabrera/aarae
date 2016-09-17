function out = mixmean_aarae(in)
% The input audio waveform is averaged across channels (dimension 2), i.e. a
% mixdown of the channels, by averaging them.
    out.audio = mean(in.audio,2);
    out.chanID = {'1'};
    out.funcallback.name = 'mixmean_aarae.m';
    out.funcallback.inarg = {};
end