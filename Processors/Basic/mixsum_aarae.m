function out = mixsum_aarae(in)
% The input audio waveform is summed across channels (dimension 2), i.e. a
% mixdown of the channels, by summing them.
    out.audio = sum(in.audio,2);
    out.chanID = {'1'};
    out.funcallback.name = 'mixsum_aarae.m';
    out.funcallback.inarg = {};
end