function out = difference_aarae(in)
% The input audio waveform is differenced, similar to differentiation
    out.audio = diff(in.audio);
    out.funcallback.name = 'difference_aarae.m';
    out.funcallback.inarg = {};
end