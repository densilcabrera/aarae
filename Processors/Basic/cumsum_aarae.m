function out = cumsum_aarae(in)
% Outputs the cumulative sum of the audio input waveform, similar to
% integration of the waveform.
    out.audio = cumsum(in.audio);
    out.funcallback.name = 'cumsum_aarae.m';
    out.funcallback.inarg = {};
end