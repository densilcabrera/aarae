function out = invert_aarae(in)
% The input audio waveform is multiplied by -1 (i.e. inverted, or phase-reversed)
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = -in.audio(indices{:});
    out.funcallback.name = 'invert_aarae.m';
    out.funcallback.inarg = {};
end