function out = abs_aarae(in)
% The input audio waveform is rectified, or if a complex waveform is input
% then its magnitude is returned.

    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = abs(in.audio(indices{:}));
    out.funcallback.name = 'abs_aarae.m';
    out.funcallback.inarg = {};