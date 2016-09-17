function out = square_aarae(in)
% The input audio waveform is squared
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = in.audio(indices{:}).^2;
    out.funcallback.name = 'square_aarae.m';
    out.funcallback.inarg = {};
end