function out = reverse_aarae(in)
% The input audio waveform is time-reversed
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = flipdim(in.audio(indices{:}),1);
    out.funcallback.name = 'reverse_aarae.m';
    out.funcallback.inarg = {};
end