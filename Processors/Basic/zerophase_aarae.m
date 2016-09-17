function out = zerophase_aarae(in)
% outputs a zero phase version of the input waveform
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = ifft(abs(fft(in.audio(indices{:}))));
    out.funcallback.name = 'zerophase_aarae.m';
    out.funcallback.inarg = {};
end