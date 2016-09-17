function out = imaginary_aarae(in)
% The imaginary part of the input audio is returned.
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = imag(in.audio(indices{:}));
    out.funcallback.name = 'imaginary_aarae.m';
    out.funcallback.inarg = {};
end