function out = real_aarae(in)
% The real part of the input audio is returned.
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = real(in.audio(indices{:}));
    out.funcallback.name = 'real_aarae.m';
    out.funcallback.inarg = {};
end