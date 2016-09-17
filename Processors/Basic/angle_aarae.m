function out = angle_aarae(in)
% The angle of the input audio is returned 
% (may be useful for complex signals).
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = angle(in.audio(indices{:}));
    out.funcallback.name = 'angle_aarae.m';
    out.funcallback.inarg = {};
end