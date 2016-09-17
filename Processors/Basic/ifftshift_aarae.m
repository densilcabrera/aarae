function out = ifftshift_aarae(in)
% Performs Matlab's ifftshift()
    out.audio = in.audio;
    indices = partial_selection(in);
    out.audio(indices{:}) = ifftshift(in.audio(indices{:}),1);
    out.funcallback.name = 'ifftshift_aarae.m';
    out.funcallback.inarg = {};
end