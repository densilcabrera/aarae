function out = removeDC_aarae(in)
% This function subtracts the average amplitude of the wave from the wave,
% and this is done independently for each channel and band (if more than
% one exists)
out.audio = in.audio;
indices = partial_selection(in);
len = size(in.audio(indices{:}),1);
out.audio(indices{:}) = in.audio(indices{:}) - repmat(mean(in.audio(indices{:})),[len,1,1,1,1,1]);
out.funcallback.name = 'removeDC_aarae.m';
out.funcallback.inarg = {};
end