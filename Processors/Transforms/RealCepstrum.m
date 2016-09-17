function out = RealCepstrum(in)
% This function returns the real cepstrum of the input audio
% using Matlab's rceps function.

out.audio = rceps(in.audio);
out.funcallback.name = 'RealCepstrum';
out.funcallback.inarg = {};