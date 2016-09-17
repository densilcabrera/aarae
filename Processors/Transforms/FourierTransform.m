function out = FourierTransform(in)

out.audio = fft(in.audio);
out.funcallback.name = 'FourierTransform.m';
out.funcallback.inarg = {};