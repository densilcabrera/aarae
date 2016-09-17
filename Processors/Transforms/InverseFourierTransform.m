function out = InverseFourierTransform(in)

out.audio = ifft(in.audio);
out.funcallback.name = 'InverseFourierTransform.m';
out.funcallback.inarg = {};