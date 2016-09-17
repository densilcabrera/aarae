function [Y y]=mpf(Y_abs)
%[Y y]=mpf(Y_abs)
%calculates milimum phase filter complex response from magnitude response
%Y_abs should be in fft form, i.e. Y_abs=abs(fft(y)) where y(t) is a signal
%This is a slightly expanded form from the Mathworks help.

Y_log=log(abs(Y_abs));
logmin=-100;
Y_log(find(Y_log<logmin))=logmin;


N=numel(Y_abs);
y = real(ifft(Y_log));
if round(N/2)==N/2
  w = [1;2*ones(N/2-1,1);ones(1-rem(N,2),1);zeros(N/2-1,1)]';
else
  error('Odd N not yet implemented, use an even fft.')
end
Y=exp(fft(w.*y));
if nargout>1
  y = real(ifft(Y));
end