function L = true_rms(sig,dbfs)
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% L = true_rms(sig,dbfs);
% calculates level (dB SPL) of input sig according to calibration of DLM
% DLM assumes that a pure tone with a maximum amplitude of "1" has a level
% of 107 dB SPL ("dB full scale")
% optional: other calibrations can be used by setting dbfs to another value
% than 107
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

sig=2*10^.5*sig;
signal=sig.^2;
summe=sum(signal);
peff=sqrt(summe/length(signal));
if peff == 0
    peff=realmin;
end
L=20*log10(peff/2e-5);

if nargin>1
    L=L-(107-dbfs);
end
end