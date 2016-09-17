function [bm,env,delay] = gammatoneFast(x,cfs,fs,align)
% Produce an array of responses from a Gammatone filter via FFT
%   
%   bm = gammatoneFast(x,cfs,fs)
%   bm = gammatoneFast(...,align)
%   [bm,env] = gammatoneFast(...)
%   [bm,env,delay] = gammatoneFast(...)
%
%   This function takes an input vector and passes it
%   through a bank of fourth-order gammatone filters, with
%   centre frequencies specified by cfs. The function
%   returns a matrix, with each row/column corresponding to
%   a filter output with a centre frequency determined by
%   the corresponding element in cfs. The orientation of the
%   output is determined by the orientation of the input: if
%   x is a row vector then the output will contain one row
%   for each filter output, and vice versa.
% 
%   Centre frequencies may be any value below the Nyquist
%   rate (determined by the sampling frequency fs).
%   Typically centre frequencies are equally spaced on the
%   ERB-rate scale and may be calculated thus:
%
%   cfs = MakeErbCFs(low_cf,high_cf,numchans)
%   
%   where low_cf is the lowest frequency in the bank,
%   high_cf is the highest, and numchans is the numbers of
%   filters in the bank.
% 
%   bm = gammatoneFast(...,align) allows phase alignment to
%   be applied. With align=false, no alignment is applied
%   (default). With align=true, fine structure and envelope
%   alignment is applied so that the impulse response peaks
%   occurs at t=0.
% 
%   [bm,env] = gammatoneFast(...) returns the instantaneous
%   envelopes env for each filter.
% 
%   [bm,env,delay] = gammatoneFast(...) returns the delay
%   (in samples) removed by the phase alignment of each
%   gammatone filter, i.e. delays(n)=0 if align=false. delay
%   is a vector the same size as cfs.
%
%   Based on code written by ZZ Jin, adapted by DLW in
%   Jan'07 and JF Woodruff in Nov'08
% 
%   See also MakeErbCFs.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

if nargin < 3
    fs = 16000; % default sampling frequency
end
if nargin < 4
    align = false; % default phase alignment
end

% check inputs
assert(isvector(x) & isnumeric(x),'x must be a vector')
assert(isvector(cfs) & isnumeric(cfs),'cfs must be a vector')
assert(isscalar(fs),'fs must be a scalar')
assert(islogical(align) & numel(align)==1,'align must be logical')

% number of frequency channels
numchans = length(cfs);

filterOrder = 4; % filter order
gL = 2^nextpow2(0.128*fs); % gammatone filter length at least 128 ms
b = 1.019.*24.7.*(4.37.*cfs./1000+1); % rate of decay or bandwidth

gt = zeros(gL,numchans);  % Initialise IR
tc = zeros(size(cfs));  % Initialise time lead
phase = 0;

tpt=(2*pi)/fs;
gain=((1.019.*b.*tpt).^filterOrder)./6; % based on integral of impulse

tmp_t = (0:gL-1)/fs;

% calculate impulse response
for i = 1:numchans
    if align
        tc(i) = (filterOrder-1)./(2*pi*b(i));
        phase = -2*pi*cfs(i)*tc(i);
    end
    gt(:,i) = gain(i)*fs^3*tmp_t.^(filterOrder-1).*exp(-2*pi*b(i)*tmp_t).*cos(2*pi*cfs(i)*tmp_t+phase);
end

% if input is row vector, transpose to column vector
rot = false;
if size(x,1)==1
    x = x';
    rot = true;
end

% gammatone filtering using FFTFILT
bm = fftfilt(gt,repmat(x,1,numchans));

% Hilbert envelope
env = abs(hilbert(bm));

% delay due to time lead
delay = round(tc.*fs);

% remove time lead
for i = 1:numchans
    bm(:,i) = [bm(delay(i)+1:end,i); zeros(delay(i),1)];
    env(:,i) = [env(delay(i)+1:end,i); zeros(delay(i),1)];
end

% transpose output if necessary
if rot
    bm = bm';
    env = env';
end

% [EOF]
