%
%	Calculate Fp2 from Fr1
%
%	function [Fp2,Fr2] = Fr1toFp2(n,b1,c1,b2,c2,frat,Fr1)
%	Toshio Irino
%	Created:  17 Nov. 2006 ( from old vesion by M.Unoki, 11 July 2002 )
%       Modified: 17 Nov. 2006 
%
%
function [Fp2, Fr2] = Fr1toFp2(n,b1,c1,b2,c2,frat,Fr1)
if nargin < 1; help Fr1toFp2; end;
if nargin < 8; SR=24000; end;
if nargin < 9; Nfft=1024*2; end;

%%%%%%% 
[dummy ERBw1] = Freq2ERB(Fr1);
Fp1 = Fr2Fpeak(n,b1,c1,Fr1);
Fr2 = frat*Fp1;
[dummy ERBw2] = Freq2ERB(Fr2);

Bw1 = b1*ERBw1;
Bw2 = b2*ERBw2;
%%%%%%%	Coef1*Fp2^3 + Coef2*Fp2^2 + Coef3*Fp2 + Coef4 = 0 %%%%%%%
Coef1 = -n;
Coef2 = c1*Bw1+c2*Bw2+n*Fr1+2*n*Fr2;
Coef3 = -2*Fr2*(c1*Bw1+n*Fr1)-n*((Bw2)^2+Fr2^2)-2*c2*Bw2*Fr1;
Coef4 = c2*Bw2*((Bw1)^2+Fr1^2)+(c1*Bw1+n*Fr1)*((Bw2)^2+Fr2^2);
p = roots([Coef1 Coef2 Coef3 Coef4]);
Fp2cand = p(imag(p)==0);
if length(Fp2cand) == 1,
    Fp2 = Fp2cand;
else
   [val ncl] = min(abs(Fp2cand - Fp1));
   Fp2 = Fp2cand(ncl); % in usual cGC range, Fp2 is close to Fp1
end;

return

%%%%% Check %%%%%%%%%%%%%%%%%%%%%%
    fs = 48000; NfrqRsl = 2048;
    [cGCresp] = CmprsGCFrsp(Fr1,fs,n,b1,c1,frat,b2,c2,NfrqRsl);
    for nn = 1:length(Fp2cand)
      [dummy nFr2(nn)] = min(abs(cGCresp.freq - Fp2cand(nn)));
    end;
    plot(cGCresp.freq, cGCresp.cGCFrsp/max(cGCresp.cGCFrsp),...
	cGCresp.freq, cGCresp.pGCFrsp)
    ax = axis;
    axis([0 max(Fp2cand)*2, ax(3:4)])
    cGCresp.cGCFrsp(nFr2)
    [val nmax] = max(cGCresp.cGCFrsp(nFr2));
    Fp2n = Fp2cand(nmax)
