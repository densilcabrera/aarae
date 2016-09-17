%
%	Calculate Fr1 from Fp2
%
%	function [Fr1,Fp1] = Fr1toFp2(n,b1,c1,b2,c2,frat,Fp2)
%	Toshio Irino
%	Created:  17 Nov. 2006 ( from old vesion by M.Unoki, 3 July 2002 )
%       Modified: 17 Nov. 2006 
%
%
function [Fr1, Fp1] = Fp2toFr1(n,b1,c1,b2,c2,frat,Fp2)
if nargin < 1; help Fp2toFr1; end;

%%%%%%% Coefficients: ERBw(Fr1)=alp1*Fr1+alp0 %%%%%%%
[dummy alp0] = Freq2ERB(0);
[dummy w1 ] =  Freq2ERB(1);
alp1 = w1 - alp0;

%%%%%%% Coefficients: fr2=bet1*fr2+bet0 %%%%%%%
bet1=frat*(1+c1*b1*alp1/n);
bet0=frat*c1*b1*alp0/n;

%%%%%%% Coefficients: ERB(fr2)=zet1*Fr1+zet0 %%%%%%%
zet1=alp1*bet1;
zet0=alp1*bet0+alp0;

%%%%%%%	Coef1*Fr1^3 + Coef2*Fr1^2 + Coef3*Fr1 + Coef4 = 0 %%%%%%%
Coef1=((b2^2*zet1^2+bet1^2)*(c1*b1*alp1+n) + (c2*b2*zet1)*(b1^2*alp1^2+1));
Coef2=((b2^2*zet1^2+bet1^2)*(c1*b1*alp0-n*Fp2) ...
    + (2*b2^2*zet1*zet0-2*bet1*(Fp2-bet0))*(c1*b1*alp1+n) ...
    + (c2*b2*zet1)*(2*b1^2*alp1*alp0-2*Fp2) + (c2*b2*zet0)*(b1^2*alp1^2+1));
Coef3=((2*b2^2*zet1*zet0-2*bet1*(Fp2-bet0))*(c1*b1*alp0-n*Fp2) ...
    + (b2^2*zet0^2+(Fp2-bet0)^2)*(c1*b1*alp1+n)...
    + (c2*b2*zet1)*(b1^2*alp0^2+Fp2^2) ...
    + (c2*b2*zet0)*(2*b1^2*alp1*alp0-2*Fp2) );
Coef4=(b2^2*zet0^2+(Fp2-bet0)^2)*(c1*b1*alp0-n*Fp2) ...
    + (c2*b2*zet0)*(b1^2*alp0^2+Fp2^2);

q=roots([Coef1 Coef2 Coef3 Coef4]);
Fr1cand=q(imag(q)==0);

if length(Fr1cand) == 1,
    Fr1 = Fr1cand;
    Fp1 = Fr2Fpeak(n,b1,c1,Fr1);
else
   Fp1cand = Fr2Fpeak(n,b1,c1,Fr1cand);
   [val ncl] = min(abs(Fp1cand - Fp2));
   Fp1 = Fp1cand(ncl); % in usual cGC range, Fp2 is close to Fp1
   Fr1 = Fr1cand(ncl);
end;

return;

%%%%% Check %%%%%%%%%%%%%%%%%%%%%%
    fs = 48000; NfrqRsl = 2048;
    for nn = 1:length(Fr1cand)
      Fr1ck = Fr1cand(nn);
      [cGCresp] = CmprsGCFrsp(Fr1ck,fs,n,b1,c1,frat,b2,c2,NfrqRsl);
      plot(cGCresp.freq, cGCresp.cGCFrsp/max(cGCresp.cGCFrsp),...
	   cGCresp.freq, cGCresp.pGCFrsp)
      hold on
    end;
    ax = axis;
    axis([0 max(Fr1cand)*2, ax(3:4)])
    drawnow
