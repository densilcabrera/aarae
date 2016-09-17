%
%	Frequency Response of Compressive GammaChirp
%	Toshio IRINO
%	Created : 9 Sept. 2004
%	Modified: 9 Sept. 2004
%	Modified: 6 July  2005 % fixing bug in NfrqRsl
%
%   function cGCresp = CmprsGCFrsp(Fr1,fs,n,b1,c1,frat,b2,c2,NfrqRsl);
%	INPUT : Fr1	: Resonance Freq. (vector)
%		fs 	: Sampling Freq.
%		n 	: Order of Gamma function t^(n-1) (vector)
%		b1      : b1 for  exp(-2*pi*b1*ERB(f))
%		c1	: c1 for  exp(j*2*pi*Fr + c1*ln(t))
%               frat    : frequency ratio. Fr2 = frat*Fp1;
%               b2      : b2 for HP-AF
%               c2      : c2 for HP-AF
%		NfrqRsl : freq. resolution 
%	OUTPUT: cGCresp : struct for cGC response
%                   pGCFrsp: passive gc frq. rsp.  (NumCh*NfrqRsl matrix)
%                   cGCFrsp: compressive gc frq. rsp.  (NumCh*NfrqRsl matrix)
%                   cGCNrmFrsp: Normalized cGCFrsp  (NumCh*NfrqRsl matrix)
%                   ACFrsp : Asym Compnstation Filter frq. rsp.
%                   AsymFunc: Asym Func
%		    freq   : frequency          (1 * NfrqRsl vector)
%                   Fp2    : peak freq.
%                   ValFp2 : peak value
%
function cGCresp = CmprsGCFrsp(Fr1,fs,n,b1,c1,frat,b2,c2,NfrqRsl);

if nargin < 1, help CmprsGCFrsp; return; end;
if nargin < 2, fs = 48000; 	end;
if length(fs) == 0, error('Specify Sampling Frequency'); end;
Fr1 = Fr1(:);
NumCh = length(Fr1);

% Setting Default 
% NOT using SetParam script for stand alone running
% Please check it with GCFBv2_SetParam.m
%
if nargin < 3, n = 4; end; 	
if length(n)  == 1, n = n*ones(NumCh,1); end;
if nargin < 4, b1 = 1.81;	         end; 	
if length(b1) == 1, b1 = b1*ones(NumCh,1); end;
if nargin < 5, c1 = -2.96;               end;	
if length(c1) == 1, c1   = c1*ones(NumCh,1); end;
if nargin < 6, frat = 1;                 end;	
if length(frat) == 1, frat = frat*ones(NumCh,1); end;
%if nargin < 7, b2 = 2.01;               end;	
if nargin < 7, b2 = 2.17;               end;	% debug 8 July 2005
if length(b2) == 1, b2   = b2*ones(NumCh,1); end;
if nargin < 8, c2 = 2.20;               end;	
if length(c2) == 1, c2   = c2*ones(NumCh,1); end;
if nargin < 9, NfrqRsl  = 1024;                end;

[pGCFrsp,freq] = GammaChirpFrsp(Fr1,fs,n,b1,c1,0,NfrqRsl);
Fp1 = Fr2Fpeak(n,b1,c1,Fr1);
Fr2 = frat.*Fp1;
[ACFFrsp,freq,AsymFunc] = AsymCmpFrspV2(Fr2,fs,b2,c2,NfrqRsl);
cGCFrsp = pGCFrsp.*AsymFunc;    %% cGCFrsp = pGCFrsp.*ACFFrsp;
[ValFp2 nchFp2] = max(cGCFrsp');
ValFp2 = ValFp2(:);
NormFactFp2 = 1./ValFp2;

%%%   function cGCresp = CmprsGCFrsp(Fr1,fs,n,b1,c1,frat,b2,c2,NfrqRsl);
cGCresp.Fr1        = Fr1;   % including original parameter as well
cGCresp.n          = n;
cGCresp.b1         = b1;
cGCresp.c1         = c1;
cGCresp.frat       = frat;
cGCresp.b2         = b2;
cGCresp.c2         = c2;
cGCresp.NfrqRsl    = NfrqRsl;
cGCresp.pGCFrsp    = pGCFrsp;
cGCresp.cGCFrsp    = cGCFrsp;
cGCresp.cGCNrmFrsp = cGCFrsp .* (NormFactFp2*ones(1,NfrqRsl));
cGCresp.ACFFrsp    = ACFFrsp;
cGCresp.AsymFunc   = AsymFunc;
cGCresp.Fp1        = Fp1;
cGCresp.Fr2        = Fr2;
cGCresp.Fp2        = freq(nchFp2)';
cGCresp.ValFp2     = ValFp2;
cGCresp.NormFctFp2 = NormFactFp2;
cGCresp.freq       = freq;

return

