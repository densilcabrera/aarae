%
%	Frequency Response of GammaChirp
%	Toshio IRINO
%	Original Version: 30 Sept 98
%	Modified:          7 June 2004 (removing transpose when NumCh==1)
%
%    function [AmpFrsp, freq, Fpeak, GrpDly, PhsFrsp] ...
%       = GammaChirpFrsp(Frs,SR,OrderG,CoefERBw,CoefC,Phase,NfrqRsl);
%	INPUT : Frs	: Resonance Freq. (vector)
%		SR 	: Sampling Freq.
%		OrderG 	: Order of Gamma function t^(OrderG-1) (vector)
%		CoefERBw: Coeficient -> exp(-2*pi*CoefERBw*ERB(f))
%		CoefC	: Coeficient -> exp(j*2*pi*Fr + CoefC*ln(t))
%		Phase	: Start Phase (0 ~ 2*pi)
%		NfrqRsl : freq. resolution 
%	OUTPUT: AmpFrsp	: abs(Response)      (NumCh*NfrqRsl matrix)
%		freq	: frequency          (1 * NfrqRsl vector)
%		Fpeak	: Peak frequency     (NumCh * 1 vector)
%               GrpDly  : Group delay        (NumCh*NfrqRsl matrix)
%               PhsFrsp : angle(Response);   (NumCh*NfrqRsl matrix)
%
function [AmpFrsp, freq, Fpeak, GrpDly, PhsFrsp ] ...
    = GammaChirpFrsp(Frs,SR,OrderG,CoefERBw,CoefC,Phase,NfrqRsl);

if nargin < 1, help GammaChirpFrsp; return; end;
if nargin < 2, SR = 48000; 	end;
if length(SR) == 0, error('Specify Sampling Frequency'); end;
Frs = Frs(:);
NumCh = length(Frs);

if nargin < 3, OrderG	= 4;                    end; 	% Default GammaTone
if length(OrderG)   == 1, OrderG = OrderG*ones(NumCh,1); end;
if nargin < 4, CoefERBw = 1.019;	        end; 	% Default GammaTone
if length(CoefERBw) == 1, CoefERBw = CoefERBw*ones(NumCh,1); end;
if nargin < 5, CoefC	= 0;	                end;	% Default GammaTone
if length(CoefC)    == 1, CoefC   = CoefC*ones(NumCh,1); end;
if nargin < 6, Phase = [];                      end;
if length(Phase)==0,  Phase = zeros(NumCh,1);	end;	% Default GammaTone
if nargin < 7, NfrqRsl  = 1024;                end;
if NfrqRsl < 256, help GammaChirpFrsp;  error('NfrqRsl < 256'); end;

[ERBrate ERBw] = Freq2ERB(Frs(:));
freq = [0:NfrqRsl-1]/NfrqRsl*SR/2;

one1 = ones(1,NfrqRsl);
bh = ( CoefERBw(:).*ERBw(:) ) * one1;
fd = ones(NumCh,1)*freq(:)' - Frs(:)*one1;
cn = ( CoefC(:)./OrderG(:) ) * one1;
n  = OrderG(:) *one1;
c  = CoefC(:)  *one1;
Phase = Phase(:) *one1;

%%% Analytic form (normalized at Fpeak) %%%
AmpFrsp = ( (1 + cn.^2)./(1 + (fd./bh).^2) ).^(n/2) ...
      .* exp(  c .* (atan(fd./bh) - atan(cn)));

if nargout > 2,
  Fpeak = Frs + CoefERBw.*ERBw.*CoefC./OrderG;
  if nargout > 3,
    GrpDly = 1/(2*pi)*(n.*bh + c.*fd)./(bh.^2 + fd.^2);
    if nargout > 4,
      PhsFrsp= -n.*atan(fd./bh) -c./2.*log((2*pi*bh).^2 +(2*pi*fd).^2 )+Phase;
    end;
  end;
end;

return



%%% No use anymore since it is somewhat confusing  7 June 2004 %%%%%%%%%%%%

if NumCh == 1, % to let the results in column vector (same as freqz)
  AmpFrsp = AmpFrsp(:); 
  freq=freq(:); 
  if nargout > 3, 
    GrpDly =  GrpDly(:);
    if nargout > 4,
      PhsFrsp =  PhsFrsp(:);
    end;
  end;
end; 

return


% Fbw   = sqrt( -1 + 2.^(1./OrderG)).*CoefERBw.*ERBw; % bandwidth when c=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OLD
% val = 2*pi*(bh + i*fd);
% arg = angle(val);
% FrspAna = 1./(abs(val)).^OrderG(nch) .* exp(CoefC(nch)* arg);
% FrspAna = FrspAna/max(FrspAna);
  

if 0
  
for nch = 1:NumCh
  bh = CoefERBw(nch)*ERBw(nch);
  fd = freq(:)' - Frs(nch);
  cn = CoefC(nch)/OrderG(nch);
  %%% Analytic form (normalized at Fpeak) %%%
  FrspAna = (sqrt( (1+cn^2)./(1+(fd/bh).^2) ) ).^OrderG(nch) ...
      .* exp(  CoefC(nch) .* (atan(fd/bh) - atan(cn)));
  %  Frsp(nch, 1:NfrqRsl) = FrspAna;
  sqrt(mean( ( Frsp(nch, 1:NfrqRsl) - FrspAna ).^2) )
  GrpDly(nch,1:NfrqRsl) = ...
      1/(2*pi)*(OrderG(nch)*bh + CoefC(nch)*fd)./(bh.^2 + fd.^2);
end;

end;


