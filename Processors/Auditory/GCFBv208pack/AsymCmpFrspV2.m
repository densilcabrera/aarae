%
%  AsymCmpFrsp Version 2
%  Amplitude Spectrum of Asymmetric Compensation IIR filter for the gammachirp 
%  corresponding to MakeAsymCmpFiltersV2.m
%
%  Toshio Irino
%  Original : 14 Apr. 99 (Version1)
%  Modified : 11 Jun 2004
%  Modified :  7 Jul 2005 % NfrqRsl
%
%  function [ACFFrsp, freq, AsymFunc]
%            = AsymCmpFrspV2(Frs,fs,b,c,NfrqRsl,NumFilt)
%
%  INPUT:    fs: Sampling frequency
%            Frs: array of the center frequencies
%            b : array or scalar of a bandwidth coefficient 
%            c : array or scalar of asymmetric parameters 
%            NfrqRsl: freq. resolution
%            NumFilt: Number of 2nd-order filters default 4
%  OUTPUT:   ACFFrsp: abs(Frsp of ACF)  (NumCh*NfrqRsl matrix)
%            freq: freq.                    (1*NfrqRsl vector)
%            AsymFunc:    Original Asymmetric Function (NumCh*NfrqRsl matrix)
%
%
function [ACFFrsp,freq,AsymFunc] = AsymCmpFrspV2(Frs,fs,b,c,NfrqRsl,NumFilt)

if nargin < 1, help AsymCmpFrspV2, end;
if nargin < 5, NfrqRsl = []; end;
if length(NfrqRsl) == 0,  NfrqRsl = 1024; end;
if nargin < 6, NumFilt = []; end;
if length(NumFilt) == 0,  NumFilt = 4; end; % default
if NumFilt ~= 4, error('NumFilter should be 4.'); end;

Frs = Frs(:);
b = b(:);
c = c(:);
NumCh = length(Frs);

SwCoef = 0;  % self consitency
%SwCoef = 1; % referece to MakeAsymCmpFiltersV2

if SwCoef == 0
% New Coefficients. NumFilter = 4; See [1]
  p0 = 2;
  p1 = 1.7818 .* (1 - 0.0791*b) .* (1 - 0.1655*abs(c));
  p2 = 0.5689 .* (1 - 0.1620*b) .* (1 - 0.0857*abs(c));
  p3 = 0.2523 .* (1 - 0.0244*b) .* (1 + 0.0574*abs(c));
  p4 = 1.0724;
else
  
  ACFcoef = MakeAsymCmpFiltersV2(fs,Frs,b,c);
end;
  
[dummy ERBw] = Freq2ERB(Frs);
freq = (0:NfrqRsl-1)/NfrqRsl*fs/2;
ACFFrsp = ones(NumCh,NfrqRsl);
freq2 = [ones(NumCh,1)*freq, Frs];


for Nfilt = 1:NumFilt

  if SwCoef == 0,
    r  = exp(-p1.*(p0./p4).^(Nfilt-1) .* 2 .* pi .*b .*ERBw /fs); 
    delfr = (p0*p4)^(Nfilt-1).*p2.*c.*b.*ERBw; 
    phi = 2*pi*max(Frs + delfr,0)/fs;
    psy = 2*pi*max(Frs - delfr,0)/fs;
    fn = Frs;    
    ap = [ones(NumCh,1), -2*r.*cos(phi),  r.^2];
    bz = [ones(NumCh,1), -2*r.*cos(psy),  r.^2]; 
  else
    ap = ACFcoef.ap(:,:,Nfilt);
    bz = ACFcoef.bz(:,:,Nfilt);
  end;

  cs1 = cos(2*pi*freq2/fs);
  cs2 = cos(4*pi*freq2/fs);

  bzz0 = (bz(:,1).^2 + bz(:,2).^2 + bz(:,3).^2 )*ones(1,NfrqRsl+1);
  bzz1 = (2* bz(:,2).*( bz(:,1) + bz(:,3)))*ones(1,NfrqRsl+1);
  bzz2 = (2* bz(:,1).*bz(:,3))*ones(1,NfrqRsl+1);
  hb = bzz0 + bzz1.* cs1 + bzz2.* cs2;

  app0 = (ap(:,1).^2 + ap(:,2).^2 + ap(:,3).^2)*ones(1,NfrqRsl+1);
  app1 = (2* ap(:,2).*( ap(:,1) + ap(:,3)))*ones(1,NfrqRsl+1);
  app2 = (2* ap(:,1).*ap(:,3))*ones(1,NfrqRsl+1);
  ha = app0 + app1.*cs1 + app2.*cs2;

  H = sqrt(hb./ha);
  Hnorm = H(:,NfrqRsl+1)*ones(1,NfrqRsl);    % Noimalization by fn value

  ACFFrsp = ACFFrsp .* H(:,1:NfrqRsl)./ Hnorm;

end;

%%% original Asymmetric Function without shift centering 
fd = (ones(NumCh,1)*freq - Frs*ones(1,NfrqRsl));
be = (b.*ERBw)*ones(1,NfrqRsl);
cc = (c.*ones(NumCh,1))*ones(1,NfrqRsl); % in case when c is scalar
AsymFunc = exp(cc.*atan2(fd,be));

return
