%
%  MakeAsymCmpFilters Version 2
%  Computes the coefficients for a bank of Asymmetric Compensation Filters 
%  This is a modified version to fix the round off problem at low freqs.
%  Use this with ACFilterBank.m
%  See also AsymCmpFrspV2.m for frequency response
%  See footnote in this m-file.
%
%  Toshio Irino 
%  Created: 18 May 2004
%  Modified: 11 Jun 2004
%  Modified:  9 Dec 2008  (removing b = b(:) etc. for speed up)
%
%  function ACFcoef=MakeAsymCmpFiltersV2(fs,fr,b,c)
%  INPUT:    fs: Sampling frequency
%            Frs: array of the center frequencies, Frs(:)
%            b : array or scalar of a bandwidth coefficient, b(:)
%            c : array or scalar of asymmetric parameters, c(:)
%
%  OUTPUT:   ACFcoef
%               fs : Sampling frequnecy
%               bz : MA coefficients  (NumCh*3*NumFilt)
%               ap : AR coefficients  (NumCh*3*NumFilt)
%
function ACFcoef = MakeAsymCmpFiltersV2(fs,Frs,b,c)

if nargin < 4, help MakeAsymCmpFiltersV2;  end;

[NumCh LenFrs] = size(Frs);
if LenFrs > 1
  error('Frs should be a column vector Frs(:).');
end;
%% Frs = Frs(:); % transpose -> column vector % It takes time!
%% b   = b(:);
%% c   = c(:);

[dummy ERBw] = Freq2ERB(Frs);
ACFcoef.fs = fs;

% New Coefficients. See [1]
  NumFilt = 4;
  p0 = 2;
  p1 = 1.7818 .* (1 - 0.0791*b) .* (1 - 0.1655*abs(c));
  p2 = 0.5689 .* (1 - 0.1620*b) .* (1 - 0.0857*abs(c));
  p3 = 0.2523 .* (1 - 0.0244*b) .* (1 + 0.0574*abs(c));
  p4 = 1.0724;
%
if NumFilt > 4, error('NumFilt > 4'); end;
ACFcoef.ap = zeros(NumCh,3,NumFilt);
ACFcoef.bz = zeros(NumCh,3,NumFilt);

for Nfilt = 1:NumFilt,

  r  = exp(-p1.*(p0./p4).^(Nfilt-1) * 2*pi .*b .*ERBw /fs); 
  delFrs = (p0*p4)^(Nfilt-1).*p2.*c.*b.*ERBw; 
  phi = 2*pi*max(Frs + delFrs,0)/fs;
  psy = 2*pi*max(Frs - delFrs,0)/fs;
  fn = Frs;    % see [2]

  %% Second order filter %%
  ap = [ones(size(r)), -2*r.*cos(phi),  r.^2];
  bz = [ones(size(r)), -2*r.*cos(psy),  r.^2];

  vwr = exp(i*2*pi*fn/fs);
  vwrs = [ones(size(vwr)), vwr, vwr.^2];
  nrm = abs( sum( (vwrs.*ap)')' ./ sum( (vwrs.*bz)')');
  bz = bz .* (nrm*ones(1,3));

  ACFcoef.ap(:,:,Nfilt) = ap;
  ACFcoef.bz(:,:,Nfilt) = bz;
end;

return

%%%%%%%%%%%%%%%%%%%%%%%%

% Note:
%
% [1] Ref for p1-p4:
% Unoki,M , Irino,T. , and Patterson, R.D. , "Improvement of an
% IIR asymmetric compensation gammachirp filter," Acost. Sci. &
% Tech. (ed. by the Acoustical Society of Japan ), 22 (6), pp. 426-430,
% Nov. 2001.
%
% [2] Conventional setting was removed.
%  fn = Frs + Nfilt* p3 .*c .*b .*ERBw/n;
%
%  This frequency fn is for normalizing GC(=GT*Hacf) filter to be unity 
%  at the peak, frequnecy. But now we use Hacf as a highpass filter as well.
%  cGC = pGC *Hacf. In this case, this normalization is useless.
%  So, it was set as the gain at Frs is unity.  (4. Jun 2004 )
%  
% [3] Removed
%  ACFcoef.fn(:,nff) = fn;
%  n : scalar of order t^(n-1) % used only in normalization 
