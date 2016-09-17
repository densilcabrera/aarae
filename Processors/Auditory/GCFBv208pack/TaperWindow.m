%
%	Taper Window Generator for signal onset/offset
%	IRINO T.
%       Created:    7 Apr. 1993
%       Modified:  29 Aug. 1996
%       Modified:   2 Jan. 2004  (overlap-adding with 'hanning' == flat )
%                                (see bottom)
%
%	function [TaperWin, TypeTaper] = ...
%		TaperWindow(LenWin,TypeTaper,LenTaper,RangeSigma,SwPlot)
%	INPUT	LenWin    : Length of Window (Number of points)
%		TypeTaper : Type of Taper (KeyWords of 3 letters)
%		  (Hamming, Hanning (=cosine^2), Blackman, Gauss, Line)
%		LenTaper  : Length of Taper  (Number of points)
%		RangeSigma: Range in Sigma (default: 3) for Gauss
%		SwPlot    : 0) Omit plotting,  1) Plot Taper
%	OUTPUT  TaperWin  : Taper Window Points (max==1);
%		TypeTaper : Type of Taper (Full Name)
%
function [TaperWin, TypeTaper] = ...
	TaperWindow(LenWin,TypeTaper,LenTaper,RangeSigma,SwPlot)

if nargin < 2,
help TaperWindow
error([ 'Specify Type of Taper : ' ...
	' Hamming, Hanning (=cosine^2), Blackman, Gauss, Line ']);
%TaperWin = ones(1,LenWin);
%return;
end;

if nargin < 3, LenTaper = fix(LenWin/2); end;
if nargin < 4, RangeSigma = 3; end;

if  (LenTaper*2+1) >= LenWin, 
	disp('Caution (TaperWindow.m) : No flat part. ');
	if LenTaper ~= fix(LenWin/2),
	disp('Caution (TaperWindow.m) : LenTaper <-- fix(LenWin/2)');
	end;
	LenTaper = fix(LenWin/2); 
end;

if nargin < 5, SwPlot = 0; end;	% changing default Swplot 29 Aug. 96

%TypeTaper = lower(TypeTaper(1:3));

% length(Taper) = LenTaper*2  --> LenTaper*2+1    changed 2 Jan 2004

if	upper(TypeTaper(1:3)) == 'HAM', 
	Taper = hamming(LenTaper*2+1)';                
	TypeTaper = 'Hamming';
elseif	upper(TypeTaper(1:3)) == 'HAN' | upper(TypeTaper(1:3)) == 'COS', 
	Taper = hanning(LenTaper*2+1)';
	TypeTaper = 'Hanning/Cosine'; 
elseif	upper(TypeTaper(1:3)) == 'BLA', 
	Taper = blackman(LenTaper*2+1)';
	TypeTaper = 'Blackman';
elseif	upper(TypeTaper(1:3)) == 'GAU',  
	if length(RangeSigma) == 0, RangeSigma = 3; end;
	% nn = -LenTaper+0.5:1:LenTaper-0.5; %    old Jun02
	nn = -LenTaper:1:LenTaper;           %    4 Jan 04
	Taper = exp(-(RangeSigma*nn/LenTaper).^2 /2);
	TypeTaper = 'Gauss';
else	%Taper = [1:LenTaper LenTaper:-1:1]/LenTaper;  % 'line',  old Jun02
        Taper = [1:LenTaper, LenTaper+1, LenTaper:-1:1]/(LenTaper+1);         
	TypeTaper = 'Line';                            %    4 Jan 04
end;

%plot(Taper)
%size(Taper);
LenTaper = fix(LenTaper);
TaperWin = [	Taper(1:LenTaper) ones(1,LenWin-LenTaper*2) ...
		Taper((LenTaper+2):(LenTaper*2+1))];   
                % skip the center value, changed 4 Jan 2004

if SwPlot == 1,

plot(TaperWin)
xlabel('Points');
ylabel('Amplitude');
title(['TypeTaper = ' TypeTaper] );

end;

return

%%%%%%%%%%%%%%%%%%%%%%

if 0 % checking the flatness for overlap-adding 

nDur = 100;
nOvrlp = 7;
TaperWin = TaperWindow(nDur+nOvrlp*2,'han',nOvrlp);

aa = TaperWin;
for nk = 1:3
aa = [aa, zeros(1,length(TaperWin)-nOvrlp)] ...
   + [zeros(1,length(aa)-nOvrlp), TaperWin];
end;

plot(aa)

end;
