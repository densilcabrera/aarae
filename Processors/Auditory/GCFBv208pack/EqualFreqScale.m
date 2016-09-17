%
%       Calculation of Equal Frequency scale on ERB/Mel/Log/Linear scale
%       Toshio IRINO
%       5 Oct 2001
%
%    function [Frs, WFvals] = EqualFreqScale(NameScale,NumCh,RangeFreq),
%      INPUT:  NameScale: 'ERB', 'mel', 'log', 'linear'
%	       NumCh:  Number of channels
%              RangeFreq: Frequency Range
%      OUTPUT: Frs:    Fr vector
%     	       WFval: Warped Freq. value
%
function [Frs, WFvals] = EqualFreqScale(NameScale,NumCh,RangeFreq),

if nargin <3, help EqualFreqScale; end;
if diff(RangeFreq) < 0, 
	help  EqualFreqScale; 
	error('RangeFreq(1) should be less than RangeFreq(2).');  
end;

if strcmp(lower(NameScale),'linear') == 1, 
	RangeWF = RangeFreq;
	dWF = diff(RangeWF)/(NumCh-1);
	WFvals = RangeWF(1):dWF:RangeWF(2)+eps*1000;
	Frs = WFvals;

elseif 	strcmp(lower(NameScale),'mel') == 1,
	RangeWF = Freq2Mel(RangeFreq);
	dWF = diff(RangeWF)/(NumCh-1);
	WFvals = RangeWF(1):dWF:RangeWF(2)+eps*1000;
	Frs = Mel2Freq(WFvals);

elseif 	strcmp(lower(NameScale),'erb') == 1,
	RangeWF = Freq2ERB(RangeFreq);
	dWF = diff(RangeWF)/(NumCh-1);
	WFvals = RangeWF(1):dWF:RangeWF(2)+eps*1000;
	Frs = ERB2Freq(WFvals);

elseif  strcmp(lower(NameScale),'log') == 1, 
	if min(RangeFreq) < 50,
	  disp([ 'min(RangeFreq) < 50. Replaced by 50.' ]);
	  RangeFreq(1) = 50;
	end;
	RangeWF = log10(RangeFreq);
	dWF = diff(RangeWF)/(NumCh-1);
	WFvals = RangeWF(1):dWF:RangeWF(2)+eps*1000;
	Frs = 10.^WFvals;

else
	help EqualFreqScale;
	error('Specify NameScale correctly');
end;

