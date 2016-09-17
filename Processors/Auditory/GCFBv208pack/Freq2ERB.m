%
%	Frequency -> ERB_N-rate and ERB_N-Bandwidth (Glasberg and Moore, 1990)
%	Toshio IRINO
%	11 Mar. 1998
%       Nodified: 26 Jul 2004 (no warning)
%       Nodified: 17 Nov 2006 (modified the comments only. ERB-> ERB_N)
%
%	function [ERBrate, ERBwidth] = Freq2ERB(cf),
%	INPUT	cf:       Center frequency
%	OUTPUT  ERBrate:  ERB_N rate
%		ERBwidth: ERB_N Bandwidth
%
%	Ref: Glasberg and Moore: Hearing Research, 47 (1990), 103-138
%            For different formulae (years), see Freq2ERBYear.m
%
function [ERBrate, ERBwidth] = Freq2ERB(cf),

if nargin < 1,  help Freq2ERB; end;

ERBrate		= 21.4.*log10(4.37*cf/1000+1);
ERBwidth	= 24.7.*(4.37*cf/1000 + 1);

return % no warning

%%% Warning for Freq. Range %%%
cfmin = 50;
cfmax = 12000;
if (min(cf) < cfmin | max(cf) > cfmax)
 disp(['Warning : Min or max frequency exceeds the proper ERB range:']);
 disp(['          ' int2str(cfmin) '(Hz) <= Fc <=  ' int2str(cfmax) '(Hz).']);
end;

%if (min(cf) < 0)
% error(['Min frequency is less than 0.']);
%end;




