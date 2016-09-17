%
%	ERB -> Frequency and ERBwidth (Glasberg and Moore, 1990)
%	Toshio IRINO
%	11 Mar. 1998
%
%	function [cf, ERBwidth] = ERB2Freq(ERBrate)
%	INPUT	ERBrate:  ERB rate
%	OUTPUT  cf:       Center frequency (Hz)
%		ERBwidth: ERB bandwidth (Hz)
%
%	Ref: Glasberg and Moore : Hearing Research, 47 (1990), 103-138
%            For different formulae (years), see ERB2FreqYear.m
%       
%
function [cf, ERBwidth] = ERB2Freq(ERBrate),

if nargin < 1,  help ERB2Freq; end;

cf = (10.^(ERBrate./21.4) - 1)./4.37 * 1000;
ERBwidth = 24.7.*(4.37*cf/1000 + 1);



