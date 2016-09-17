function y=ErbRateToHz(x)
% Convert ERB rate to Hz.
% 
%   y = ErbRateToHz(x)
% 
%   y = ErbRateToHz(x) converts the ERB number x to the
%   eqivalent frequency y (in Hz).
% 
% See also HZTOERBRATE.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

y=(10.^(x/21.4)-1)/4.37e-3;

% [EOF]
