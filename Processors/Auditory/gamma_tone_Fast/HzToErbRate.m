function y=HzToErbRate(x)
% Convert Hz to ERB rate
%
%   y=HzToErbRate(x)
%
%   y = HzToErbRate(x) converts the frequency X (in Hz) to
%   the eqivalent ERB number.
% 
%   See also ERBRATETOHZ, MAKEERBCFS.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

y=(21.4*log10(4.37e-3*x+1));

% [EOF]
