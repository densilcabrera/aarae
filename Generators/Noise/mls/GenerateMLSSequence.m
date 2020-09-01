function [sequence mls] = GenerateMLSSequence(repetitions, N, impalign)

% GenerateMLSSequence
%
% Generates an impulse followed by a burst of MLS
%
%   [sequence mls] = GenerateMLSSequence([repetitions],[N],[impalign])
%
%   Inputs:
%       repetitions:    Number of repetitions of each amplitude (default 4)
%       N:              MLS order, length 2^N-1 (default 14)
%       impalign:       1 to include alignment impulse, 0 to exclude
%   Outputs:
%       sequence:       The generated sequence.
%       mls:            One MLS burst
%
% References:

%**************************************************************************
% Author:           M. R. P. Thomas 
% Date:             21 Feb 2007
% Last modified:    09 May 2007
%**************************************************************************

narginchk(0,3);
if (nargin < 3)
    impalign = 1;
end
if (nargin <2 )
    N=14; 
end
if (nargin < 1)
    repetitions=4;
end

if (repetitions <2)
   error('Repetitions must be at least 2'); 
end

sequence = [];
if(impalign)
    starttoimpulse = 10000;
    impulsetoburst = 2^N;

    impulse = [1 zeros(1,impulsetoburst-1)];
    sequence = [zeros(1,starttoimpulse) impulse];
end

mls = GenerateMLS(N,1); % for same MLS every time
%mls = GenerateMLS(N,0); % for random seed

for j=1:1:repetitions
   sequence = [sequence mls];
end

sequence = [sequence zeros(1,2^N-1)];

sequence = sequence(:); % Make column vector