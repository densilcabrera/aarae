function [tspSeq tsp] = GenerateTSPSequence(reps, N, impalign, lincirc)

%   GenerateTSPSequence
%
%   [tspSeq tsp] = GenerateTSPSequence([reps], [N], [impalign])
%
%   Generates a TSP sequence of order N with optional alignment impulse
%
%   tsp = GenerateTSPSequence(reps, N, impalign)
%
%   Inputs:
%       reps:           (Optional) No. of Repetitions (default 4)
%       N:              (Optional) TSP order, length 2^N (defualt 14)
%       impalign:       (Optional) 1 to include alignment impulse, 0 to
%                       exclude
%       lincirc:        (Optional) 'lin' uses linear convolution (with 2^N
%                       pause between bursts, 'circ' uses circular
%                       convolution (with no burst but 2^N zeros appended
%                       at end).
%   Outputs:
%       tspSeq:         The generated TSP sequence
%       tsp:            One TSP burst
%
%   References:
%       Y. Suzuki, F. Asani, H.-Y. Kim and T. Sone, "An optimum
%       computer-generated pulse signal suitable for the measurement of
%       very long impulse responses," J. Acoust. Soc. Am. Vol. 97(2), pp.
%       1119-1123, 1995.
%
%**************************************************************************
% Author:           M. R. P. Thomas
% Date:             17 Aug 2009
% Last modified:    17 Aug 2009
%**************************************************************************

error(nargchk(0,4,nargin));
if (nargin < 4)
    lincirc = 'lin';
end
if (nargin < 3)
    impalign = 1;
end
if (nargin <2 )
    N=14; 
end
if (nargin < 1)
    reps=4;
end

if ~(strcmp(lincirc,'lin') || strcmp(lincirc,'circ'))
    error('Argument strcmp must be either ''lin'' or ''circ''');
end

tspSeq = [];
if(impalign)
    starttoimpulse = 10000;
    impulsetoburst = 2^N;

    impulse = [1; zeros(impulsetoburst-1,1)];
    tspSeq = [zeros(starttoimpulse,1); impulse];
end

tsp = GenerateTSP(N,1);         % Normalize on.

for j=1:1:reps
   tspSeq = [tspSeq; tsp];
   if(strcmp(lincirc,'lin'))
        tspSeq = [tspSeq; zeros(2^N,1)];
   end
end

if(strcmp(lincirc,'circ'))
    tspSeq = [tspSeq; zeros(2^N,1)];
end