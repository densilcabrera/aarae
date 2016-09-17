function clrMap = nicolor(N,periodic)
%
%NICOLOR    Color maps optimised for B&W printing and color blindness.
%
%  clrMap = nicolor(N,periodic)
% 
%  Input: - N: the number of colors.
%         - (optional) periodic: if set to true, the output is a colormap
%         designed for plotting periodic functions (e.g. the phase of a
%         complex field).
% 
%  Output: - clrMap is the Nx3 array of the colormap's RGB values.
%
%  See also JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.
%
%  N. Epain, 2011
%
%  Reference:
%   - "Rainbow Color Map (Still) Considered Harmful", D. Borland and 
%   R.M. Taylor, Computer Graphics and Applications (IEEE), March-April 
%   2007, 27(2), pp. 14-17.

% Default input arguments
if nargin < 2
   periodic = false ;
end
if (nargin < 1) || isempty(N)
   N = size(get(gcf,'colormap'),1) ;
end

% Base colors
if periodic == true
    clrMap = [ 0.00 0.01 0.68
               0.00 0.18 0.70
               0.00 0.31 0.61
               0.00 0.40 0.53
               0.00 0.49 0.49
               0.00 0.57 0.43
               0.00 0.65 0.32
               0.00 0.74 0.19
               0.00 0.82 0.00
               0.23 0.90 0.00
               0.48 0.96 0.00
               0.75 1.00 0.00
               1.00 1.00 0.92
               1.00 0.92 0.67
               1.00 0.82 0.63
               1.00 0.71 0.61
               1.00 0.59 0.59
               1.00 0.45 0.59
               1.00 0.26 0.63
               0.94 0.00 0.71
               0.77 0.00 0.77
               0.61 0.00 0.81
               0.42 0.00 0.83
               0.20 0.00 0.80
               0.00 0.00 0.68 ] ;
           
else
    clrMap = [ 0.00 0.13 0.07
               0.00 0.18 0.11
               0.00 0.23 0.20
               0.11 0.26 0.39
               0.21 0.28 0.56
               0.38 0.27 0.71
               0.52 0.26 0.83
               0.68 0.26 0.85
               0.80 0.30 0.84
               0.93 0.36 0.69
               0.98 0.47 0.59
               1.00 0.57 0.55
               1.00 0.67 0.50
               1.00 0.77 0.42
               1.00 0.86 0.40
               1.00 0.95 0.40
               1.00 1.00 1.00 ] ;
end

% Interpolate to obtain the required number of colors
clrMap = interp1(linspace(0,1,size(clrMap,1)),clrMap,linspace(0,1,N)) ;


