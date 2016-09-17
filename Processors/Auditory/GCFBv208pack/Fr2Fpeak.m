%
%    Estimate fr from fpeak 
%    Toshio IRINO
%    10 June 98
%
% function [fpeak, ERBw] = Fr2Fpeak(n,b,c,fr)
%    INPUT:  n,b,c : gammachirp param.
%            fr    : fr
%    OUTPUT: fpeak : peak freq.
%            ERBw  : ERBw(fr)
%
function [fpeak, ERBw] = Fr2Fpeak(n,b,c,fr)

if nargin < 4, help Fr2Fpeak; end;

n = n(:);
b = b(:);
c = c(:);
fr = fr(:);

[dummy ERBw] = Freq2ERB(fr);
fpeak = fr + c.*ERBw.*b./n;
