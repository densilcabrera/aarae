function fun=StruveL1(z)
%
% StruveL0 calculates the function StruveL1 for complex argument z
%
% Author : T.P. Theodoulidis
% Date   : 28 June 2012
%
% Arguments
% z : can be scalar, vector, matrix
% 
% External routines called: StruveH1
%
fun=-StruveH1(1i*z);
%