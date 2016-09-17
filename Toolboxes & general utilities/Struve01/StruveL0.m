function fun=StruveL0(z)
%
% StruveL0 calculates the function StruveL0 for complex argument z
%
% Author : T.P. Theodoulidis
% Date   : 28 June 2012
%
% Arguments
% z : can be scalar, vector, matrix
% 
% External routines called: StruveH0
%
fun=-1i*StruveH0(1i*z);
%