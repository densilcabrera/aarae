function eval = cheval(Ctype,An,xv)
%
% cheval evaluates any one of the four types of Chebyshev series.
% It is a Matlab translation of the Fortran function EVAL
% found in page 20 of:
% Y.L. Luke, Algorithms for the computation of mathematical functions
% Academic Press, 1977, p.20
%
% Author : T.P. Theodoulidis
% Date   : 11 June 2012
%
% Ctype (string): type of Chebyshev polynomial
%                 'regular' T_n(x)
%                 'shifted' T*_n(x)
%                 'even'    T_2n(x)
%                 'odd'     T_2n+1(x)
% An            : vector of Chebyshev coefficients
% z             : argument, can be scalar, vector, matrix
%
switch Ctype
    case {'regular'}
        log1=1; log2=1;
    case {'shifted'}
        log1=1; log2=0;
    case {'even'}
        log1=0; log2=1;
    case {'odd'}
        log1=0; log2=0;
    otherwise
        return;
end
%
x=xv(:);
%
xfac=2*(2*x-1);
if log1 && log2, xfac=2*x; end
if ~log1, xfac=2*(2*x.*x-1); end
n=length(An);
n1=n+1;
Bn=zeros(length(x),n1+1);
%
for j=1:n
    Bn(:,n1-j)=xfac.*Bn(:,n1+1-j)-Bn(:,n1+2-j)+An(n1-j);
end
eval=Bn(:,1)-xfac.*Bn(:,2)/2;
if ~log1 && ~log2, eval=x.*(Bn(:,1)-Bn(:,2)); end
eval=reshape(eval,size(xv));
%
