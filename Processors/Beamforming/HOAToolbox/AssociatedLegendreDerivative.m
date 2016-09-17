function PlDer = AssociatedLegendreDerivative(l,x)
%ASSOCIATEDLEGENDREDERIVATIVE Associated Legendre function derivative.
%
%  PlDer = AssociatedLegendreDerivative(l,x) computes the values of the
%  Schmidt-seminormalized associated Legendre function derivatives Plm'(x)
%  with m = 0, 1, ... l.
%
%  The output format is the same as in MATLAB's legendre.m function.
%
%  Based on Frederik J. Simons' "libbrecht.m" routine (thanks Frederik!)
%
%  References: Masters and Richards-Dinger, 1998 
%              Libbrecht, 1985.
%
%  See also LEGENDRE
%

%  N.Epain, 2010.

% Check input arguments
if nargin < 2
    error('Not enough input arguments.')
end

if numel(l) > 1 || ~isreal(l) || l < 0 || l ~= round(l)
    error('l must be a positive scalar integer.');
end

if ~isreal(x) || max(abs(x)) > 1 || ischar(x)
    error('x must be real and in the range (-1,1).')
end

% Very simple case : l = 0
if l == 0
    PlDer = zeros(1,length(x)) ;
    return
end

% Initializations
x     = x(:)' ;
Pl    = zeros(l+1,length(x)) ;
PlDer = zeros(l+1,length(x)) ;
sint  = sin(acos(x)) ;

% Comment by F.J. Simons :
% Here is the Libbrecht algorithm. Compute starting prefactor
% sqrt(1/factorial(2*N))*factorial(2*N)/(2^N)/factorial(N)
% by computing (1/2)*(3/4)*(5/6)*...*((2l-1)/(2l)):
f1 = 1 ;
for I = 1 : l
    f1 = f1 * (2*I-1) / (2*I) ;
end
f1 = sqrt(f1) ;

% Initial values for m = l
Pl(l+1,:)    = f1 ;
PlDer(l+1,:) = 0 ;

% Comment by F.J. Simons :
% Compute prefactors
m = 1 : l ; 
f2 = sqrt((l+m).*(l-m+1)) ;

% Comment by F.J. Simons :
% For all m downgoing (MRD (1998) Eq. (5))
% Dividing out f2 progressively yields sqrt(1/factorial(2*l))
% Note that you're switching sign here, too; which you need for the
% algorithm but need to undo for the Schmidt harmonics
for m = l : -1 : 1
    Pl(m,:)    = -(sint.*PlDer(m+1,:) + 2*m*x.*Pl(m+1,:)) /f2(m) ;
    PlDer(m,:) = sint.*Pl(m+1,:) * f2(m) ;
end

% Comment by F.J. Simons :
% Now convert back to ordinary spherical harmonics
f3 = 1 ;
for m = 2 : (l+1)
    PlDer(m,:) = (sint.*PlDer(m,:)+(m-1)*x.*Pl(m,:)) .* f3 ;
    f3 = f3 .* sint ;
end

% Normalization
PlDer(2:end,:) = PlDer(2:end,:)*sqrt(2) ;
for m = l : -2 : 1
    PlDer(m,:) = -PlDer(m,:) ;
end

