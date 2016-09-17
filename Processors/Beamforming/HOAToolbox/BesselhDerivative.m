function HDer = BesselhDerivative(nu,k,z)
%BESSELHDERIVATIVE Derivative of the Bessel function of the third kind
%(Hankel function).
%
%  HDer = BesselhDerivative(nu,k,z), for k = 1 or 2, computes the
%  derivative of the Hankel functions H1_nu'(z) or H2_nu'(z) for each
%  element of the complex array Z.
%
%  See also BESSELH

%  N.Epain, 2010

HDer = 1/2 * ( besselh(nu-1,k,z) - besselh(nu+1,k,z) ) ;
                           
end