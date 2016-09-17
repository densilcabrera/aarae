function JDer = BesseljDerivative(nu,z)
%BESSELJDERIVATIVE Derivative of the Bessel function of the first kind
%
%  JDer = BesseljDerivative(nu,k,z), for k = 1 or 2, computes the
%  derivative of the Bessel functions J_nu'(z) for each element of the
%  complex array Z.
%
%  See also BESSELJ
%

%  N.Epain, 2010.

JDer = 1/2 * ( besselj(nu-1,z) - besselj(nu+1,z) ) ;

JDer(z==0) = 1/3 * (nu(z==0)==1) ;
                           
end