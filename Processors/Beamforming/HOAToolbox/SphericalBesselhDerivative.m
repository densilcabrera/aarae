function y = SphericalBesselhDerivative(nu,k,x)

if (isvector(nu)&&isvector(x)) && (size(nu,1)~=size(x,1))
    
    y = 1/2 * ( SphericalBesselh(nu-1,k,x) - SphericalBesselh(nu+1,k,x) ...
          - repmat(1./x(:),1,length(nu))  .* SphericalBesselh(nu,k,x) ) ;
      
else
    
    y = 1/2 * ( SphericalBesselh(nu-1,k,x) - SphericalBesselh(nu+1,k,x) ...
                                 - (1./x) .* SphericalBesselh(nu,k,x) ) ;

end

end