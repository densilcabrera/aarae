function y = SphericalBesseljDerivative(nu,x)

if (isvector(nu)&&isvector(x)) && (size(nu,1)~=size(x,1))
    
    y = 1/2 .* ( SphericalBesselj(nu-1,x) - SphericalBesselj(nu+1,x) ...
          - repmat(1./x(:),1,length(nu)) .* SphericalBesselj(nu,x) ) ;
       
else
   
    y = 1/2 .* ( SphericalBesselj(nu-1,x) - SphericalBesselj(nu+1,x) ...
                                - (1./x) .* SphericalBesselj(nu,x) ) ;
    
end

y(x==0) = 1/3 * (nu(x==0)==1) ;
                           
end