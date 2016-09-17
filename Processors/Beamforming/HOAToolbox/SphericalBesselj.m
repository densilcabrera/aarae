function y = SphericalBesselj(nu,x)

if (isvector(nu)&&isvector(x)) && (size(nu,1)~=size(x,1))
    
    y = repmat(sqrt(pi/2./x(:)),1,length(nu)) .* besselj(nu(:)'+1/2,x(:)) ;

else
    
    y = sqrt(pi/2./x) .* besselj(nu+1/2,x) ;
    
end

y(x==0) = (nu(x==0)==0) ;

end