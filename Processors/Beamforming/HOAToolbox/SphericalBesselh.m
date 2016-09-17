function y = SphericalBesselh(nu,k,x)

if (isvector(nu)&&isvector(x)) && (size(nu,1)~=size(x,1))
    
    y = repmat(sqrt(pi/2./x(:)),1,length(nu)) .* besselh(nu(:)'+1/2,k,x(:)) ;
    
else
    
    y = sqrt(pi/2./x) .* besselh(nu+1/2,k,x) ;    
    
end

end