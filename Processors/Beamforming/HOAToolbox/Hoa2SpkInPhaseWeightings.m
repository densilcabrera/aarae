function g  = Hoa2SpkInPhaseWeightings(order,dimOpt)

if nargin < 2
    dimOpt = '3d' ;
end

switch dimOpt
    
    case '3d'
        
        g = gamma(order+1)*gamma(order+2) ...
            ./ gamma(order+2+(0:order).') ...
            ./ gamma(order+1-(0:order).') ;
        
    case '2d'
        
        g = gamma(order+1)^2 ...
            ./ gamma(order+1+(0:order).') ...
            ./ gamma(order+1-(0:order).') ;

end

