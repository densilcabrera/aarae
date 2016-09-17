function g  = Hoa2SpkMaxrEWeightings(order,dimOpt)

if nargin < 2
    dimOpt = '3d' ;
end

switch order
    case 0
        
        g = 1 ;
        
    otherwise
        
        switch dimOpt
            
            case '3d'
                
                nbPts = 3*order + 10 ;
                xMin = .5 ;
                xMax = 1 ;
                x = (1:nbPts) * (xMax-xMin)/nbPts + xMin ;
                p = legendre(order+1,x) ;
                dx = 1 ;
                while dx > 1e-7
                    ind = find((p(1,1:end-1)<=0)&(p(1,2:end)>0),1,'last') ;
                    dx = x(ind+1) - x(ind) ;
                    dp = p(1,ind+1) - p(1,ind) ;
                    xc = ( x(ind)*p(1,ind+1) - x(ind+1)*p(1,ind) ) / dp ;
                    dx2 = max(xc-x(ind),x(ind+1)-xc) ;
                    xMin = xc - dx2 ;
                    xMax = xc + dx2 ;
                    x = (1:nbPts) * (xMax-xMin)/nbPts + xMin ;
                    p = legendre(order+1,x) ;
                end
                rE = interp1(p(1,:),x,0,'cubic') ;
                
                g = zeros(order+1,1) ;
                g(1) = 1 ;
                g(2) = rE ;
                for m = 1 : order-1
                    g(m+2) = ( (2*m+1)*rE*g(m+1) - m*g(m) ) / (m+1) ;
                end
                
            case '2d'
                
                g = cos( (0:order)'*pi / (2*order+2) ) ;
                
        end
        
end
