function PlotRobinsonProject(SphCoo,Dat,res,axsDat,phaOpt)
%PLOTROBINSONPROJECT  2d plot of data defined on a sphere (worldmap style).
%
%   PlotRobinsonProject(SphCoo,Dat,[res],[axsDat])
%
%   Inputs:
%      - SphCoo [Nx2] is the array of the spherical coordinate values
%        ( SphCoo = [ azm elv ])
%      - Dat [NxM] is the data to be plotted
%      - res (optional) is the resolution of the plot in degrees. The
%        default value is 3.
%      - axsDat (optional) if set to 'both' (default), both the axes and 
%        data are plotted. If set to axes, only the axes are plotted. If
%        set to 'data', only the data are plotted.
%      - phaOpt (optional) if true, the data is considered as a phase and
%        the interpolation is done in a different way. First, exp(1i*Dat)
%        is interpolated, then the phase of the resulting complex numbers
%        is plotted.
%
%   Nicolas Epain, 2011

% Check if the data corresponds to the spherical coordinate values
if size(SphCoo,1) ~= size(Dat,1)
    error([ 'The number of rows in ''SphCoo'' must be equal to ' ...
        'the number of rows in ''Dat''']) ; 
end

% Default "phase" option: false
if ( nargin < 5 ) || isempty(phaOpt)
    phaOpt = false ;
end

% Default "axes/data/both" option: both
if ( nargin < 4 ) || isempty(axsDat)
    axsDat = 'both' ;
end

% Default "resolution" option: 3 deg
if ( nargin < 3 ) || isempty(res)
    res = 3 ;
end

% Plot the axes/data/both ?
switch lower(axsDat)
    case 'both'
        pltAxs = true  ; pltDat = true ;
    case 'axes'  
        pltAxs = true  ; pltDat = false ; 
        res = 30 ;
    case 'data'
        pltAxs = false ; pltDat = true ;
end

% HRIR data
azmInp = SphCoo(:,1) ;
elvInp = SphCoo(:,2) ;

% Number of plots (number of columns in Dat)
nmbPlt = size(Dat,2) ;


%%% ROBINSON PROJECTION MAPPING DATA

% Elevation values
elvMap = (-90:5:90)'*pi/180 ;

% Corresponding parallel lengths / height relative to the equator
% (Source: Wikipedia - Robinson projection)
parMap = [ ...
    0.5322 -1.0000
    0.5722 -0.9761
    0.6213 -0.9394
    0.6732 -0.8936
    0.7186 -0.8435
    0.7597 -0.7903
    0.7986 -0.7346
    0.8350 -0.6769
    0.8679 -0.6176
    0.8962 -0.5571
    0.9216 -0.4958
    0.9427 -0.4340
    0.9600 -0.3720
    0.9730 -0.3100
    0.9822 -0.2480
    0.9900 -0.1860
    0.9954 -0.1240
    0.9986 -0.0620
    1.0000 +0.0000
    0.9986 +0.0620
    0.9954 +0.1240
    0.9900 +0.1860
    0.9822 +0.2480
    0.9730 +0.3100
    0.9600 +0.3720
    0.9427 +0.4340
    0.9216 +0.4958
    0.8962 +0.5571
    0.8679 +0.6176
    0.8350 +0.6769
    0.7986 +0.7346
    0.7597 +0.7903
    0.7186 +0.8435
    0.6732 +0.8936
    0.6213 +0.9394
    0.5722 +0.9761
    0.5322 +1.0000 ] ;
parMap(:,2) = parMap(:,2) * 0.5072/2 ;


%%% INTERPOLATE THE DATA

% Azimuth/elevation values for the plot
azmPlt = (-180:res:180)'*pi/180 ;
elvPlt = (-090:res:090)'*pi/180 ;
[azmPlt,elvPlt] = meshgrid(azmPlt,elvPlt) ;

% Corresponding x/y values on the map
wdtPlt = interp1(elvMap,parMap(:,1),elvPlt) ;
yyyPlt = interp1(elvMap,parMap(:,2),elvPlt) ;
xxxPlt = azmPlt/pi.*wdtPlt/2 ;

% Interpolate the data
if ~phaOpt
    DatInt = SphericalSplineInterp(azmInp,elvInp,Dat,azmPlt(:),elvPlt(:)) ;
else
    DatInt = SphericalSplineInterp(azmInp,elvInp, ...
        exp(1i*Dat),azmPlt(:),elvPlt(:)) ;
    DatInt = angle(DatInt) ;
end


%%% PLOT THE DATA

for I = 1 : nmbPlt
    
    figure('color','white')
    
    % Surf
    if pltDat == true
        surf(xxxPlt,yyyPlt,reshape(DatInt(:,I),size(xxxPlt))) ;
        view(2) ; shading interp ; axis equal ; axis off
    end
    hold on
    
    % Height at which the axes will be plotted
    zzz = max(DatInt(:,I))+.01*max(DatInt(:,I)) ;
    
    % Plot the axes ticks etc
    if pltAxs == true
        
        % Plot the axes
        PlotMy2dAxes(elvMap,parMap,zzz) ;
                
    end
     
end
    

% Sub routine: plot the axes in the 2d case
function PlotMy2dAxes(elvMap,parMap,zzz)

% Plot the parallels and meridians
elvPar = (-090:30:090)'*pi/180 ;
yyyPar = interp1(elvMap,parMap(:,2),elvPar) ;
lngPar = interp1(elvMap,parMap(:,1),elvPar) ;
azmMer = (-180:30:180)'*pi/180 ;
elvMer = (-090:03:090)'*pi/180 ;
yyyMer = interp1(elvMap,parMap(:,2),elvMer) ;
xxxMer = interp1(elvMap,parMap(:,1),elvMer) ;
for J = 2 : length(elvPar)-1
    if yyyPar(J)~=0
        plot3([-.5 lngPar(J)/2], ...
            [yyyPar(J) yyyPar(J)],[zzz zzz], ...
            ':','color',[64 64 64]/255,'linewidth',1) ;
    else
        plot3([-.5 lngPar(J)/2], ...
            [yyyPar(J) yyyPar(J)],[zzz zzz],'k','linewidth',1) ;
    end
end
for J = 2 : length(azmMer)-1
    if rem(azmMer(J),pi/2)~=0
        plot3(xxxMer*azmMer(J)/pi/2,yyyMer, ...
            zzz*ones(size(yyyMer)), ...
            ':','color',[64 64 64]/255,'linewidth',1) ;
    else
        plot3(xxxMer*azmMer(J)/pi/2,yyyMer, ...
            zzz*ones(size(yyyMer)),'k','linewidth',1) ;
    end
end
plot3([-.5 lngPar(1)/2],[yyyPar(1) yyyPar(1)], ...
    [zzz zzz],'k','linewidth',1) ;
plot3([-.5 lngPar(end)/2], ...
    [yyyPar(end) yyyPar(end)],[zzz zzz],'k','linewidth',1) ;
plot3(xxxMer*azmMer(1)/pi/2,yyyMer, ...
    zzz*ones(size(yyyMer)),'k','linewidth',1) ;
plot3(xxxMer*azmMer(end)/pi/2,yyyMer, ...
    zzz*ones(size(yyyMer)),'k','linewidth',1) ;

% Azimuth and elevation labels
azmLab = (-180:90:180)'*pi/180 ;
xxxLab = interp1(elvMap,parMap(:,1),-pi/2) ;
for J = 1 : length(azmLab)
    text(xxxLab*azmLab(J)/pi/2,-.28,num2str(azmLab(J)*180/pi), ...
        'horizontalalignment','center','FontSize',16) ;
end
text(0,-.33,'Azimuth [deg]', ...
    'horizontalalignment','center','FontSize',16) ;
elvLab = (-90:30:90)'*pi/180 ;
yyyLab = interp1(elvMap,parMap(:,2),elvLab) ;
xxxLab = -.51*ones(size(yyyLab)) ;
for J = 1 : length(elvLab)
    text(xxxLab(J),yyyLab(J),num2str(elvLab(J)*180/pi), ...
        'horizontalalignment','right','FontSize',16) ;
end
text(-.61,0,'Elevation [deg]','Rotation',90, ...
    'horizontalalignment','center','FontSize',16) ;

% Set the axes etc
xlim([-.51 .51]), axis equal off
set(gcf,'color',[1 1 1])


