function [directCorr] = ...
     Mic2Hoa2SpkDirectiveBeamsCorrelation(mic2HoaCfg,Dec,frq,varargin)
 
% Default option values
plotOption  = true ;
dimMode     = '3d' ;
freqScale  = 'log' ;
angularRes  = 2 ;
[souSphCoord(:,1),souSphCoord(:,2)] = RegularPolyhedron(642) ;
angularResOr2dModePassed = false ;

% Number of optional arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% New option values are checked and assigned
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        case 'plotOption'
            if islogical(varargin{I+1})
                plotOption = varargin{I+1} ;
            else
                error('plotOption must be a logical value') ;
            end
        case 'freqScale'
            switch varargin{I+1}
                case {'linear','log'}
                    freqScale = varargin{I+1} ;
                otherwise
                    error(['Frequency scale available options ' ...
                        'are ''linear'' and ''log''']) ;
            end
        case 'dimMode'
            if strcmp(varargin{I+1},'2d') || strcmp(varargin{I+1},'3d')
                dimMode = varargin{I+1} ;
                if strcmp(dimMode,'2d')
                    angularResOr2dModePassed = true ;
                end
            else
                error('Unknown plotting mode') ;
            end
        case 'angularRes'
            if isreal(varargin{I+1}) && isscalar(varargin{I+1})
                angularRes = varargin{I+1} ;
                angularResOr2dModePassed = true ;
            else
                error('Angular resolution must be a real value') ;
            end
        case 'souSphCoord'
            if isreal(varargin{I+1}) && (size(varargin{I+1},2)==2)
                souSphCoord = varargin{I+1} ;
            else
                error(['Source spherical coordinates must be passed ' ...
                    'in a Nx2 array of real values']) ;
            end
        otherwise
            error('Unknown option') ;
    end
end

% Source positions in the '2d' and/or specified angular res. case
if angularResOr2dModePassed
    if strcmp(dimMode,'2d')
        souSphCoord = (-180:angularRes:180).'*pi/180 ;
        souSphCoord = [souSphCoord zeros(size(souSphCoord))] ;
    else
        souAzm = (-180:angularRes:180)'*pi/180 ;
        souElv = (-090:angularRes:090)'*pi/180 ;
        [souAzm,souElv] = meshgrid(souAzm,souElv) ;
        souSphCoord = [ souAzm(:) souElv(:) ] ;
    end
end

% Directive beams
[achBeams,refBeams] = Mic2Hoa2SpkDirectiveBeams(mic2HoaCfg,Dec,frq, ...
    'souSphCoord',souSphCoord,'plotOption',false) ;

% Spatial correlation
directCorr = zeros(length(frq),1) ;
for I = 1 : length(frq)
    correlation = zeros(size(Dec,1),1) ;
    for J = 1 : size(Dec,1)
        correlation(J) = sum(abs(achBeams(:,J,I)).*abs(refBeams(:,J))) ...
            / (norm(abs(achBeams(:,J,I)))*norm(abs(refBeams(:,J)))) ;
    end
    directCorr(I) = mean(correlation) ;
end


    

% Plot of the directivity correlation
if plotOption
    figure ;
    plot([frq(1) frq(end)],[100 100],'lineStyle',':','color',[.8 .8 .8]);
    hold on
    plotFreq = frq + 1*(frq==0) ;
    plot(plotFreq,100*directCorr) ;
    set(gca,'Xscale',freqScale,'layer','top') ;
    set(gca,'YTick',0:10:100,'YTickLabel', ...
        {'0','10','20','30','40','50','60','70','80','90','100'}) ;
    axis([max(plotFreq(1),20) plotFreq(end) 0 110]) ;
    xlabel('Frequency [Hz]') ;
    ylabel('Directivity correlation [%]') ;
end