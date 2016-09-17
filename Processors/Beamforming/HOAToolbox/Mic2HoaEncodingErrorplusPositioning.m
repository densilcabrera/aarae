function [encodingError,freq] = Mic2HoaEncodingErrorplusPositioning(mic2HoaCfg,micFmtIdeal,varargin)

% Default option values
noiseLvl   = -40 ;
plotOption = true ;
mode       = 'harmonics' ;
freqScale  = 'log' ;
colorMap   = 'co' ;

% Number of optional arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% New option values are checked and assigned
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        case 'noiseLvl'
            if isscalar(varargin{I+1}) && isreal(varargin{I+1})
                noiseLvl = varargin{I+1} ;
            else
                error('noiseLvl must be a real value') ;
            end
        case 'plotOption'
            if islogical(varargin{I+1})
                plotOption = varargin{I+1} ;
            else
                error('plotOption must be a logical value') ;
            end
        case 'mode'
            switch varargin{I+1}
                case {'harmonics','orders','global'}
                    mode = varargin{I+1} ;
                otherwise
                    error(['Available modes ' ...
                        'are ''harmonics'', ''orders'' and ''global''']) ;
            end
        case 'freqScale'
            switch varargin{I+1}
                case {'linear','log'}
                    freqScale = varargin{I+1} ;
                otherwise
                    error(['Frequency scale available options ' ...
                        'are ''linear'' and ''log''']) ;
            end
        case 'colorMap'
            switch varargin{I+1}
                case {'co','bw'}
                    colorMap = varargin{I+1} ;
                otherwise
                    error(['Colormap available options are ' ...
                        '''co'' and ''bw''']) ;
            end
        otherwise
            error('Unknown option') ;
    end
end

% Shortcuts
hoaFmt     = mic2HoaCfg.hoaFmt ;
micFmt     = mic2HoaCfg.micFmt ;
mic2HoaOpt = mic2HoaCfg.mic2HoaOpt ;
filters    = mic2HoaCfg.filters ;

% Encoding error is only computed for the highest sampling frequency 
% value available
[fs,indFs] = max(mic2HoaOpt.sampFreq) ;
nbFft = mic2HoaOpt.filterLength(indFs) ;

% Frequency vector (and corresponding wave number)
freq = [1, fs/nbFft : fs/nbFft : fs/2]' ;

% HoaFmt required for an accurate simulation
ordPhys = Hoa2MicRequiredOrder(micFmtIdeal,max(freq)) ;
ordPhys = max(max(hoaFmt.index(:,1)),ordPhys) ;
hoaFmtPhys = GenerateHoaFmt('res2d',ordPhys,'res3d',ordPhys, ...
    'type',hoaFmt.type,'conv',hoaFmt.conv) ;

% hoaFmtPhys -> hoaFmt conversion matrix
CnvMat = Hoa2HoaConversionMatrix(hoaFmtPhys,hoaFmt) ;

% Encoding filter responses
switch mic2HoaOpt.filterType
    case 'postEqFir'
        % Equalization filter Fourier Transforms are computed
        postEqFir = fft(fftshift(filters(indFs).postEqFir,1),[],1) ;
        postEqFir = postEqFir(:,hoaFmt.index(:,1)+1) ;
        % Global encoding filter responses (including gain matrix)
        filterResponses = repmat(postEqFir,[1 1 micFmt.nbMic]) ...
            .* repmat(permute(filters(indFs).gainMatrix,[3 1 2]), ...
            [size(postEqFir,1) 1 1]) ; 
    case 'firMatrix'
        % Encoding filter Fourier Transforms are computed
        filterResponses = fft(fftshift(filters(indFs).firMatrix,1),[],1) ;
end

% HOA to mic transfer matrix
T = Hoa2MicTransferMatrix(hoaFmtPhys,micFmtIdeal,freq) ;

% Encoding error due to aliasing / amplification
encodingError = zeros(length(freq),hoaFmt.nbComp) ;
for I = 1 : length(freq)
    encodingError(I,:) = sum(abs( ...
        permute(filterResponses(I,:,:),[2 3 1]) ...
        * T(:,:,I) - CnvMat ).^2 ,2) ;
end

% Adding the error due to measurement noise
encodingError = encodingError ...
    + sum(abs(filterResponses(1:nbFft/2+1,:,:)).^2,3) ...
    * 10^(noiseLvl/10) ;

% Format the error depending on the "mode"
switch mode
    
    case 'harmonics'
        
        % One error value for each harmonic
        encodingError = sqrt(encodingError) ;
        
    case 'orders'
        
        % List of the HOA orders encoded
        hoaOrder = unique(hoaFmt.index(:,1)) ; 
        
        % Encoding error for each order
        err = zeros(length(freq),length(hoaOrder)) ;
        for I = 1 : length(hoaOrder)
            ordIdx = hoaFmt.index(:,1) == hoaOrder(I) ;
            err(:,I) = sqrt(mean(encodingError(:,ordIdx),2)) ;
        end
        encodingError = err ;
        
    case 'global'
        
        % Average value on every harmonic
        encodingError = sqrt(mean(encodingError,2)) ;
end

% Figure
if plotOption

    switch mode

        case 'harmonics'
            
            % colormap is chosen
            switch colorMap
                case 'co'
                    couleurMap = hot(64) ;
                    couleurDeg = [.7 .7 .7] ;
                case 'bw'
                    couleurMap = repmat((.2:.8/63:1)',1,3) ;
                    couleurDeg = [0 0 0] ;
            end
            
            % The harmonics are sorted by increasing order
            [hoaFmt.index,prm] = sortrows(hoaFmt.index,1) ;
            encodingError = encodingError(:,prm) ;
            
            % "surf" or "imagesc" of the encoding error
            switch freqScale
                case 'log'
                    % Surfs the encoding error
                    err = 20*log10(encodingError) ;
                    err(end+1,end+1) = 0 ;
                    figure
                    hold on ;
                    surf([freq;freq(end)],0:hoaFmt.nbComp,err') ;
                    shading flat ;
                    set(gca,'Xscale','log')
                case 'linear'
                    % Uses 'imagesc' to plot the encoding error
                    figure
                    hold on ;
                    imagesc([0;freq(2:end)],0.5:hoaFmt.nbComp-.5, ...
                        20*log10(encodingError)') ;

            end
            
            % list of the orders and the number of harmonics in each order
            orders = unique(hoaFmt.index(:,1)) ;
            nmbOrd = size(orders,1) ;
            for I = 1 : nmbOrd
                orders(I,2) = length(find(hoaFmt.index(:,1)==orders(I))) ;
            end
            
            % Compute the y-axis ticks and corresponding labels
            J = 0 ;
            YTick = zeros(nmbOrd,1) ;
            for I = 1 : nmbOrd
                plot3([freq(1);freq(end)],(J+orders(I,2))*ones(2,1), ...
                    max(max(encodingError))*ones(2,1),':', ...
                    'color',couleurDeg) ;
                YTick(I) = J + orders(I,2)/2 ;
                J = J + orders(I,2) ;
            end
            set(gca,'YTick',YTick,'YTicklabel',orders(:,1));
            
            % Other axis and tick settings
            axis([20 0.998*freq(end) 0 hoaFmt.nbComp]);
            grid off ;
            box on ;
            set(gca,'layer','top');
            xlabel('frequency [Hz]');
            ylabel('order');
            caxis([-40 0]) ;
            colormap(couleurMap) ;
            colorbar('Ytick',[-40 -30 -20 -10 0], ...
                'Yticklabel',{'-40','-30','-20','-10','0'});
            XTick = get(gca,'XTick') ;
            for I = 1 : length(XTick)
                plot3(XTick(I)*ones(2,1),[0;hoaFmt.nbComp], ...
                    max(max(encodingError))*ones(2,1),':', ...
                    'color',couleurDeg);
            end
            hold off;

        case {'orders','global'}
            
            figure
            plot(freq,20*log10(encodingError))
            if strcmp(freqScale,'log')
                set(gca,'Xscale','log')
            end
            axis([ 20 max(freq) ...
                   20*floor(min(min(log10(encodingError)))) ...
                   20*ceil(max(max(log10(encodingError)))) ]) ;
            xlabel('Frequency [Hz]') ;
            ylabel('Encoding error [dB]')
    end
end
