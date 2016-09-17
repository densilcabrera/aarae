function [snr,freq] = Mic2HoaSignalToNoiseRatio(mic2HoaCfg,varargin)

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

% Calculates the SNR (inverse of the encoding error)
[err,freq] = Mic2HoaEncodingError(mic2HoaCfg,'noiseLvl',noiseLvl, ...
    'mode',mode,'freqScale',freqScale,'colorMap',colorMap, ...
    'plotOption',false) ;
snr = 1./err ;

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
            snr = snr(:,prm) ;
            
            % "surf" or "imagesc" of the encoding error
            switch freqScale
                case 'log'
                    % Surfs the encoding error
                    snrTmp = 20*log10(snr) ;
                    snrTmp(end+1,end+1) = 0 ;
                    figure
                    hold on ;
                    surf([freq;freq(end)],0:hoaFmt.nbComp,snrTmp') ;
                    shading flat ;
                    set(gca,'Xscale','log')
                case 'linear'
                    % Uses 'imagesc' to plot the encoding error
                    figure
                    hold on ;
                    imagesc([0;freq(2:end)],0.5:hoaFmt.nbComp-.5, ...
                        20*log10(snr)') ;

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
                    max(max(snr))*ones(2,1),':','color',couleurDeg);
                YTick(I) = J + orders(I,2)/2 ;
                J = J + orders(I,2) ;
            end
            set(gca,'YTick',YTick,'YTicklabel',orders(:,1));
            
            % Other axis and tick settings
            axis([50 0.998*freq(end) 0 hoaFmt.nbComp]);
            grid off ;
            box on ;
            set(gca,'layer','top');
            xlabel('frequency [Hz]');
            ylabel('order');
            caxis([0 40]) ;
            colormap(couleurMap) ;
            colorbar('Ytick',[0 10 20 30 40], ...
                'Yticklabel',{'0','10','20','30','40'});
            XTick = get(gca,'XTick') ;
            for I = 1 : length(XTick)
                plot3(XTick(I)*ones(2,1),[0;hoaFmt.nbComp], ...
                    max(max(snr))*ones(2,1),':','color',couleurDeg) ;
            end
            hold off;

        case {'orders','global'}
            
            figure
            plot(freq,20*log10(snr))
            if strcmp(freqScale,'log')
                set(gca,'Xscale','log')
            end
            axis([ 20 max(freq) ...
                   20*ceil(max(max(log10(snr))))-80 ...
                   20*ceil(max(max(log10(snr)))) ]) ;
            xlabel('Frequency [Hz]') ;
            ylabel('Signal to noise ratio [dB]')
    end
end
