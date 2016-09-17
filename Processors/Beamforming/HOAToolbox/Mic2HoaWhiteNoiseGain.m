function [wng,frq] = Mic2HoaWhiteNoiseGain(mic2HoaCfg,varargin)

% Default option values
plotOption = true ;
freqScale  = 'log' ;

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
        otherwise
            error('Unknown option') ;
    end
end

% Shortcuts
hoaFmt     = mic2HoaCfg.hoaFmt ;
micFmt     = mic2HoaCfg.micFmt ;
mic2HoaOpt = mic2HoaCfg.mic2HoaOpt ;
filters    = mic2HoaCfg.filters ;

% The W.N.G is only computed for the highest sampling frequency 
% value available
[fs,indFs] = max(mic2HoaOpt.sampFreq) ;
nbFft = mic2HoaOpt.filterLength(indFs) ;

% Frequency vector (and corresponding wave number)
frq = [1, fs/nbFft : fs/nbFft : fs/2]' ;

% Hoa to mic transfer functions
HoaMic = Hoa2MicTransferMatrix(hoaFmt,micFmt,frq) ;

% Encoding filter responses
switch mic2HoaOpt.filterType
    case 'postEqFir'
        % Equalization filter Fourier Transforms are computed
        postEqFir = fft(fftshift(filters(indFs).postEqFir,1),[],1) ;
        postEqFir = postEqFir(:,hoaFmt.index(:,1)+1) ;
        % Global encoding filter responses (including gain matrix)
        FltRsp = repmat(postEqFir,[1 1 micFmt.nbMic]) ...
            .* repmat(permute(filters(indFs).gainMatrix,[3 1 2]), ...
            [size(postEqFir,1) 1 1]) ; 
    case 'firMatrix'
        % Encoding filter Fourier Transforms are computed
        FltRsp = fft(fftshift(filters(indFs).firMatrix,1),[],1) ;
end
FltRsp = permute(FltRsp,[2 3 1]) ;

% White noise gain
wng = zeros(length(frq),hoaFmt.nbComp) ;
for I = 1 : length(frq)
    for J = 1 : hoaFmt.nbComp
        wng(I,J) = abs(FltRsp(J,:,I)*HoaMic(:,J,I)).^2 / norm(FltRsp(J,:,I)).^2 ;
    end
end


% Plot the White noise gain
if plotOption
    figure
    plot(frq,10*log10(wng))
    maxLvl = round(max(max(log10(wng))))*10 ;
    axis([50,max(frq),maxLvl-70,maxLvl+10]) ;
    xlabel('Frequency [Hz]') ;
    ylabel('White Noise Gain [dB]')
    if strcmp(freqScale,'log')
        set(gca,'xscale','log') ;
    end
end
