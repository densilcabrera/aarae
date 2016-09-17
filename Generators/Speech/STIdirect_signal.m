function OUT = STIdirect_signal(duration, fs, octavebandlevel, silenceduration, STImethod)
% Generates a signal for STI direct measurement (direct method) based on
% IEC 60268-16 (2011). One of the following four approaches can be taken:
% STImethod = 0:
%   generate seven STIPA-like signals that populate the entire 98-value MTF
%   matrix. This is one of the faster methods.
% STImethod = 1:
%   generate seven pairs of synchronously double-modulated signals - like
%   STIPA, but with the same modulation frequency pairs in each octave
%   band. This takes the same amount of time as the first method, but is
%   probably easier to understand and provides a direct experience of what
%   is meant by modulation transfer function.
% STImethod = 2:
%   generate fourteen moduated signals that populate the entire 98-value
%   MTF matrix one modulation frequency at a time. This is the slower
%   method, but it is easier to understand, and probably closer to the
%   intention of IEC 60268-16. It provides a direct experience of what is
%   meant by modulation transfer function. This method is faster than the
%   full direct method (see below) because all seven octave bands are
%   measured simultaneously.
% STImethod = 3:
%   generate 98 modulated signals one octave band frequency and modulation
%   frequency at a time. This is impractically slow for typical measurement
%   situations (a single measurement could even take 45 minutes), but can
%   be used as a reference direct method to compare with the other more
%   efficient methods. Note that instead of using half-octave bands of
%   noise, this method uses full octave bands of noise, which could be a
%   cause of deviations from the other methods. If you wish to generate
%   this signal, first make sure that you have sufficient time (and RAM) to
%   work with it.
%
% The duration of each signal should normally be 20 s or longer (longer
% durations normally yield more accurate results), and the default duration
% is 30 s. Because the direct method uses random noise, the noise
% introduces minor statistical errors (deviations from ideal MTF values) in
% the analysis - but these errors become smaller as the duration of the
% signal is increased.
%
% The seven or 14 (or 98) signals are separated by a gap of 1 s by default,
% but a longer gap may be helpful in reverberant conditions.
%
% Equalization (of the signal for the source) and calibration (of source
% and receiver) are required to use this properly (in the same way as for
% STIPA measurements). Simple equalization can be done by correctly setting
% the octave band levels (in this functions input dialog box) for the
% source being used.
%
% Analysis of the generated and re-recorded signals should be done using
% AARAE's STI_direct analyser (in the Speech folder). This analyser
% automatically detects which type of signal was generated (from the
% 'properties' data) and analyses it accordingly. Note that you may use the
% original signal as a reference signal as part of the analysis, and this
% should reduce (or ideally remove) the random influence of the noise in
% the analysis (although this process does not strictly follow IEC
% 60268-16).
%
% Code by Densil Cabrera
% version 1.01 (30 March 2015)

if nargin < 5
    
    param = inputdlg({'STIPA-like dual modulation method [0], Synchronous dual modulation method [1], Single modulation method (slower) [2], or Single modulation one octave band at a time (very slow reference method) [3]';...
        'Duration of each of signal in the sequence [s]';...
        'Duration of gap between each signal [s]';...
        '125 Hz octave band level (dB)';...
        '250 Hz octave band level (dB)';...
        '500 Hz octave band level (dB)';...
        '1000 Hz octave band level (dB)';...
        '2000 Hz octave band level (dB)';...
        '4000 Hz octave band level (dB)';...
        '8000 Hz octave band level (dB)';...
        'Broadband gain (dB)'; ...
        'Sampling frequency [samples/s]'}, ...
        'Noise input parameters',1,{'0';'30'; '1'; ...
        '2.9';'2.9';'-0.8';'-6.8';'-12.8';'-18.8';'-24.8'; ...
        '0';'48000'});
    param = str2num(char(param));
    if length(param) < 12, param = []; end
    if ~isempty(param)
        STImethod = param(1);
        duration = param(2);
        silenceduration = param(3);
        octavebandlevel(1) = param(4);
        octavebandlevel(2) = param(5);
        octavebandlevel(3) = param(6);
        octavebandlevel(4) = param(7);
        octavebandlevel(5) = param(8);
        octavebandlevel(6) = param(9);
        octavebandlevel(7) = param(10);
        Lgain = param(11);
        fs = param(12);
    end
else
    param = [];
    Lgain = 0;
end
if exist('fs','var') && fs<44100, fs = 44100; end

if ~isempty(param) || nargin ~= 0
    gaplen = round(silenceduration * fs);
    if duration < 2
        duration = 2;
        warndlg('STIdirect_signal minimum duration is 2 s, but 20 s or longer is recommended!')
    end
    
    % Generate pink noise
    % number of samples of the output wave
    nsamples = round(duration * fs);
    % even number of samples
    if rem(nsamples,2) == 1, nsamples = nsamples + 1; end
    
    fexponent = -1; % spectral slope exponent
    
    
    
    bandnumbers = 21:3:39;
    fc = 10.^(bandnumbers/10);
    if STImethod == 3
        bandwidth = 1; % full-octave bandwidths
    else
        bandwidth = 0.5; % half-octave bandwidths
    end
    
    noisebands = zeros(nsamples,7); % preallocate
    
    for k = 1:7
        
        % magnitude slope function (for half spectrum, not including DC and
        % Nyquist)
        magslope = ((1:nsamples/2-1)./(nsamples/4)).^(fexponent*0.5)';
        
        flowpass=fc(k)/10^(0.15*bandwidth); % high cut-off frequency in Hz
        fhighpass=fc(k)*10^(0.15*bandwidth); % low cut-off frequency in Hz
        
        lowcut = floor(flowpass * nsamples / fs);
        magslope(1:lowcut) = 0;
        
        highcut = ceil(fhighpass * nsamples / fs);
        magslope(highcut:end) = 0;
        
        % generate noise in the frequency domain, by random phase
        noisyslope = magslope .* exp(1i*2*pi.*rand(nsamples/2-1,1));
        clear magslope;
        
        % transform to time domain
        noisebands(:,k) = ifft([0;noisyslope;0;flipud(conj(noisyslope))]);
        clear noisyslope;
    end
    
    switch STImethod
        
        
        case 3
            Fm = [0.63,0.8,1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5];
            %m = 1; % modulation depth for simple STI
            t = ((1:nsamples)'-1)./fs;
            y = zeros(ceil(98*nsamples + 98*gaplen),1);
            startindex = zeros(98,1);
            
            sigcount = 1;
            for section = 1:14
                envelope = (1 + cos(2 * pi * Fm(section) .* t)).^0.5;
                noisebands1 = noisebands .* repmat(envelope,[1,7]);
                % adjust octave band level
                gainchange = 10.^(octavebandlevel./20) ./ mean(noisebands1.^2).^0.5;
                noisebands1 = noisebands1 .* repmat(gainchange,[nsamples,1]);
                
                for oct = 1:7
                    y1 = noisebands1(:,oct) * 10^(Lgain/20) / 20;
                    startindex(sigcount) = 1+nsamples*(sigcount-1)+gaplen*(sigcount-1); 
                    y(startindex(sigcount):startindex(sigcount)+nsamples-1) = y1;
                    sigcount = sigcount+1;
                end
            end
            Fm = repmat(Fm,[1,7]);
            
        case 2
            Fm = [0.63,0.8,1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5];
            %m = 1; % modulation depth for simple STI
            t = ((1:nsamples)'-1)./fs;
            y = zeros(ceil(14*nsamples + 13*gaplen),1);
            startindex = zeros(14,1);
            
            for section = 1:14
                envelope = (1 + cos(2 * pi * Fm(section) .* t)).^0.5;
                noisebands1 = noisebands .* repmat(envelope,[1,7]);
                % adjust octave band level
                gainchange = 10.^(octavebandlevel./20) ./ mean(noisebands1.^2).^0.5;
                noisebands1 = noisebands1 .* repmat(gainchange,[nsamples,1]);
                y1 = sum(noisebands1,2) * 10^(Lgain/20) / 20;
                startindex(section) = 1+nsamples*(section-1)+gaplen*(section-1);
                y(startindex(section):startindex(section)+nsamples-1) = y1;
            end
            
            
        case 1
            Fm = cat(3,[0.63, 3.15],[0.8, 4],[1, 5],...
                [1.25, 6.25],[1.6, 8],[2, 10],[2.5, 12.5]);
            Fm = repmat(Fm,[7,1,1]);
            m = 0.55; % Modulation depth for STIPA signal
            
            noisebands1 = zeros(size(noisebands));
            t = ((1:nsamples)'-1)./fs;
            y = zeros(ceil(7*nsamples+6*gaplen),1);
            startindex = zeros(7,1);
            for section = 1:7
                for k = 1:7
                    envelope = (0.5*(1+m*cos(2*pi*Fm(k,1,section)*t)-m*cos(2*pi*Fm(k,2,section)*t))).^0.5;
                    noisebands1(:,k) = noisebands(:,k) .* envelope;
                    
                    % adjust octave band level
                    gainchange = 10.^(octavebandlevel(k)/20) / mean(noisebands1(:,k).^2).^0.5;
                    noisebands1(:,k) = noisebands1(:,k) * gainchange;
                end
                % dividing by 20 provides a reasonable level signal (prior to Lgain
                % adjustment)
                y1 = sum(noisebands1,2) * 10^(Lgain/20) / 20;
                startindex(section) = 1+nsamples*(section-1)+gaplen*(section-1);
                y(startindex(section):startindex(section)+nsamples-1) = y1;
            end
            
            
        case 0
            % PART 1 is the same as a STIPA signal
            Fm = [1.6, 8;...
                1, 5;...
                0.63, 3.15;...
                2, 10;...
                1.25, 6.25;...
                0.8, 4;...
                2.5, 12.5];
            Fm = repmat(Fm,[1,1,7]);
            for section = 2:7
                Fm(:,:,section) = circshift(Fm(:,:,section-1),1);
            end
            m = 0.55; % Modulation depth for STIPA signal
            
            noisebands1 = zeros(size(noisebands));
            t = ((1:nsamples)'-1)./fs;
            y = zeros(ceil(7*nsamples+6*gaplen),1);
            startindex = zeros(7,1);
            for section = 1:7
                for k = 1:7
                    envelope = (0.5*(1+m*cos(2*pi*Fm(k,1,section)*t)-m*cos(2*pi*Fm(k,2,section)*t))).^0.5;
                    noisebands1(:,k) = noisebands(:,k) .* envelope;
                    
                    % adjust octave band level
                    gainchange = 10.^(octavebandlevel(k)/20) / mean(noisebands1(:,k).^2).^0.5;
                    noisebands1(:,k) = noisebands1(:,k) * gainchange;
                end
                % dividing by 20 provides a reasonable level signal (prior to Lgain
                % adjustment)
                y1 = sum(noisebands1,2) * 10^(Lgain/20) / 20;
                startindex(section) = 1+nsamples*(section-1)+gaplen*(section-1);
                y(startindex(section):startindex(section)+nsamples-1) = y1;
            end
    end
    
    
    
    
    tag = ['STIdirect' num2str(duration) 's' num2str(STImethod)];
    
    OUT.audio = y;
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 'STIdirect_signal.m';
    OUT.funcallback.inarg = {duration, fs, octavebandlevel, silenceduration, STImethod};
    OUT.properties.startindex=startindex;
    OUT.properties.Fm = Fm;
    OUT.properties.gap = gaplen;
else
    OUT = [];
end

end % End of function