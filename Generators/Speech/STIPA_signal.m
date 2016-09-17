function OUT = STIPA_signal(duration, fs, octavebandlevel)
% Generates a signal for STIPA measurement (direct method) according to IEC
% 60268-16 (2011)
%
% Equalization and calibration is required to use this properly.
%
% Code by Densil Cabrera
% version 1.01 (27 December 2013)

if nargin < 3
    
    param = inputdlg({'Duration of the wave [s]';...
        '125 Hz octave band level (dB)';...
        '250 Hz octave band level (dB)';...
        '500 Hz octave band level (dB)';...
        '1000 Hz octave band level (dB)';...
        '2000 Hz octave band level (dB)';...
        '4000 Hz octave band level (dB)';...
        '8000 Hz octave band level (dB)';...
        'Broadband gain (dB)'; ...
        'Sampling frequency [samples/s]'}, ...
        'Noise input parameters',1,{'30'; ...
        '2.9';'2.9';'-0.8';'-6.8';'-12.8';'-18.8';'-24.8'; ...
        '0';'48000'});
    param = str2num(char(param));
    if length(param) < 3, param = []; end
    if ~isempty(param)
        duration = param(1);
        octavebandlevel(1) = param(2);
        octavebandlevel(2) = param(3);
        octavebandlevel(3) = param(4);
        octavebandlevel(4) = param(5);
        octavebandlevel(5) = param(6);
        octavebandlevel(6) = param(7);
        octavebandlevel(7) = param(8);
        Lgain = param(9);
        fs = param(10);
    end
else
    param = [];
Lgain = 0;
end
if exist('fs','var') && fs<44100, fs = 44100; end

if ~isempty(param) || nargin ~= 0
    
    % Generate pink noise
    % number of samples of the output wave
    nsamples = duration * fs;
    % even number of samples
    if rem(nsamples,2) == 1, nsamples = nsamples + 1; end
    
    fexponent = -1; % spectral slope exponent
    

    
    bandnumbers = 21:3:39;
    fc = 10.^(bandnumbers/10);
    bandwidth = 0.5; % half-octave bandwidths
    
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
    
    % Modulation frequencies in Hz
%     Fm = 10.^([2, 9; ...
%         0, 7;...
%         -2, 5;...
%         3, 10;...
%         1, 8;...
%         -1, 6;...
%         4, 11]./10);

    Fm = [1.6, 8;...
        1, 5;...
        0.63, 3.15;...
        2, 10;...
        1.25, 6.25;...
        0.8, 4;...
        2.5, 12.5];
     
    m = 0.55; % Modulation depth for STIPA signal
    
    
    t = ((1:nsamples)'-1)./fs;
    for k = 1:7
        
envelope = (0.5*(1+m*cos(2*pi*Fm(k,1)*t)-m*cos(2*pi*Fm(k,2)*t))).^0.5;
noisebands(:,k) = noisebands(:,k) .* envelope;
        
        % adjust octave band level
        gainchange = 10.^(octavebandlevel(k)/20) / mean(noisebands(:,k).^2).^0.5;
        noisebands(:,k) = noisebands(:,k) * gainchange;
    end
    
    % dividing by 20 provides a reasonable level signal (prior to Lgain
    % adjustment)
    y = sum(noisebands,2) * 10^(Lgain/20) / 20;
    
    tag = ['STIPA' num2str(duration)];
    
    OUT.audio = y;
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 'STIPA_signal.m';
    OUT.funcallback.inarg = {duration, fs, octavebandlevel};
else
    OUT = [];
end

end % End of function