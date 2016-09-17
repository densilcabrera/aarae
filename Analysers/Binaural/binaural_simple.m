function out = binaural_simple(in, fs, f_low, f_high, t_start, t_end)
% This function calculates parameters from binaural impulse response
% measurements.
%
% You can use it to analyse the direct sound of binaural impulse
% responses if you want to know the ITD and ILD over particular frequency
% ranges. These should indicate what direction the sound comes from, in
% terms of lateralization.
%
% When you analyse longer portions of binaural impulse responses (e.g., 50
% ms or 80 ms), the interaural cross correlation function is likely to be
% of interest (including the IACC value, which is used to indicate the
% 'apparent source width' or spaciousness of the sound). A value close to 1
% indicates that the sound mainly comes from 1 direction, whereas a value
% close to 0 indicates that the sound is spread out in space.
%
% You can filter the binaural impulse response by specifying values of
% f_low and f_high (these must be between 0 Hz and the Nyquist frequency).
% If you omit these input arguments, then no filtering is done. When they
% are included, filtering is performed in the frequency domain by zeroing
% the spectrum components outside of the specified range. Filtering affects
% most of the output arguments (IACC, ITD, IACF, ILD, Level) and the IACF plot.
% It does not affect the plots of the wave, envelope and spectrum nor the
% Levelspectrum and frequencies output arguments.
%
% This function does not aim to replicate standard IACC methods, but
% instead is intended as an exploratory tool. It is not suitable for
% analysing binaural impulse responses filtered into multiple bands.
%
% Code by Densil Cabrera
% version 1.0 (17 October 2013)


% INPUT ARGUMENTS
% wave is two channels of audio (a binaural impulse response). This
% function will crash if you input 1 channel!
% f_low is the lower cut-off frequency of a bandpass filter in Hz.
% f_high is the upper cut-off frequency of a bandpass filter in Hz.
% f_high must be greater than f_low, and you should consider the frequency
% resolution of the fft in choosing these values.
% fs is the audio sampling rate in Hz.
%
% OUTPUT ARGUMENTS
% IACC is the interaural cross correlation coefficient of the filtered
% spectrum.
% ITD is the interaural time difference of the filtered spectrum in ms.
% IACF is the interaural cross correlation function (sampling rate is fs).
% ILD is the interaural level difference of the filtered spectrum, in dB.
% Level is the level in dB of each channel (after filtering). Note that
% this is an energy level, rather than a power level - i.e. it will
% increase with the duration of the input wave.
% Levelspectrum is the magnitude spectrum of each channel expressed in dB.
% frequencies is a vector with the frequencies pertaining to Levelspectrum

if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1); 
    wave = in.audio;
    fs = in.fs;
else
    wave = audio;
    % default values of the bandpass filter are 0 Hz and Nyquist frequency
    % (i.e., no filtering)
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

S = size(wave); % size of the audio
ndim = length(S); % number of dimensions
switch ndim
    case 1
        len = S(1); % number of samples in audio
        chans = 1; % number of channels
        %bands = 1; % number of bands
    case 2
        len = S(1); % number of samples in audio
        chans = S(2); % number of channels
        %bands = 1; % number of bands
    case 3
        len = S(1); % number of samples in audio
        chans = S(2); % number of channels
        %bands = S(3); % number of bands
        wave = mean(wave,3); % mixdown bands
end

% only allow analysis of 2 channels
if chans == 1
    wave = [wave, wave];
    chans = 2;
elseif chans > 2
    wave = wave(:,1:2);
    chans = 2;
end

if nargin < 6, t_end = num2str(floor(1000*(len-1)/fs)); end
if nargin < 5, t_start = 0; end
if nargin < 4, f_high = fs/2; end
if nargin < 3, f_low = 0; end

%if isstruct(in)
%dialog box for settings
if nargin < 2
    prompt = {'Upper cutoff frequency (Hz)', ...
        'Lower cutoff frequency (Hz)', ...
        'Start time (ms)', ...
        'End time (ms)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {num2str(fs/2), '0', '0', num2str(floor(1000*(len-1)/fs))};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        f_high = str2num(answer{1,1});
        f_low = str2num(answer{2,1});
        t_start = str2num(answer{3,1});
        t_end = str2num(answer{4,1});
    end
end
    
if ~isempty(wave) && ~isempty(fs) && ~isempty(f_low) && ~isempty(f_high) && ~isempty(t_start) && ~isempty(t_end)
    % check crop points
    startsample = round(t_start * fs/1000) +1;
    if startsample < 1, startsample = 1; end
    endsample = round(t_end * fs/1000) +1;
    if endsample > len, endsample = len; end
    if endsample <= startsample
        startsample = 1;
        endsample = len;
    end
    % crop the wave
    wave = wave(startsample:endsample,:);
    len = length(wave);

    % plot of the input data and its envelope
    figure('Name', 'Binaural Input data')
    subplot(2,1,1)
    plot(((1:len)-1)./fs *1000, wave)
    xlabel('Time (ms)')
    ylabel('Amplitude')
    title('Wave')
    subplot(2,1,2)
    plot(((1:len)-1)./fs *1000, abs(hilbert(wave)))
    xlabel('Time (ms)')
    ylabel('Amplitude')
    title('Envelope')

    % minimum length of the fft (zero-pads the wave if the wave is shorter)
    % 10 ms gives 100 Hz resolution
    if len < ceil(fs *0.01), len = ceil(fs *0.01); end
    % the wave will be zero-padded in the fft by its length - this allows us to
    % do linear (rather than circular) cross correlation
    len = len * 2;

    % find the spectrum (of both channels)
    % len determines the length of the fft in samples
    spectrum = fft(wave, len);

    % spectrum in dB below the Nyquist frequency
    Levelspectrum = 10*log10(abs(spectrum(1:round(len/2),:)).^2);
    frequencies = fs .* ((1:round(len/2))-1) ./ len;

    % plot the magnitude spectrum
    figure('Name', 'Binaural Magnitude Spectrum')
    plot(frequencies, Levelspectrum)
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude spectrum')

    % The bandpass filter is implemented in the frequency domain by zero-ing
    % components outside of the desired range
    % lowpass filter: this zeros the spectrum from the cutoff frequency to its
    % image above the Nyquist frequency.
    cutoff_high = ceil(f_high * len /fs + 1);
    spectrum(cutoff_high:(len - cutoff_high +2),:)=0;

    % highpass filter: this zeros the spectrum from the cutoff frequency down
    % to 0 Hz, and does the same for the image frequencies above the Nyquist
    % frequency.
    cutoff_low = floor(f_low * len /fs + 1);
    spectrum(1:cutoff_low,:)=0;
    spectrum(end-cutoff_low+2:end,:)=0;

    % IACF is the interaural cross correlation function.
    % We do cross correlation by multiplication of spectra (i.e. the frequency
    % domain equivalent to time-domain cross correlation). Complex conjugate is
    % used for one of the spectra (which is the equivalent of time-reversal).
    IACF = ifft(spectrum(:,1) .* conj(spectrum(:,2)));

    % ACF is the autocorrelation function (of each channel).
    % This is used for normalising the IACF.
    ACF = ifft(spectrum .* conj(spectrum));
    % The peak of the ACF is at zero time lag (by definition).
    % We will use this for normalization.
    ACF_peak = ACF(1,:);

    % Normalise the IACF so that it ranges between -1 and 1.
    % This is done by dividing by the square root of the product of the two
    % autocorrelation function peaks.
    IACF = ifftshift(IACF / (ACF_peak(1)*ACF_peak(2))^0.5);

    % We are only interested in lag times of +/- 1 ms because this is the
    % natural range of ITD values (rounded up).
    timerange = round(len/2)+round(-fs/1000):round(len/2)+round(fs/1000);

    % IACC is the peak value of the IACF within the +/- 1 ms range. Usually
    % the absolute value of the IACF is taken before finding the maximum, but I
    % think the merits of that are debatable.
    IACC = max(IACF(timerange));

    % The interaural time difference (ITD) is the time at which the peak of the
    % IACF occurs. This is likely to be useful for direct sound, but may be
    % useless if you include room reflections and reverberation in the input
    % wave.
    ITD = find(IACF(timerange)==IACC); % find the peak index
    ITD = (ITD-1)./fs*1000-1; % convert to milliseconds

    % Sound level of each ear in dB
    Level = 10*log10(sum(abs(spectrum).^2));

    % Interaural level difference
    ILD = Level(1)-Level(2);

    % plot the IACF
    figure('Name', 'Binaural Interaural Cross Correlation')
    plot((timerange-round(len/2))./fs*1000, IACF(timerange),'r')
    hold on
    % plot the IACC datapoint
    scatter(ITD,IACC,'o')
    xlabel('Lag (ms)')
    ylabel('Coefficient')
    title('Interaural cross correlation function')
    ylim([-1 1])
    grid on
    hold off
    doresultleaf(IACF(timerange),'Coefficient',{'Lag'},...
                 'Lag',      (timerange-round(len/2))./fs*1000, 'ms',           true,...
                 'Function', {'Unique'},                        'categorical',  [],...
                 'name','IACF');

    % Group delay
    GD = diff(unwrap(angle(((spectrum(1:round(len/2),1)).*conj(spectrum(1:round(len/2),2))) ...
        ./(conj(spectrum(1:round(len/2),1)).*spectrum(1:round(len/2),1)))))...
        .*len/(fs*2*pi);
    GD = GD.*1000; % change to milliseconds
    % Plot group delay
    figure('Name', 'Binaural Group Delay')
    semilogx(frequencies(3:end), GD(2:end),'r')
    xlabel('Frequency (Hz)')
    xlim([frequencies(2) frequencies(end)])
    ylabel('Group Delay (ms)')
    grid on
    %median(GD)
    if isstruct(in)
        doresultleaf(GD(2:end),'Group Delay (ms)',{'Frequency'},...
                     'Frequency', frequencies(3:end), 'Hz',           true,...
                     'Function',  {'Unique'},         'categorical',  [],...
                     'name','Group_delay');
    end
    % organise output structure
    out.IACC = IACC;
    out.ITD = ITD;
    out.IACF = IACF;
    out.ILD = ILD;
    out.Level = Level;
    out.Levelspectrum = Levelspectrum;
    out.frequencies = frequencies;
    out.funcallback.name = 'binaural_simple.m';
    out.funcallback.inarg = {fs,f_low,f_high,t_start,t_end};
end