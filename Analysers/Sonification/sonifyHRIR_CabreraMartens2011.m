function OUT = sonifyHRIR_CabreraMartens2011(HRIR,fs)
% A function to sonify head-related impulse responses, as described by
% Densil Cabrera and William L. Martens,
% "Sonifying head-related transfer functions",
% in Y. Suzuki, D. Brungart, Y. Iwaya, K. Iida, D.Cabrera and H. Kato (eds),
% "Principles and Applications of Spatial Hearing",
% New Jersey, World Scientific, 2011.
%
% Code by Densil Cabrera
% version 1.1 (9 October 2013)

% INPUT ARGUMENTS
%
% HRIR is a head-related impulse response. It can either be input as a
%   matrix (2-channels of audio usually, although 1 is ok too) or it can be
%   a structure. If it is a structure it must contain the following leaves:
%   HRIR.audio - the audio data (2-channels usually, although 1 is ok);
%   HRIR.fs - the audio sampling rate.
%
% fs is the audio sampling rate (only for when HRIR is input as a matrix)
%
% OUTPUT STRUCTURE
% OUT.audio is the sonification as an audio wave (2 channels)
% OUT.fs is the audio sampling rate of the sonification

% Interpret input
if isstruct(HRIR)
    HRIR = choose_from_higher_dimensions(HRIR,2,1); 
    OUT = HRIR;
    data = squeeze(HRIR.audio(:,:,1)); % discard 3rd dimension if present
    fs = HRIR.fs;
else
    data = squeeze(HRIR(:,:,1)); % discard 3rd dimension if present
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(data) && ~isempty(fs)
    [len, chans] = size(data);

    % if 1 channel, then duplicate it
    if chans == 1
        data = repmat(data,[1 2]);
    end

    % limit the length of the impulse response to 5 ms
    if len>ceil(fs*0.005)

        % find channel with peak value, and its index
        maxchan = (max(abs(data(:,1))) < max(abs(data(:,2))))+1;
        maxindex = find(abs(data(:,maxchan)) == max(abs(data(:,maxchan))),1);

        % determine the start time (1 ms pre peak, if available)
        startindex = maxindex - ceil(fs * 0.001);
        if startindex < 1, startindex = 1; end

        % crop data
        data = data(startindex:startindex+floor(fs*0.005)-1,:);
        len = length(data);
    end

    % sonification parameters
    hisslevel = -40; % relative level of hiss in dB
    e1	= 1;         % exponent for untransposed magnitude spectrum
    e2	= 3;         % exponent for transposed magnitude spectrum
    e3	= 1.5;       % exponent for envelope
    s1	= 10;        % time-stretch for carrier
    s2	= 1000;      % time-stretch for envelope

    % data length and amplitude
    data	= data(1:2*floor(end/2),1:2); % even length, 2 channels
    outlength	= length(data) * s2;      % length of output wave
    rms	= mean(data.^2) .^ 0.5;  % root-mean-square of each HRIR channel

    % generate hiss (steady state rendering of HRIR)
    hiss	= steadystate(data, outlength, e1, rms) .* 10 .^ (hisslevel/20);

    % generate carrier
    carrier	= steadystate(resample(data,s1,1), outlength,e2,rms);

    % generate envelope (normalized)
    envelope	= abs(hilbert(resample(data,s2,1))) .^ e3;
    envelope	= envelope ./ repmat(max(abs(envelope)), length(envelope), 1);

    % combine elements to make the sonification
    y	= (hiss+carrier.*envelope) ./ (1+10^(hisslevel/20));
    
    if isstruct(HRIR)
        OUT.audio = y;
        OUT.fs = fs;
        OUT.funcallback.name = '';
        OUT.funcallback.inarg = {};
    else
        OUT = y;
    end
else
    OUT = [];
end


% Plot sonification
figure('name', 'HRIR visualisation')

ILDthreshold = -20;

subplot(3,1,1)
L = 10*log10(abs(hilbert(data)).^2);
t = 1000*((1:len)'-1)/fs;
plot(repmat(t,[1,2]),L-max(max(L)))
ylim([-30 0])
xlabel('Time (ms)')
ylabel('Magnitude (dB)')
title('HRIR envelope and ILD')
grid on
hold on
ILDthresholdindex = (L(:,1)-max(max(L))>ILDthreshold) |(L(:,2)-max(max(L))>ILDthreshold);
ILD = L(:,1)-L(:,2);
ILD(~ILDthresholdindex) = NaN;
plot(t,ILD-15,'Color', [0.8 0.4 0.4])
plot([t(1), t(end)*0.9], [-15, -15], '--', 'Color', [0.6 0.4 0.4])
text(t(end)*0.91,-15, 'ILD 0 dB', 'Color', [0.6 0.4 0.4])
hold off

subplot(3,1,2)
% minimum length of the fft (zero-pads the wave if the wave is shorter)
% 10 ms gives 100 Hz resolution
if len < ceil(fs *0.01), len = ceil(fs *0.01); end
% the wave will be zero-padded in the fft by its length - this allows us to
% do linear (rather than circular) cross correlation
len = len * 2;

% find the spectrum (of both channels)
% len determines the length of the fft in samples
spectrum = fft(data, len);

% spectrum in dB below the Nyquist frequency
Levelspectrum = 10*log10(abs(spectrum(1:round(len/2),:)).^2);
frequencies = fs .* ((1:round(len/2))-1) ./ len;

% plot the magnitude spectrum
plot(frequencies, Levelspectrum - max(max(Levelspectrum)))
ylim([-30 0])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('HRTF spectrum and ILD')
grid on
hold on
ILDthresholdindex =  ...
    (Levelspectrum(:,1) - max(max(Levelspectrum))>ILDthreshold) ...
    | (Levelspectrum(:,2) - max(max(Levelspectrum))>ILDthreshold);
ILD = Levelspectrum(:,1)-Levelspectrum(:,2);
ILD(~ILDthresholdindex) = NaN;
plot(frequencies, ILD-15,'Color', [0.8 0.4 0.4])
plot([frequencies(1), frequencies(end)*0.9], [-15, -15], '--', 'Color', [0.6 0.4 0.4])
text(frequencies(end)*0.91,-15, 'ILD 0 dB', 'Color', [0.6 0.4 0.4])
hold off

subplot(3,1,3)


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
plot((timerange-round(len/2))./fs*1000, IACF(timerange),'r')
hold on
% plot the IACC datapoint
scatter(ITD,IACC,'ko')
ytextpos = IACC;
if ytextpos > 0.8, ytextpos = 0.8; end
text(ITD+0.05,ytextpos,['tau = ', num2str(round(ITD*1000)/1000), ' ms'])
xlabel('Lag (ms)')
ylabel('Coefficient')
title('Interaural cross correlation function')
ylim([-1 1])
grid on
hold off




function y	= steadystate(HRIR, outlength, exponent, rms)
% function to generate steady-state sound using a magnitude spectrum
% and a spectrum exponent, at a specified rms amplitude.

% zero-padded fast Fourier transform
spectrum	= fft(HRIR, outlength);

% apply exponent to magnitude spectrum
magnitude	= abs(spectrum).^ exponent;

% aquire phase data - we will preserve interaural phase difference
phase	= angle(spectrum);

% random phase for noise generation
randphase	= exp(1i*2*pi.*rand(outlength /2-1,1));
noise	= [0; randphase; 0; flipud(conj(randphase))];

% generate signal (note that real() is used in case of rounding errors)
y	= real(ifft(magnitude.*exp(1i*phase).*[noise,noise]));

% normalise to rms
y	= repmat(rms, outlength,1) .* y ./ repmat(mean(y.^2).^0.5, outlength,1);
