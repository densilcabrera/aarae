function OUT = autoconvolve_spectexponent(in,fs,exponent,normalize,hiF,loF,audioplay)
% Performs n-th power autoconvolution by raising the signal's spectrum to
% an exponent. 
%
% Use an exponent of 2 for first order autoconvolution.
% Higher integers yield more extreme results (and longer duration output).
% -1 inverts the spectrum.
% 0 yields an impulse.
% 1 yields the original (notwithstanding rounding errors).
% Note that fractional exponents will yield complex-valued
% signals.
%
% The output is normalized by default (but this can be bypassed).
% An optional lowpass and/or highpass filter can be applied.
%
% Code by Densil Cabrera & Daniel Jimenez.
% version 2.0 (20 March 2014)
%
% INPUT AND OUTPUT ARGUMENTS
% The input argument is a structure with the following fields required:
% in.audio: contains the audio data (dim1 is time, dim2 is chans, dim3 is
%   bands in the case of pre-filtered signals).
% in.fs: is the audio sampling rate in Hz.
%
% The output waveform retains the same dim2 and dim3 configuration as the
% input, and has the same sampling rate, but will probably be a different 
% length.
if nargin < 7, audioplay = 0; end
if nargin < 6, loF = 0; end
if nargin < 5, hiF = 48000; end
if nargin < 4, normalize = 1; end
if nargin < 3
    exponent = 2;
    if isstruct(in)
        % fft must be long enough to avoid circular autoconvolution
        audio = in.audio;
        len = length(audio);
        fs = in.fs;
        Nyquist = fs/2;
    else
        audio = in;
        len = length(audio);
        if nargin < 2
            fs = inputdlg({'Sampling frequency [samples/s]'},...
                               'Fs',1,{'48000'});
            fs = str2num(fs);
            Nyquist = fs/2;
        else
            Nyquist = fs/2;
        end
    end
    % dialog box to get user settings
    prompt = {'Exponent','Normalize (0 | 1)', 'High cutoff frequency (Hz)', ...
        'Low cutoff frequency'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'2','1',num2str(Nyquist),'0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer)
        exponent = str2num(answer{1,1});
        normalize = str2num(answer{2,1});
        hiF = str2num(answer{3,1});
        loF = str2num(answer{4,1});
    else
        OUT = [];
        return
    end
end
if isstruct(in)
    % fft must be long enough to avoid circular autoconvolution
    audio = in.audio;
    len = length(audio);
    fs = in.fs;
    Nyquist = fs/2;
else
    audio = in;
    len = length(audio);
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(fs);
        Nyquist = fs/2;
    else
        Nyquist = fs/2;
    end
end
fftlen = len * abs(exponent);
if fftlen < len, fftlen = len; end
fftlen = 2*(ceil(fftlen/2)); % even length to simplify filtering

% derive spectrum and apply exponent
spectrum = fft(audio,fftlen).^exponent;

% lowpass filter
if (hiF < Nyquist) && (hiF > loF)  && (hiF > 0)
    hicomponent = ceil(hiF / (in.fs / fftlen)) + 1;
    spectrum(hicomponent:fftlen - hicomponent+2,:,:,:,:,:) = 0;
end

% hipass filter
if (loF > 0) && (loF < hiF) && (loF < Nyquist)
    locomponent = floor(loF / (fs / fftlen)) + 1;
    spectrum(1:locomponent,:,:,:,:,:) = 0;
    spectrum(fftlen-locomponent+2:fftlen,:,:,:,:,:) = 0;
end

% return to time domain
out = ifft(spectrum);

% normalization
if normalize
    out = out / max(max(max(max(max(max(abs(abs(out))))))));
end

if isstruct(in)
    OUT.audio = out;
    OUT.funcallback.name = 'autoconvolve_spectexponent.m';
    OUT.funcallback.inarg = {fs,exponent,normalize,hiF,loF,audioplay};
else
    OUT = out;
end

if audioplay
    if ~isreal(out)    
        disp('Autoconvolution output is complex.')
    else
        wavout = out;
    end

    sound(sum(sum(sum(sum(wavout,3),4),5),6)...
        ./max(max(abs(sum(sum(sum(sum(wavout,3),4),5),6)))), in.fs)

    % Loop for replaying, saving and finishing
    choice = 'x'; % create a string
    
    % loop until the user presses the 'Done' button
    while ~strcmp(choice,'Done')
        choice = questdlg('What next?', ...
            'Autoconvolution', ...
            'Play again', 'Save audio as .wav', 'Done','Done');
        switch choice
            case 'Play again'
                sound(sum(wavout,3)./max(max(abs(sum(wavout,3)))), in.fs)
            case 'Save audio as .wav'
                [filename, pathname] = uiputfile({'*.wav'},'Save as');
                if ~filename == 0 && size(waveout,3)<3
                    audiowrite([pathname,filename],wavout,in.fs);
                end
                if size(waveout,3)>2
                    disp('unable to save multiband audio as wav file')
                end
        end % switch
    end % while
end
