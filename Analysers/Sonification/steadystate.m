function OUT = steadystate(in,fs,duration)
% Steady state rendering of audio (which can be useful for sonification).
%
% This is achieved by adding random phase offsets to the spectrum 
% (interchannel phase difference is preserved).
%
% Code by Densil Cabrera
% version 1.0 (12 October 2013)
%
% INPUT ARGUMENTS
% in is audio data, or a structure containing audio data (.audio) and 
% sampling rate (.fs) . If it is a structure, then the other input
% arguments are not used.
% fs is sampling rate in Hz
% duration is the duration of the sonification in seconds

if nargin < 3
    % dialog box for settings
    prompt = {'Output audio duration (s):'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if isempty(answer)
        OUT = [];
        return
    else
        duration = str2num(answer{1,1});
    end
end
if isstruct(in)
    in = choose_from_higher_dimensions(in,2,1); 
    data = in.audio;
    fs = in.fs;
else
    data = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(data) && ~isempty(fs) && ~isempty(duration)
    % discard 3rd dimension if it exists
    data = squeeze(data(:,:,1));

    % check number of channels
    [~, chan] = size(data);

    % ensure fft length is positive and even
    fftlen = 2*ceil(abs(duration*fs)/2); 

    % do fft of BIR with zero-padding
    spectrum = fft(data,fftlen); 

    % operate from DC to Nyquist frequency
    halfspec = spectrum(1:fftlen/2+1,:); 

    % random phase offsets, evenly distributed between –pi and pi
    % (don't change phase of DC or Nyquist!)
    randphases = rand(length(halfspec)-2,1).*2*pi;
    % same random phase offsets on ch 2 preserves interchannel phase difference
    if chan == 2, randphases(:,2) = randphases(:,1); end
    % add random phase offsets
    halfspec(2:end-1,:) = abs(halfspec(2:end-1,:)).*...
        exp(1i.*(randphases + angle(halfspec(2:end-1,:))));

    % return to time domain
    y = ifft([halfspec; conj(flipud(halfspec(2:end-1,:)))]);

    % 10 ms fade-in & fade-out
    windfunc = tukeywin(fftlen, 0.01/duration); 
    if chan == 2, windfunc(:,2) = windfunc(:,1); end
    y = y .* windfunc;

    if isstruct(in)
        OUT = in;
        OUT.audio = y;
        OUT.funcallback.name = 'steadystate.m';
        OUT.funcallback.inarg = {fs,duration};
    else
        OUT = y;
    end
else
    OUT = [];
end