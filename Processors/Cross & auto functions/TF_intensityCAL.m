function OUT = TF_intensityCAL(IN)
% This function calculates the transfer function between two recordings of
% swept sinusoids. It is assumed that:
% * the two recordings have been made using identical parameters
% * a multicycle measurement has been made, including a silent cycle and
%   at least one sweep cycle
% * the recordings have not already been converted to impulse responses


OUT = [];
if ~(isfield(IN,'audio') && ...
        isfield(IN,'audio2') && ...
        isfield(IN,'fs') && ...
        isfield(IN,'properties'))
    return
end
if ~(isfield(IN.properties,'startflag') && ...
        isfield(IN.properties,'relgain'))
    return
end

% Get second audio input
IN2 = choose_audio; % call AARAE's choose_audio function
if isempty(IN2)
    return
end

if ~(isfield(IN2,'audio') && ...
        isfield(IN2,'audio2') && ...
        isfield(IN2,'fs') && ...
        isfield(IN2,'properties'))
    return
end
if ~(isfield(IN2.properties,'startflag') && ...
        isfield(IN2.properties,'relgain'))
    return
end

if IN2.fs ~= fs
    disp('Sampling rate mis-match')
    return
end

if length(IN.audio2) ~= length(IN2.audio2)
    return
end

for n = 1:length(IN.properties.startflag)
    if ~isinf(IN.properties.startflag)
        
    end
end

        
[len, chans, bands, cycles, dim5, dim6] = size(IN.audio); 
[len2, chans2, bands2, cycles2, dim5_2, dim6_2] = size(IN2.audio); 
fftlen = max([IN.startflag(2)-1 IN2.startflag(2)-1]) + length(IN.audio2);



X = fft(IN.audio,fftlen);
Y = fft(IN2.audio,fftlen);
Z = fft(IN.audio2, fftlen); % inverse filter of sweep

% calculate TF both ways.
% multiply by 0.5, since we want the TF to the midpoint between the mics
TF1 = 0.5 * X.*conj(Y) ./ (conj(Y) .* Y);
TF2 = 0.5 * X.*conj(X) ./ (conj(X) .* Y);

TF = (TF1 + TF2) ./2; % average the two TFs

% Here we could find the weak spectrum components of the inverse filter,
% and attenuate them further
Zthreshold = max(abs(Z))./1000;
Zlow = Z(abs(Z)<Zthreshold);


% apply the inverse filter of the sweep
OUT.audio = TF1 .* Z;


% return to the time domain
OUT.audio = ifft(OUT.audio);

% time domain clean-up


