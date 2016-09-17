%
%  Making / Showing Calibration Tone 
%  Meddis(1986) == Schroeder and Hall (1974)
%    30 dB SPL <--> RMS level of 1.0
%    90 dB SPL <--> RMS level of 1000.0
%  Slaney's Gammatone filter Level (described in MeddisHairCell.m)
%    0  dB SPL <--> RMS level of 1.0
%    60 dB SPL <--> RMS level of 1000.0
%
%  Irino, T.
%  Original : 10 May 2004
%  Updated  :  9 Jun 2004
%  Updated  : 14 Aug 2004  (Renamed MkCalibTone.m -> MkCalibTone4AFB.m)
%
%  INPUT:  fs: sampling frequency (default 16 kHz);
%          fsig: signal frequnecy (default 1000 Hz);
%          SignalLeveldB: (default 60 dB); 
%          Tsnd : Duration of sound in ms (default 500 ms);
%          SwRef: Reference 'Meddis'(Default)  or 'Slaney' 
%
%  OUTPUT: SndCalib: Calibration Sound
%
function SndCalib = MkCalibTone4AFB(fs,fsig,SignalLeveldB,Tsnd,SwRef)

if nargin < 1, fs = []; end;
if length(fs) == 0, fs = 16000; end;
if nargin < 2, fsig = []; end;
if length(fsig) == 0, fsig = 1000; end;
if nargin < 3, SignalLeveldB = []; end;
if length(SignalLeveldB) == 0, SignalLeveldB = 60; end;
if nargin < 4, Tsnd = []; end;
if length(Tsnd) == 0, Tsnd = 500; end;
if nargin < 5, SwRef = []; end;
if length(SwRef) == 0, SwRef = 'Meddis'; end;

LenSnd = Tsnd*fs/1000;
disp(['MkCalibTone4AFB (' SwRef ') : ' ...
num2str(fsig) '-Hz Tone, Sampling freq. = ' num2str(fs) ' (Hz)']);

ValdB = SignalLeveldB;
if strcmp(upper(SwRef(1:3)),'MED') == 1, 
 ValdB = SignalLeveldB - 30; 
end;

Amp = sqrt(2)*10^(ValdB/20);  % RMS value of 1
t = (0:LenSnd-1)/fs;
Snd = Amp*sin(2*pi*fsig*t);
plot(t*1000,Snd);
xlabel('Time (ms)');
ylabel('Amplitude');
rmsVal = sqrt(mean(Snd.^2));
title([ SwRef ' : 1 kHz Tone: ' int2str(SignalLeveldB)  ...
        ' dB SPL = RMS Level of ' num2str(rmsVal) ...
        ' =  Max Level of ' num2str(max(abs(Snd)))]);
axis([0 50 [-1.2 1.2]*Amp]);

SndCalib  = Snd.*TaperWindow(LenSnd,'han',10*fs/1000);

%SndPlay1 = SndCalib/(2^15-1); % normalization by 16bit as an example
%sound(SndPlay1,fs);

return


%%%% Comment in Meddis (1986) and Shroeder and Hall (1976) %%%%%
%"s(t) is normalized such that ave(s^2(t)) = 1 corresponds to a sound level
%of approximately 30 dB re 0.0002 dyn/cm^2." (= 20 uPa. 1 dyne = 10^(-5) N)

%%%% Comment from MeddisHairCell.m %%%%%

h = 50000;      % This parameter scales the discharge rate. 
                % Adjust as necessary.
                % In combination with the gammatone filterbank (ERBFilterBank),
                % h=50000 will produce a steady-state average discharge
                % probability of about 135 spikes/s within the 1kHz channel,
                % for an input consisting of a 1 kHz sinewave at 60 dB SPL
                % (0 dB SPL corresponds to an RMS level of 1.0 at the
                % input of the gammatone filter).  Scaling and constant
                % courtesy of Alain de Cheveigne'

