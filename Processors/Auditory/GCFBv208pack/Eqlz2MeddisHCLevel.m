%    Equalizing Signal RMS Level to the Level for MeddisHairCell
%    Irino, T.
%    Created:   9 Jun. 2004
%    Modified:  9 Jun. 2004
%    Modified:  23 Sept 2005 (remove fs. fs has been remained for some reason.)
%  Eqlz2MeddisHCLevel(Snd,fs,OutLeveldB) --> Eqlz2MeddisHCLevel(Snd,OutLeveldB)
%
%    function [SndEqM, AmpdB] = Eqlz2MeddisHCLevel(Snd,OutLeveldB);
%    INPUT  Snd: input sound
%           OutLeveldB : Output level (No default value,  RMS level)
%
%    OUTPUT SndEqM: Equalized Sound (rms value of 1 is 30 dB SPL)
%           AmpdB: 3 values in dB 
%                  [OutputLevel_dB, CompensationValue_dB, SourceLevel_dB]
%
% Ref: Meddis (1986), JASA, 79(3),pp.702-711.
%
% rms(s(t)) == sqrt(mean(s.^2)) == 1   --> 30 dB SPL
% rms(s(t)) == sqrt(mean(s.^2)) == 10  --> 50 dB SPL
% rms(s(t)) == sqrt(mean(s.^2)) == 100 --> 70 dB SPL
%
function [SndEqM, AmpdB] = Eqlz2MeddisHCLevel(Snd,OutLeveldB,dummy);

if nargin < 2, help Eqlz2MeddisHCLevel; end;
% if nargin < 3, OutLeveldB = []; end;  
% if length(OutLeveldB) == 0, OutLeveldB = 50; end; % for speech
% No default value! 
if nargin > 2 | OutLeveldB > 120 % for checking inconsistency
  disp('Eqlz2MeddisHCLevel was modified to take 2 input arguments. Sept2005.')
    error('function [SndEqM, AmpdB] = Eqlz2MeddisHCLevel(Snd,OutLeveldB)');
end;

SourceLevel = sqrt(mean(Snd.^2))*10^(30/20); % level in terms of Meddis Level

Amp = (10^(OutLeveldB/20))/SourceLevel;
SndEqM = Amp * Snd; 

AmpdB = [OutLeveldB  20*log10([Amp, SourceLevel])];
