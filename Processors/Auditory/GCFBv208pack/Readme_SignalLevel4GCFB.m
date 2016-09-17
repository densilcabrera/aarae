%
%  Introductory script for signal level normalization for the GCFB
%  Irino, T.
%  Original : 5 Sep 2006
%  Updated  : 5 Sep 2006
%
%  Usage : execute this file in matlab
%          matlab> Readme_SignalLevel4GCFB
%
clf

disp(' ');
disp('Readme_SignalLevel4GCFB:');
disp('This script describes a method to equalize the signal level ');
disp('for the gammachirp filterbank GCFBv2.xx.');
disp(' ');
disp('The default signal level for the GCFBv2.xx is the same as ');
disp('that defined by Meddis(1986) and Schroeder and Hall (1974). ');
disp('>>>  30 dB SPL <--> RMS level of 1.0 ');
disp('>>>  90 dB SPL <--> RMS level of 1000.0');
disp(' ');
disp(['Any arbitrary signal can be equalized using "Eqlz2MeddisHCLevel.m"']);
disp(' ');

disp('Return to continue: ');
pause
     disp(' ');
     disp(' ');

disp(['*** Example 1 (sin wave) ***']);
     echo on
     echo on
     fs = 48000;         % (Hz)  sampling rate
     fsig = 1000;        % (Hz)  frequency
     SignalLeveldB = 60; % (dB SPL) re. Meddis HC
     Tsnd = 3;           % (sec) duration
     amp  = 0.9;         % max amplitude 
     Snd1 = amp*sin(2*pi*fsig*(0:Tsnd*fs-1)/fs);  % sinusoid
     [SndEqM, AmpdB] = Eqlz2MeddisHCLevel(Snd1,SignalLeveldB);
     nn = 1:0.01*fs;
     plot((nn-1)/fs*1000, SndEqM(nn));
     xlabel('Time (ms)');
     ylabel('Amplitude');
     RmsLeveldB = 20*log10(sqrt(mean(SndEqM.^2)));
     disp(['RMS Level (re. Meddis HC) = ' num2str(RmsLeveldB+30) ' dB SPL.']);
     echo off

     disp(' ');
     disp('The values in the second argument "AmpdB": ')
     disp(sprintf('   %7.3f    %7.3f      %7.3f', AmpdB))
     disp(['      |          |            | ']);
     disp(['OutputLevel_dB   |            | ']);
     disp(['      CompensationValue_dB    | ']);
     disp(['                        SourceLevel_dB (re. Meddis HC)']);
     disp(' ');
     disp('You may use AmpdB(2) or "CompensationValue_dB" for amplification');
     disp('of other signals which are normalized in advance.');
     disp('For example, the values in .wav files are restricted ');
     disp('within the range between -1 and 1. The common scaling');
     disp('factor may be defined by using "CompensationValue_dB".');
     disp(' ');
     disp(' ');

disp('Return to continue: ');
pause
     disp(' ');
     disp(' ');

disp(['*** Example 2 (white noise) ***']);
     echo on
     echo on
     Snd2 = randn(1,Tsnd*fs);
     [SndEqM2, AmpdB] = Eqlz2MeddisHCLevel(Snd2,SignalLeveldB);
     echo off
     RmsLeveldB4 = 20*log10(sqrt(mean(SndEqM2.^2)));
     disp(['RMS Level (re. Meddis HC) = ' num2str(RmsLeveldB4+30) ' dB SPL.']);
     disp('The values in the second argument "AmpdB": ')
     disp(sprintf('   %7.3f    %7.3f      %7.3f', AmpdB))
     disp(['      |          |            | ']);
     disp(['OutputLevel_dB   |            | ']);
     disp(['     CompensationValue_dB     | ']);
     disp(['                        SourceLevel_dB (re. Meddis HC)']);
     disp(' ');

