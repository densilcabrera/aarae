%
%      Test GCFBv207 (for GCFBpack) 
%      Toshio IRINO
%      Created:   15 Sep 2005 (for v205)
%      Modified:   7 Apr 2006 (v206, Compensation of Group delay OutMidCrct)
%      Modified:  23 Dec 2006 (v207, 'dynamic' rather than 'time-varying')
%      Modified:  18 Mar 2007 (check. modifying title)
%      Modified:  25 Jul 2007 (GCresp)
%      Modified:   5 Aug 2007 (Elapsed time)
%      Modified:  26 Jun 2013
%
%
%clear

%%%% Stimuli : a simple pulse train %%%%
fs = 48000;
Tp = 10; % (ms) 100 Hz pulse train
Snd1 = [1, zeros(1,Tp*fs/1000-1)];
Snd = [];

for nnn = 1:10,
%for nnn = 1:3,
    Snd = [Snd1 Snd];
end;
Tsnd = length(Snd)/fs;
disp(['Duration of sound = ' num2str(Tsnd*1000) ' (ms)']);


SigSPLlist = [40:20:80]

for nfig = 1:2
figure(nfig);

for nnn = 1:3
SigSPL = SigSPLlist(nnn);
Snd =  Eqlz2MeddisHCLevel(Snd,SigSPL);

%%%% GCFB %%%%
GCparam.fs     = fs;
GCparam.NumCh  = 100;
GCparam.FRange = [100, 6000];

GCparam.OutMidCrct = 'No';
%GCparam.OutMidCrct = 'ELC';

if nfig == 1, GCparam.Ctrl = 'dynamic'; % used to be 'time-varying'
else    GCparam.Ctrl = 'static'; % or 'fixed' 
end;

tic
[cGCout, pGCout, GCparam, GCresp] = GCFBv208(Snd,GCparam);
tm = toc;
disp(['Elapsed time is ' num2str(tm,4) ' (sec) = ' ...
      num2str(tm/Tsnd,4) ' times RealTime.']);

subplot(3,1,nnn)
imagesc(max(cGCout,0));
set(gca,'YDir','normal');
title(['GCFB control = "' GCparam.Ctrl '";  Signal Level = ' ...
       int2str(SigSPL) ' dB SPL']);

GCparam 
GCresp
drawnow
end;
end;
