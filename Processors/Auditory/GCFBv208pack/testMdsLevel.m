%
%  test Level set for Meddis IHC
%  Irino T.
%  9 June 2004
%
%clear

% 50 dB sinusoid for Meddis IHC
fs = 48000;
fsig = 1000;
Snd50dBSPL = MkCalibTone(fs,fsig,50,1000);

OutLeveldB = 50;
[SndEq AmpdB ] = Eqlz2MeddisHCLevel(Snd50dBSPL,fs,OutLeveldB);

tms = (0:length(SndEq)-1)/fs*1000;
plot(tms,SndEq,tms,SndEq+10);

AmpdB
disp('The level of Snd50dBSPL is little smaller than the intended one');
disp('because the effect of the taper window.');
return

%%%%%%%%%%
%v = 0.01:1:1000;

v0 = 10*sqrt(2);
a = v0^0.8;

bb = v0/log10(v0*10);
CmprsLog = bb*max(0,log10(max(v*10,0.1)));

loglog(v, a*v.^0.2, v, CmprsLog);

IndB = 0:100;
OutdB = BMCmprsFuncGM(IndB);

a = 10*sqrt(2)*10.^((IndB-50)/20);
b = BMCmprsFuncGCFB(a);
plot(IndB,OutdB,20*log10(a)+27,20*log10(b)+54.5,'--');

