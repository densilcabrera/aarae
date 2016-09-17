%
%       Dynamic Compressive Gammachirp Filterbank 
%       Toshio IRINO
%       Created:   6 Sep 2003
%       Modified:  7 Jun 2004
%       Modified:  12 Jul 2004  (PpgcEstShiftERB)
%       Modified:  14 Jul 2004  (LinPpgc)
%       Modified:  4  Aug 2004  (introducing GCresp)
%       Modified:  16 Aug 2004  (ExpDecayVal)
%       Modified:  31 Aug 2004  (introducing GCFBv2_SetParam)
%       Modified:  8  Sep 2004 (TTS. tidy up the names. 2.00 -> 2.01)
%       Modified:  10 Sep 2004 (Normalization at Level estimation path)
%       Modified:  7 Oct 2004   (c2val is level dependent 2.02)
%       Modified:  22 Oct 2004  (level estimation  2.03) 
%       Modified:  8 Nov 2004   (error detection of SndIn)
%       Modified:  30 Nov 2004  (c2val control)
%       Modified:  23 May 2005  (v205. Pc == average of two input, RMS2dBSPL,
%   			 Fast filtering when 'fix' : under construction)
%       Modified:  24 May 2005  (v205 Mod in LinLvl1 =..., LvldB= ...)
%       Modified:   3 Jun 2005  (v205)
%       Modified:   1 Jun 2005  (v205, GCparam.GainCmpnstdB)
%       Modified:  14 Jul 2005  (v205, GCparam.LvlEst.RefdB, Pwr, Weight)
%       Modified:  15 Sep 2005  (v205, rename GCparam.LvlRefdB --> GainRefdB)
%       Modified:   7 Apr 2006  (v206, Compensation of Group delay OutMidCrct)
%       Modified:  16 Apr 2006  (v206, Minimum phase OutMidCrct: NoGD-cmpnst)
%       Modified:  27 Jun 2006  (v206, GCresp.GainFactor)
%       Modified:  22 Dec 2006  (v207, speed up for 'fix' condition)
%       Modified:   7 Jan 2007  (v207, output GCresp.Fp2 for 'fix' condition)
%       Modified:  19 Jan 2007  (v207, GCparam.Ctrl: 'static'=='fixed')
%       Modified:   5 Aug 2007  (v207, GCresp.GainFactor --> vector)
%       Modified:  19 Dec 2011  (v208, non-sample-by-sample coefficient update)
%       Modified:  18 Dec 2012  (v208, no AF-HP update in Level-estimation path)
%       Modified:  19 Dec 2012  (v208, clean up level-estimation path)
%
% function [cGCout, pGCout, GCparam, GCresp] = GCFB2(Snd,GCparam)
%      INPUT:  Snd:    Input Sound
%              GCparam:  Gammachirp parameters 
%                  GCparam.fs:     Sampling rate          (48000)
%                  GCparam.NumCh:  Number of Channels     (75)
%                  GCparam.FRange: Frequency Range of GCFB [100 6000]
%                           specifying asymptotic freq. of passive GC (Fr1)
%        
%      OUTPUT: cGCout:  Compressive GammaChirp Filter Output
%              pGCout:  Passive GammaChirp Filter Output
%              Ppgc:    power at the output of passive GC
%              GCparam: GCparam values
%              GCresp : GC response result
%
% Note
%   1)  This version is completely different from GCFB v.1.04 (obsolete).
%       We introduced the "compressive gammachirp" to accomodate both the 
%       psychoacoustical simultaneous masking and the compressive 
%       characteristics (Irino and Patterson, 2001). The parameters were 
%       determined from large dataset (See Patterson, Unoki, and Irino, 2003.)
%
%
% References: 
%  Irino,T. and Unoki,M.:  IEEE ICASSP98, pp.3653-3656, May 1998.
%  Irino,T. and Patterson,R.D. :  JASA, Vol.101, pp.412-419, 1997.
%  Irino,T. and Patterson,R.D. :  JASA, Vol.109, pp.2008-2022, 2001.
%  Patterson,R.D., Unoki,M. and Irino,T. :  JASA, Vol.114,pp.1529-1542,2003.
%  Irino,T. and and Patterson,R.D. : IEEE Trans.ASLP, Vol.14, Nov. 2006.
%  
%
function [cGCout, pGCout, GCparam, GCresp] = GCFBv208(SndIn,GCparam)

%%%% Handling Input Parameters %%%%%
if nargin < 2,         help GCFBv208;           end;    
Tstart0 = clock;

[nc, LenSnd] = size(SndIn);
if nc ~= 1, 
error('Check SndIn. It should be 1 ch (Monaural) and  a single row vector.' ); 
end;

GCparam = GCFBv208_SetParam(GCparam);

%%%%%%%%%%%%%
fs   = GCparam.fs;
NumCh = GCparam.NumCh;
if length(GCparam.b1) == 1 & length(GCparam.c1) == 1
  b1 = GCparam.b1(1);  % Freq. independent 
  c1 = GCparam.c1(1);  % Freq. independent 
else
  error('Not prepared yet: Freq. dependent b1, c1');
end;

[Fr1, ERBrate1]  = EqualFreqScale('ERB',NumCh,GCparam.FRange);
Fr1 = Fr1(:);
ERBspace1 = mean(diff(ERBrate1));

GCresp.Fr1 = Fr1;
Fp1 = Fr2Fpeak(GCparam.n,b1,c1,Fr1);
GCresp.Fp1 = Fp1;

[ERBrate ERBw] = Freq2ERB(Fr1);
[ERBrate1kHz ERBw1kHz] = Freq2ERB(1000);
Ef = ERBrate/ERBrate1kHz - 1;
GCresp.Ef = Ef;

b2val = GCparam.b2(1,1)*ones(NumCh,1) + GCparam.b2(1,2)*Ef(:);
c2val = GCparam.c2(1,1)*ones(NumCh,1) + GCparam.c2(1,2)*Ef(:); 
GCresp.b2val = b2val;
GCresp.c2val = c2val;

%%%%
%  fratVal = frat(1,1)      + frat(1,2)*Ef + 
%            frat(2,1)*Ppgc + frat(2,2)*Ef*Ppgc;
%%%

%%%%% Outer-Mid Ear Compensation %%%%
if length(GCparam.OutMidCrct) > 2
  disp(['*** Outer/Middle Ear correction (minimum phase) : ' ...
        GCparam.OutMidCrct ' ***']);
  CmpnOutMid = OutMidCrctFilt(GCparam.OutMidCrct,fs,0,2); % 2) minimum phase
  % 1kHz: -4 dB, 2kHz: -1 dB, 4kHz: +4 dB (ELC)
  % Now we use Minimum phase version of OutMidCrctFilt (modified 16 Apr. 2006).
  % No compensation is necessary.  16 Apr. 2006

  Snd = filter(CmpnOutMid,1,SndIn);
else 
  disp('*** No Outer/Middle Ear correction ***');
  Snd = SndIn;
end;

% for inverse filer,  use OutMidCrctFilt('ELC',fs,0,1);

%%%%% Gammachirp  %%%
disp('*** Gammmachirp Calculation ***');
if 0,  disp(GCparam), end;
tic;

  
SwFastPrcs = 1; % ON:  Fast Processing for static filter 
if SwFastPrcs ~= 1, error('SwFastPrcs should be 1.'); end;
if SwFastPrcs == 1  &  strncmp(GCparam.Ctrl,'sta',3) == 1  
   % error('Fast processing for linear cGC:  Under construction.');
   %%%%%%%%%%% for HP-AF %%%

   LvldB = GCparam.GainRefdB;
   fratVal = GCparam.frat(1,1) + GCparam.frat(1,2)*Ef(:) + ...
            (GCparam.frat(2,1) + GCparam.frat(2,2)*Ef(:))*LvldB;
   Fr2val = Fp1(:).*fratVal;
   GCresp.Fr2 = Fr2val;
   [ACFcoefFastPrcs] = MakeAsymCmpFiltersV2(fs,Fr2val,b2val,c2val);   

else % HP-AF for dynamic-GC level estimation path. 18 Dec 2012 Checked
   Fr2LvlEst = GCparam.LvlEst.frat * Fp1(:);
   [ACFcoefLvlEst] = ...
   MakeAsymCmpFiltersV2(fs,Fr2LvlEst,GCparam.LvlEst.b2, GCparam.LvlEst.c2);

end;


%%%% Start calculation %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Passive Gammachirp & Level estimation filtering  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Tstart  = clock;
  cGCout  = zeros(NumCh, LenSnd);
  pGCout  = zeros(NumCh, LenSnd);
  Ppgc    = zeros(NumCh, LenSnd); 
  cGCoutLvlEst = zeros(NumCh, LenSnd);
  disp('--- Channel-by-channel processing ---');

  for nch=1:NumCh

    % passive gammachirp
    pgc = GammaChirp(Fr1(nch),fs,GCparam.n,b1,c1,0,'','peak'); % pGC
    pGCout(nch,1:LenSnd)=fftfilt(pgc,Snd);       % fast fft based filtering 
                                                 
    %%% Fast processing for fixed cGC: checked. %%%
    if SwFastPrcs == 1  &  strncmp(GCparam.Ctrl,'sta',3) == 1   
      StrGC = 'Static (Fixed) Compressive-Gammachirp';
      GCout1 = pGCout(nch,:);
      for Nfilt = 1:4;
        GCout1 = filter(ACFcoefFastPrcs.bz(nch,:,Nfilt), ...
                        ACFcoefFastPrcs.ap(nch,:,Nfilt), GCout1);
      end;
      cGCout(nch,:) = GCout1;
      GCresp.Fp2(nch) = Fr1toFp2(GCparam.n,b1,c1,b2val(nch),c2val(nch),...
                                 fratVal(nch),GCresp.Fr1(nch));
      if nch == NumCh, GCresp.Fp2 = GCresp.Fp2(:); end;
      
    else  %  Level estimation pass for Dynamic.  18 Dec 2012 Checked
      StrGC = 'Passive-Gammachirp & Level estimation filter';
      GCout1 = pGCout(nch,:);
      for Nfilt = 1:4;
        GCout1 = filter(ACFcoefLvlEst.bz(nch,:,Nfilt), ...
                        ACFcoefLvlEst.ap(nch,:,Nfilt), GCout1);
      end;
      cGCoutLvlEst(nch,:) = GCout1;    

    end;

    if nch == 1 | rem(nch,20)==0  
         disp([StrGC ' ch #' num2str(nch) ' / #' num2str(NumCh) ...
    '.    elapsed time = ' num2str(fix(etime(clock,Tstart)*10)/10) ' (sec)']);
    end;  

  end;

%
% Passive filter -->  jump to Gain Normalization
%

% Dynamic filter here
if strncmp(GCparam.Ctrl,'dyn',3) == 1  % Dynamic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dynamic Compressive Gammachirp filtering  %%%%%
%%%%% Sample-by-sample processing               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Initial settings %%%%%%%%%%%%%%%
nDisp          = 20*fs/1000; % display every 20 ms
cGCout         = zeros(NumCh,LenSnd);
GCresp.Fr2     = zeros(NumCh,LenSnd);
GCresp.fratVal = zeros(NumCh,LenSnd);
GCresp.Fp2     = []; % No output

LvldB          = zeros(NumCh,LenSnd);
LvlLinPrev     = zeros(NumCh,2);
ExpDecayVal    = exp(-1/(GCparam.LvlEst.DecayHL*fs/1000)*log(2)); % decay exp.
NchShift       = round(GCparam.LvlEst.LctERB/ERBspace1);
NchLvlEst      = min(max(1, (1:NumCh)'+NchShift),NumCh);  % shift in NumCh [1:NumCh]
LvlLinMinLim   = 10^(-GCparam.LvlEst.RMStoSPLdB/20); % minimum should be 0 dBSPL
LvlLinRef      = 10.^(( GCparam.LvlEst.RefdB - GCparam.LvlEst.RMStoSPLdB)/20); 

%%%%% Sample-by-sample processing %%%%%%%
disp('--- Sample-by-sample processing ---');
Tstart = clock;
 for nsmpl=1:LenSnd

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Level estimation circuit %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%% Modified:  24 May 05 
      LvlLin(1:NumCh,1) = ...
         max([max(pGCout(NchLvlEst,nsmpl),0), LvlLinPrev(:,1)*ExpDecayVal]')';
      LvlLin(1:NumCh,2) = ...
         max([max(cGCoutLvlEst(NchLvlEst,nsmpl),0), LvlLinPrev(:,2)*ExpDecayVal]')';
      LvlLinPrev = LvlLin;
      
      %%%%% Modified: 14 July 05
      LvlLinTtl = GCparam.LvlEst.Weight * ...
          LvlLinRef.*(LvlLin(:,1)/LvlLinRef).^GCparam.LvlEst.Pwr(1) ...
       + ( 1 - GCparam.LvlEst.Weight ) * ...
          LvlLinRef.*(LvlLin(:,2)/LvlLinRef).^GCparam.LvlEst.Pwr(2);

      LvldB(:,nsmpl) = 20*log10( max(LvlLinTtl,LvlLinMinLim) ) ...
                                   + GCparam.LvlEst.RMStoSPLdB; 

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%% Signal path %%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Filtering High-Pass Asym. Comp. Filter
     fratVal = GCparam.frat(1,1) + GCparam.frat(1,2)*Ef(:) + ...
              (GCparam.frat(2,1) + GCparam.frat(2,2)*Ef(:)).*LvldB(:,nsmpl);
     Fr2Val = Fp1(:).*fratVal;
 
     if rem(nsmpl-1, GCparam.NumUpdateAsymCmp) == 0 % update periodically
        [ACFcoef] = MakeAsymCmpFiltersV2(fs,Fr2Val,b2val,c2val);  
     end;

     if nsmpl == 1, 
       [dummy,ACFstatus] =  ACFilterBank(ACFcoef,[]);  % initiallization
     end;

     [SigOut,ACFstatus] = ACFilterBank(ACFcoef,ACFstatus,pGCout(:,nsmpl));
     cGCout(:,nsmpl) = SigOut;
     GCresp.Fr2(:,nsmpl) = Fr2Val;
     GCresp.fratVal(:,nsmpl) = fratVal;
     % Derivation of GCresp.Fp2 is too time consuming.
     % please use CalFp2GCFB.m

    if nsmpl==1 | rem(nsmpl,nDisp)==0,
     %%% [  20*log10([max(LvlLin(:,1)) max(LvlLin(:,2)) max(LvlLinTtl) ])...
     %%%  + GCparam.LvlEst.RMStoSPLdB      max(LvldB(:,nsmpl))]
     disp(['Dynamic Compressive-Gammachirp: Time ' int2str(nsmpl/fs*1000) ...
	    '(ms) / ' int2str(LenSnd/fs*1000) '(ms).  elapsed time = ' ...
	    num2str(fix(etime(clock,Tstart)*10)/10) ' (sec)']);
    end;
 end;

end; % if strncmp(GCparam.Ctrl,'dyn',3) == 1  % Dynamic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Signal path Gain Normalization at Reference Level (GainRefdB) %%%
%%%%  for static & dynamic filters                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fratRef = GCparam.frat(1,1) + GCparam.frat(1,2)*Ef(:) + ...
          (GCparam.frat(2,1) + GCparam.frat(2,2)*Ef(:)).*GCparam.GainRefdB;

cGCRef = CmprsGCFrsp(Fr1,fs,GCparam.n,b1,c1,fratRef,b2val,c2val); 
GCresp.cGCRef = cGCRef;
GCresp.LvldB  = LvldB;

GCresp.GainFactor = 10^(GCparam.GainCmpnstdB/20)*cGCRef.NormFctFp2;
cGCout = (GCresp.GainFactor*ones(1,LenSnd)).*cGCout;

%%%%%%%%%%


return









%%%% Trush %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 19 Dec 2012  Checked OK
%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%% Level estimation path %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Fr2LvlEst = GCparam.LvlEst.frat * Fp1LvlEst;
      % no need to update for every sampling point 75->54,   18 Dec 2012
      %if rem(nsmpl-1, GCparam.NumUpdateAsymCmp) == 0 % update periodically
      %    [ACFcoefLvlEst] = ...
      %    MakeAsymCmpFiltersV2(fs,Fr2LvlEst,GCparam.LvlEst.b2, GCparam.LvlEst.c2);
      % end;

%tlap(1) = etime(clock,Tstart1);
      
%      if nsmpl == 1,       %%  initialization 
%        [ACFcoefLvlEst] = ...
%         MakeAsymCmpFiltersV2(fs,Fr2LvlEst,GCparam.LvlEst.b2, GCparam.LvlEst.c2);
%        [dummy,ACFstatusLvlEst] = ACFilterBank(ACFcoefLvlEst,[]); 
%      end;
%tlap(2) = etime(clock,Tstart1);
%
%      [cGCLvlEst,ACFstatusLvlEst] =...
%	 ACFilterBank(ACFcoefLvlEst,ACFstatusLvlEst,pGCout(NchLvlEst,nsmpl));



  
% if SwFastPrcs ~= 1  |  strncmp(GCparam.Ctrl,'sta',3) ~= 1  % Dynamic

%   if strncmp(GCparam.Ctrl,'sta',3) == 1  
%      LvldB(:,nsmpl) = GCparam.GainRefdB*ones(NumCh,1); % fixed value

%   elseif strncmp(GCparam.Ctrl,'dyn',3) == 1,  % when dynamic (time-varying)

%   else 
%    error([ 'GCparam.Ctrl should be "fix[ed]" or "dyn[amic]/tim[e-varying]"'])
%   end;
% end; % if SwFastPrcs ~= 1  |  strncmp(GCparam.Ctrl,'sta',3) ~= 1  % Dynamic


%%%%

% Obsolete, 5 Aug 07
%GCresp.GainFactor ...
%    = 10^(GCparam.GainCmpnstdB/20)*(cGCRef.NormFctFp2 * ones(1,LenSnd));
%cGCout = GCresp.GainFactor.*cGCout;



  % ----- obsolete -----
  % CmpnOutMid = OutMidCrctFilt(GCparam.OutMidCrct,fs,0);
  % You need to be careful about the group delay caused by this process.
  % about 6.27 ms (= helf of 12.5 ms total length of linear phase filter)
  % --> introducing Group delay compensation on 7 Apr. 2006
  %
  %     
  % SwGDcmpnst = 1;  % change here if you need no compensation.
  % if SwGDcmpnst == 1, % Group-delay compensation ON
  %  LenGD = fix(length(CmpnOutMid)/2);
  %  SndIn = [SndIn, zeros(1,LenGD)];
  %  Snd = filter(CmpnOutMid,1,SndIn);
  %  Snd = Snd((LenGD+1):end);
  %  % [LenGD (length(SndIn)-LenGD) length(Snd)]
  % disp(['*** Outer/Middle Ear correction: Group delay compensation ON ***']);
  % else
  % --------------------

  
