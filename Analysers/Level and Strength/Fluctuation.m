function [OUT, varargout] = Fluctuation(IN,soundfield)
% This function attempts to implement a novel fluctuation strength model
% for modulated and unmodulated signals in free and diffuse sound fields as
% detailed in Tingting Zhou et al. "A Model for Calculating
% Psychoacoustical Fluctuation Strength", JAES, September 2015.
%
% The model uses 75 ERB channels to calculate generalized modulation depth
% (GMD) values, which are converted into specific fluctuation strengths.
% Test signal: 60dB 1kHz 100% mod 4Hz = 1 vacil
%
% INPUT ARGUMENTS:
% IN - audio signal (1 channel only)
% soundfield - free field = 0; diffuse field = 1
%
% OUTPUTS:
% For a single frame signals (up to 10 seconds long):
% - Total Fluctuation strength
% - Specific fluctuation strength (vacil/ERB channel)
%
% For multiple frame signals (more than 10 seconds long):
% - Total Fluctuation strength
% - Fluctuation strength statistics
% - Time-averaged specific Fluctuation strength
% - Time-varying Fluctuation strength
%
% MATLAB code by Ella Manor for AARAE (2017)
% *************************************************************************
if isstruct(IN)

   IN = choose_from_higher_dimensions(IN,1,1); % only 1 channel analysis at present
   audio = IN.audio;
   fs = IN.fs;     

   if isfield(IN,'cal')
       cal = IN.cal;
   else
       warndlg('Calibration data missing - please calibrate prior to calling this function.','AARAE info','modal');
       IN = cal_aarae(IN);
       cal = IN.cal;
   end
   if isfield(IN,'name')
       name = IN.name;
   else
       name = '';
   end
else
   audio = IN(:,1,1,1,1,1);
   cal = cal(1);
   name = '';
end
% *************************************************************************
if nargin == 1

   param = inputdlg({'Sound field: Free (0); Diffuse (1)';...
                     'Sound signal: Unmodulated (0); Modulated (1)'},...
                     'Parameters',...
                     [1 50],...
                     {'0';'1'});

   if length(param) < 2, param = []; end


   if ~isempty(param)

       soundfield = str2double(param(1));
       signaltype = str2double(param(2));
   else
       % get out of here if the user presses 'cancel'
       OUT = [];
       return
   end

else
   param = [];
end
% *********************************************************************
if ~isempty(audio) && ~isempty(fs)

   calconstant = 0;
   cal = cal-calconstant;
   if length(cal) == 1
       audio = audio .* 10.^(cal/20);
   else
       audio(:,1) = audio(:,1) .* 10.^(cal(1)/20);
       audio(:,2) = audio(:,2) .* 10.^(cal(2)/20);
   end
   disp(['rms level of the entire wave ', num2str(10*log10(mean(audio.^2)+10e-99)+calconstant), ' dB'])
[samples,nchan] =  size(audio);
   if nchan > 1
       audio = audio(:,1);
       disp('Only one channel signals are analysed by FluctStrength.')
       nchan = 1;
   end
if fs ~= 48000
   gcd_fs = gcd(48000,fs); % greatest common denominator
   audio = resample(audio,48000/gcd_fs,fs/gcd_fs);
   fs = 48000;
   samples = size(audio,1);
end   
% *********************************************************************
% Signal preparation
n = 10; % 10 seconds window
N = n*fs; % convert to samples
N = 2^nextpow2(N); % for faster processing
if length(audio) < N
   audio = repmat(audio,ceil(N/length(audio)),1);
   audio = audio(1:N);
else
    nn = ceil(samples/N);
    zeropad = (N*nn)-samples;
    audio = [audio;zeros(zeropad,1)];   
end
N_b = buffer(audio,N);
nFrames = size(N_b,2);
window = blackman(N);
% *********************************************************************
% Filters design and weighting functions
%
filt_bandpass = [0 0.;
   (0.001:0.001:0.099)' zeros(99,1);
                   0.1 0.01;
                   0.2 0.1;
                   0.3 0.13;
                   0.4 0.2;
                   0.5 0.25;
                   0.6 0.35;
                   0.7 0.36;
                   0.8 0.38;
                   0.9 0.43;
                   1. 0.45;
                   2 0.71;
                   3 0.83;
                   4 1.;
                   5 0.9;
                   6 0.75;
                   7 0.63;
                   8 0.48;
                   9 0.42;
                   10 0.38;
                   20 0.17;
                   30 0.09;
                   40 0.05;
                   50 0.02;
                   60 0.01;
                   70 0.;
                   80 0.;
                   90 0.;
                   100 0.;
                   125 0.;
                   160 0.;
                   200 0.;
                   250 0.;
                   315 0.;
                   400 0.;
                   500 0.];
mag = interp1(filt_bandpass(:,1),filt_bandpass(:,2),0:0.1:500); % interpolate from 0 Hz to Nyquist
mag = minphasefreqdomain([mag flip(mag(2:end-1))]',120,1); % make minmum phase
bpFilt = ifft(mag); % transform to time domain
bpFilt = bpFilt(1:end/2); % optionally discard the filter's tail
wf = hann(2*length(bpFilt)); % optionally fade-out using half-Hann
wf = wf(end/2+1:end);
bpFilt = bpFilt .* wf;
%  
% figure;
% subplot(2,1,1)
% plot((0:length(bpFilt)-1)./1000, bpFilt);
% subplot(2,1,2)
% TF = fft(bpFilt);
% f = 1000*(0:round(length(TF)/2))/length(TF);
% semilogx(f, abs(TF(1:length(f))))

% Calculate centre frequencies based on erb number scale
ERBchn = 75;
ERBnum = linspace((21.366*log10(4.368e-3*50+1)),(21.366*log10(4.368e-3*15000+1)),ERBchn);
cfs = (10.^(ERBnum/21.366)-1)/4.368e-3;
% Weighting coefficients
g0(:,1) = cfs;
g0(:,2)=[0.85 0.91 0.95 1.03 1.1 1.16 1.18 1.2 1.223 1.228,...
   1.233 1.238 1.243 1.248 1.25 1.26 1.27 1.275 1.3 1.32,...
   1.34 1.33 1.34 1.343 1.348 1.35 1.34 1.343 1.35 1.345,...
   1.344 1.343 1.342 1.34 1.338 1.333 1.33 1.32 1.31 1.3,...
   1.29 1.285 1.28 1.27 1.26 1.25 1.235 1.23 1.225 1.21,...
   1.205 1.2 1.17 1.16 1.14 1.11 1.09 1.07 1.05 1.04,...
   1.035 1. 0.99 0.95 0.93 0.92 0.91 0.86 0.84 0.82 0.79 0.76 0.75 0.74 0.71];

% Transfer functions for the outer and middle ear.
s = 4095;
if soundfield == 0 % free field
% Free field to eardrum transfer function (frontal incinece). First column
% = frequnecy band in Hz. Second column = level at the eardrum - level
% measured in free field in dB [ANSI S3.4-2007]
filt_outer = [0 0.;
               20 0.;
               25 0.;
               31.5 0.;
               40 0.;
               50 0.;
               63 0.;
               80 0.;
               100 0.;
               125 0.1;
               160 0.3;
               200 0.5;
               250 0.9;
               315 1.4;
               400 1.6;
               500 1.7;
               630 2.5;
               750 2.7;
               800 2.6;
               1000 2.6;
               1250 3.2;
               1500 5.2;
               1600 6.6;
               2000 12.;
               2500 16.8;
               3000 15.3;
               3150 15.2;
               4000 14.2;
               5000 10.7;
               6000 7.1;
               6300 6.4;
               8000 1.8;
               9000 -0.9;
               10000 -1.6;
               11200 1.9;
               12500 4.9;
               14000 2.0;
               15000 -2.0;
               16000 2.5;
               20000 2.5];
b_outer = fir2(s,[filt_outer(:,1)./(0.5*fs);1],[db2mag(filt_outer(:,2));0]);
elseif soundfield == 1 % diffuse field
% Diffuse field to eardrum transfer function. First column
% = frequnecy band in Hz. Second column = level at the eardrum - level
% measured in diffuse field in dB [ANSI S3.4-2007]
filt_outer = [0 0.;
               20 0.;
               25 0.;
               31.5 0.;
               40 0.;
               50 0.;
               63 0.;
               80 0.;
               100 0.;
               125 0.1;
               160 0.3;
               200 0.4;
               250 0.5;
               315 1.;
               400 1.6;
               500 1.7;
               630 2.2;
               750 2.7;
               800 2.9;
               1000 3.8;
               1250 5.3;
               1500 6.8;
               1600 7.2;
               2000 10.2;
               2500 14.9;
               3000 14.5;
               3150 14.4;
               4000 12.7;
               5000 10.8;
               6000 8.9;
               6300 8.7;
               8000 8.5;
               9000 6.2;
               10000 5.;
               11200 4.5;
               12500 4.0;
               14000 3.3;
               15000 2.6;
               16000 2.;
               20000 2.];
b_outer = fir2(s,[filt_outer(:,1)./(0.5*fs);1],[db2mag(filt_outer(:,2));0]);
end
% Transfer function of middle ear. First column
% = frequnecy band in Hz. Second column = effective level at the cochlea
% relative to the level at the eardrum in dB. Values of 14kHz-20kHz have
% not been validated.
filt_middle = [0 -40;
           20 -39.6;
           25 -32.;
           31.5 -25.85;
           40 -21.4;
           50 -18.5;
           63 -15.9;
           80 -14.1;
           100 -12.4;
           125 -11.;
           160 -9.6;
           200 -8.3;
           250 -7.4;
           315 -6.2;
           400 -4.8;
           500 -3.8;
           630 -3.3;
           750 -2.9;
           800 -2.6;
           1000 -2.6;
           1250 -4.5;
           1500 -5.4;
           1600 -6.1;
           2000 -8.5;
           2500 -10.4;
           3000 -7.3;
           3150 -7.;
           4000 -6.6;
           5000 -7.;
           6000 -9.2;
           6300 -10.2;
           8000 -12.2;
           9000 -10.8;
           10000 -10.1;
           11200 -12.7;
           12500 -15.;
           14000 -18.2;
           15000 -23.8;
           16000 -32.3;
           18000 -45.5
           20000 -50.];           

b_middle = fir2(s,[filt_middle(:,1)./(0.5*fs);1],[db2mag(filt_middle(:,2));0]);

HTres= [0 130;
           20 70;
           25 65;
           31.5 58;
           40 48;
           50 40;
           63 35;
           80 30;
           100 25;
           125 22;
           160 18;
           200 12;
           250 10;
           315 9;
           400 8;
           500 5;
           630 3;
           750 3;
           800 3;
           1000 3;
           1250 3;
           1500 2;
           1600 1;
           2000 0;
           2500 -3;
           3000 -4;
           3150 -4;
           4000 -5;
           5000 -1;
           6000 6;
           6300 7;
           8000 15;
           9000 17;
           10000 15;
           11200 13;
           12500 10;
           14000 13;
           15000 19;
           16000 48;
           18000 60
           20000 130]; 

f=fs*(0:N/2)./N;
MinExcdB = interp1(HTres(:,1),HTres(:,2),f);

% *********************************************************************
% Main analyser
Ftotal = zeros(1,nFrames);
FI = zeros(nFrames,ERBchn);
TimePoints = zeros(1,nFrames);

for currentFrame = 1:nFrames

   dataIn = fft((N_b(:,currentFrame) .* window),N); % apply window to current frame
   filteredOuter = filtfilt(b_outer,1,filtfilt(b_middle,1,dataIn)); % apply outer and middle ear filters  
   filteredOuterSpect = real(ifft(filteredOuter));
   filteredOuterNyq = filteredOuterSpect(1:N/2+1);
   filteredOuterNyq(2:end-1) = 2*filteredOuterNyq(2:end-1);
   dataOut = filteredOuterNyq'; 
   
   % ignore components below hearing threshold 
   spectrumdata = abs(dataOut);
   LdB = mag2db(spectrumdata);
   whichL  = find(LdB>MinExcdB); 
   dataLeft = dataOut(1,whichL(1,:));
   
   % calculate the envelope for each ERB channel
   [bm,env] = gammatoneFast(dataLeft,cfs,fs); % includes ERB overlapped channels filtering  
   [chno,bmlength] = size(bm);
   bmidx = 1:1:bmlength;

    if signaltype == 0
        env = zeros(chno,bmlength);
        for k = 1:chno
            [pks,locs] = findpeaks(bm(k,:));
            env(k,:) = interp1(locs,pks,bmidx,'spline','extrap');
        end
    end
    
   % decimate the signal to fs=1kHz after calculating the envelopes 
   envDecimated = [];
   for k = 1:chno
       tmp1 = decimate(env(k,:),4,12);
       tmp2 = decimate(tmp1,6,12);
       envDecimated(k,:) = decimate(tmp2,2,12);
   end

   % calculate modulation depth
   [h0,hBPrms,mdept] = deal(zeros(1,chno));
  
   hBPi   = zeros(chno,length(envDecimated)+length(bpFilt)-1); % time-domain filtering    

   for k = 1:chno
       etmp      = envDecimated(k,:);
       h0(k)     = mean(etmp);
       etmpp     = etmp-h0(k);
       hBPi(k,:) = conv(bpFilt,etmpp);
       hBPrms(k) = rms(hBPi(k,:));

       if h0(k) > 0 
           mdept(k) = hBPrms(k)/h0(k);
           if mdept(k) > 1
               mdept(k) = 1.;
           end
       else
           mdept(k) = 0;
       end   
   end
   
%    calculate correlation coefficients
   ki = zeros(1,chno-2);
   for k=1:chno-2
       cfac = cov(hBPi(k,:),hBPi(k+2,:));
       den  = diag(cfac);
       den  = sqrt(den*den');
       if den(2,1) > 0
           ki(k) = cfac(2,1)/den(2,1);
       else
           ki(k) = 0;
       end
   end
   
   % calculate specific fluctuation
   fi = zeros(1,chno);
   for k = 1:chno
       Gzi = g0(k,2);
       md  = mdept(k); % delta_fi
       if k == chno-1 || k == chno
           Ki = 1;
       else
           Ki = ki(k);
       end
       if k == 1 || k == 2
           Ki2 = 1;
       else
           Ki2 = ki(k-2);
       end
       fi(k) = (Gzi*md*(Ki*Ki2)^2);
   end

   % calculate total fluctuation
    Cal = 0.05; % constant used so that 1.00 vacil = 1kHz AM at 4Hz modfreq and 100% modulation rate calibrated to 60dB
    Ftotal(currentFrame) = Cal * sum(fi);
    TimePoints(currentFrame) = currentFrame*(n*fs);
    FI(currentFrame,:) = fi;
end
TimePointsT = TimePoints./fs;
mean_Ftotal = mean(Ftotal);
Ftotal = [Ftotal 0];
TimePointsT = [0 TimePointsT];

%%
% *************************************************************************
% Data Presentation
% *************************************************************************
if isstruct(IN)
   % ********* TABLES *********
   LineWidth = 1.5;
    
   if nFrames > 1
   % Fluctuation statistics
   Fmean = mean(Ftotal);
   Fstd = std(Ftotal);
   Fmax = max(Ftotal);
   F1 = prctile(Ftotal,99);
   F2 = prctile(Ftotal,98);
   F3 = prctile(Ftotal,97);
   F4 = prctile(Ftotal,96);
   F5 = prctile(Ftotal,95);
   F10 = prctile(Ftotal,90);
   F20 = prctile(Ftotal,80);
   F30 = prctile(Ftotal,70);
   F40 = prctile(Ftotal,60);
   F50 = median(Ftotal);
   F60 = prctile(Ftotal,40);
   F70 = prctile(Ftotal,30);
   F80 = prctile(Ftotal,20);
   F90 = prctile(Ftotal,10);
   Fmin = min(Ftotal);

   dataF = [Fmean;Fstd;Fmax;F1;F2;F3;F4;F5;F10;F20;F30;F40;F50;F60;F70;F80;F90;Fmin];

   % generate tables of results

   fig1 = figure('Name','Time-varying Fluctuation Strength Statistics');
   table1 = uitable('Data',dataF,...
       'ColumnName',{'Fluctuation'},...
       'RowName',{'Mean','Standard deviation','Maximum',...
       'F1','F2','F3','F4',...
       'F5','F10','F20','F30','F40','F50 (median)','F60',...
       'F70','F80','F90','Minimum'});

   [~,tables] = disptables(fig1,table1); % AARAE function

   OUT.tables = tables;

   % ********* CHARTS *********
   % Figure for charts
   figure('Name',['Fluctuation Strength of ',name])

   subplot(2,1,1)
    
   % Time-varying Fluctuation
   stairs(TimePointsT,Ftotal,'r-','LineWidth',LineWidth);
   title ('Time-Varying Total Fluctuation');
   xlabel('Time (10 seconds frames)');
   ylabel('Fluctuation (vacil)');

   % Specific Fluctuation spectrogram
   subplot(2,1,2)
   imagesc(TimePointsT, 1:1:ERBchn,FI');
       cH = colorbar('southoutside');
       cH.Label.String = 'Fluctuation (vacil)';
   set(gca,'YDir','normal');
   ax=gca;
   axis tight;
   ax.Title.String = 'Specific Fluctuation';
   ax.XLabel.String = 'Time (10 seconds frames)';
   ax.YLabel.String = 'ERB channel';
   hold off;
   
   % ******** AARAE RESULTS LEAVES **********

   % Time-varying Fluctuation results leaf
   doresultleaf(Ftotal','Fluctuation [vacil]',{'time'},...
       'Time (10 seconds frames)',TimePointsT','s',true,...
       'fluctuationtype', {'Fluctuation over time'}, 'categorical',[],...
       'name','Time_varying_fluctuation');

   % Time-averaged Specific Fluctuation results leaf
   doresultleaf(FI','Specific fluctuation [vacil/erb]',{'Equivalent rectangular bandwidth'},...
       'Equivalent rectangular bandwidth',1:ERBchn,'erb',true,...
       'fluctuationtype', {'Fluctuation over erb'}, 'categorical', [],...
       'name','Time_averaged_specific_fluctuation');

   % Specific fluctuation spectrogram results leaf
   doresultleaf(FI','Specific Fluctuation [vacil/erb]',{'time','Equivalent rectangular bandwidth'},...
       'Equivalent rectangular bandwidth',1:ERBchn,'erb',true,...
       'Time',TimePointsT','s',true,...
       'fluctuationtype', {'Fluctuation over erb'}, 'categorical', [],...
       'name','Time_varying_specific_fluctuation');
   OUT.funcallback.name = 'FluctStrength.m';
   OUT.funcallback.inarg = {fs,cal};

   else
       % generate tables of results

   fig1 = figure('Name','Fluctuation Strength');
   table1 = uitable('Data',Ftotal(1),...
       'ColumnName',{'Fluctuation'},...
       'RowName',{'Overall Fluctuation'});

   [~,tables] = disptables(fig1,table1); % AARAE function

   OUT.tables = tables;

   % Time-averaged Fluctuation as a function of erb number
   % figure
   figure('Name',['Fluctuation Strength of ',name])
   plot((1:ERBchn), FI,'r-','LineWidth',LineWidth);
   ax=gca;
   ax.Title.String = 'A single frame Fluctuation';
   ax.XLabel.String = 'ERB channel';
   ax.YLabel.String = 'Specific Fluctuation (vacil/erb)';
   hold off;

   % Time-averaged Specific Fluctuation results leaf
   doresultleaf(FI','Specific fluctuation [vacil/erb]',{'Equivalent rectangular bandwidth'},...
       'Equivalent rectangular bandwidth',1:ERBchn,'erb',true,...
       'fluctuationtype', {'Fluctuation over erb'}, 'categorical', [],...
       'name','Time_averaged_specific_fluctuation');

   end

else
   % output numbers only
   OUT = mean_Ftotal; % fluctuation
   varargout{1} = Ftotal; % time-varying specific roughness
   varargout{2} = FI; % mean specific roughness
   varargout{3} = TimePointsT; % time
  varargout{4} = 1:ERBchn; % critical band rate (for specific roughness)
end
end
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.