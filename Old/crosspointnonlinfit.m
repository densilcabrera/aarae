function [crosspoint,Tlate,C] = crosspointnonlinfit(IR2, fs, fc, noisemultiplier)
% This function attempts to find the point at which a reverberant decay
% is affected by a steady state noise floor by non-linear curve fitting of
% the smoothed IR.
%
% IR2 is the squared IR (please square prior to calling this function!)
%
% fs is audio sampling rate in Hz
%
% fc is the band centre frequencies in Hz
%
% noisemultiplier is used to multiply the noise power (in the non-linear
% model) to find the crosspoint. If it is 1, then the crosspoint is where
% the exponential decay is equal to the noise. If it is 2, then the
% crosspoint is earlier, where the decay is equal to the noise power x 2
% (i.e. noise + 3 dB).

if ~exist('noisemultiplier','var'), noisemultiplier = 10; end % x10 corresponds to 10 dB margin

[len,chans,bands,dim4,dim5,dim6] = size(IR2);
%crosspoint = deal(round(0.9*len)*ones(1,chans,bands,dim4,dim5,dim6));
crosspoint = repmat(len,[1,chans,bands,dim4,dim5,dim6]);
% derive a smoothed envelope function for each band
% lopassfreq = 4; % smoothing filter cutoff frequency in Hz
% halforder = 1; % smoothing filter order
% Nyquist = fs/2;
% [num, den] = butter(halforder, lopassfreq/Nyquist, 'low');
% %envelopes = 10*log10(filtfilt(num, den, RIRoct .^2));
% envelopes = filtfilt(num, den, 10*log10(IR2));
% % increased filter order for lower frequency bands
% if bands>2
%     for n = 2:bands-1
%         envelopes(:,:,1:end-n) = filtfilt(num, den, envelopes(:,:,1:end-n));
%     end
% end

longestwindur = 50; % in ms (50 ms is used in Lundeby)
shortestwindur = 10; %in ms (10 ms is used in Lundeby)
if ~exist('fc','var') && bands > 1
    winlen = round(0.001*fs*linspace(longestwindur,shortestwindur,bands));
elseif exist('fc','var') && bands > 1
    winlen = zeros(1,bands);
    fc_low = fc<=125;
    winlen(fc_low) = round(0.001*fs*longestwindur);
    fc_hi = fc>=8000;
    winlen(fc_hi) = round(0.001*fs*shortestwindur);
    fc_mid = ~(fc_low | fc_hi);
    if sum(fc_mid) >1
        winlen(fc_mid) = round(0.001*fs*linspace(longestwindur,shortestwindur,sum(fc_mid)));
    elseif sum(fc_mid) ==1
        winlen(fc_mid) = round(0.001*fs*mean([shortestwindur longestwindur]));
    end
else
    winlen(1:bands) = round(0.001*fs*25);
end
winlen = repmat(permute(winlen,[1,3,2]),[1,chans,1]);
Tlate = zeros(1,chans,bands,dim4,dim5,dim6);
for d4 = 1:dim4
    for d5 = 1:dim5
        for d6 = 1:dim6
            % 1. AVERAGE SQUARED IR IN LOCAL TIME INTERVALS
            IR2smooth = zeros(len,chans,bands); %just for preallocation
            for b = 1:bands
                IR2smooth(:,:,b) = fftfilt(ones(winlen(1,1,b),1)./winlen(1,1,b),IR2(:,:,b,d4,d5,d6));
            end
            IR2smooth(IR2smooth<=0)=1e-300;
            IR2smoothdB =  10*log10(IR2smooth);
            maxIR2smoothdB = max(IR2smoothdB);
            maxind = ones(1,chans,bands);
            for ch = 1:chans
                for b = 1:bands
                    maxind(1,ch,b) = find(IR2smoothdB(:,ch,b) == maxIR2smoothdB(1,ch,b),1,'first');
                    IR2smoothdB(:,ch,b) = IR2smoothdB(:,ch,b) - maxIR2smoothdB(1,ch,b);
                end
            end
            
            
            
            %IR2smoothdB = IR2smoothdB - repmat(max(IR2smoothdB),[len,1,1]); % make max = 0 dB
            % preallocate
            maxsample = zeros(1, chans, bands);
            Tstart = zeros(1, chans, bands);
            a = zeros(1, bands);
            b = zeros(1, bands);
            alltimes = (0:length(IR2smoothdB)-1)./fs;
            Tstartoffset = -10;
            for ch = 1:chans
                for band = 1:bands
                    maxsample(ch,band) = find(IR2smoothdB(:,ch,band) == 0, 1, 'last');
                    Tstart(ch,band) = find(IR2smoothdB(maxsample(ch,band):end,ch,band) <= Tstartoffset, 1, 'first');
                    Tstart(ch,band) = Tstart(ch,band)+maxsample(ch,band);
                    times = alltimes(Tstart(ch,band):end);
                    s = fitoptions('Method','NonlinearLeastSquares',...
                        'Lower',[-100,0],...
                        'Upper',[0,0.1],...
                        'Startpoint',[1 1]);
%                    s = fitoptions('Method','NonlinearLeastSquares');
                    f = fittype('10*log10(10^(a*x/10)+b)','options',s);
                    [c,~] = fit(times',IR2smoothdB(Tstart(ch,band):end,ch,band),f);
                    a(ch,band) = c.a;
                    b(ch,band) = c.b;
                    
                    % find where decay intersects noise * noisemultiplier
                    x = find(10.^(a(ch,band).*times./10) <= b(ch,band)*noisemultiplier,1,'first');
                    if ~isempty(x)
                        crosspoint(1,ch,band,d4,d5,d6) = x+Tstart(ch,band)-1;
                    else
                        x = len;
                    end
                    
                    % calculate the constant C (in ISO3382-1) from the
                    % decay starting 10 dB above the crosspoint to the
                    % crosspoint
                    
                    % find 10 dB above crosspoint (multiplying noise by 10
                    % give 10 dB increase
                    % we do this on the regression to avoid problems with
                    % short term irregularities in the smoothed IR
                    x10 = find(10.^(a(ch,band).*times./10) <= b(ch,band)*noisemultiplier*10,1,'first');
                    
                    % linear fit to crosspoint (subtracting noise floor)
                    subtractnoise=10.^(IR2smoothdB(:,ch,band)./10)-b(ch,band);
                    subtractnoise(subtractnoise<0)=0;
                    IR2smoothdB(:,ch,band) = 10*log10(subtractnoise);
                    if ~isempty(x10)
                    q = polyfit((times(x10:x))', ...
                    IR2smoothdB(x10:x,ch,band),1)'; % linear regression with noise subtraction
                
%                 figure; 
%                 plot(q(1)*times)
%                 hold on
%                 plot([zeros(Tstart(ch,band),1); 10*log10(10.^(a(ch,band)*times'./10)+b(ch,band))],'r')
%                 plot(IR2smoothdB(:,ch,band),'g')
                
                    % determine late reverb time
                    Tlate(1,ch,band) = -60/q(1);
                    end
                    
                    % sum energy from crosspoint to end
                    C = sum(10.^(a(ch,band).*times(crosspoint(1,ch,band,d4,d5,d6):end)./10))./fs;
                    %C = sum(10.^(q(1)*times(crosspoint(1,ch,band,d4,d5,d6):end)+q(2))./10)./fs;
                    
                    
                end
            end
        end
    end
end