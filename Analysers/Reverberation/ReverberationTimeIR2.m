function OUT = ReverberationTimeIR2(IN,startthresh,bpo,f_low,f_hi,filterstrength,phasemode,noisecomp,autotrunc,customstart,customend,doplot)
% This function calculates reverberation time and associated parameters
% from an impulse response.
%
% The function accepts up to six dimensional data, and dimensions are
% treated in the following ways:
% * Dimensions 2 (input channels), 5 (output channels) and 6 are
% power-averagerd prior to calculating reverberation time
% * Dimension 4 (cycles) is synchronously averaged
% * Dimension 3 (bands) is normally singleton as input so that filtering
% can be done. If not, filtering is not done.
%
% This function is intended to replace ReverberationTime_IR1, which was
% originally written by Grant Cuthbert and Densil Cabrera. However, the
% numerous revisions of and extensions to to the older function over some
% years made a re-organisation, simplification and re-write needed. Several
% changes to the calculation have been made, and it is highly likely that
% further changes will be made in the future (following a performance
% evaluation that was conducted in late 2015 - early 2016). 
% 



if isstruct(IN)
    
    % deal with higher dimensions
    [len,chans,bands,dim4,dim5,dim6] = size(IN.audio);
    
    % synchronously average dimension 4,
    % excluding the silent cycle if it exists
    if dim4 > 1
        if isfield(IN,'properties')
            if isfield(IN.properties,'relgain')
                IN.audio = mean(IN.audio(:,:,:,~isinf(IN.properties.relgain),:,:),4);
            else
                IN.audio = mean(IN.audio,4);
            end
        else
            IN.audio = mean(IN.audio,4);
        end
    end
    
    % if the input is multiband, then we will not filter, and we will try
    % to discover what the band frequencies are.
    if bands > 1
        dofilter = false;
        if isfield(IN,'bandID')
            bandfc = IN.bandID;
        else
            % estimate the band centre frequency from power centroid of mixdown
            % of first 1 s of the recording
            tempaudio = abs(fft(permute(mean(mean(mean(IN.audio,6),5),2),[1,3,2,4,5,6]),fs)).^2;
            f = (0:round(fc/2))'; % frequencies are at 1 Hz intervals due to fft length of fs
            bandfc = sum(tempaudio(1:length(f),:).*repmat(f,[1,bands]))./...
                sum(tempaudio(1:length(f),:));
            clear tempaudio f
        end
    else
        dofilter = true;
    end
    
    
    ir = IN.audio;
    fs = IN.fs;
else
    warndlg('This function requires an AARAE audio structure as input')
    return
end





% Dialog box
if nargin == 1
    prompt = {'Threshold for IR start detection', ...
        'Bands per octave (0 | 1 | 3)', ...
        'Lowest centre frequency (Hz)', ...
        'Highest centre frequency (Hz)', ...
        'Filter strength (use 1 for 12th order 72 dB/oct; use 2 for 24th order 144 dB/oct)', ...
        'Zero phase (0), Maximum phase (-1) or Minimum phase (1) filters',...
        'Noise compensation: None (0), Subtract final 10% Chu (1)', ...
        'Automatic end truncation: None (0), Crosspoint from nonlinear fitting (1), Crosspoint from nonlinear fitting with truncated decay energy adjustment (C in ISO3382-1) (2)',...
        'Custom evaluation range start (dB)',...
        'Custom evaluation range end (dB)',...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = [1 60];
    def = {'-20','1','125','8000',...
        '1','-1','1','2','-5','-15','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if length(answer) < 1, answer = []; end
    if ~isempty(answer)
        startthresh = str2num(answer{1,1});
        bpo = str2num(answer{2,1});
        f_low = str2num(answer{3,1});
        f_hi = str2num(answer{4,1});
        filterstrength = str2num(answer{5,1});
        phasemode = str2num(answer{6,1});
        noisecomp = str2num(answer{7,1});
        autotrunc = str2num(answer{8,1});
        customstart = str2num(answer{9,1});
        customend = str2num(answer{10,1});
        doplot = str2num(answer{11,1});
    else
        OUT = [];
        return
    end
end






% *************************************************************************
% Truncation at the start of the IR
% *************************************************************************

% Get last 10 % for Chu noise compensation if set, and to return the
% Peak to noise ratio (this needs to be done before the circular shift
% is done as part of startpoint detection/allignment)
ir_end10 = ir(round(0.9*len):end,:,:,:,:,:);

% if noisecomp == 2;
%     autotrunc = 1;
% end


m = zeros(1,chans,bands,1,dim5,dim6); % maximum value of the IR
startpoint = zeros(1,chans,bands,1,dim5,dim6); % the auto-detected start time of the IR
for d2 = 1:chans
    for d3 = 1:bands
        for d5 = 1:dim5
            for d6 = 1:dim6
                m(1,d2,d3,1,d5,d6) = max(ir(:,d2,d3,1,d5,d6).^2); % maximum value of the IR
                startpoint(1,d2,d3,1,d5,d6) = ...
                    find(ir(:,d2,d3,1,d5,d6).^2 >= m(1,d2,d3,1,d5,d6)./ ...
                    (10^(abs(startthresh)/10)),1,'first'); % Define start point
                startpoint(isempty(startpoint))=1; %for cases where startpoint was not found
                if startpoint(1,d2,d3,1,d5,d6) >1
                    % zero the data before the startpoint
                    ir(1:startpoint(1,d2,d3,1,d5,d6)-1,d2,d3,1,d5,d6) = 0;
                    % rotate the zeros to the end (to keep a constant data length)
                    ir(:,d2,d3,1,d5,d6) = circshift(ir(:,d2,d3,1,d5,d6),-(startpoint(1,d2,d3,1,d5,d6)-1));
                end
            end
        end
    end
end


try
    early50 = [zeros(fs,chans,bands,1,dim5,dim6);ir(1:1+floor(fs*0.05),:,:,:,:,:);zeros(fs,chans,bands,1,dim5,dim6)]; % Truncate Early50
    late50 = [zeros(fs,chans,bands,1,dim5,dim6);ir(ceil(fs*0.05):end,:,:,:,:,:);zeros(fs,chans,bands,1,dim5,dim6)]; % Truncate Late50
    do50 = true;
catch
    [early50, late50] = deal(NaN(1,chans,bands,1,dim5,dim6));
    do50 = false;
end

try
    early80 = [zeros(fs,chans,bands,1,dim5,dim6);ir(1:1+floor(fs*0.08),:,:,:,:,:);zeros(fs,chans,bands,1,dim5,dim6)]; % Truncate Early80
    late80 = [zeros(fs,chans,bands,1,dim5,dim6);ir(ceil(fs*0.08):end,:,:,:,:,:);zeros(fs,chans,bands,1,dim5,dim6)]; % Truncate Late80
    do80 = true;
catch
    [early80, late80] = deal(NaN(1,chans,bands,1,dim5,dim6));
    do80 = false;
end






% *************************************************************************
% Filtering (octave and 1/3-octave band)
% *************************************************************************


if dofilter == true && ((bpo ==1) || (bpo == 3))
    noctaves = log2(f_hi/f_low);
    if bpo == 1
        f_low = 1000*2.^round(log2((f_low/1000))); % make sure it is oct
        fc = f_low .* 2.^(0:round(noctaves)); % approx freqencies
    elseif bpo == 3
        fc = f_low .* 2.^(0:1/3:round(3*noctaves)/3);% approx freqencies
    end
    bandfc = exact2nom_oct(fc); % nominal frequencies
    
    if bpo == 1
        order = [12,12]*filterstrength;
        iroct = octbandfilter_viaFFT(ir,fs,bandfc,order,0,1000,0,phasemode);
        if do50
            early50 = octbandfilter_viaFFT(early50,fs,bandfc,order,0,1000,0,phasemode);
            early80 = octbandfilter_viaFFT(early80,fs,bandfc,order,0,1000,0,phasemode);
        end
        if do80
            late50 = octbandfilter_viaFFT(late50,fs,bandfc,order,0,1000,0,phasemode);
            late80 = octbandfilter_viaFFT(late80,fs,bandfc,order,0,1000,0,phasemode);
        end
        ir_end10 = octbandfilter_viaFFT(ir_end10,fs,bandfc,order,0,1000,0,phasemode);
    else % if bpo == 3
        order = [36,24] * filterstrength;
        iroct = thirdoctbandfilter_viaFFT(ir,fs,bandfc,order,0,1000,0,phasemode);
        if do50
            early50 = thirdoctbandfilter_viaFFT(early50,fs,bandfc,order,0,1000,0,phasemode);
            early80 = thirdoctbandfilter_viaFFT(early80,fs,bandfc,order,0,1000,0,phasemode);
        end
        if do80
            late50 = thirdoctbandfilter_viaFFT(late50,fs,bandfc,order,0,1000,0,phasemode);
            late80 = thirdoctbandfilter_viaFFT(late80,fs,bandfc,order,0,1000,0,phasemode);
        end
        ir_end10 = thirdoctbandfilter_viaFFT(ir_end10,fs,bandfc,order,0,1000,0,phasemode);
    end
    
    bands = size(iroct,3);
    bandfc = bandfc(1:bands);
    
else
    iroct = ir;
end
crosspoint = repmat(len-startpoint,[1,1,bands,1,1,1]);







% *************************************************************************
% Find preliminary crosspoint and truncate the end (for zero and maximum
% phase filter types only). The following constants are also used for
% finding the final crosspoint.
% *************************************************************************
% window length of smoothing filter for approximate IR envelope
C = zeros(1,chans,bands,dim4,dim5,dim6);
if ~exist('bandfc','var'), bandfc = 1000; end % used for soft-cropping
if autotrunc == 1 || autotrunc == 2
    longestwindur = 50; % in ms (50 ms is used in Lundeby)
    shortestwindur = 10; %in ms (10 ms is used in Lundeby)
    if bands > 1
        winlen = zeros(1,bands);
        fc_low = bandfc<=125;
        winlen(fc_low) = round(0.001*fs*longestwindur);
        fc_hi = bandfc>=8000;
        winlen(fc_hi) = round(0.001*fs*shortestwindur);
        fc_mid = ~(fc_low | fc_hi);
        if sum(fc_mid) >1
            winlen(fc_mid) = round(0.001*fs*linspace(longestwindur,shortestwindur,sum(fc_mid)));
        elseif sum(fc_mid) ==1
            winlen(fc_mid) = round(0.001*fs*mean([shortestwindur longestwindur]));
        end
    else
        winlen = round(0.001*fs*50); % 50 ms window
    end
    winlen = repmat(permute(winlen,[1,3,2]),[1,chans,1]);
    Tlate = zeros(1,chans,bands,dim4,dim5,dim6);
    alltimes = (0:len-1)./fs;
    Tstartoffset = -10; % level in dB relative to max from which reverberation will be fitted
    if dofilter == true && (phasemode == 0 || phasemode == -1)
        if noisecomp == 1
            noisemultiplier = 0.125; % lower the modelled noise floor by 9 dB when finding the preliminary crosspoint
        else
            noisemultiplier = 0.25; % -6 dB
        end
        % make short half-Hann window fade-outs to reduce problems with
        % discontinuities when preliminary trunation is done (soft-cropping over 4 periods)
        winperiods = 4;
        fadewin = zeros(ceil(fs/1000 + winperiods*fs/bandfc(1)),bands);
        for b = 1:bands
            winx = hann(round(2*(fs/1000 + winperiods*fs/bandfc(b))));
            winx = winx(round(end/2):end);
            fadewin(1:length(winx),b) = winx;
        end
        fadewinlen = length(fadewin);
        
        for d5 = 1:dim5
            for d6 = 1:dim6
                IR2smooth = zeros(len,chans,bands); %just for preallocation
                for b = 1:bands
                    IR2smooth(:,:,b) = fftfilt(ones(winlen(1,1,b),1)./winlen(1,1,b),iroct(:,:,b,1,d5,d6).^2);
                end
                IR2smooth(IR2smooth<=0)=1e-300;
                IR2smoothdB =  10*log10(IR2smooth);
                maxIR2smoothdB = max(IR2smoothdB);
                maxind = ones(1,chans,bands);
                for ch = 1:chans
                    for b = 1:bands
                        try
                            maxind(1,ch,b) = find(IR2smoothdB(:,ch,b) == maxIR2smoothdB(1,ch,b),1,'first');
                            IR2smoothdB(:,ch,b) = IR2smoothdB(:,ch,b) - maxIR2smoothdB(1,ch,b);
                            Tstart = [];
                            tweak = 0;
                            while isempty(Tstart) && (tweak < -Tstartoffset)
                                Tstart = find(IR2smoothdB(maxind(1,ch,b):crosspoint(1,ch,b,1,d5,d6),ch,b) <= Tstartoffset+tweak, 1, 'first');
                                Tstart = Tstart+maxind(1,ch,b)-1;
                                tweak = tweak+1;
                            end
                            times = alltimes(Tstart:crosspoint(1,ch,b,1,d5,d6));
                            levels = IR2smoothdB(Tstart:crosspoint(1,ch,b,1,d5,d6),ch,b);
                            s = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower',[-100,0],...
                                'Upper',[0,0.1],...
                                'Startpoint',[1 1]);
                            f = fittype('10*log10(10^(a*x/10)+b)','options',s);
                            [c,~] = fit(times',levels,f);
                            % find where decay intersects noise * noisemultiplier
                            x = find(10.^(c.a.*times./10) <= c.b*noisemultiplier,1,'first');
                            if ~isempty(x)
                                crosspoint(1,ch,b,1,d5,d6) = x+Tstart-1;
                            end
                            iroct(:,ch,b,1,d5,d6) = ir(:,ch,1,1,d5,d6); % replace filtered with unfiltered
                            if crosspoint(1,ch,b,1,d5,d6)+ fadewinlen <= len
                                iroct(crosspoint(1,ch,b,1,d5,d6):crosspoint(1,ch,b,1,d5,d6)+fadewinlen-1,ch,b) =...
                                    iroct(crosspoint(1,ch,b,1,d5,d6):crosspoint(1,ch,b,1,d5,d6)+fadewinlen-1,ch,b,1,d5,d6).*fadewin(:,b);
                                iroct(crosspoint(1,ch,b,1,d5,d6)+fadewinlen-1:end,ch,b,1,d5,d6) = 0;
                            else
                                win1 = fadewin(1:len-crosspoint(1,ch,b,1,d5,d6)+1,b);
                                iroct(crosspoint(1,ch,b,1,d5,d6):end,ch,b,1,d5,d6) = iroct(crosspoint(1,ch,b,1,d5,d6):end,ch,b,1,d5,d6).*win1;
                            end
                        catch
                        end
                    end
                end
            end
        end
        % Re-apply the filter-bank now that excess noise in the tail has been
        % removed
        if bpo == 1
            for b = 1:bands
                iroct(:,:,b,:,:,:) = octbandfilter_viaFFT(iroct(:,:,b,:,:,:),fs,bandfc(b),order,0,1000,0,phasemode);
            end
        elseif bpo == 3
            for b = 1:bands
                iroct(:,:,b,:,:,:) = thirdoctbandfilter_viaFFT(iroct(:,:,b,:,:,:),fs,bandfc(b),order,0,1000,0,phasemode);
            end
        end
        for ch = 1:chans
            for b = 1:bands
                for d5 = 1:dim5
                    for d6 = 1:dim6
                        iroct(crosspoint(1,ch,b,1,d5,d6):end,ch,b,1,d5,d6) = 0;
                    end
                end
            end
        end
    end
    
    
    
    
    
    % *************************************************************************
    % Find final crosspoint and truncate the end.
    % *************************************************************************
    
    if noisecomp == 1
        noisemultiplier = 0.5; % lower the modelled noise floor by 3 dB to anticipate noise power subtraction
    else
        noisemultiplier = 1; % 0 dB
    end
    for d5 = 1:dim5
        for d6 = 1:dim6
            for b = 1:bands
                IR2smooth(:,:,b) = fftfilt(ones(winlen(1,1,b),1)./winlen(1,1,b),iroct(:,:,b,1,d5,d6).^2);
            end
            IR2smooth(IR2smooth<=0)=1e-300;
            IR2smoothdB =  10*log10(IR2smooth);
            maxIR2smoothdB = max(IR2smoothdB);
            maxind = ones(1,chans,bands);
            for ch = 1:chans
                for b = 1:bands
                    try
                        maxind(1,ch,b) = find(IR2smoothdB(:,ch,b) == maxIR2smoothdB(1,ch,b),1,'first');
                        IR2smoothdB(:,ch,b) = IR2smoothdB(:,ch,b) - maxIR2smoothdB(1,ch,b);
                        Tstart = [];
                        tweak = 0;
                        while isempty(Tstart) && (tweak < -Tstartoffset)
                            Tstart = find(IR2smoothdB(maxind(1,ch,b):crosspoint(1,ch,b,1,d5,d6),ch,b) <= Tstartoffset+tweak, 1, 'first');
                            Tstart = Tstart+maxind(1,ch,b)-1;
                            tweak = tweak+1;
                        end
                        adjustedcrosspoint = round(0.95*crosspoint(1,ch,b,1,d5,d6));
                        times = alltimes(Tstart:adjustedcrosspoint);
                        levels = IR2smoothdB(Tstart:adjustedcrosspoint,ch,b);
                        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',[-100,0],...
                            'Upper',[0,0.1],...
                            'Startpoint',[1 1]);
                        f = fittype('10*log10(10^(a*x/10)+b)','options',s);
                        [c,~] = fit(times',levels,f);
                        % find where decay intersects noise * noisemultiplier
                        x = find(10.^(c.a.*times./10) <= c.b*noisemultiplier,1,'first');
                        if ~isempty(x)
                            crosspoint(1,ch,b,1,d5,d6) = x+Tstart-winlen(b);
                        end
                        iroct(crosspoint(1,ch,b,1,d5,d6):end,ch,b,1,d5,d6) = 0;
                        if autotrunc == 2
%                             % calculate the constant C (in ISO3382-1) from the
%                             % decay starting 10 dB above the crosspoint to the
%                             % crosspoint
%                             
%                             % find 10 dB above crosspoint (multiplying noise by 10
%                             % give 10 dB increase
%                             % we do this on the regression to avoid problems with
%                             % short term irregularities in the smoothed IR
%                             x10 = find(10.^(c.a.*times./10) <= c.b*noisemultiplier*10,1,'first');
%                             % linear fit to crosspoint (subtracting noise floor)
%                             subtractnoise=10.^(IR2smoothdB(:,ch,b)./10)-c.b;
%                             subtractnoise(subtractnoise<0)=0;
%                             IR2smoothdB(:,ch,b) = 10*log10(subtractnoise);
%                             if ~isempty(x10)
%                                 q = polyfit((times(x10:x))', ...
%                                     IR2smoothdB(x10:x,ch,b),1)'; % linear regression with noise subtraction
%                                 
%                                 %                 figure;
%                                 %                 plot(q(1)*times)
%                                 %                 hold on
%                                 %                 plot([zeros(Tstart(ch,band),1); 10*log10(10.^(a(ch,band)*times'./10)+b(ch,band))],'r')
%                                 %                 plot(IR2smoothdB(:,ch,band),'g')
%                                 
%                                 % determine late reverb time
%                                 Tlate(1,ch,b) = -60/q(1);
%                             end
                            
                            % sum energy from crosspoint to end - this
                            % could be done analytically instead
                            C(1,ch,b,1,d5,d6) = sum(10.^(c.a.*times(crosspoint(1,ch,b,1,d5,d6):end)./10))./fs;
                            %C = sum(10.^(q(1)*times(crosspoint(1,ch,band,d4,d5,d6):end)+q(2))./10)./fs;
                        end
                    catch
                    end
                end
            end
        end
    end
end


% *************************************************************************
% Energy ratio parameters
% *************************************************************************
alloct = permute(sum(sum(sum(sum(iroct.^2),2),5),6),[1,3,2,4,5,6]);
if do50
    early50 = permute(sum(sum(sum(sum(early50.^2),2),5),6),[1,3,2,4,5,6]);
    late50 = permute(sum(sum(sum(sum(late50.^2),2),5),6),[1,3,2,4,5,6]);
    C50 = 10*log10(early50 ./ late50); % C50
    %D50 = (early50 ./ alloct); % D50
    D50 = early50 ./ (early50+late50); % D50
else
    [C50,D50] = deal(NaN(1,bands));
end

if do80
    early80 = permute(sum(sum(sum(sum(early80.^2),2),5),6),[1,3,2,4,5,6]);
    late80 = permute(sum(sum(sum(sum(late80.^2),2),5),6),[1,3,2,4,5,6]);
    C80 = 10*log10(early80 ./ late80); % C80
    %D80 = (early80 ./ alloct); % D80
    D80 = early80 ./ (early80+late80); % D50
else
    [C80,D80] = deal(NaN(1,bands));
end

Ts = (sum(permute(sum(sum(sum(iroct.^2,2),5),6),[1,3,2,4,5,6]) .* ...
    repmat((0:len-1)'./fs,[1,bands])))./alloct; % Ts






% *************************************************************************
% Noise and truncation corrections, reverse integration & level decay
% *************************************************************************

% mean square of last 10%
ir_end10 = mean(ir_end10.^2);
% dynamic range, or peak-to-noise ratio (using last 10% as 'noise')
PNR = permute(10*log10(mean(mean(mean(max(iroct.^2)./ir_end10,2),5),6)),[1,3,2,4,5,6]);
if noisecomp ~= 1
    ir_end10 = zeros(1,chans,bands,1,dim5,dim6);
end
noisepowerfactor = 1; % this is just to facilitate experimentation
% power units, with noise subtraction
iroct = iroct.^2 - noisepowerfactor*repmat(ir_end10,[len,1,1,1,1,1]);
for ch = 1:chans
    for b = 1:bands
        for d5 = 1:dim5
            for d6 = 1:dim6
                iroct(crosspoint(1,ch,b,1,d5,d6),ch,b,1,d5,d6) =...
                    iroct(crosspoint(1,ch,b,1,d5,d6),ch,b,1,d5,d6)+...
                    C(1,ch,b,1,d5,d6); % add C to the final sample
            end
        end
    end
end
iroct(iroct<0) = 0;

% Reverse integrate squared IR,
iroct = flip(cumsum(flip(iroct,1)),1);

% normalise so max = 1
iroct = iroct./repmat(iroct(1,:,:,:,:,:),[len,1,1,1,1,1]);

% average dimensions 2, 5 and 6, and convert to decibels
% also put bands into dimension 2 (so that dimensions 3 and above are
% singleton)
levdecay = permute(10*log10(mean(mean(mean(iroct,5),6),2)),[1,3,2,4,5,6]);







% *************************************************************************
% Calculate reverberation time parameters
% *************************************************************************


% preallocate
[o,p,q,x] = deal(zeros(2,bands));
[EDT,T20,T30,Tcustom,EDTr2,T20r2,T30r2,Tcustomr2,irstart,tstart,edtend,...
    t20end,t30end,customstartind,customendind] = deal(zeros(1,bands));

for b = 1:bands
    % find the indices for the relevant start and end samples
    irstart(1,b) = find(levdecay(:,b) <= 0, 1, 'first'); % 0 dB
    tstart(1,b) = find(levdecay(:,b) <= -5, 1, 'first'); % -5 dB
    edtend(1,b) = find(levdecay(:,b) <= -10, 1, 'first'); % -10 dB
    t20end(1,b) = find(levdecay(:,b) <= -25, 1, 'first'); % -25 dB
    t30end(1,b) = find(levdecay(:,b) <= -35, 1, 'first'); % -35 dB
    customstartind(1,b) = find(levdecay(:,b) <= -abs(customstart), 1, 'first');
    customendind(1,b) = find(levdecay(:,b) <= -abs(customend), 1, 'first');
    
    % linear regression for EDT
    o(:,b) = polyfit((irstart(1,b):edtend(1,b))', ...
        levdecay(irstart(1,b):edtend(1,b),b),1)';
    EDT(1,b) = 6*((o(2,b)-10)/o(1,b)-(o(2,b)-0)/o(1,b))/fs;
    EDTr2(1,b) = corr(levdecay(irstart(1,b):edtend(1,b),b), ...
        (irstart(1,b):edtend(1,b))' * o(1,b) + o(2,b)).^2; % correlation coefficient
    
    % linear regression for T20
    p(:,b) = polyfit((tstart(1,b):t20end(1,b))', ...
        levdecay(tstart(1,b):t20end(1,b),b),1)';
    T20(1,b) = 3*((p(2,b)-25)/p(1,b)-(p(2,b)-5)/p(1,b))/fs; 
    T20r2(1,b) = corr(levdecay(tstart(1,b):t20end(1,b),b), ...
        (tstart(1,b):t20end(1,b))' * p(1,b) + p(2,b)).^2; % correlation coefficient
    
    % linear regression for T30
    q(:,b) = polyfit((tstart(1,b):t30end(1,b))', ...
        levdecay(tstart(1,b):t30end(1,b),b),1)';
    T30(1,b) = 2*((q(2,b)-35)/q(1,b)-(q(2,b)-5)/q(1,b))/fs; 
    T30r2(1,b) = corr(levdecay(tstart(1,b):t30end(1,b),b), ...
        (tstart(1,b):t30end(1,b))' * q(1,b) + q(2,b)).^2; % correlation coefficient
    
    % linear regression for Tcustom
    x(:,b) = polyfit((customstartind(1,b):customendind(1,b))', ...
        levdecay(customstartind(1,b):customendind(1,b),b),1)';
    Tcustom(1,b) = (60/(abs(customend)-abs(customstart)))...
                    *((x(2,b)-abs(customend))/x(1,b) ...
                    -(x(2,b)-abs(customstart))/x(1,b))/fs;
    Tcustomr2(1,b) = corr(levdecay(customstartind(1,b):customendind(1,b),b), ...
        (customstartind(1,b):customendind(1,b))' * x(1,b) + x(2,b)).^2; % correlation coefficient
end


% percentage ratio of T20 to T30
T20T30ratio = (T20./T30)*100;




% *************************************************************************
% Calculate low, mid and high frequency values for octave band
% reverberation time
% *************************************************************************



if  bpo == 1
    fc500 = find(bandfc == 500);
    fc1000 = find(bandfc == 1000);
    if ~isempty(fc500) && ~isempty(fc1000) && fc1000-fc500 == 1
        T30mid = mean([T30(1,fc500); T30(1,fc1000)]);
        T20mid = mean([T20(1,fc500); T20(1,fc1000)]);
        EDTmid = mean([EDT(1,fc500); EDT(1,fc1000)]);
    else
        [EDTmid,T20mid,T30mid] = deal(nan);
    end
    
    fc125 = find(bandfc == 125);
    fc250 = find(bandfc == 250);
    if ~isempty(fc125) && ~isempty(fc250) && fc250-fc125 == 1
        T30low = mean([T30(1,fc125); T30(1,fc250)]);
        T20low = mean([T20(1,fc125); T20(1,fc250)]);
        EDTlow = mean([EDT(1,fc125); EDT(1,fc250)]);
    else
        [EDTlow,T20low,T30low] = deal(nan);
    end
    
    fc2000 = find(bandfc == 2000);
    fc4000 = find(bandfc == 4000);
    if ~isempty(fc2000) && ~isempty(fc4000) && fc4000-fc2000 == 1
        T30high = mean([T30(1,fc2000); T30(1,fc4000)]);
        T20high = mean([T20(1,fc2000); T20(1,fc4000)]);
        EDThigh = mean([EDT(1,fc2000); EDT(1,fc4000)]);
    else
        [EDThigh,T20high,T30high] = deal(nan);
    end
    
    % Bass ratio
    if ~isempty(T30mid) && ~isempty(T30low)
        BR_T30 = T30low ./ T30mid;
        BR_T20 = T20low ./ T20mid;
        BR_EDT = EDTlow ./ EDTmid;
    else
        [BR_EDT,BR_T20,BR_T30] = deal(nan);
    end
    
    if ~isempty(T30mid) && ~isempty(T30high)
        TR_T30 = T30high ./ T30mid;
        TR_T20 = T20high ./ T20mid;
        TR_EDT = EDThigh ./ EDTmid;
    else
        [TR_EDT,TR_T20,TR_T30] = deal(nan);
    end
    
end




%--------------------------------------------------------------------------
% AARAE TABLE
%--------------------------------------------------------------------------
OUT.tables = [];
f = figure('Name','Reverberation Parameters', ...
    'Position',[200 200 620 360]);
dat1 = [EDT;T20;T30;Tcustom;...
    C50;C80;D50;D80;Ts; ...
    EDTr2;T20r2;T30r2;Tcustomr2;...
    T20T30ratio;PNR];
if bpo == 1 || bpo == 3
    cnames1 = num2cell(bandfc);
elseif bpo == 0 && bands == 1
    cnames1 = {'Broadband'};
elseif length(bandfc) == bands
    cnames1 = num2cell(bandfc);
end
rnames1 = {'Early decay time (s)',...
    'Reverberation time T20 (s)',...
    'Reverberation time T30 (s)',...
    ['Reverberation time ' num2str(-abs(customstart)) ' to ' num2str(-abs(customend)) ' dB (s)'],...
    'Clarity index C50 (dB)',...
    'Clarity index C80 (dB)',...
    'Definition D50',...
    'Definition D80',...
    'Centre time Ts (s)',...
    'Correlation coefficient EDT r^2',...
    'Correlation coefficient T20 r^2',...
    'Correlation coefficient T30 r^2',...
    'Correlation coefficient Tcustom r^2',...
    'Ratio of T20 to T30 %%',...
    'Dynamic range: Peak to last 10%% (dB)'};
t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
set(t1,'ColumnWidth',{60});

if bpo == 1
    dat2 = [EDTlow,EDTmid,EDThigh,BR_EDT,TR_EDT;...
        T20low,T20mid,T20high,BR_T20,TR_T20;...
        T30low,T30mid,T30high,BR_T30,TR_T30];
    cnames2 = {'Low Freq','Mid Freq','High Freq','Bass Ratio','Treble Ratio'};
    rnames2 = {'EDT', 'T20', 'T30'};
    t2 =uitable('Data',dat2,'ColumnName',cnames2,'RowName',rnames2);
    set(t2,'ColumnWidth',{90});
    [~,tables] = disptables(f,[t1 t2],{'Reverb parameters','Reverb parameters summary'});
else
    [~,tables] = disptables(f,t1,{'Reverb parameters'});
end
OUT.tables = [OUT.tables tables];




if doplot
figure('Name','Level Decay and Regression Lines')
[r,c] = subplotpositions(bands+1,0.4);
levdecayend = len./fs;
for b = 1:bands
    subplot(r,c,b)
    hold on
    
    % plot the level decay(s) on a single subplot
    plot(((1:len)-1)./fs, levdecay(:,b),'Color',[0.2 0.2 0.2], ...
        'LineStyle',':','DisplayName','Level Decay')
    if (noisecomp == 0) || (noisecomp == 1) || (noisecomp == 2)
        % linear regression for T30
        plot(((tstart(1,b):t30end(1,b))./fs), ...
            (tstart(1,b):t30end(1,b)).* ...
            q(1,b)+q(2,b), ...
            'Color',[0 0 0.6],'DisplayName', 'T30')
        
        % linear regression for T20
        plot(((tstart(1,b):t20end(1,b))./fs), ...
            (tstart(1,b):t20end(1,b)).* ...
            p(1,b)+p(2,b), ...
            'Color',[0 0.6 0],'DisplayName','T20')
        
        % linear regression for EDT
        plot(((irstart(1,b):edtend(1,b))./fs), ...
            (irstart(1,b):edtend(1,b)).* ...
            o(1,b)+o(2,b), ...
            'Color',[0.9 0 0],'DisplayName','EDT')
    end
    

    
    % x axis label (only on the bottom row of subplots)
    if b > (c*r - c)
        xlabel('Time (s)')
    end
    
    % y axis label (only on the left column of subplots)
    if mod(b-1, c) == 0
        ylabel('Level (dB)')
    end
    
    xlim([0 levdecayend])
    ylim([-65 0])
    
    % text on subplots
    text(levdecayend*0.45,-5,...
        ['EDT ',num2str(0.01*round(100*EDT(1,b)))],'Color',[0.9 0 0])
    
    text(levdecayend*0.45,-10,...
        ['T20 ',num2str(0.01*round(100*T20(1,b)))],'Color',[0 0.6 0])
    
    text(levdecayend*0.45,-15,...
        ['T30 ',num2str(0.01*round(100*T30(1,b)))],'Color',[0 0 0.6])
    
    
    % subplot title
    if bpo == 1 || bpo == 3
        title([num2str(bandfc(b)),' Hz'])
    elseif bpo == 0 && bands == 1
        title('Broadband')
    elseif length(bandfc) == bands
        title([num2str(bandfc(b)),' Hz'])
    end
    
end % for band

% DIY legend
subplot(r,c,b+1)

plot([0.1,0.4], [0.8,0.8],'Color',[0.2 0.2 0.2], ...
    'LineStyle',':','DisplayName','Level Decay')
xlim([0,1])
ylim([0,1])
hold on
text(0.5,0.8,'Decay');
plot([0.1,0.4], [0.6,0.6],'Color',[0.9 0 0],'DisplayName','EDT')
text(0.5,0.6,'EDT');
plot([0.1,0.4], [0.4,0.4],'Color',[0 0.6 0],'DisplayName','T20')
text(0.5,0.4,'T20');
plot([0.1,0.4], [0.2,0.2],'Color',[0 0 0.6],'DisplayName', 'T30')
text(0.5,0.2,'T30');
set(gca,'YTickLabel','',...
    'YTick',zeros(1,0),...
    'XTickLabel','',...
    'XTick',zeros(1,0))
hold off
end


doresultleaf(levdecay,'Level [dB]',{'Time'},...
                'Time',      ((1:len)-1)./fs,  's',           true,...
                'Frequency', num2cell(bandfc), 'Hz',          false,...
                'name','Reverb_time');

            
OUT.funcallback.name = 'ReverberationTimeIR2.m';
OUT.funcallback.inarg = {startthresh,bpo,f_low,f_hi,filterstrength,phasemode,noisecomp,autotrunc,customstart,customend,doplot};
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2015,2016, Densil Cabrera
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
