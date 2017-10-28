function OUT = THD_via_ESS_v2(IN,octsmooth,thresholddB,winlen,winpow1,winpow2,subtractnoise)
% This function analyses harmonic distortion based on a multicycle
% exponential sweep recording or impulse response
%
% Use an exponential sweep as the test signal, potentially at different
% gains, and include the silent cycle. A 60 s sweep is recommended.
%
% It is strongly advised to obtain a high signal to noise ratio, so that
% distortion can be effectively separated from noise.
%
% This function accepts either the recording itself, or an IR stacked in
% dimension 4.
%
% THD is calculated from the first 6 harmonics only


if ~isstruct(IN)
    warndlg('This function currently only works for AARAE structures');
    OUT = [];
    return
end

% currently this function only works on dimensions 1,2 and 4
IN = choose_from_higher_dimensions(IN,4,1);
if size(IN.audio,3) > 1
    IN.audio = sum(IN.audio,3);
    disp('Multiband audio has been summed.');
end

% convert to IR if it has not been done already
if isfield(IN,'audio2')
    IN = convolveaudiowithaudio2(IN,1,0,1);
end
IR = IN.audio;
fs = IN.fs;
%sweeplen = round(IN.properties.sig_len * fs);
%gaplen = round(IN.properties.gapdur * fs);

[~,chans,~,dim4,dim5,dim6] = size(IR);
if dim4 == 1 && chans > 1
    IR = permute(IR,[1,4,3,2]);
    disp('Assuming IRs are stacked in dimension 2 instead of dimension 4');
end
if isfield(IN,'chanID')
    chanID = IN.chanID;
end
if isfield(IN,'properties') % check for required special fields
    if isfield(IN.properties,'freq') ...
            && isfield(IN.properties,'dur')
        T=IN.properties.dur;% Extract the length of the sweep.
        freqs=IN.properties.freq; %Extract the highest and lowest frequencies in the sweep
    else
        warndlg('The input audio structure is missing required ''properties'' fields (freq and dur) for THD analysis','AARAE info','modal');
        %disp('The input audio structure is missing required fields for THD analysis')
        OUT = [];
        return
    end
else
    warndlg('The input audio structure is missing required ''properties'' fields (freq and dur) for THD analysis','AARAE info','modal');
    %disp('The input audio structure is missing required fields for THD analysis')
    OUT = [];
    return
end

if isfield(IN,'properties') % check for relgain
    if isfield(IN.properties,'relgain')
        relgain=IN.properties.relgain; % Extract the relative gain of the sweeps
    elseif size(IN.audio,4) ==1
        relgain = 0; % we don't need relgain if there is only one cycle
    else
        warndlg('The input audio structure is missing required ''properties'' field ''relgain'' for multicycle THD analysis','AARAE info','modal');
        OUT = [];
        return
    end
elseif size(IN.audio,4) ==1
    relgain = 0;
else
    warndlg('The input audio structure is missing required ''properties'' field ''relgain'' for multicycle THD analysis','AARAE info','modal');
    OUT = [];
    return
end
    

% dialog box
if nargin == 1
    param = inputdlg({'Fractional octave smoothing (use 0 for no smoothing, 3 for 1/3-octave smoothing, 1 for octave smoothing)';...
        'Threshold above background noise to display as distortion [dB]';
        'Analysis window duration [s]';
        'Window fade-in power';
        'Window fade-out power';
        'Subtract noise [0 | 1]'},...
        'Silence Sweep Analysis',...
        [1 60],...
        {'0','0','1','4','1','0'});
    
    param = str2num(char(param));
    
    if length(param) < 6, param = []; end
    if ~isempty(param)
        octsmooth = param(1);
        thresholddB = param(2);
        winlen = round(param(3) * fs);
        winpow1 = param(4);
        winpow2 = param(5);
        subtractnoise = param(6);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end

% Find IR peak index for each channel and cycle
if isinf(relgain(1))
    peakind = zeros(size(IR,2),size(IR,4)-1);
    for ch = 1:size(IR,2)
        for d4 = 2:size(IR,4)
            peakind(ch,d4-1) = find(IR(:,ch,1,d4).^2==max(IR(:,ch,1,d4).^2));
        end
    end
else
    peakind = zeros(size(IR,2),size(IR,4));
    for ch = 1:size(IR,2)
        for d4 = 1:size(IR,4)
            peakind(ch,d4) = find(max(IR(:,ch,1,d4).^2)==IR(:,ch,1,d4).^2);
        end
    end
end


if dim5 > 1
    if isfield(IN,'dim5ID')
        dim5ID = IN.dim5ID;
    else
        dim5ID = makechanID(dim5,10);
    end
else
    dim5ID = {''};
end

if dim6 > 1
    if isfield(IN,'dim6ID')
        dim6ID = IN.dim6ID;
    else
        dim6ID = makechanID(dim6,20);
    end
else
    dim6ID = {''};
end

[~,chans,~,cycles,dim5,dim6] = size(IR);

IRstartindex = round(mean(mean(peakind))-winlen/2); % improve this, e.g. by cross correlation with harmonic offsets
wf = window(@blackmanharris,winlen);
wf(1:round(end/2)) = wf(1:round(end/2)).^winpow1;
wf((round(end/2)+1):end) = wf((round(end/2)+1):end).^winpow2;
fftlen = max([winlen fs]);
% ANALYSE IR FROM SWEEP, INCLUDING HARMONIC DISTORTION
harmonicoffsets=round((log10(1:6))./...
    (1./IN.properties.sig_len...
    .*log10(IN.properties.freq(2)./IN.properties.freq(1))));
sweepIR1 = IR(IRstartindex:IRstartindex+winlen-1,:,:,:,:,:);
sweepIR2 = IR(IRstartindex-harmonicoffsets(2):IRstartindex-harmonicoffsets(2)+winlen-1,:,:,:,:,:);
sweepIR3 = IR(IRstartindex-harmonicoffsets(3):IRstartindex-harmonicoffsets(3)+winlen-1,:,:,:,:,:);
sweepIR4 = IR(IRstartindex-harmonicoffsets(4):IRstartindex-harmonicoffsets(4)+winlen-1,:,:,:,:,:);
sweepIR5 = IR(IRstartindex-harmonicoffsets(5):IRstartindex-harmonicoffsets(5)+winlen-1,:,:,:,:,:);
sweepIR6 = IR(IRstartindex-harmonicoffsets(6):IRstartindex-harmonicoffsets(6)+winlen-1,:,:,:,:,:);
if numel(sweepIR1)<1e6
    sweepspect1 = abs(fft(sweepIR1.*repmat(wf,[1,chans,1,cycles,dim5,dim6]),fftlen)).^2;
    sweepspect2 = abs(fft(sweepIR2.*repmat(wf,[1,chans,1,cycles,dim5,dim6]),fftlen)).^2;
    sweepspect3 = abs(fft(sweepIR3.*repmat(wf,[1,chans,1,cycles,dim5,dim6]),fftlen)).^2;
    sweepspect4 = abs(fft(sweepIR4.*repmat(wf,[1,chans,1,cycles,dim5,dim6]),fftlen)).^2;
    sweepspect5 = abs(fft(sweepIR5.*repmat(wf,[1,chans,1,cycles,dim5,dim6]),fftlen)).^2;
    sweepspect6 = abs(fft(sweepIR6.*repmat(wf,[1,chans,1,cycles,dim5,dim6]),fftlen)).^2;
else
    % nested for loops to reduce the chance of memory blow-out for
    % multidimensional measurements
    [sweepspect1,sweepspect2,sweepspect3,sweepspect4,sweepspect5] = ...
        deal(zeros(size(sweepIR1)));
    for ch = 1:chans
        for d4 = 1:cycles
            for d5 = 1:dim5
                for d6 = 1:dim6
                    sweepspect1 = abs(fft(sweepIR1(:,ch,1,d4,d5,d6).*wf,fftlen)).^2;
                    sweepspect2 = abs(fft(sweepIR2(:,ch,1,d4,d5,d6).*wf,fftlen)).^2;
                    sweepspect3 = abs(fft(sweepIR3(:,ch,1,d4,d5,d6).*wf,fftlen)).^2;
                    sweepspect4 = abs(fft(sweepIR4(:,ch,1,d4,d5,d6).*wf,fftlen)).^2;
                    sweepspect5 = abs(fft(sweepIR5(:,ch,1,d4,d5,d6).*wf,fftlen)).^2;
                    sweepspect6 = abs(fft(sweepIR6(:,ch,1,d4,d5,d6).*wf,fftlen)).^2;
                end
            end
        end
    end
end



t = (0:winlen-1)'./fs;
f = fs*(0:round(fftlen/2))'./fftlen;

if octsmooth > 0
    for ch = 1:chans
        for d4 = 1:cycles
            for d5 = 1:dim5
                for d6 = 1:dim6
                    sweepspect1(:,ch,1,d4,d5,d6) = octavesmoothing(sweepspect1(:,ch,1,d4,d5,d6),octsmooth, fs);
                    sweepspect2(:,ch,1,d4,d5,d6) = octavesmoothing(sweepspect2(:,ch,1,d4,d5,d6),octsmooth, fs);
                    sweepspect3(:,ch,1,d4,d5,d6) = octavesmoothing(sweepspect3(:,ch,1,d4,d5,d6),octsmooth, fs);
                    sweepspect4(:,ch,1,d4,d5,d6) = octavesmoothing(sweepspect4(:,ch,1,d4,d5,d6),octsmooth, fs);
                    sweepspect5(:,ch,1,d4,d5,d6) = octavesmoothing(sweepspect5(:,ch,1,d4,d5,d6),octsmooth, fs);
                    sweepspect6(:,ch,1,d4,d5,d6) = octavesmoothing(sweepspect6(:,ch,1,d4,d5,d6),octsmooth, fs);
                end
            end
        end
    end
    lowlimit = 128/(winlen/fs); % avoid very low freq hump error.
    % we zero (rather than delete) the LF components because this makes the
    % interpolation of the frequency scale for harmonic distortion simpler
    % to implement.
    sweepspect1(f<lowlimit,:,:,:,:,:)=0;
    sweepspect2(f<lowlimit,:,:,:,:,:)=0;
    sweepspect3(f<lowlimit,:,:,:,:,:)=0;
    sweepspect4(f<lowlimit,:,:,:,:,:)=0;
    sweepspect5(f<lowlimit,:,:,:,:,:)=0;
    sweepspect6(f<lowlimit,:,:,:,:,:)=0;
    
    sweepspect1(sweepspect1<0)=0;
    sweepspect2(sweepspect2<0)=0;
    sweepspect3(sweepspect3<0)=0;
    sweepspect4(sweepspect4<0)=0;
    sweepspect5(sweepspect5<0)=0;
    sweepspect6(sweepspect6<0)=0;
end

if isinf(relgain(1))
    relgainindices = 2:length(relgain);
    includenoise = true;
    if subtractnoise
        sweepspect1(:,:,1,relgainindices,:,:) = sweepspect1(:,:,1,relgainindices,:,:) - repmat(sweepspect1(:,:,1,1,:,:),[1,1,1,length(relgainindices),1,1]);
        sweepspect1(sweepspect1<0)=0;
        sweepspect2(:,:,1,relgainindices,:,:) = sweepspect2(:,:,1,relgainindices,:,:) - repmat(sweepspect2(:,:,1,1,:,:),[1,1,1,length(relgainindices),1,1]);
        sweepspect2(sweepspect2<0)=0;
        sweepspect3(:,:,1,relgainindices,:,:) = sweepspect3(:,:,1,relgainindices,:,:) - repmat(sweepspect3(:,:,1,1,:,:),[1,1,1,length(relgainindices),1,1]);
        sweepspect3(sweepspect3<0)=0;
        sweepspect4(:,:,1,relgainindices,:,:) = sweepspect4(:,:,1,relgainindices,:,:) - repmat(sweepspect4(:,:,1,1,:,:),[1,1,1,length(relgainindices),1,1]);
        sweepspect4(sweepspect4<0)=0;
        sweepspect5(:,:,1,relgainindices,:,:) = sweepspect5(:,:,1,relgainindices,:,:) - repmat(sweepspect5(:,:,1,1,:,:),[1,1,1,length(relgainindices),1,1]);
        sweepspect5(sweepspect5<0)=0;
        sweepspect6(:,:,1,relgainindices,:,:) = sweepspect6(:,:,1,relgainindices,:,:) - repmat(sweepspect6(:,:,1,1,:,:),[1,1,1,length(relgainindices),1,1]);
        sweepspect6(sweepspect6<0)=0;
    end
else
    relgainindices = 1:length(relgain);
    includenoise = false;
end



linecolor = HSVplotcolours2(cycles, 1, 1); % H,S,V

% Impulse response figure(s)
for d6 = 1:dim6
    for d5 = 1:dim5
        for ch = 1:chans
            compplot=figure('Name',[char(chanID(ch)), ' ', char(dim5ID(d5)), ' ' char(dim6ID(d6))]);
            subplot(3,2,1)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                plot(t,10*log10(sweepIR1(:,ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            plot(t,10*log10(wf.^2),...
                'color',[0.5 0.5 0.5],...
                'DisplayName','wf');
            title('Linear IR')
            ylabel('Level [dB]')
            grid on
            
            subplot(3,2,2)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                plot(t,10*log10(sweepIR2(:,ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            plot(t,10*log10(wf.^2),...
                'color',[0.5 0.5 0.5],...
                'DisplayName','wf');
            title('2nd harmonic IR')
            grid on
            
            subplot(3,2,3)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                plot(t,10*log10(sweepIR3(:,ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            plot(t,10*log10(wf.^2),...
                'color',[0.5 0.5 0.5],...
                'DisplayName','wf');
            title('3rd harmonic IR')
            ylabel('Level [dB]')
            grid on
            
            subplot(3,2,4)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                plot(t,10*log10(sweepIR4(:,ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            plot(t,10*log10(wf.^2),...
                'color',[0.5 0.5 0.5],...
                'DisplayName','wf');
            title('4th harmonic IR')
            grid on
            
            subplot(3,2,5)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                plot(t,10*log10(sweepIR5(:,ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            plot(t,10*log10(wf.^2),...
                'color',[0.5 0.5 0.5],...
                'DisplayName','wf');
            title('5th harmonic IR')
            ylabel('Level [dB]')
            xlabel('Time [s]')
            grid on
            
            subplot(3,2,6)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                plot(t,10*log10(sweepIR6(:,ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            plot(t,10*log10(wf.^2),...
                'color',[0.5 0.5 0.5],...
                'DisplayName','wf');
            title('6th harmonic IR')
            xlabel('Time [s]')
            grid on
            
            iplots = get(compplot,'Children');
            xlims = cell2mat(get(iplots,'Xlim'));
            set(iplots,'Xlim',[min(xlims(:,1)) max(xlims(:,2))])
            ylims = cell2mat(get(iplots,'Ylim'));
            set(iplots,'Ylim',[min(ylims(:,1)) max(ylims(:,2))])
            uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
                'Position', [0 0 65 30],...
                'Callback', 'setaxeslimits');
        end
    end
end

% Spectrum figure(s)
for d6 = 1:dim6
    for d5 = 1:dim5
        for ch = 1:chans
            compplot=figure('Name',[char(chanID(ch)), ' ', char(dim5ID(d5)), ' ' char(dim6ID(d6))]);
            subplot(3,2,1)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                semilogx(f,10*log10(sweepspect1(1:length(f),ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            title('Linear IR')
            ylabel('Level [dB]')
            grid on
            
            subplot(3,2,2)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                semilogx(f,10*log10(sweepspect2(1:length(f),ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            title('2nd harmonic IR')
            grid on
            
            subplot(3,2,3)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                semilogx(f,10*log10(sweepspect3(1:length(f),ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            ylabel('Level [dB]')
            grid on
            
            subplot(3,2,4)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                semilogx(f,10*log10(sweepspect4(1:length(f),ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            title('4th harmonic IR')
            grid on
            
            subplot(3,2,5)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                semilogx(f,10*log10(sweepspect5(1:length(f),ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            title('5th harmonic IR')
            ylabel('Level [dB]')
            xlabel('Frequency [Hz]')
            grid on
            
            subplot(3,2,6)
            for cycleindex = cycles:-1:1
                labelstring = [num2str(relgain(cycleindex)),' dB '];
                colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                semilogx(f,10*log10(sweepspect6(1:length(f),ch,1,cycleindex,d5,d6).^2),...
                    'color',colr,...
                    'DisplayName',labelstring);
                hold on
            end
            title('6th harmonic IR')
            xlabel('Frequency [Hz]')
            grid on
            
            iplots = get(compplot,'Children');
            set(iplots,'Xlim',[20 20000])
        ylims = cell2mat(get(iplots,'Ylim'));
        set(iplots,'Ylim',[min(ylims(:,1)) max(ylims(:,2))])
        uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
            'Position', [0 0 65 30],...
            'Callback', 'setaxeslimits');
        end
    end
end





% Harmonic distortion figure - one fig per non-silent cycle
for d6 = 1:dim6
    for d4 = relgainindices
        compplot = figure('Name',['Harmonic distortion relative to linear response, ' num2str(relgain(d4)) ' dB']);
        numberofsubplots = chans * dim5;
        maxsubplots = 100;
        if numberofsubplots > maxsubplots, numberofsubplots = maxsubplots; end
        [r, c] = subplotpositions(numberofsubplots,0.7);
        plotnum = 1;
        for d5 = 1:dim5
            for ch = 1:chans
                if plotnum <= maxsubplots
                    subplot(r,c,plotnum)
                    h = sweepspect2(1:length(f),ch,1,d4,d5,d6);
                    if includenoise
                        h(h<=sweepspect2(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h0 = interp(sweepspect1((1:length(f)),ch,1,d4,d5,d6),2);
                    h = 10*log10(h ./ h0(1:length(f)));
                    hindex = isreal(h) | ~isnan(h);
                    h(isinf(h))=nan;
                    semilogx(f(hindex)./2,real(h(hindex)),...
                        'color',[1,0,0],...
                        'DisplayName','2nd harmonic');
                    title([char(chanID(ch)) ' ' char(dim5ID(d5))])
                    hold on
                    
                    h = sweepspect3(1:length(f),ch,1,d4,d5,d6);
                    if includenoise
                    h(h<=sweepspect3(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h0 = interp(sweepspect1((1:length(f)),ch,1,d4,d5,d6),3);
                    h = 10*log10(h ./ h0(1:length(f)));
                    hindex = isreal(h) | ~isnan(h);
                    h(isinf(h))=nan;
                    semilogx(f(hindex)./3,real(h(hindex)),...
                        'color',[0.8,0.8,0],...
                        'DisplayName','3rd harmonic');
                    
                    h = sweepspect4(1:length(f),ch,1,d4,d5,d6);
                    if includenoise
                    h(h<=sweepspect4(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h0 = interp(sweepspect1((1:length(f)),ch,1,d4,d5,d6),4);
                    h = 10*log10(h ./ h0(1:length(f)));
                    hindex = isreal(h) | ~isnan(h);
                    h(isinf(h))=nan;
                    semilogx(f(hindex)./4,real(h(hindex)),...
                        'color',[0,0.8,0],...
                        'DisplayName','4th harmonic');
                    
                    h = sweepspect5(1:length(f),ch,1,d4,d5,d6);
                    if includenoise
                    h(h<=sweepspect5(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h0 = interp(sweepspect1((1:length(f)),ch,1,d4,d5,d6),5);
                    h = 10*log10(h ./ h0(1:length(f)));
                    hindex = isreal(h) | ~isnan(h);
                    h(isinf(h))=nan;
                    semilogx(f(hindex)./5,real(h(hindex)),...
                        'color',[0,0,1],...
                        'DisplayName','5th harmonic');
                    
                     h = sweepspect6(1:length(f),ch,1,d4,d5,d6);
                    if includenoise
                    h(h<=sweepspect6(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h0 = interp(sweepspect1((1:length(f)),ch,1,d4,d5,d6),6);
                    h = 10*log10(h ./ h0(1:length(f)));
                    hindex = isreal(h) | ~isnan(h);
                    h(isinf(h))=nan;
                    semilogx(f(hindex)./6,real(h(hindex)),...
                        'color',[0.5,0,0.5],...
                        'DisplayName','6th harmonic');
                    
                    semilogx([20 20000], [0 0],...
                        'color',[0.6 0.6 0.6]);
                    xlim([20 10000])
                    ylim([-120 20])
                    xlabel('Excitation Freq [Hz]')
                    ylabel('Relative Level [dB]')
                    grid on
                end
                plotnum = plotnum+1;
            end
        end
        iplots = get(compplot,'Children');
        if length(iplots) > 1
            %xlims = cell2mat(get(iplots,'Xlim'));
            %set(iplots,'Xlim',[min(xlims(:,1)) max(xlims(:,2))])
            set(iplots,'Xlim',[20 20000])
            ylims = cell2mat(get(iplots,'Ylim'));
            set(iplots,'Ylim',[min(ylims(:,1)) max(ylims(:,2))])
            uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
                'Position', [0 0 65 30],...
                'Callback', 'setaxeslimits');
        end
    end
end



% THD figure
for d6 = 1:dim6
    compplot = figure('Name',['THD, ' IN.name char(dim6ID(d6))]);
    numberofsubplots = chans * dim5;
    maxsubplots = 100;
    if numberofsubplots > maxsubplots, numberofsubplots = maxsubplots; end
    [r, c] = subplotpositions(numberofsubplots,0.7);
    plotnum = 1;
    for d5 = 1:dim5
        for ch = 1:chans
            if plotnum <= maxsubplots
                subplot(r,c,plotnum)
                for cycleindex = flip(relgainindices)
                    h2 = sweepspect2(1:length(f),ch,1,cycleindex,d5,d6);
                    if includenoise
                        h2(h2<=sweepspect2(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h2 = decimate(h2,2);
                    h3 = sweepspect3(1:length(f),ch,1,cycleindex,d5,d6);
                    if includenoise
                        h3(h3<=sweepspect3(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h3 = decimate(h3,3);
                    h4 = sweepspect4(1:length(f),ch,1,cycleindex,d5,d6);
                    if includenoise
                        h4(h4<=sweepspect4(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h4 = decimate(h4,4);
                    h5 = sweepspect5(1:length(f),ch,1,cycleindex,d5,d6);
                    if includenoise
                        h5(h5<=sweepspect5(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h5 = decimate(h5,5);
                    h6 = sweepspect6(1:length(f),ch,1,cycleindex,d5,d6);
                    if includenoise
                        h6(h6<=sweepspect5(1:length(f),ch,1,1,d5,d6)*db2pow(thresholddB)) = 0;
                    end
                    h6 = decimate(h6,6);
                    h = h2;
                    h(1:length(h3)) = h(1:length(h3))+h3;
                    h(1:length(h4)) = h(1:length(h4))+h4;
                    h(1:length(h5)) = h(1:length(h5))+h5;
                    h(1:length(h6)) = h(1:length(h6))+h6;
                    h0 = sweepspect1((1:length(f)),ch,1,cycleindex,d5,d6);
                    h = 10*log10(h ./ h0(1:length(h)));
                    hindex = isreal(h) | ~isnan(h);
                    h(isinf(h))=nan;
                    colr = permute(linecolor(cycleindex,1,1,:),[1,4,2,3]);
                    semilogx(f(hindex),real(h(hindex)),...
                        'color',colr,...
                        'DisplayName',[num2str(relgain(cycleindex)) ' dB']);
                    title([char(chanID(ch)) ' ' char(dim5ID(d5))])
                    hold on
                end
                xlim([20 20000])
                ylim([-120 20])
                xlabel('Excitation Freq [Hz]')
                ylabel('Relative Level [dB]')
                grid on
            end
            plotnum = plotnum+1;
        end
    end
    iplots = get(compplot,'Children');
    if length(iplots) > 1
        %xlims = cell2mat(get(iplots,'Xlim'));
        %set(iplots,'Xlim',[min(xlims(:,1)) max(xlims(:,2))])
        set(iplots,'Xlim',[20 20000])
        ylims = cell2mat(get(iplots,'Ylim'));
        set(iplots,'Ylim',[min(ylims(:,1)) max(ylims(:,2))])
        uicontrol('Style', 'pushbutton', 'String', 'Axes limits',...
            'Position', [0 0 65 30],...
            'Callback', 'setaxeslimits');
    end
    
end





OUT.funcallback.name = 'THD_via_ESS_v2.m';
OUT.funcallback.inarg = {octsmooth,thresholddB,winlen,winpow1,winpow2,subtractnoise};

%**************************************************************************
% Copyright (c) 2017, Densil Cabrera
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
%  * Neither the name of the The University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
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
%**************************************************************************