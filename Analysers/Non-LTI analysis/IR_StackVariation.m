function OUT = IR_StackVariation(IN, autocropthresh, calval, windur,winhop,winfun,scaling, domain, fs, cal)
% This function is designed to analyse an impulse response stack - either
% in dimension 4 or dimension 2. The purpose of such an analysis would
% probably be to assess time variance in a system.
%
% You can create an IR stack in the following steps:
% * Using AARAE's Generator GUI, choose measurement signal such as a swept
%   sinusoid, MLS or even an impulse (impuse1).
% * Before generating the signal, choose the number of cycles (>1), and
%   include a silent cycle if you wish to also assess background noise.
% * Generate the test signal
% * Play and record the test signal through the system
% * For signals such as swept sinusoids and the impulse1, create and IR
%   stack by using the '*' (convolve audio with audio2) button. If the
%   audio is single channel, you can choose to stack in dimension 2 or in
%   dimension 4; if the audio is multichannel you can only stack in
%   dimension 4.
% * For MLS signals and some others, you will need to use the appropriate
%   processor to derive an impulse response. Not all processors support the
%   generation of IR stacks.
%
% This function outputs a results leaf with statistics for each time sample
% (or spectrum component), as well as figures with results over time
% periods.
%
% Code by Densil Cabrera
% Version 1 beta (17 September 2015)



if nargin ==1 

    param = inputdlg({'Autocrop start threshold [-ve dB, or 0 to omit autocrop]';... 
                      'Reference calibration level [dB]';...
                      'Window duration for time-variance level evolution (s)';...
                      'Hop size between windows (s)';...
                      'Window function for time-variance level evolution: Rectangular (0), Hann (1)';...
                      'Result leaf stats for: Amplitude [1],Power [2],Level [3]';...
                      'Result leaf stats for: Time domain wave [0] Time domain envelope [1] or Frequency domain [2]'},...
                      'IR stack variation analysis',... 
                      [1 90],... 
                      {'-20';'0';'0.02';'0.01';'0';'1';'0'}); 

    param = str2num(char(param)); 

    if length(param) < 7, param = []; end 
    if ~isempty(param) 
        autocropthresh = param(1);
        calval = param(2);
        windur = param(3);
        winhop = param(4);
        winfun = param(5);
        scaling = param(6);
        domain = param(7);
    else
        OUT=[];
        return
    end
else
    param = [];
end


% *************************************************************************
if isstruct(IN) 
    IN = choose_from_higher_dimensions(IN,4,1); 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    if isfield(IN,'cal') % Get the calibration offset if it exists
        cal = IN.cal;
    else
        cal = 0;
    end
    
    % chanID is a cell array of strings describing each channel
    if isfield(IN,'chanID') % Get the channel ID if it exists
        chanID = IN.chanID;
    elseif size(audio,4) > 1
        chanID = makechanID(size(audio,2),0);
    else
        chanID = makechanID(1,0);
    end
    % in this analyser, dim2 might be being used for the IRstack, rather
    % than for channels - in which case it probably won't have a chanID
    
    % bandID is a vector, usually listing the centre frequencies of the
    % bands
    if isfield(IN,'bandID') % Get the band ID if it exists
        bandID = IN.bandID;
    else
        % asssign ordinal band numbers if bandID does not exist
        bandID = 1:size(audio,3);
    end
    

    % *********************************************************************
    
    
elseif ~isempty(param) || nargin > 1
                      
    audio = IN;
   
    fs = input_3;
    cal = input_4;
end
% *************************************************************************






if ~isempty(audio) && ~isempty(fs) && ~isempty(cal)

    [len,chans,bands,dim4] = size(audio);
    if dim4 > 1
        % assume IR stack is in dimension 4
        stackdim = 4;
        if isstruct(IN)
            if isfield(IN.properties,'relgain')
                if isinf(IN.properties.relgain(1))
                    % remove silent cycle if it exists
                    audio = audio(:,:,:,2:end);
                    IN.properties.relgain = IN.properties.relgain(2:end);
                    dim4 = dim4-1;
                end
            end
        end
        ncycles = dim4;
    elseif chans > 1
        % assume IR stack is in dimension 2
        stackdim = 2;
        ncycles = chans;
    else
        disp('IR stack not found in dimensions 2 or 4 - unable to analyse')
        OUT = [];
        return
    end
    
    if exist('cal','var')
        % typically use calval of 0 dB, or 94 dB for Pa if appropriately calibrated
        audio = cal_reset_aarae(audio,calval,cal);
    end
    
    if autocropthresh ~=0
        audio = autocropstart_aarae(audio,autocropthresh,2);
        len = size(audio,1);
    end
    
    % power vs synchronous averaging of the audio
    audiosynchMS = mean(audio,stackdim).^2;
    audiopowermean = mean(audio.^2,stackdim);
    
    % overall TV of the whole audio
    TV = permute(10*log10(mean(audiopowermean)./mean(audiosynchMS)),[2,3,4,1])...
        ./(10*log10(ncycles));
    % chan in dim1, bands in dim 2, cycles (must be singleton) in dim 3
    
    % variation in audio gain (coefficient of variation, in dB)
    if stackdim == 4
        gainvar = permute(10*log10(mean(mean(audio.^2,1) ./ repmat(mean(audiopowermean,1),[1,1,1,dim4]),4)),[2,3,4,1]);
    else
        gainvar = permute(10*log10(mean(mean(audio.^2,1) ./ repmat(mean(audiopowermean,1),[1,chans,1,1]),2)),[2,3,4,1]);
    end
    
    WindowLength = round(windur*fs);
    Offset = round(winhop*fs); % hopefully no need for rounding!
    nwin = round((len-WindowLength)/Offset); % number of windows
    if nwin < 1
        warndlg('Audio signal is shorter than the window length - unable to process with IR_StackVariation')
        OUT = [];
        return
    end
    [TVtime, gainvartime, Lwindow] = deal(zeros(nwin,size(TV,1),size(TV,2)));
    if winfun == 1
        w = repmat(hann(WindowLength),size(TV,1),size(TV,2));
    end
    for n = 1:nwin
        start = round((n-1)*Offset + 1);
        finish = start + WindowLength - 1;
        if winfun == 1
            TVtime(n,:,:) = 10*log10(mean(w.^2.*audiopowermean(start:finish,:,:))./...
                mean(w.^2.*audiosynchMS(start:finish,:,:)))./(10*log10(ncycles));
%             TVtime(n,:,:) = ncycles * (mean(w.^2.*audiosynchMS(start:finish,:,:))./ ...
%                  mean(w.^2.*audiopowermean(start:finish,:,:)) - (ncycles-1)/ncycles);
            Lwindow(n,:,:) = 10*log10(mean(w.^2.*audiopowermean(start:finish,:,:)));
            if stackdim == 4
                gainvartime(n,:,:) = 10*log10(mean(mean(repmat(w,[1,1,1,dim4])...
                    .*audio(start:finish,:,:,:).^2,1)./...
                    repmat(mean(w.^2.*audiopowermean(start:finish,:,:),1),[1,1,1,dim4]),4));
            else
                gainvartime(n,:,:) = 10*log10(mean(mean(repmat(w,[1,chans,1])...
                    .*audio(start:finish,:,:,:).^2,1)./...
                    repmat(mean(w.^2.*audiopowermean(start:finish,:,:),1),[1,1,1,dim4]),2));
            end
        else
            TVtime(n,:,:) = 10*log10(mean(audiopowermean(start:finish,:,:))./...
                mean(audiosynchMS(start:finish,:,:)))./(10*log10(ncycles));
%               TVtime(n,:,:) = ncycles * (mean(audiosynchMS(start:finish,:,:))./ ...
%                  mean(audiopowermean(start:finish,:,:)) - (ncycles-1)/ncycles);
            Lwindow(n,:,:) = 10*log10(mean(audiopowermean(start:finish,:,:)));
            if stackdim == 4
                gainvartime(n,:,:) = 10*log10(mean(mean(...
                    audio(start:finish,:,:,:).^2,1)./...
                    repmat(mean(audiopowermean(start:finish,:,:),1),[1,1,1,dim4]),4));
            else
                gainvartime(n,:,:) = 10*log10(mean(mean(...
                    audio(start:finish,:,:,:).^2,1)./...
                    repmat(mean(audiopowermean(start:finish,:,:),1),[1,1,1,dim4]),2));
            end
        end
    end
    
    % calculate coherence (using mean spectrum as reference) - THIS
    % IS USELESS
%     Gy = fft(audio);
%     if stackdim == 4
%         Gx = repmat(mean(Gy,4),[1,1,1,ncycles]);
%     else
%         Gx = repmat(mean(Gy,2),[1,ncycles,1]);
%     end
%     coherence = abs(conj(Gx).*Gy).^2 ./ ...
%         ((Gx.*conj(Gx)).*(Gy.*conj(Gy)));
%     clear Gx Gy
    
    figure('Name','IR stack variation')
    t = winhop*((1:nwin)-1)'; 
    f = fs*((1:len)'-1)./len;
    M = HSVplotcolours2(size(TVtime,3), size(TVtime,2), 1);
    subplot(3,1,1:2)
    for ch = 1:size(TVtime,2)
        for b = 1:size(TVtime,3)
            plot(t,TVtime(:,ch,b),'Color',permute(M(b,ch,:),[1,3,2]),...
                'DisplayName',[chanID{ch} ', ' num2str(bandID(b)) ' Hz'])
            hold on
        end
    end
    title(['Time Variance Index using ' num2str(ncycles) ' cycles'])
    xlabel('Time (s)')
    ylabel('TVI')
    

    
    subplot(3,1,3)
    for ch = 1:size(Lwindow,2)
        for b = 1:size(Lwindow,3)
            plot(t,Lwindow(:,ch,b),'Color',permute(M(b,ch,:),[1,3,2]),...
                'DisplayName',[chanID{ch} ', ' num2str(bandID(b)) ' Hz'])
            hold on
        end
    end
    %title('Level in each window')
    xlabel('Time (s)')
    ylabel('Level (dB)')
    
    
%         subplot(6,1,4:5)
%     for ch = 1:size(coherence,2)
%         for b = 1:size(coherence,3)
%             plot(f(1:round(len/2)),coherence(1:round(len/2),ch,b),'Color',permute(M(b,ch,:),[1,3,2]),...
%                 'DisplayName',[chanID{ch} ', ' num2str(bandID(b)) ' Hz'])
%             hold on
%         end
%     end
%     %title('Variation in level')
%     xlabel('Frequency (Hz)')
%     ylabel('Coherence')
%     
%         subplot(6,1,6)
%         spectrum = mean(abs(fft(audio)).^2,4);
%     for ch = 1:size(gainvartime,2)
%         for b = 1:size(gainvartime,3)
%             plot(f(1:round(len/2)),spectrum(1:round(len/2),:,:),'Color',permute(M(b,ch,:),[1,3,2]),...
%                 'DisplayName',[chanID{ch} ', ' num2str(bandID(b)) ' Hz'])
%             hold on
%         end
%     end
%     xlabel('Frequency (Hz)')
%     ylabel('Level (dB)')
    
    switch domain 
        case 2 % frequency domain
            audio = abs(fft(audio));
            audio = audio(1:ceil(end/2),:,:,:);
            xval = fs * ((1:length(audio))-1)./len; % frequencies
            xstring = 'Frequency';
            xunit = 'Hz';
        case 1 % time domain, absolute values
            audio = abs(audio);
            xval = ((1:len)-1)./fs; % times
            xstring = 'Time';
            xunit = 's';
        otherwise
            xval = ((1:len)-1)./fs; % times
            xstring = 'Time';
            xunit = 's';
    end
    
    if scaling ~=1
        if scaling < 3
            audio = audio.^scaling;
        else
            audio = 10*log10(abs(audio).^2);
        end
    end
    
    % main calculations
    % written for quick editing of the values and their order (rather than
    % for optimised speed)
    Datastack = [];
    ValNames = {};
    
    val = mean(audio,stackdim);
    ValNames = [ValNames, 'mean'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end
    
    val = std(audio,[],stackdim);
    ValNames = [ValNames, 'std'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end

    
    val = skewness(audio,[],stackdim);
    ValNames = [ValNames, 'skewness'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end

    
    val = kurtosis(audio,[],stackdim);
    ValNames = [ValNames, 'kurtosis'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end    
    
    val = max(audio,[],stackdim);
    ValNames = [ValNames, 'max'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end
    
    val = min(audio,[],stackdim);
    ValNames = [ValNames, 'min'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end
          
    val = std(audio,[],stackdim)./mean(audio,stackdim);
    ValNames = [ValNames, 'coefvar'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end
    
    val = rms(audio,stackdim);
    ValNames = [ValNames, 'rms'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end    
    
    
    val = std(audio,[],stackdim)./rms(audio,stackdim);
    ValNames = [ValNames, 'stdonrms'];
    if stackdim == 4
        Datastack = cat(4,Datastack,val);
    else
        Datastack = cat(2,Datastack,val);
    end
    
    if stackdim == 4
        doresultleaf(Datastack, 'Value', {xstring},...
            xstring,        xval,       xunit,         true,...
            'channels',     chanID,     'categorical', [],...
            'bands',        num2cell(bandID), 'Hz',    false,...
            'statistic',    ValNames,      'categorical',        [],...
            'name','IR_Stack_Stats');
    else
        doresultleaf(Datastack, 'Value', {xstring},...
            xstring,        xval,       xunit,         true,...
            'statistic',    ValNames,      'categorical',        [],...
            'bands',        num2cell(bandID), 'Hz',    false,...
            'name','IR_Stack_Stats');
    end
    
   

    if isstruct(IN)
        fig1 = figure('Name','Time variance');
        table1 = uitable('Data',TV,...
                'ColumnName',cellstr(num2str((bandID(:)))),...
                'RowName',chanID);
        OUT.tables = disptables(fig1,table1);
       
        OUT.funcallback.name = 'IR_StackVariation.m';
        OUT.funcallback.inarg = {autocropthresh, calval, windur,winhop,winfun,scaling, domain, fs, cal}; 
       
    else
       
        OUT = Datastack;
    end

else
    
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2014-15, Densil Cabrera
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
%  * Neither the name of the University of Sydney nor the names of its contributors
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