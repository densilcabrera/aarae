function OUT = ModeReverbTime(IN,fs,TUstart,TUend,dononlinear,orderout,cutoff,maxRTestimate,autotrunc)
% This function is intended for detailed analysis of the reverberation
% times in the low frequency modal part of a room's response (i.e., below
% the Schroeder frequency).
%
% The input audio should be a room impulse response with a high
% signal-to-noise ratio. The impulse response can have multiple input
% channels (e.g. microphones, in dimension 2) and multiple output channels
% (e.g. loudspeakers, in dimension 5). An asynchronous-output
% multiple-input impulse response recorded by AARAE will be in this format.
% Also, AARAE's calculator rectroommodesynthesiser creates pseudo-IRs in
% this format.
%
% First, peaks are identified (below the specified cutoff frequency). A
% narrow-band maximum phase filter is applied at each peak, and the
% reverberation time is calculated from this.
%
% STILL UNDER DEVELOPMENT!!!

if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
    if isfield(IN,'chanID')
        [chandata,chformat] = readchanID(IN.chanID);
    else
        chformat = 0;
        chandata = (1:size(audio,2))';
    end
    if isfield(IN,'properties')
        if isfield(IN.properties,'dim5ID')
            [outchandata,outchformat] = readchanID(IN.properties.dim5ID);
        else
            outchformat = 0;
            outchandata = (1:size(audio,5))';
        end
    else
        outchformat = 0;
        outchandata = (1:size(audio,5))';
    end
else
    audio = IN;
    chformat = 0;
    chandata = (1:size(audio,2))';
    outchformat = 0;
    outchandata = (1:size(audio,5))';
end

[len,chans,bands,dim4,outchans,dim6] = size(audio);

% we should only have 1, 2 and 5 as nonsingleton dimensions
if bands > 1 || dim4 > 1 || dim6 > 1
    warndlg('Analysis abandonded because this function requires input audio dimensions 3 (bands) and 4 (cycles) and 6 to be singleton')
    OUT = [];
    return
end
% even length
if rem(len,2)==1
    audio = [audio; zeros(1,chans,1,1,outchans)];
    len = len+1;
end

% user settings
if ~exist('cutoff','var'), cutoff = 200; end
if isempty(cutoff), cutoff = 200; end
if ~exist('zeropaddur','var'), zeropaddur = len/fs; end
if isempty(zeropaddur), zeropaddur = len/fs; end
if ~exist('TUstart','var'), TUstart = -5; end
if isempty(TUstart), TUstart = -5; end
if ~exist('TUend','var'), TUend = -15; end
if isempty(TUend), TUend = -15; end

if nargin == 1
    prompt = {'High cutoff frequency of the analysis (Hz)', ...
        'User-specified T evaluation range start level (dB)',...
        'User-specified T evaluation range end level (dB)',...
        'Do non-linear fitting [0 | 1] (NOT WORKING YET)',...
        'Filter order',...
        'If you wish to do exponential growth function pre-processing, what is the maximum estimated reverberation time (s)? Otherwise use 0',...
        'Lundeby automatic end truncation? [0 | 1]'};
    dlg_title = 'Mode Reverberation Time Settings';
    num_lines = 1;
    def = {num2str(cutoff),num2str(TUstart),num2str(TUend),'0','24','0','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        cutoff = str2double(answer{1,1});
        TUstart = abs(str2double(answer{2,1}));
        TUend = abs(str2double(answer{3,1}));
        dononlinear = abs(str2double(answer{4,1}));
        orderout = str2double(answer{5,1});
        maxRTestimate = str2double(answer{6,1});
        autotrunc = str2double(answer{7,1});
    end
end

cutoffmag = db2mag(-3);
if maxRTestimate == 0
    adjustmentexponent = 0;
else
    adjustmentexponent = -3*log(10)/maxRTestimate;
end

TUstart = -abs(TUstart); TUend = -abs(TUend);
if TUend >= TUstart
    TUstart = -5;
    TUend = -15;
end
% startthreshold = -20;
% % find ensemble start threshold
% startindex = find();



% optionally multiply the IR by an exponential growth function, to remove
% BT limit (if there is one!)
if adjustmentexponent ~= 0
    t = (0:(len-1))'./fs;
    audio = audio .* exp(adjustmentexponent*repmat(t,[1,chans,1,1,outchans]));
end


zeropadlen = round(zeropaddur*fs);
if rem(len+zeropadlen,2)==1, zeropadlen = zeropadlen+1; end % even length fft
%fftlen = len+zeropadlen;

% Identify the spectrum peaks in all channels up to the cutoff frequency

spectrapadded = fft([zeros(zeropadlen,chans,1,1,outchans); audio]);
spectra = fft(audio); % spectra for filter design
f_all = fs*((0:round(len/2))'-1)./len;
f = f_all(f_all<=cutoff);
nf = length(f);
foundpeaksval = zeros(nf,chans,outchans);
for ch = 1:chans
    for och = 1:outchans
        [val,ind] = findpeaks(abs(spectra(1:nf,ch,1,1,och)));
        foundpeaksval(ind,ch,och) = val;
    end
end
foundpeaksdB = 10*log10(foundpeaksval.^2./max(max(max(foundpeaksval.^2))));
% (this could be used to get rid of peaks X dB below max)
foundpeaks = foundpeaksval;
foundpeaks(foundpeaks>0) = 1;

% identify the channel combinations that have the highest amplitude for
% each peak
% [chmaxval, chmaxind] = max(foundpeaksval,[],2);
% [ochmaxval, ochmaxind] = max(foundpeaksval,[],5);


% consider using a sparse matrix instead of the following pre-allocation
fpeaks = sum(sum(foundpeaks,2),3); % number of peaks at each frequency
numberofdecays = sum(fpeaks);
chanID = cell(numberofdecays,1);
[fc, sourceind,recind,Level] = deal(zeros(numberofdecays,1));
if chformat == 0
    receiver = zeros(numberofdecays,1);
else
    receiver = zeros(numberofdecays,3);
end
if outchformat == 0
    source = zeros(numberofdecays,1);
else
    source = zeros(numberofdecays,3);
end
bandfiltered = zeros(2*len, numberofdecays);
filtercoefficients = zeros(len, numberofdecays);
% create filterbank and reorganise the IR for filtering
filtercount = 1;
columncount = 1;
for findex = 2:nf-1
    if fpeaks(findex) > 0
        mag = zeros(len,1);
        mag(findex) = 1; % 0 dB
        if findex > 2
            mag(2:findex-1) = ...
                (f_all(2:findex-1)./ f_all(findex-1) ).^(orderout) .* cutoffmag;
        end
        mag(findex+1:len/2+1) = ...
            (f_all(findex+1:len/2+1) ./ f_all(findex+1)).^(-orderout) .* cutoffmag;
        mag(len/2+2:end) = flipud(mag(2:len/2));
        % use zero phase filters for pre-filtering
        %mag = conj(minphasefreqdomain(mag,120,1)); % maximum phase filter coefficients
        selectedspectra = zeros(size(spectrapadded,1),fpeaks(findex));
        selectedspectracount = 0;
        for ch = 1:chans
            for och = 1:outchans
                if foundpeaks(findex,ch,och)==1
                    selectedspectracount = selectedspectracount+1;
                    selectedspectra(:,selectedspectracount) = ...
                        spectrapadded(:,ch,1,1,och);
                    chanID{columncount,1} = ...
                        [num2str(f(findex)), ' Hz. Source: ',  ...
                        num2str(chandata(ch)), '. Receiver: ', ...
                        num2str(outchandata(och)), '.'];
                    fc(columncount) = f(findex);
                    if chformat == 0
                        receiver(columncount) = chandata(ch);
                    else
                        receiver(columncount,:) = chandata(ch,:);
                    end
                    if outchformat == 0
                        source(columncount) = outchandata(och);
                    else
                        source(columncount,:) = outchandata(och,:);
                    end
                    Level(columncount) = foundpeaksdB(findex,ch,och);
                    sourceind(columncount) = och;
                    recind(columncount) = ch;
                    columncount = columncount+1;
                end
            end
        end
        bandfiltered(:,filtercount:filtercount+fpeaks(findex)-1) = ...
            selectedspectra;
        filtercoefficients(:,filtercount:filtercount+fpeaks(findex)-1) = ...
            repmat(mag,[1,fpeaks(findex)]);
        filtercount = filtercount + fpeaks(findex);
    end
end

filtercoefficients = conj(minphasefreqdomain(filtercoefficients,120,1)); % maximum phase filter coefficients
prefiltered = ifft(bandfiltered.*fft([ifft(filtercoefficients); zeros(zeropadlen,numberofdecays)]));


% calculate reverberation time from filtered IRs

% truncation:
% * remove maxphase pre-ring
prefiltered = prefiltered(1:len,:);

% undo exponential decay adjustment (now that filtering has been done)
if adjustmentexponent ~= 0
    t = (0:(size(prefiltered,1)-1))'./fs;
    prefiltered = prefiltered .* exp(-adjustmentexponent*repmat(t,[1,numberofdecays]));
end


% * Lundeby autotruncation
if autotrunc
    [crosspoint, Tlate, okcrosspoint] =...
        lundebycrosspoint(prefiltered.^2,fs,fc);
else
    crosspoint = len*ones(numberofdecays,1);
end
clear prefiltered

% final filtering - replace audio with autotruncated version
for n = 1:numberofdecays
    if okcrosspoint(n) && crosspoint(n)<len
        bandfiltered(:,n) = fft([zeros(zeropadlen,1);...
            audio(1:crosspoint(n),recind(n),1,1,sourceind(n));...
            zeros(len-crosspoint(n),1)]);
    end
end

bandfiltered = ifft(bandfiltered.*fft([ifft(filtercoefficients);zeros(len,numberofdecays)])); % apply filterbank & return to time domain

% undo exponential decay adjustment (now that filtering has been done)
if adjustmentexponent ~= 0
    t = (0:(size(bandfiltered,1)-1))'./fs;
    bandfiltered = bandfiltered .* exp(-adjustmentexponent*repmat(t,[1,numberofdecays]));
end

% remove zeropadding and square
%bandfiltered = bandfiltered(zeropadlen+1:end,:).^2;
bandfiltered = bandfiltered(1:len,:).^2;

% reverse-integration & convert to decibels
bandfiltered = 10*log10(flip(cumsum(flip(bandfiltered,1),1),1));

[o,p,q,r] = deal(zeros(2, numberofdecays));
s = deal(zeros(3, numberofdecays));
[EDT,T20,T30,EDTr2,T20r2,T30r2,Tuser,Tuserr2,EDTnonlin,EDTr2nonlin] =...
    deal(zeros(1, numberofdecays));
levdecaymodelEDT = zeros(size(bandfiltered));

for n = 1:numberofdecays
    bandfiltered(:,n) = bandfiltered(:,n) - bandfiltered(1,n); % start at 0 dB
    irstart = find(bandfiltered(:,n) <= 0, 1, 'first'); % 0 dB
    tstart = find(bandfiltered(:,n) <= -5, 1, 'first'); % -5 dB
    edtend = find(bandfiltered(:,n) <= -10, 1, 'first'); % -10 dB
    t20end = find(bandfiltered(:,n) <= -25, 1, 'first'); % -25 dB
    t30end = find(bandfiltered(:,n) <= -35, 1, 'first'); % -35 dB
    tuserstart = find(bandfiltered(:,n) <= TUstart, 1, 'first'); 
    tuserend = find(bandfiltered(:,n) <= TUend, 1, 'first'); 
    
    o(:,n) = polyfit((irstart:edtend)', ...
        bandfiltered(irstart:edtend,n),1)'; % EDT linear regression
    EDT(1,n) = 6*((o(2,n)-10)/o(1,n) ...
        -(o(2,n)-0)/o(1,n))/fs; % EDT value
    EDTr2(1,n) = corr(bandfiltered(irstart:edtend,n), ...
        (irstart:edtend)' * o(1,n) + ...
        o(2,n)).^2; % correlation coefficient, EDT
    
    p(:,n) = polyfit((tstart:t20end)', ...
        bandfiltered(tstart:t20end,n),1)'; % T20 regression
    T20(1,n) = 3*((p(2,n)-25)/p(1,n) ...
        -(p(2,n)-5)/ ...
        p(1,n))/fs; % reverberation time, T20
    T20r2(1,n) = corr(bandfiltered(tstart:t20end,n), ...
        (tstart:t20end)'*p(1,n) ...
        + p(2,n)).^2; % correlation coefficient, T20
    
    q(:,n) = polyfit((tstart:t30end)', ...
        bandfiltered(tstart:t30end,n),1)'; % T30 regression
    T30(1,n) = 2*((q(2,n)-35)/q(1,n) ...
        -(q(2,n)-5)/ ...
        q(1,n))/fs; % reverberation time, T30
    T30r2(1,n) = corr(bandfiltered(tstart:t30end,n), ...
        (tstart:t30end)'*q(1,n) ...
        + q(2,n)).^2; % correlation coefficient, T30
    
        r(:,n) = polyfit((tuserstart:tuserend)', ...
        bandfiltered(tuserstart:tuserend,n),1)'; % Tuser regression
    Tuser(1,n) = 60/(TUstart-TUend) *((r(2,n)+TUend)/r(1,n) ...
        -(r(2,n)+TUstart)/ ...
        r(1,n))/fs; % reverberation time, Tuser
    Tuserr2(1,n) = corr(bandfiltered(tuserstart:tuserend,n), ...
        (tuserstart:tuserend)'*r(1,n) ...
        + r(2,n)).^2; % correlation coefficient, Tuser
    
    if dononlinear == 1
        % Non-linear regression of level decay
        dataend1 = find(isinf(bandfiltered(irstart:end,n)),1,'first');
        if isempty(dataend1)
            dataend = length(bandfiltered(irstart:end,n));
        else
            dataend = dataend1-1;
        end
        
        y = bandfiltered(irstart:dataend,n);
        datalen = 1+dataend-irstart;
        x = ((irstart:dataend)'-irstart) ./ fs; % discrete time
        modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
            + abs(B(3)).*(datalen/fs-x)));
        % derive coefficients
        % nonlinear curve fitting seed coefficients
        b = [1*rand;-20*rand;1e-7*rand];
        maxloopcount = 10; % maximum number of fitting attempts for each parameter
        jthreshold = 0.9; % value less than 1, the higher the value, the stricter the threshold
        fitted = 0;
        loopcount = 0;
        while ~fitted
            loopcount = loopcount+1;
            [s(:,n),~,j] = nlinfit(x,y,modelfun,b);
            jflag = min(sum(j.^2));
            if jflag > jthreshold || loopcount == maxloopcount,
                fitted = 1;
                b = s(:,n)';
            else
                b = [100*rand;-100*rand;0.1*rand];
                disp('*')
            end
            
        end
        s(1,n) = abs(s(1,n));
        s(2,n) = -abs(s(2,n));
        s(3,n) = abs(s(3,n));
        %                 disp(['EDT nonlinear parameters: ',num2str(o(1,n)),...
        %                     ' ',num2str(o(2,n)),' ',num2str(o(3,n))])
        % create acoustic parameter-specific model using coefficients
        levdecaymodelEDT(1:crosspoint(n)-irstart+1,n) = ...
            10*log10(s(1,n).*exp(-abs(s(2,n)).*...
            ((0:crosspoint(n)-irstart)'./fs)));
        levdecaymodelEDT(:,n) = ...
            levdecaymodelEDT(:,n)-max(levdecaymodelEDT(:,n));
        
        startmodelEDT = ...
            find(levdecaymodelEDT(:,n) <= 0, 1, 'first'); % 0 dB
        endmodelEDT = ...
            find(levdecaymodelEDT(:,n) <= -10, 1, 'first'); % -10 dB
        
        % length of specific model
        edtmodlen = endmodelEDT-startmodelEDT+1;
        
        EDTnonlin(1,n) = 6.*(edtmodlen)./fs; % EDT in seconds
        
        EDTr2nonlin(1,n) = corr(y(irstart:edtend), ...
            (irstart:edtend)' * s(1,n) + ...
            s(2,n)).^2; % correlation coefficient, EDT
    end
end

% quality criteria:
% level (from initial spectrum analysis)
% distance between peaks, taking power into account
% R^2
% double slope decay detection
% modulation detection
% peak to noise floor ratio


% write output table
fig1 = figure('Name', 'Mode Reverberation Analysis');

if dononlinear == 1
    tabledata = [fc,source,receiver,Level,EDT',EDTr2',T20',T20r2',T30',T30r2',Tuser',Tuserr2',EDTnonlin',EDTr2nonlin'];
else
    tabledata = [fc,source,receiver,Level,EDT',EDTr2',T20',T20r2',T30',T30r2',Tuser',Tuserr2'];
end

if size(source,2)==1 && size(receiver,2)==1
    columnhead = {'Frequency (Hz)','Source','Receiver','Level (dB)',...
        'EDT (s)','EDT r^2','T20 (s)','T20 r^2','T30 (s)','T30 r^2',...
        'T user','T user r^2'};
elseif size(source,2)==1
    columnhead = {'Frequency (Hz)','Source',...
        'Receiver x','Receiver y','Receiver z','Level (dB)',...
        'EDT (s)','EDT r^2','T20 (s)','T20 r^2','T30 (s)','T30 r^2',...
        'T user','T user r^2'};
elseif size(receiver,2)==1
    columnhead = {'Frequency (Hz)','Source x','Source y','Source z',...
        'Receiver','Level (dB)',...
        'EDT (s)','EDT r^2','T20 (s)','T20 r^2','T30 (s)','T30 r^2',...
        'T user','T user r^2'};
else
    columnhead = {'Frequency (Hz)','Source x','Source y','Source z',...
        'Receiver x','Receiver y','Receiver z','Level (dB)',...
        'EDT (s)','EDT r^2','T20 (s)','T20 r^2','T30 (s)','T30 r^2',...
        'T user','T user r^2'};
end
if dononlinear == 1
    columnhead = [columnhead, {'Nonlinear','r^2 Nonlinear'}];
end
table1 = uitable('Data',tabledata,...
    'ColumnName',columnhead,...
    'RowName',{num2str((1:numberofdecays)')},...
    'ColumnWidth', {60});
[~,table] = disptables(fig1,table1);
OUT.tables = table;

% make plots


figure; plot(bandfiltered)

figure
plotspectra = 10*log10(abs(spectra(1:nf,:,:,:,:)).^2./max(max(max(foundpeaksval.^2))));

M = HSVplotcolours1(outchans, chans, [0.4 1]);
subplot(5,1,1)
plot(fc,EDT','LineStyle',':','Color',[0.5 0.5 0.5])
hold on
for ch = 1:chans
    for och = 1:outchans
        indices = (recind==ch) & (sourceind==och);
        plot(fc(indices),EDT(indices)','LineStyle','none','Marker','x',...
            'MarkerEdgeColor',permute(M(ch,och,:),[1,3,2]));
    end
end
ylabel('EDT (s)')
xlim([f(1) f(end)]);

subplot(5,1,2)
plot(fc,T20','LineStyle',':','Color',[0.5 0.5 0.5])
hold on
for ch = 1:chans
    for och = 1:outchans
        indices = (recind==ch) & (sourceind==och);
        plot(fc(indices),T20(indices)','LineStyle','none','Marker','x',...
            'MarkerEdgeColor',permute(M(ch,och,:),[1,3,2]));
    end
end
ylabel('T20 (s)')
xlim([f(1) f(end)]);


subplot(5,1,3)
plot(fc,T30','LineStyle',':','Color',[0.5 0.5 0.5])
hold on
for ch = 1:chans
    for och = 1:outchans
        indices = (recind==ch) & (sourceind==och);
        plot(fc(indices),T30(indices)','LineStyle','none','Marker','x',...
            'MarkerEdgeColor',permute(M(ch,och,:),[1,3,2]));
    end
end
ylabel('T30 (s)')
xlim([f(1) f(end)]);


subplot(5,1,4)
plot(fc,Tuser','LineStyle',':','Color',[0.5 0.5 0.5])
hold on
for ch = 1:chans
    for och = 1:outchans
        indices = (recind==ch) & (sourceind==och);
        plot(fc(indices),Tuser(indices)','LineStyle','none','Marker','x',...
            'MarkerEdgeColor',permute(M(ch,och,:),[1,3,2]));
    end
end
ylabel('Tuser (s)')
xlim([f(1) f(end)]);


subplot(5,1,5)
hold on
for ch = 1:chans
    for och = 1:outchans
        plot(f,plotspectra(:,ch,1,1,och),...
            'Color',permute(M(ch,och,:),[1,3,2]));
        indices = (recind==ch) & (sourceind==och);
        plot(fc(indices),Level(indices)','LineStyle','none','Marker','o',...
            'MarkerEdgeColor',permute(M(ch,och,:),[1,3,2]));
    end
end
ylabel('Level (dB)')
xlabel('Frequency (Hz)')
xlim([f(1) f(end)]);

OUT.funcallback.name = 'ModeReverbTime.m';
OUT.funcallback.inarg = {fs,TUstart,TUend,dononlinear,orderout,cutoff,maxRTestimate,autotrunc};

%**************************************************************************
% Copyright (c) 2015, Densil Cabrera
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