function [OUT] = Speech_Based_STI_via_ER_AARAE(IN,Received_Signal,wd,ovlp,fs)
% This function estimates the Speech Transmission Index (STI) using
% recordings of actual human speech. It uses the Envelope Regression Method
% (ER) to compare a speech signal recorded at the source location with a
% speech signal concurrently recorded at a receiver location.
% 
% To use this function select the reference audio recording (made at the
% source position). Next, select this function and press ?Analyse?. Choose
% your preferred analysis window length (it is a good idea for this to be
% at least twice as long as the reverb time), and the amount you wish the
% windows to overlap. Finally select the audio recording made at the
% received position.
% 
% This function is implemented based on the method described in Payton et
% al (2013).
% 
% Opsata et al (2014) found that when using short analysis windows (less
% than 40 s) the 90th percentile of the ER-STI results is a better
% indicator of the intelligibility than the mean. Once the analysis window
% exceeds 100 s the mean is shown to be a better indicator of speech
% intelligibility.

% Adam Opsata (2014)

if nargin ==1 
    
    param = inputdlg({'Analysis Window Size (s)';... % input box titles 
        'Window Overlap'},...% inputdlg window.
        'User Input Parameters',... % dialog window title.
        [1 30],... %  define the number of rows 
        ...        % 
        {'5';'0.75'}); % the preset answers for the dialog.
    
    param = str2num(char(param)); % turn strings into numbers.
    
    if length(param) < 2, param = []; end % check that the user
    % has input all the required
    % fields.
    if ~isempty(param) % assign the dialog's
        % inputs to your function's input parameters.
        wd = param(1); %window duration: window size in s
        ovlp = param(2); %Overlap
    end
    
    if isstruct(IN) % if the reference signal is a structure
    IN = choose_from_higher_dimensions(IN,1,1); 
    Reference_Signal = IN.audio; % Extract the audio data
    fs = IN.fs;
    else
        Reference_Signal=IN;
        fs_from_dlg = inputdlg({'Sampling Rate (Hz)'},...% inputdlg window.
        'User Input Parameters',... % This is the dialog window title.
        [1 30],...
        {'48000'});
        fs=fs_from_dlg(1);
    end
    
    
    if true % Use a menu & dialog box to select a received signal file within AARAE
        
        selection = choose_audio; % call AARAE's choose_audio function
        if ~isempty(selection)
            Received_Signal = selection.audio; % additional audio data
            fs2 = selection.fs; % sampling rate
            if ~(fs2 == fs)
                % match sampling rates if desired
                gcd_fs = gcd(fs,fs2); % greatest common denominator
                Received_Signal = resample(Received_Signal,fs/gcd_fs,fs2/gcd_fs);
                
            end
            %[len2, chans2, bands2] = size(Received_Signal); % new wave dimensions
        end
    end
    
end
    
    %these values add dimensions, which will be helpful when this function
    %is implemented to process multichannel signals:
    wdv=1;
    m=1;
    t=1;
    

    
    
    %%%%%%%% design the filter
    BandsPerOctave = 1;
    N = 6;           % Filter Order
    F0 = 1000;       % Center Frequency (Hz)
    f = fdesign.octave(BandsPerOctave,'Class 1','N,F0',N,F0,fs);
    %
    % Get all the valid center frequencies in the audio range to design the
    % filter bank:
    
    F0 = validfrequencies(f);
    F0= F0(3:9);
    
    Nfc = length(F0);
    for i=1:Nfc,
        f.F0 = F0(i);
        Hd(i) = design(f,'butter');
    end
    
    %Defining colors
    colour = [127, 0, 255; % violet
        0, 0, 0; ... % black
        255, 0, 0; ... % red
        255, 128, 0; ... % orange
        204, 204, 0; ... % dark yellow
        0, 204, 0; ... % mid green
        0, 204, 204; ... % dark cyan
        0, 0, 255]; ... % blue
        
    colour = colour / 255; % rescale to 0-1 range
    
    
    %time aligning the reference and recorded signals
    dlya=finddelay(Reference_Signal,Received_Signal); %finding the delay between the signals
    cpa=zeros((abs(dlya)),1); %a matrix of zeros to be used as delay compensation
    
    if dlya>0 % if the delay is greater than zero, add the delay to the ref signal
        
        Reference_Signal=[cpa;Reference_Signal];
        Received_Signal=[Received_Signal;cpa];
        
    end
    if dlya<0 %if the delay is less than zero, add the delay to the rec signal
        
        Reference_Signal=[Reference_Signal;cpa];
        Received_Signal=[cpa;Received_Signal];
    end
    lenRef=floor((length(Reference_Signal)));
    lenRec=floor((length(Received_Signal)));

    lenSml=min([lenRef,lenRec]);
    
    Reference_Signal=Reference_Signal(1:lenSml); %Making the signals the same length
    Received_Signal=Received_Signal(1:lenSml);
    
    audio=[Reference_Signal,Received_Signal];% the audio is now contained in a 4D matrix. three of them are in use, the fourth will be channels.
    
    one=length(audio); %signal
    two=2; % channels
    thr=7; %octave bands
    fou=1; %takes, or segments
    
    sze=[one,two,thr,fou]; %common size of the matrix
    
    ints=zeros(sze); %the envelope of the signals in each octave band
    for i=1:7;
        for j=1:2;
            ints(:,j,i) = (filter(Hd(i),audio(:,j))).^2;
        end
    end
    
    hamlen=ceil(fs*0.01); %convolving with a hamming window to smooth
    HAM=hamming(hamlen);
    con_length=hamlen+one; %convolution length
    smoooth=fft(HAM, con_length);
    
    nfs=408; % sample rate after downsampling
    dec=ceil(con_length*nfs/fs); % length of new down-sampled matrix
    
    rsmp= zeros(dec, two, thr); %the smoothed envelope of the bands, down-sampled.
    for b=1:two
        for c=1:thr
            rsmp(:,b,c)=resample((ifft((fft(ints(:,b,c), con_length).*smoooth)./con_length)),408,fs);
        end
    end
    
    wl=wd*nfs; % window length in samples
    
    lenfa=floor(length(rsmp(:,1,1,1)));
    lenca=floor(length(rsmp(:,2,1,1)));
    
    lens=[lenfa;lenca];
    
    lenm=min(lens); % the smallest signal length, could be user defined if a portion of the signal is desired
    wn=floor((lenm-wl)./(wl-(wl*ovlp))); %number of windows that can fit into the smallest signal
    SSTI=zeros(wn,1);
    
    weight=[ 0.085, 0.127, 0.23,  0.233, 0.309, 0.224, 0.173]; %weighting
    bta=   [ 0.085, 0.078, 0.065, 0.011, 0.047, 0.095]; % B weightings
    
    strt=1;
    
    %Calculating the modulation metric as described in Payton 2013
    for w=1:wn;
        stp=strt+wl;
        M=zeros(1,7); %the modulation metric
        
        for i=1:7
            x=rsmp(strt:stp,1,i,1);
            y=rsmp(strt:stp,2,i,1);
            
            N=length(x);
            ux=mean(x);
            uy=mean(y);
            
            norm=ux/uy;
            numer=(1/N).*(sum(x.*y))-(ux.*uy);
            denom=(1/N).*(sum(x.^2))-(ux.^2);
            
            M(:,i)=norm*numer/denom;
            % the following lines are very similar to taking the abs of the ssti result
            % if M(:,i)>1; % I added these 3 lines out of logic not instruction from a
            %     M(:,i)=.99; % paper. Without them the sSTI value is complex.
            % end
        end
        
        M(M>=1)=0.999999999;
        M(M<0)=0;
        %%%%% Clip
        aSNR=10.*log10(M./(1-M));
        for b=1:7 % this was also found to be necessary to avoid complex results
            if aSNR(b)<-15
                aSNR(b)=-15;
            end
            if aSNR(b)>15
                aSNR(b)=15;
            end
        end
        TI=(aSNR+15)./30;
        
        a_weighted=sum(weight.*TI);
        
        b_weights=zeros(1,6);
        for i=1:6
            b_weights(1,i)= bta(1,i).*sqrt(TI(1,i).*TI(1,(i+1)));
        end
        
        b_weighted=sum(b_weights);
        
        
        SSTI(w)=real(a_weighted-b_weighted); % real values to avoid complex
        %results (which are somewhat rare)
        
        strt=ceil(stp-(wl*ovlp));
        
    end
    
    sstia_mean(m,t)=mean(SSTI(SSTI>0.001));
    Mean_STI=sstia_mean(m,t);
    
    p=90;
    HUP(m,t,wdv) = prctile(SSTI,p); %the 90th percentile STI Value
    L10_STI=HUP(m,t,wdv);
    
    %%%%%%%% PLOTTING AND COMPARING %%%%%%%%%%%
    drtn=floor(lenm./nfs); %Duration
    vec=(linspace(0,drtn,wn)');
    set(figure, 'Position', [100, 100, 1400, 700]);
    plot(vec,SSTI,...
        'Color',colour(8,:),...
        'DisplayName',['Speech-Based STI']);
    xlabel('Time (s)','FontSize',16,...
        'FontWeight','normal','Color','k')
    ylabel('STI','FontSize',16,...
        'FontWeight','normal','Color','k')
    title(['Estimate of STI using the Envelope Regression Method, Window Size: ', num2str(wd),' s  '],...
        'FontSize',20,'FontWeight','normal','Color','k')
%     title(['Envelope Regression STI, Measurement ', num2str(m), ', Position ', num2str(t),', Window Size: ', num2str(wd),' s  '],...
%         'FontSize',20,'FontWeight','normal','Color','k')
    ylim([0 1]);
    xlim([1 drtn]);
    hold on
    plot(HUP(m,t,wdv),...
        'MarkerSize',8,'LineStyle','none',...
        'MarkerFaceColor',[0 1 1],'Marker','o','Color',[ 0 1 1],...
        'DisplayName',[strcat('90th percentile of ER results =  ', num2str(HUP(m,t,wdv)))]);
    legend('Show','Location','EastOutside')
    
    hold off
    
    OUT = [Mean_STI,L10_STI];
    
    
    %**************************************************************************
% Copyright (c) <2014>, <Adam Opsata>
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
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
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