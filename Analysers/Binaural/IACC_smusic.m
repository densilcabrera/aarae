function OUT = IACC_smusic(IN)
% This function calculates IACC from running signals (2 channel) in octave
% and 1/3-octave bands, returning the mean and maximum IACC values.
%
% For further information, see:
%
% E. Manor, W.L. Martens & D. Cabrera, "Preferred spatial post-processing
% of popular stereophonic music for headphone reproduction," 133rd Audio
% Engineering Convention, October 2012, paper 8779.
%
% This function uses ERBFilterBank and MakeERBFilters by Malcolm Slaney, 
% which are in Dan Ellis' gammatonegram functions.


if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
end


if ~isempty(audio) && ~isempty(fs)
    
    
    X = audio;
    [nlen,chans,bands] = size(X);
    
    if chans ~= 2
        OUT = [];
        warndlg('IACC_smusic is only suitable for 2-channel audio input','AARAE info','modal')
        return
    end
    
    if bands > 1
        OUT = [];
        warndlg('IACC_smusic is only suitable for single band audio input','AARAE info','modal')
        return
    end
    
    % Octave and 1/3-octave spaced centre frequencies
    freqs = [31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000];
    freqsThird = [25,31.5,40,50,63,80,99,125,158,198,250,315,396,500,630,793,1000,1260,1587,2000,2520,3174,4000,5040,6349,8000,10080,12699,16000,20159];
    
    %%%%For 1/1 octave bandwidth
    
    % make and apply auditory filter
    % yL and yR are sets of intensity level in dB for each critical band
    fcoefs = MakeERBFilters(fs, freqs);
    XL = X(:,1);
    XR = X(:,2);
    yL = ERBFilterBank(XL, fcoefs);
    yL = yL';
    yR = ERBFilterBank(XR, fcoefs);
    yR = yR';
    
    % combine yL and yR
    yALL = zeros(nlen,20);
    a = 1;
    for i=1:10
        yALL(1:nlen,a) = yL(1:nlen,i);
        a = a+1;
        yALL(1:nlen,a) = yR(1:nlen,i);
        a = a+1;
    end
    
    % compute correlation coefficient
    overlaps = ceil(fs / 1000);  % +/- 1ms
    
    % create 3 matrices: IACC(max), IACC(min), IACC(med)
    
    IACC_max = zeros(1,10);
    IACC_mean = zeros(1,10);
    
    % compute all IACC values between L and R channels, for each critical band
    
    b = 1;
    for k=1:10
        yC(1:nlen,1) = yALL(1:nlen,b);
        b = b+1;
        yC(1:nlen,2) = yALL(1:nlen,b);
        b = b+1;
        [cc, lags] = xcorr(yC, overlaps, 'coeff');
        z = cc(:, 2);
        
        [zmx] = max(z);
        [zmn] = mean(z);
        
        IACC_max(k) = zmx;
        IACC_mean(k) = zmn;
        
    end
    
    IACC_max_val = mean(IACC_max);
    IACC_mean_val = mean(IACC_max);
    
    %%%%%%For 1/3 octave bandwidth
    
    % make and apply auditory filter
    % yL and yR are sets of intensity level in dB for each critical band
    fcoefs = MakeERBFilters(fs, freqsThird);
    XL = X(:,1);
    XR = X(:,2);
    yL = ERBFilterBank(XL, fcoefs);
    yL = yL';
    yR = ERBFilterBank(XR, fcoefs);
    yR = yR';
    
    % combine yL and yR
    yALL = zeros(nlen,60);
    a = 1;
    for i=1:30
        yALL(1:nlen,a) = yL(1:nlen,i);
        a = a+1;
        yALL(1:nlen,a) = yR(1:nlen,i);
        a = a+1;
    end
    
    % compute correlation coefficient
    overlaps = ceil(fs / 1000);  % +/- 1ms
    
    % create 3 matrices: IACC(max), IACC(mean), IACC(med)
    
    IACC_maxT = zeros(1,30);
    IACC_meanT = zeros(1,30);
    
    % compute all IACC values between L and R channels, for each critical band
    
    b = 1;
    for k=1:30
        yC(1:nlen,1) = yALL(1:nlen,b);
        b = b+1;
        yC(1:nlen,2) = yALL(1:nlen,b);
        b = b+1;
        [cc, lags] = xcorr(yC,overlaps, 'coeff');
        z = cc(:, 2);
        
        [zp] = max(z);
        [za] = mean(z);
        
        IACC_maxT(k) = zp;
        IACC_meanT(k) = za;
        
    end
    
    IACC_maxT_val = mean(IACC_maxT);
    IACC_meanT_val = mean(IACC_meanT);
    
    %%%%%Plot IACC amx and mean response for 1/1 octave bandwidth
    
    figure (1), clf;
    
    h = semilogx(freqs,IACC_max(1,1:10), 'b-');
    set(h, 'LineWidth', 2.);
    set(gca, 'FontSize', 18.);
    grid on;
    h = title('IACC response 1/1 octave');
    set(h, 'FontWeight', 'bold');
    set(h, 'FontSize', 16);
    set(gca, 'XScale', 'log');
    set(gca, 'LineWidth', 2.);
    h=ylabel('IACC');
    set(h, 'FontSize', 16.);
    h=xlabel('Frequency (Hz)');
    set(h, 'FontSize', 16.);
    set(gca, 'Ylim', [-1. 1.2]);
    
    hold on;
    
    h = semilogx(freqs,IACC_mean(1,1:10),'r-');
    set(h, 'LineWidth', 2.);
    
    hleg = legend('IACCmax','IACCmean');
    set(hleg,'Location','SouthWest');
    
    %print -dpng -r140 IACC1-1octave.png
    if isstruct(IN)
        IACC_oct = cat(2,IACC_max(1,1:10)',IACC_mean(1,1:10)');
        doresultleaf(IACC_oct,'Coefficient',{'Frequency'},...
                     'Frequency', num2cell(freqs),    'Hz',          true,...
                     'IACC',      {'Maximum','Mean'}, 'categorical', [],...
                     'name','IACC_oct');
    end
    %%%%%Plot IACC amx and mean response for 1/3 octave bandwidth
    
    figure (2), clf;
    
    h = semilogx(freqsThird,IACC_maxT(1,1:30), 'b-');
    set(h, 'LineWidth', 2.);
    set(gca, 'FontSize', 18.);
    grid on;
    h = title('IACC response 1/3 octave');
    set(h, 'FontWeight', 'bold');
    set(h, 'FontSize', 16);
    set(gca, 'XScale', 'log');
    set(gca, 'LineWidth', 2.);
    h=ylabel('IACC');
    set(h, 'FontSize', 16.);
    h=xlabel('Frequency (Hz)');
    set(h, 'FontSize', 16.);
    set(gca, 'Ylim', [-1. 1.2]);
    
    hold on;
    
    h = semilogx(freqsThird,IACC_meanT(1,1:30),'r-');
    set(h, 'LineWidth', 2.);
    
    hleg = legend('IACCmax','IACCmean');
    set(hleg,'Location','SouthWest');
    
    %print -dpng -r140 IACC1-3octave.png
    if isstruct(IN)
        IACC_thirdoct = cat(2,IACC_maxT(1,1:30)',IACC_meanT(1,1:30)');
        doresultleaf(IACC_thirdoct,'Coefficient',{'Frequency'},...
                     'Frequency', num2cell(freqsThird), 'Hz',          true,...
                     'IACC',      {'Maximum','Mean'},   'categorical', [],...
                     'name','IACC_thirdoct');
    end
    if isstruct(IN)
        
        OUT.IACC_max_val = IACC_max_val;
        OUT.IACC_mean_val = IACC_mean_val;
        OUT.IACC_maxT_val = IACC_maxT_val;
        OUT.IACC_meanT_val = IACC_meanT_val;
        OUT.funcallback.name = 'IACC_smusic.m';
        OUT.funcallback.inarg = {};
    end
    
else
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2013, Ella Manor
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