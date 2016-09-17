function out = ReverberationTime_IR1(data,fs,startthresh,bpo,doplot,filterstrength,phasemode,noisecomp,autotrunc,f_low,f_hi,nonlin,customstart,customend)
% This function calculates reverberation time and related ISO 3382-1
% sound energy parameters (C50, D50, etc.) from a room impulse response.
%
% The impulse response can be multichannel, and also multiband (although it
% is preferable to do bandpass filtering within this function instead of
% prior to it).
%
% Code by Grant Cuthbert & Densil Cabrera
% Version 1.08 (December 2015)

%--------------------------------------------------------------------------
% INPUT VARIABLES
%--------------------------------------------------------------------------
%
% data = impulse response
%
% fs = Sampling frequency
%
% startthresh = Defines beginning of the IR as the first startthresh sample
%               as the maximum value in the IR (e.g if startthresh = -20,
%               the new IR start is the first sample that is >= 20 dB below
%               the maximum
%
% bpo = frequency scale to analyse (bands per octave)
%       (1 = octave bands (default); 3 = 1/3 octave bands)
%
% doplot = Output plots (1 = yes; 0 = no)
%
% filterstrength = a value that increases (or reduces) filter selectivity
% in the octave and 1/3-octave band filterbanks. A value of 2 doubles the
% effective filter order from the default.
%
% phasemode = -1 maximum, 0 zero, 1 minimum phase filters.
%
% noisecomp = noise compensation method:
% 0: none
% 1: Chu
% 2: Extrapolate late decay from crosspoint
%
% 3: Nonlinear fitting of entire reverse-integrated decay, following Xiang
% 4: Nonlinear fitting of entire reverse-integrated decay in dB, similar to
% Xiang (not recommended)
% 5: Nonlinear fitting of EDT, T20 and T30 evaluation ranges in dB
%
% SEE: N. Xiang  (1995) 'Evaluation of reverberation times using a
% nonlinear regression approach,' Journal of the Acoustical Society of
% America, 98(4), 2112-2121
% Nonlinear fitting assumes that background noise is steady state and the
% reverberation decay envelope is exponential.
%
% autotrunc = automatic truncation (0 = none, 1 = Lundeby)
%
% SEE: A. Lundeby, T.E. Vigran, H. Bietz and M. Vorlaender, "Uncertainties
% of Measurements in Room Acoustics", Acustica 81, 344-355 (1995)
%
%
% f_low, f_hi = approximate lowest and highest band frequencies to be
% filtered in the octave or 1/3 octave band filterbank.
%
%--------------------------------------------------------------------------
% OUTPUT VARIABLES
%--------------------------------------------------------------------------
%
% EDT = Decay curve time interval from 0 dB to -10 dB, in seconds,
%       multiplied by 6
% T20 = Decay curve time interval from -5 dB to -25 dB, in seconds,
%       multiplied by 3
% T30 = Decay curve time interval from -5 dB to -35 dB, in seconds,
%       multiplied by 2
%
% squared correlation coefficients:
% EDTr2 = EDT reverse-integrated decay curve to EDT linear regression line
% T20r2 = T20 reverse-integrated decay curve to T20 linear regression line
% T30r2 = T30 reverse-integrated decay curve to T30 linear regression line
%
% T20T30r = T20 to T30 ratio, as a percentage
% T30mid = T30 average of 500 Hz and 1000 Hz octave bands
%
% C50 = Early (<50 ms) to late energy ratio, in decibels
% C80 = Early (<80 ms) to late energy ratio, in decibels
% D50 = Early (<50 ms) to total energy ratio, as a percentage
% D80 = Early (<80 ms) to total energy ratio, as a percentage
% Ts = Time of the centre of gravity of the squared IR, in seconds
%
% BR = bass ratio
% TR = treble ratio
%
%--------------------------------------------------------------------------

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2013,2014, Grant Cuthbert and Densil Cabrera
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



if isstruct(data)
    % AARAE's check dimensions utility to make a partial selection of data
    % if it has more dimensions than this analyser can deal with. The
    % second argument is the maximum number of dimensions that the analyser
    % can handle, and the third argument determines whether a dialog box is
    % used for selection. If the number of dimensions in the input is ok
    % then no change is made.
    data = choose_from_higher_dimensions(data,3,1);
    ir = data.audio;
    fs = data.fs;
    if size(ir,3)>1
        bpo=0;
        f_low = 0;
        f_hi = 0;
    else
        if ~exist('bpo','var')
            bpo=1;
            f_low = 125;
            f_hi = 8000;
        end
    end
    if isfield(data,'chanID')
        chanID = data.chanID;
    else
        chanID = cellstr([repmat('Chan',size(data.audio,2),1) num2str((1:size(data.audio,2))')]);
    end
else
    ir = data;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end
if nargin < 10, f_low = 125; f_hi = 8000; end
if nargin < 9, autotrunc = 2; end
if nargin < 8, noisecomp = 1; end
if nargin < 7, phasemode = -1; end
if nargin < 6, filterstrength = 2; end
if nargin < 5, doplot = 1; end
if nargin < 4, bpo = 1; end
if nargin < 3
    startthresh = -20;
    % Dialog box for settings
    prompt = {'Threshold for IR start detection', ...
        'Bands per octave (0 | 1 | 3)', ...
        'Lowest centre frequency (Hz)', ...
        'Highest centre frequency (Hz)', ...
        'Filter strength (use 1 for 12th order 72 dB/oct; use 2 for 24th order 144 dB/oct)', ...
        'Zero phase (0), Maximum phase (-1) or Minimum phase (1) filters',...
        'Noise compensation: None (0), Subtract final 10% Chu (1), Extrapolate Late Decay (2)', ...
        'Automatic end truncation: None (0), Lundeby (1), Crosspoint from nonlinear fitting (2)',...
        'Also do non-linear curve fitting: No (0); fit the reverse integrated decay plus noise (1); fit in dB (2); or fit using specified exponent of power (e.g. 0.5 for amplitude)',...
        'Custom evaluation range start (dB)',...
        'Custom evaluation range end (dB)',...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = [1 60];
    def = {'-20',num2str(bpo),num2str(f_low),num2str(f_hi),...
        '2','-1','1','2','0','-5','-15','1'};
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
        nonlin = str2num(answer{9,1});
        customstart = str2num(answer{10,1});
        customend = str2num(answer{11,1});
        doplot = str2num(answer{12,1});
    else
        out = [];
        return
    end
end

if ~isempty(ir) && ~isempty(fs) && ~isempty(startthresh) && ~isempty(bpo) && ~isempty(doplot) && ~isempty(filterstrength) && ~isempty(phasemode) && ~isempty(noisecomp) && ~isempty(autotrunc) && ~isempty(f_low) && ~isempty(f_hi)
    %--------------------------------------------------------------------------
    % START TRUNCATION
    %--------------------------------------------------------------------------
    
    [len,chans,bands] = size(ir);
    if (bands>1)
        multibandIR = 1;
    else
        multibandIR = 0;
    end
    
    % Get last 10 % for Chu noise compensation if set, and to return the
    % Peak to noise ratio (this needs to be done before the circular shift
    % is done as part of startpoint detection/allignment)
    ir_end10 = ir(round(0.9*len):end,:,:);

    if noisecomp == 2;
        autotrunc = 1;
    end
    
    % Preallocate
    m = zeros(1,chans,bands); % maximum value of the IR
    startpoint = zeros(1,chans,bands); % the auto-detected start time of the IR
    
    for dim2 = 1:chans
        for dim3 = 1:bands
            m(1,dim2,dim3) = max(ir(:,dim2,dim3).^2); % maximum value of the IR
            startpoint(1,dim2,dim3) = find(ir(:,dim2,dim3).^2 >= m(1,dim2,dim3)./ ...
                (10^(abs(startthresh)/10)),1,'first'); % Define start point
            
            %for cases where startpoint was not found
            startpoint(isempty(startpoint))=1;
            
            %startpoint = min(startpoint,[],3);
            if startpoint(1,dim2,dim3) >1
                
                % zero the data before the startpoint
                ir(1:startpoint(1,dim2,dim3)-1,dim2,dim3) = 0;
                
                % rotate the zeros to the end (to keep a constant data length)
                ir(:,dim2,dim3) = circshift(ir(:,dim2,dim3),-(startpoint(1,dim2,dim3)-1));
                
            end % if startpoint
        end
    end % for dim2
    
    
    try
        early50 = ir(1:1+floor(fs*0.05),:,:); % Truncate Early50
        late50 = ir(ceil(fs*0.05):end,:,:); % Truncate Late50
        do50 = true;
    catch
        [early50, late50] = deal(NaN(1,chans,bands));
        do50 = false;
    end
    
    try
        early80 = ir(1:1+floor(fs*0.08),:,:); % Truncate Early80
        late80 = ir(ceil(fs*0.08):end,:,:); % Truncate Late80
        do80 = true;
    catch
        [early80, late80] = deal(NaN(1,chans,bands));
        do80 = false;
    end
    
    
    %--------------------------------------------------------------------------
    % FILTERING
    %--------------------------------------------------------------------------
    
    if (multibandIR == 0) && ((bpo ==1) || (bpo == 3))
        noctaves = log2(f_hi/f_low);
        if bpo == 1
            f_low = 1000*2.^round(log2((f_low/1000))); % make sure it is oct
            fc = f_low .* 2.^(0:round(noctaves)); % approx freqencies
        elseif bpo == 3
            fc = f_low .* 2.^(0:1/3:round(3*noctaves)/3);% approx freqencies
        end
        bandfc = exact2nom_oct(fc); % nominal frequencies
        bands = length(bandfc);
    else
        if isfield(data,'bandID');
            bandfc = data.bandID;
            if find(bandfc==1000) - find(bandfc==500) == 1
                bpo = 1;
            end
        else
            bandfc = 1:bands;
        end
        if length(bandfc) ~= bands
            bandfc = 1:bands;
        end
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (multibandIR == 0) && ((bpo ==1) || (bpo == 3))
        
        % use AARAE's fft-based filters
        
        if bpo == 1
            order = [12,12]*filterstrength;
            iroct = octbandfilter_viaFFT(ir,fs,bandfc,order,0,1000,0,phasemode);
            if do50
                early50oct = octbandfilter_viaFFT(early50,fs,bandfc,order,0,1000,0,phasemode);
                early80oct = octbandfilter_viaFFT(early80,fs,bandfc,order,0,1000,0,phasemode);
            end
            if do80
                late50oct = octbandfilter_viaFFT(late50,fs,bandfc,order,0,1000,0,phasemode);
                late80oct = octbandfilter_viaFFT(late80,fs,bandfc,order,0,1000,0,phasemode);
            end
            ir_end10oct = octbandfilter_viaFFT(ir_end10,fs,bandfc,order,0,1000,0,phasemode);
            
            
        else
            order = [36,24] * filterstrength;
            iroct = thirdoctbandfilter_viaFFT(ir,fs,bandfc,order,0,1000,0,phasemode);
            if do50
                early50oct = thirdoctbandfilter_viaFFT(early50,fs,bandfc,order,0,1000,0,phasemode);
                early80oct = thirdoctbandfilter_viaFFT(early80,fs,bandfc,order,0,1000,0,phasemode);
            end
            if do80
                late50oct = thirdoctbandfilter_viaFFT(late50,fs,bandfc,order,0,1000,0,phasemode);
                late80oct = thirdoctbandfilter_viaFFT(late80,fs,bandfc,order,0,1000,0,phasemode);
            end
            ir_end10oct = thirdoctbandfilter_viaFFT(ir_end10,fs,bandfc,order,0,1000,0,phasemode);

        end
        
        % check the number of bands again, in case some were not ok to filter
        bands = size(iroct,3);
        fc = fc(1:bands);
        bandfc = bandfc(1:bands);
        
        %----------------------------------------------------------------------
        % AUTO-TRUNCATION OF END
        %----------------------------------------------------------------------
        
        
        C = zeros(1,chans,bands); % decay energy that has been chopped off by truncation (set initially to zero)
        refilter = false;
        if autotrunc == 1 || autotrunc == 2
            switch autotrunc
                case 1
            % Lundeby crosspoint
            [crosspoint, Tlate, okcrosspoint] =...
                lundebycrosspoint(iroct(1:(end-max(max(max(startpoint)))),:,:).^2, fs,fc);
            C=exp(-crosspoint*6*log(10)./(Tlate*fs));
                case 2

                    if noisecomp == 1
                        noisemultiplier = 0.5; % 3 dB lower if noise is being subtracted
                    else
                        noisemultiplier = 1; 
                    end
                    [crosspoint,Tlate,C] = crosspointnonlinfit(iroct(1:(end-max(max(max(startpoint)))),:,:).^2, fs, fc,noisemultiplier);
                    %C=exp(-crosspoint*6*log(10)./(Tlate*fs));
                    %const = log(1e6)./Tlate;
                    %C = const.*exp(-const.*crosspoint/fs);
                    okcrosspoint=ones(1,chans,bands);
                case 3
                    % use initial preliminary crosspoint from Lundeby
                    [~, Tlate, okcrosspoint,crosspoint] =...
                lundebycrosspoint(iroct(1:(end-max(max(max(startpoint)))),:,:).^2, fs,fc);
            end
            % Consider redoing the filtering here to minimise ringing from
            % background noise contaminating the IR. This could be slow
            % because each band has a different crosspoint, and so needs to
            % be filtered individually.
            if bpo == 1 || bpo == 3
                refilter = true;
            else
                refilter = false;
            end
            if refilter
                iroct = repmat(ir,[1,1,bands]); % this overwrites the filtered audio!
            end
            softcrop = 1;
            if softcrop == 1
                % short half-Hann fade-out to reduce problems with
                % discontinuities when cropping is done
                winperiods = 4;
                win = zeros(ceil(fs/1000 + winperiods*fs/bandfc(1)),bands);
                for b = 1:bands
                    winx = hann(round(2*(fs/1000 + winperiods*fs/bandfc(b))));
                    winx = winx(round(end/2):end);
                    win(1:length(winx),b) = winx;
                end 
                winlen = length(win);
            end
            % autotruncation
            for ch = 1:chans
                for b = 1:bands
                    if crosspoint(1,ch,b) < len
                        if softcrop == 0
                            iroct(crosspoint(1,ch,b):end,ch,b) = 0;
                        % soft crop
                        else
                            if crosspoint(1,ch,b)+ winlen <= len
                                iroct(crosspoint(1,ch,b):crosspoint(1,ch,b)+winlen-1,ch,b) =...
                                    iroct(crosspoint(1,ch,b):crosspoint(1,ch,b)+winlen-1,ch,b).*win(:,b);
                                iroct(crosspoint(1,ch,b)+winlen-1:end,ch,b) = 0;
                            else
                                win1 = win(1:len-crosspoint(1,ch,b)+1,b);
                                iroct(crosspoint(1,ch,b):end,ch,b) = iroct(crosspoint(1,ch,b):end,ch,b).*win1;
                            end
                        end
                        if crosspoint(1,ch,b) > fs*0.05
                            if do50
                                late50oct(round(crosspoint(1,ch,b)-fs*0.05+1):end,ch,b)=0;
                            end
                        end
                        if crosspoint(1,ch,b) > fs*0.08
                            if do80
                                late80oct(round(crosspoint(1,ch,b)-fs*0.08+1):end,ch,b)=0;
                            end
                        end
                    end
                end
            end
        end
        
        
        if refilter
            % redo the filtering, now that the trash has been emptied
            % this is quite inefficient to do one band at a time, but
            % necessary if they all have different crosspoints
            for b = 1:bands
                if bpo == 1
                    iroct(:,:,b) = octbandfilter_viaFFT(iroct(:,:,b),fs,bandfc(b),order,0,1000,0,phasemode);
                else
                    iroct(:,:,b) = thirdoctbandfilter_viaFFT(iroct(:,:,b),fs,bandfc(b),order,0,1000,0,phasemode);
                end
            end
        end
        
        % Zero the data that should be zero from startpoint detection (but
        % isn't due to filter ringing).
        % This is especially important for non-linear curve fitting.
        % This is done regardless of the above end-truncation in case a
        % crosspoint was not found.
        for ch = 1:chans
            iroct(end-startpoint(1,ch,1):end,ch,:) = 0;
        end
        
        
        %----------------------------------------------------------------------
        % CALCULATE ENERGY PARAMETERS
        %----------------------------------------------------------------------
        alloct = sum(iroct.^2);
        if do50
            early50oct = sum(early50oct.^2);
            late50oct = sum(late50oct.^2);
            C50 = 10*log10(early50oct ./ late50oct); % C50
            D50 = (early50oct ./ alloct); % D50
        else
            [C50,D50] = deal(NaN(1,chans,bands));
        end
        
        if do80
            
            early80oct = sum(early80oct.^2);
            late80oct = sum(late80oct.^2);
            C80 = 10*log10(early80oct ./ late80oct); % C80
            D80 = (early80oct ./ alloct); % D80
        else
            [C80,D80] = deal(NaN(1,chans,bands));
        end
        % time values of IR in seconds
        tstimes = (0:(length(iroct)-1))' ./ fs;
        
        Ts = (sum(iroct.^2 .* ...
            repmat(tstimes,[1,chans,bands])))./alloct; % Ts
        %     end
        
        
        
        
    end % if multibandIR == 0
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    if (multibandIR == 1) || (bpo == 0)
        bands = size(ir,3);
        iroct = ir;
        
        %----------------------------------------------------------------------
        % END AUTO-TRUNCATION
        %----------------------------------------------------------------------
        
        if autotrunc == 1
            % Lundeby crosspoint
            [crosspoint, Tlate, okcrosspoint] =...
                lundebycrosspoint(iroct(1:(end-max(max(max(startpoint)))),:,:).^2, fs,bandfc);
        elseif autotrunc == 2
            crosspoint = crosspointnonlinfit(iroct(1:(end-max(max(max(startpoint)))),:,:).^2, fs, bandfc);
            okcrosspoint=ones(1,chans,bands);
        end
        if autotrunc == 1 || autotrunc == 2
            % autotruncation
            softcrop = 1;
            if softcrop == 1
                % short half-Hann fade-out to avoid problems with
                % discontinuities when cropping is done
                win = hann(round(fs/1000)); 
                win = win(round(end/2:end));
                winlen = length(win);
            end
            for ch = 1:chans
                for b = 1:bands
                    if crosspoint(1,ch,b) < len
                        % hard crop - can introduce discontinuities
                        if softcrop == 0
                            iroct(crosspoint(1,ch,b):end,ch,b) = 0;
                        % soft crop
                        else
                            if crosspoint(1,ch,b)+ winlen <= len
                                iroct(crosspoint(1,ch,b):crosspoint(1,ch,b)+winlen-1,ch,b) =...
                                    iroct(crosspoint(1,ch,b):crosspoint(1,ch,b)+winlen-1,ch,b).*win;
                                iroct(crosspoint(1,ch,b)+winlen-1:end,ch,b) = 0;
                            else
                                win1 = win(1:len-crosspoint);
                                iroct(crosspoint(1,ch,b):end,ch,b) = iroct(crosspoint(1,ch,b):end,ch,b).*win1;
                            end
                        end
                        if crosspoint(1,ch,b) > fs*0.05
                            late50(round(crosspoint(1,ch,b)-fs*0.05+1):end,ch,b)=0;
                        end
                        if crosspoint(1,ch,b) > fs*0.08
                            late80(round(crosspoint(1,ch,b)-fs*0.08+1):end,ch,b)=0;
                        end
                    end
                end
            end
        end
        
        
        %----------------------------------------------------------------------
        % CALCULATE ENERGY PARAMETERS
        %----------------------------------------------------------------------
        if multibandIR == 1
            disp('Ideally energy ratio parameters should use an unfiltered impulse response as input.')
            disp('Filtering should be done after the early and late parts of the IR have been determined.')
        end
        early50oct = sum(early50.^2);
        early80oct = sum(early80.^2);
        late50oct = sum(late50.^2);
        late80oct = sum(late80.^2);
        alloct = sum(ir.^2);
        
        
        C50 = 10*log10(early50oct ./ late50oct); % C50
        C80 = 10*log10(early80oct ./ late80oct); % C80
        D50 = (early50oct ./ alloct); % D50
        D80 = (early80oct ./ alloct); % D80
        
        % time values of IR in seconds
        tstimes = (0:(length(iroct)-1))' ./ fs;
        
        Ts = (sum(iroct.^2 .* ...
            repmat(tstimes,[1,chans,bands])))./alloct; % Ts
        
    end % if multibandIR == 1
    
    % mean square of last 10%
    if (multibandIR == 0) && ((bpo==1) || (bpo==3))
        ir_end10oct = mean(ir_end10oct.^2);
    else
        ir_end10oct = mean(ir_end10.^2);
    end
    % dynamic range, or peak-to-noise ratio (using last 10% as 'noise')
    PNR = 10*log10(max(iroct.^2)./ir_end10oct);
    if noisecomp ~= 1
        ir_end10oct = zeros(1,chans,bands);
    end
    
    %--------------------------------------------------------------------------
    % REVERBERATION TIME (Reverse Integration, Linear Regression)
    %--------------------------------------------------------------------------
    
    %***************************************
    % Derive the reverse-integrated decay curve(s)
    
    % square and correct for background noise (Chu method) if set
    % and correct for exponential decay truncation (C)
    noisepowerfactor = 1; % arguably this could be greater than 1 to compensate for in-phase addition of noise with IR
    iroct2 = iroct.^2 - noisepowerfactor*repmat(ir_end10oct,[len,1,1]);
    for ch = 1:chans
        for b = 1:bands
            iroct2(crosspoint(1,ch,b),:,:) = iroct2(crosspoint(1,ch,b),:,:)+C; % add C to final sample
        end
    end
    
    iroct2(iroct2<0) = 0;
    
    if noisecomp == 2
        % extrapolate decay beyond crosspoint
        tau = 0.8*Tlate./log(1e6); % 0.8 is more conservative than 1
        for ch = 1:chans
            for b = 1:bands
                if okcrosspoint(1,ch,b)
                    endmean = mean(iroct2(crosspoint(1,ch,b)-50:crosspoint(1,ch,b)-1,ch,b));
                    iroct2((crosspoint(1,ch,b)):end,ch,b) = ...
                        endmean .* exp(((((crosspoint(1,ch,b)):len)...
                        -(crosspoint(1,ch,b)+25))./fs)./-tau(1,ch,b));
                end
            end
        end
    end
    
    % Reverse integrate squared IR, and express the result in decibels
    levdecay = 10*log10(flip(cumsum(flip(iroct2,1)),1));
    
    for dim2 = 1:chans
        for dim3 = 1:bands
            % Adjust so that the IR starts at 0 dB
            levdecay(:,dim2,dim3) = levdecay(:,dim2,dim3) - levdecay(1,dim2,dim3);
        end
    end
    
    %if (noisecomp == 0) || (noisecomp == 1) || (noisecomp == 2)
        %***************************************
        % Derive decay times, via linear regressions over the appropriate
        % evaluation ranges
        
        % preallocate
        [o,p,q,x] = deal(zeros(2, chans, bands));
        [EDT,T20,T30,Tcustom,EDTr2,T20r2,T30r2,Tcustomr2] = deal(zeros(1, chans, bands));
        
        for dim2 = 1:chans
            for dim3 = 1:bands
                
                % find the indices for the relevant start and end samples
                irstart = find(levdecay(:,dim2,dim3) <= 0, 1, 'first'); % 0 dB
                tstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                edtend = find(levdecay(:,dim2,dim3) <= -10, 1, 'first'); % -10 dB
                t20end = find(levdecay(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
                t30end = find(levdecay(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
                customstartind = find(levdecay(:,dim2,dim3) <= -abs(customstart), 1, 'first');
                customendind = find(levdecay(:,dim2,dim3) <= -abs(customend), 1, 'first');
                
                %******************************************************************
                % linear regression for EDT
                
                o(:,dim2,dim3) = polyfit((irstart:edtend)', ...
                    levdecay(irstart:edtend,dim2,dim3),1)';
                EDT(1,dim2,dim3) = 6*((o(2,dim2,dim3)-10)/o(1,dim2,dim3) ...
                    -(o(2,dim2,dim3)-0)/o(1,dim2,dim3))/fs; % EDT
                EDTr2(1,dim2,dim3) = corr(levdecay(irstart:edtend,dim2,dim3), ...
                    (irstart:edtend)' * o(1,dim2,dim3) + ...
                    o(2,dim2,dim3)).^2; % correlation coefficient, EDT
                
                
                %******************************************************************
                % linear regression for T20
                
                p(:,dim2,dim3) = polyfit((tstart:t20end)', ...
                    levdecay(tstart:t20end,dim2,dim3),1)';
                T20(1,dim2, dim3) = 3*((p(2,dim2,dim3)-25)/p(1,dim2,dim3) ...
                    -(p(2,dim2,dim3)-5)/ ...
                    p(1,dim2,dim3))/fs; % reverberation time, T20
                T20r2(1,dim2, dim3) = corr(levdecay(tstart:t20end,dim2,dim3), ...
                    (tstart:t20end)'*p(1,dim2,dim3) ...
                    + p(2,dim2,dim3)).^2; % correlation coefficient, T20
                
                
                
                %******************************************************************
                % linear regression for T30
                
                q(:,dim2,dim3) = polyfit((tstart:t30end)', ...
                    levdecay(tstart:t30end,dim2,dim3),1)'; % linear regression
                T30(1,dim2, dim3) = 2*((q(2,dim2,dim3)-35)/q(1,dim2,dim3) ...
                    -(q(2,dim2,dim3)-5)/ ...
                    q(1,dim2,dim3))/fs; % reverberation time, T30
                T30r2(1,dim2, dim3) = corr(levdecay(tstart:t30end,dim2,dim3), ...
                    (tstart:t30end)'*q(1,dim2,dim3) ...
                    + q(2,dim2,dim3)).^2; % correlation coefficient, T30
                
                %******************************************************************
                % linear regression for custom (user-specified) evaluation range
                
                x(:,dim2,dim3) = polyfit((customstartind:customendind)', ...
                    levdecay(customstartind:customendind,dim2,dim3),1)'; % linear regression
                Tcustom(1,dim2, dim3) = (60/(abs(customend)-abs(customstart)))...
                    *((x(2,dim2,dim3)-abs(customend))/x(1,dim2,dim3) ...
                    -(x(2,dim2,dim3)-abs(customstart))/ ...
                    x(1,dim2,dim3))/fs; % reverberation time, custom
                Tcustomr2(1,dim2, dim3) = corr(levdecay(customstartind:customendind,dim2,dim3), ...
                    (customstartind:customendind)'*x(1,dim2,dim3) ...
                    + x(2,dim2,dim3)).^2; % correlation coefficient, custom
                
                
            end % dim3
        end % dim2
        
   % elseif (noisecomp == 3) || (noisecomp == 4) || (noisecomp == 5)
        % nonlinear curve fit following Xiang 1995
        
        if nonlin >0
        % nonlinear curve fitting seed coefficients
        b1 = 1*rand;
        b2 = -20*rand;
        b3 = 1e-7*rand;
        b = [b1;b2;b3];
        maxloopcount = 10; % maximum number of fitting attempts for each parameter
        jthreshold = 0.9; % value less than 1, the higher the value, the stricter the threshold
        
        % preallocate
        t = zeros(len, chans, bands);
%         [o,p,q,x] = deal(zeros(3, chans, bands));
%         levdecaymodelEDT = zeros(len, chans, bands);
%         levdecaymodelT20 = zeros(len, chans, bands);
%         levdecaymodelT30 = zeros(len, chans, bands);
%         edtmodlen = zeros(1, chans, bands);
%         t20modlen = zeros(1, chans, bands);
%         t30modlen = zeros(1, chans, bands);
%         EDT = zeros(1, chans, bands);
%         T20 = zeros(1, chans, bands);
%         T30 = zeros(1, chans, bands);
%         EDTr2 = zeros(1, chans, bands);
%         T20r2 = zeros(1, chans, bands);
%         T30r2 = zeros(1, chans, bands);
%         dataend = zeros(1, chans, bands);
        xx = zeros(3,chans,bands);
        levdecaymodelALL = zeros(len, chans, bands);
        [edtmodlen,Tnonlin,Tnonlinr2,dataend] = deal(zeros(1, chans, bands));
        
        
        for dim2 = 1:chans
            for dim3 = 1:bands
                % find the end of the sound (find first -inf)
                dataend1 = find(isinf(levdecay(:,dim2,dim3)),1,'first');
                if isempty(dataend1)
                    dataend(1,dim2,dim3) = length(levdecay);
                else
                    dataend(1,dim2,dim3) = dataend1-1;
                end
            end
        end
        
        
        for dim2 = 1:chans
            for dim3 = 1:bands
                irstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); 
                Tend = find(levdecay(:,dim2,dim3) <= -25, 1, 'first'); % just used for correlation coef
                t(1:dataend(:,dim2,dim3),dim2,dim3) = ((1:dataend(:,dim2,dim3))-1)./fs;
                
                irlen = (t(dataend(:,dim2,dim3),dim2,dim3)); % duration of ir
                
                
                % decay curve model
                if nonlin == 1
                        % model directly the reverse integrated power decay
                        y = 10.^(levdecay((irstart:dataend(1,dim2,dim3)),dim2,dim3)./10);
                        x = t((irstart:dataend(1,dim2,dim3)),dim2,dim3); % discrete time
                        modelfun = @(B,x)(abs(B(1)).*exp(-abs(B(2)).*x)...
                            + abs(B(3)).*(irlen-x));
                elseif nonlin == 2
                        % curve-fit in decibels
                        y = levdecay((irstart:dataend(1,dim2,dim3)),dim2,dim3);
                        x = t((irstart:dataend(1,dim2,dim3)),dim2,dim3); % discrete time
                        modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
                            + abs(B(3)).*(irlen-x)));
%                 elseif nonlin == 3
%                         % curve-fit exp (almost certainly of no use!)
%                         % this excessively weights the start
%                         y = exp(10.^(levdecay((irstart:dataend(1,dim2,dim3)),dim2,dim3)./10));
%                         x = t((irstart:dataend(1,dim2,dim3)),dim2,dim3); % discrete time
%                         modelfun = @(B,x)exp((abs(B(1)).*exp(-abs(B(2)).*x)...
%                             + abs(B(3)).*(irlen-x)));
                else
                        % curve-fit using exponent specified by nonlin
                        % a value of 0.5 gives pressure weighting, which
                        % places less weight on the start than power. A
                        % smaller value (e.g. 0.1) places even more
                        % emphasis on getting the tail modelled well.
                        y = 10.^(levdecay((irstart:dataend(1,dim2,dim3)),dim2,dim3)./10).^nonlin;
                        x = t((irstart:dataend(1,dim2,dim3)),dim2,dim3); % discrete time
                        modelfun = @(B,x)((abs(B(1)).*exp(-abs(B(2)).*x)...
                            + abs(B(3)).*(irlen-x))).^nonlin;
                end
%                 elseif noisecomp == 5
%                     % use edt length for method 5
%                     y = levdecay((irstart:edtend),dim2,dim3);
%                     x = t((irstart:edtend),dim2,dim3); % discrete time
%                     modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
%                         + abs(B(3)).*(irlen-x)));
%                 end
                
                % derive coefficients
                try
                fitted = 0;
                loopcount = 0;
                while ~fitted
                    loopcount = loopcount+1;
                    [xx(:,dim2,dim3),~,j] = nlinfit(x,y,modelfun,b);
                    jflag = min(sum(j.^2));
                    if jflag > jthreshold || loopcount == maxloopcount,
                        fitted = 1;
                        b = xx(:,dim2,dim3)';
                    else
                        b = [100*rand;-100*rand;0.1*rand];
                        %disp('*')
                    end
                    
                end
                xx(1,dim2,dim3) = abs(xx(1,dim2,dim3));
                xx(2,dim2,dim3) = -abs(xx(2,dim2,dim3));
                xx(3,dim2,dim3) = abs(xx(3,dim2,dim3));
%                 disp(['Nonlinear parameters: ',num2str(xx(1,dim2,dim3)),...
%                     ' ',num2str(xx(2,dim2,dim3)),' ',num2str(xx(3,dim2,dim3))])
                % create acoustic parameter-specific model using coefficients
                levdecaymodelALL(:,dim2,dim3) = ...
                    10*log10(xx(1,dim2,dim3).*exp(-abs(xx(2,dim2,dim3)).*t(:,dim2,dim3)));
                levdecaymodelALL(:,dim2,dim3) = ...
                    levdecaymodelALL(:,dim2,dim3)-max(levdecaymodelALL(:,dim2,dim3));
                
                startmodelEDT = ...
                    find(levdecaymodelALL(:,dim2,dim3) <= -5, 1, 'first'); 
                endmodelEDT = ...
                    find(levdecaymodelALL(:,dim2,dim3) <= -35, 1, 'first'); 
                
                % length of specific model
                edtmodlen(:,dim2,dim3) = length(t(startmodelEDT:endmodelEDT,dim2,dim3));
                
                Tnonlin(1,dim2,dim3) = 2.*(edtmodlen(1,dim2,dim3))./fs; 
                
                Tnonlinr2(1,dim2,dim3) = corr(levdecay(irstart:Tend,dim2,dim3), ...
                    (irstart:Tend)' * xx(1,dim2,dim3) + ...
                    xx(2,dim2,dim3)).^2; % correlation coefficient
                
                catch
                    % if nonlinear fitting returns an error
                     Tnonlin(1,dim2,dim3) = nan;
                     Tnonlinr2(1,dim2,dim3) = nan;
                end
            end % dim3
        end % dim2
        else
            [Tnonlin,Tnonlinr2] = deal(nan(1, chans, bands));
        end
%         % no need for further fitting etc for noisecomp 3 and 4
%         if (noisecomp == 3) || (noisecomp ==4)
%             [p,q] = deal(o);
%             %[T20,T30] = deal(EDT);
%             disp('Noise compensation methods 3 and 4 use a non-linear regression over the entire reverse-integrated decay.')
%             disp('Hence the evaluation ranges of EDT, T20 and T30 are not used in their full sense,')
%             disp('and their values are likely to be very similar to each other.')
%         end
%         
%         % T20 nonlinear model
%         for dim2 = 1:chans
%             for dim3 = 1:bands
%                 tstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
%                 t20end = find(levdecay(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
%                 if (isempty(t20end)) || (isempty(tstart))|| (t20end > dataend(:,dim2,dim3))
%                     T20(1,dim2,dim3) = NaN;
%                     T20r2(1,dim2,dim3) = NaN;
%                     p(:,dim2,dim3) = [0;0;0];
%                     %disp('insufficient dynamic range for T20')
%                 else
%                     irlen = (t(dataend(:,dim2,dim3),dim2,dim3)); % finite length of IR (ULI)
%                     
%                     % decay curve model
%                     
%                     if noisecomp == 5
%                         % vector of response (dependent variable) values
%                         x = t((tstart:t20end),dim2,dim3); % discrete time
%                         y = levdecay((tstart:t20end),dim2,dim3);
%                         modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
%                             + abs(B(3)).*(irlen-x)));
%                         
%                         % derive coefficients
%                         fitted = 0;
%                         loopcount = 0;
%                         while ~fitted
%                             loopcount = loopcount+1;
%                             [p(:,dim2,dim3),~,j] = nlinfit(x,y,modelfun,b);
%                             jflag = min(sum(j.^2));
%                             if jflag > jthreshold || loopcount == maxloopcount,
%                                 fitted = 1;
%                                 b = p(:,dim2,dim3)';
%                             else
%                                 b = [100*rand;100*rand;rand];
%                                 disp('*')
%                             end
%                             
%                         end
%                         p(1,dim2,dim3) = abs(p(1,dim2,dim3));
%                         p(2,dim2,dim3) = -abs(p(2,dim2,dim3));
%                         p(3,dim2,dim3) = abs(p(3,dim2,dim3));
%                         
%                         disp(['T20 nonlinear parameters: ',num2str(p(1,dim2,dim3)),...
%                             ' ',num2str(p(2,dim2,dim3)),' ',num2str(p(3,dim2,dim3))])
%                     end
%                     % create acoustic parameter-specific model using coefficients
%                     levdecaymodelT20(:,dim2,dim3) = ...
%                         10*log10(p(1,dim2,dim3).*exp(-abs(p(2,dim2,dim3)).*t(:,dim2,dim3)));
%                     levdecaymodelT20(:,dim2,dim3) = ...
%                         levdecaymodelT20(:,dim2,dim3)-max(levdecaymodelT20(:,dim2,dim3));
%                     
%                     startmodelT20 = ...
%                         find(levdecaymodelT20(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
%                     endmodelT20 = ...
%                         find(levdecaymodelT20(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
%                     t20modlen(:,dim2,dim3) = length(t(startmodelT20:endmodelT20,dim2,dim3));
%                     
%                     T20(1,dim2,dim3) = 3.*(t20modlen(1,dim2,dim3))./fs; % T20 in seconds
%                     
%                     T20r2(1,dim2,dim3) = corr(levdecay(tstart:t20end,dim2,dim3), ...
%                         (tstart:t20end)'*p(1,dim2,dim3) ...
%                         + p(2,dim2,dim3)).^2; % correlation coefficient, T20
%                 end
%             end % dim 3
%         end % dim 2
%         
%         
%         % T30 nonlinear model
%         for dim2 = 1:chans
%             for dim3 = 1:bands
%                 %tstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
%                 t30end = find(levdecay(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
%                 if (isempty(t30end)) || (isempty(tstart)) || (t30end > dataend(:,dim2,dim3))
%                     T30(1,dim2,dim3) = NaN;
%                     T30r2(1,dim2,dim3) = NaN;
%                     q(:,dim2,dim3) = [0;0;0];
%                     %disp('insufficient dynamic range for T30')
%                 else
%                     
%                     irlen = (t(dataend(:,dim2,dim3),dim2,dim3)); % finite length of IR (ULI)
%                     
%                     
%                     % decay curve model
%                     
%                     if noisecomp == 5
%                         % vector of response (dependent variable) values
%                         x = t((tstart:t30end),dim2,dim3); % discrete time
%                         y = levdecay((tstart:t30end),dim2,dim3);
%                         modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
%                             + abs(B(3)).*(irlen-x)));
%                         
%                         
%                         % derive coefficients
%                         fitted = 0;
%                         loopcount = 0;
%                         while ~fitted
%                             loopcount = loopcount+1;
%                             [q(:,dim2,dim3),~,j] = nlinfit(x,y,modelfun,b);
%                             jflag = min(sum(j.^2));
%                             if jflag > jthreshold || loopcount == maxloopcount,
%                                 fitted = 1;
%                                 b = q(:,dim2,dim3)';
%                             else
%                                 b = [100*rand;100*rand;rand];
%                                 disp('*')
%                             end
%                             
%                         end
%                         q(1,dim2,dim3) = abs(q(1,dim2,dim3));
%                         q(2,dim2,dim3) = -abs(q(2,dim2,dim3));
%                         q(3,dim2,dim3) = abs(q(3,dim2,dim3));
%                         disp(['T30 nonlinear parameters: ',num2str(q(1,dim2,dim3)),...
%                             ' ',num2str(q(2,dim2,dim3)),' ',num2str(q(3,dim2,dim3))])
%                         
%                     end
%                     % create acoustic parameter-specific model using coefficients
%                     
%                     levdecaymodelT30(:,dim2,dim3) = ...
%                         10*log10(q(1,dim2,dim3).*exp(-abs(q(2,dim2,dim3)).*t(:,dim2,dim3)));
%                     
%                     levdecaymodelT30(:,dim2,dim3) = ...
%                         levdecaymodelT30(:,dim2,dim3)-max(levdecaymodelT30(:,dim2,dim3));
%                     
%                     startmodelT30 = ...
%                         find(levdecaymodelT30(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
%                     endmodelT30 = ...
%                         find(levdecaymodelT30(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
%                     t30modlen(:,dim2,dim3) = length(t(startmodelT30:endmodelT30,dim2,dim3));
%                     
%                     T30(1,dim2,dim3) = 2.*(t30modlen(1,dim2,dim3))./fs; % T30 in seconds
%                     
%                     T30r2(1,dim2,dim3) = corr(levdecay(tstart:t30end,dim2,dim3), ...
%                         (tstart:t30end)'*q(1,dim2,dim3) ...
%                         + q(2,dim2,dim3)).^2; % correlation coefficient, T30
%                 end
%             end % dim 3
%         end % dim 2
   % end
    
    %--------------------------------------------------------------------------
    
    
    if  bpo == 1
        fc500 = find(bandfc == 500);
        fc1000 = find(bandfc == 1000);
        if ~isempty(fc500) && ~isempty(fc1000) && fc1000-fc500 == 1
            T30mid = mean([T30(1,:,fc500); T30(1,:,fc1000)]);
            T20mid = mean([T20(1,:,fc500); T20(1,:,fc1000)]);
            EDTmid = mean([EDT(1,:,fc500); EDT(1,:,fc1000)]);
        else
            [EDTmid,T20mid,T30mid] = deal(nan(1,chans));
        end
        
        fc125 = find(bandfc == 125);
        fc250 = find(bandfc == 250);
        if ~isempty(fc125) && ~isempty(fc250) && fc250-fc125 == 1
            T30low = mean([T30(1,:,fc125); T30(1,:,fc250)]);
            T20low = mean([T20(1,:,fc125); T20(1,:,fc250)]);
            EDTlow = mean([EDT(1,:,fc125); EDT(1,:,fc250)]);
        else
            [EDTlow,T20low,T30low] = deal(nan(1,chans));
        end
        
        fc2000 = find(bandfc == 2000);
        fc4000 = find(bandfc == 4000);
        if ~isempty(fc2000) && ~isempty(fc4000) && fc4000-fc2000 == 1
            T30high = mean([T30(1,:,fc2000); T30(1,:,fc4000)]);
            T20high = mean([T20(1,:,fc2000); T20(1,:,fc4000)]);
            EDThigh = mean([EDT(1,:,fc2000); EDT(1,:,fc4000)]);
        else
            [EDThigh,T20high,T30high] = deal(nan(1,chans));
        end
        
        % Bass ratio
        if ~isempty(T30mid) && ~isempty(T30low)
            BR_T30 = T30low ./ T30mid;
            BR_T20 = T20low ./ T20mid;
            BR_EDT = EDTlow ./ EDTmid;
        else
            [BR_EDT,BR_T20,BR_T30] = deal(nan(1,chans));
        end
        
        if ~isempty(T30mid) && ~isempty(T30high)
            TR_T30 = T30high ./ T30mid;
            TR_T20 = T20high ./ T20mid;
            TR_EDT = EDThigh ./ EDTmid;
        else
            [TR_EDT,TR_T20,TR_T30] = deal(nan(1,chans));
        end
        
    end
    
    % percentage ratio of T20 to T30
    T20T30r = (T20./T30)*100;
    
    %--------------------------------------------------------------------------
    % OUTPUT
    %--------------------------------------------------------------------------
    
    
    % Create output structure
    
    out.bandfc = bandfc;
    
    
    out.EDT = permute(EDT,[2,3,1]);
    out.T20 = permute(T20,[2,3,1]);
    out.T30 = permute(T30,[2,3,1]);
    out.Tcustom = permute(Tcustom,[2,3,1]);
    out.Tnonlin = permute(Tnonlin,[2,3,1]);
    out.C50 = permute(C50,[2,3,1]);
    out.C80 = permute(C80,[2,3,1]);
    out.D50 = permute(D50,[2,3,1]);
    out.D80 = permute(D80,[2,3,1]);
    out.Ts = permute(Ts,[2,3,1]);
    out.EDTr2 = permute(EDTr2,[2,3,1]);
    out.T20r2 = permute(T20r2,[2,3,1]);
    out.T30r2 = permute(T30r2,[2,3,1]);
    out.Tcustomr2 = permute(Tcustomr2,[2,3,1]);
    out.Tnonlinr2 = permute(Tnonlinr2,[2,3,1]);
    out.T20T30ratio = permute(T20T30r,[2,3,1]);
    out.PNR = permute(PNR,[2,3,1]);
    if exist('EDTmid','var')
        out.EDTmid = permute(EDTmid,[2,3,1]);
        out.T20mid = permute(T20mid,[2,3,1]);
        out.T30mid = permute(T30mid,[2,3,1]);
    end
    if exist('EDTlow','var')
        out.EDTlow = permute(EDTlow,[2,3,1]);
        out.T20low = permute(T20low,[2,3,1]);
        out.T30low = permute(T30low,[2,3,1]);
    end
    if exist('EDThigh','var')
        out.EDThigh = permute(EDThigh,[2,3,1]);
        out.T20high = permute(T20high,[2,3,1]);
        out.T30high = permute(T30high,[2,3,1]);
    end
    if exist('BR_EDT','var')
        out.BR_EDT = permute(BR_EDT,[2,3,1]);
        out.BR_T20 = permute(BR_T20,[2,3,1]);
        out.BR_T30 = permute(BR_T30,[2,3,1]);
    end
    if exist('TR_EDT','var')
        out.TR_EDT = permute(TR_EDT,[2,3,1]);
        out.TR_T20 = permute(TR_T20,[2,3,1]);
        out.TR_T30 = permute(TR_T30,[2,3,1]);
    end
    
    % if chans == 1
    %     disp(out)
    % else
    %disp(['bandfc:' ,num2str(out.bandfc), ' Hz'])
    %disp('Early Decay Time (s):')
    %disp(out.EDT)
    %disp('Reverberation Time T20 (s):')
    %disp(out.T20)
    %disp('Reverberation Time T30 (s):')
    %disp(out.T30)
    
    %disp('Clarity Index C50 (dB):')
    %disp(out.C50)
    %disp('Clarity Index C80 (dB):')
    %disp(out.C80)
    %disp('Definition D50:')
    %disp(out.D50)
    %disp('Definition D80:')
    %disp(out.D80)
    %disp('Centre Time (s):')
    %disp(out.Ts)
    
    %disp('EDT squared correlation coefficient:')
    %disp(out.EDTr2)
    %disp('T20 squared correlation coefficient:')
    %disp(out.T20r2)
    %disp('T30 squared correlation coefficient:')
    %disp(out.T30r2)
    %disp('Ratio of T20 to T30 (%):')
    %disp(out.T20T30ratio)
    %if exist('T20low','var')
    %    disp('Low-frequency EDT (s):')
    %    disp(out.EDTlow)
    %    disp('Low-frequency T20 (s):')
    %    disp(out.T20low)
    %    disp('Low-frequency T30 (s):')
    %    disp(out.T30low)
    %end
    %if exist('T20mid','var')
    %    disp('Mid-frequency EDT (s):')
    %    disp(out.EDTmid)
    %    disp('Mid-frequency T20 (s):')
    %    disp(out.T20mid)
    %    disp('Mid-frequency T30 (s):')
    %    disp(out.T30mid)
    %end
    %if exist('T20high','var')
    %    disp('High-frequency EDT (s):')
    %    disp(out.EDThigh)
    %    disp('High-frequency T20 (s):')
    %    disp(out.T20high)
    %    disp('High-frequency T30 (s):')
    %    disp(out.T30high)
    %end
    %if exist('BR_T20','var')
    %    disp('Bass ratio EDT:')
    %    disp(out.BR_EDT)
    %    disp('Bass ratio T20:')
    %    disp(out.BR_T20)
    %    disp('Bass ratio T30:')
    %    disp(out.BR_T30)
    %end
    %if exist('TR_T20','var')
    %    disp('Treble ratio EDT:')
    %    disp(out.TR_EDT)
    %    disp('Treble ratio T20:')
    %    disp(out.TR_T20)
    %    disp('Treble ratio T30:')
    %    disp(out.TR_T30)
    %end
    % end
    
    %--------------------------------------------------------------------------
    % AARAE TABLE
    %--------------------------------------------------------------------------
    out.tables = [];
    if isstruct(data)
        for ch = 1:chans
            f = figure('Name','Reverberation Parameters', ...
                'Position',[200 200 620 360]);
            dat1 = [out.EDT(ch,:);out.T20(ch,:);out.T30(ch,:);out.Tcustom(ch,:);...
                out.Tnonlin(ch,:);out.C50(ch,:);out.C80(ch,:);out.D50(ch,:); ...
                out.D80(ch,:);out.Ts(ch,:); ...
                out.EDTr2(ch,:);out.T20r2(ch,:);out.T30r2(ch,:);out.Tcustomr2(ch,:);...
                out.Tnonlinr2(ch,:);out.T20T30ratio(ch,:);out.PNR(ch,:)];
            cnames1 = num2cell(bandfc);
            rnames1 = {'Early decay time (s)',...
                'Reverberation time T20 (s)',...
                'Reverberation time T30 (s)',...
                ['Reverberation time ' num2str(-abs(customstart)) ' to ' num2str(-abs(customend)) ' dB (s)'],...
                'Reverberation time nonlinear (Xiang) (s)',...
                'Clarity index C50 (dB)',...
                'Clarity index C80 (dB)',...
                'Definition D50',...
                'Definition D80',...
                'Centre time Ts (s)',...
                'Correlation coefficient EDT r^2',...
                'Correlation coefficient T20 r^2',...
                'Correlation coefficient T30 r^2',...
                'Correlation coefficient Tcustom r^2',...
                'Correlation coefficient Tnonlinear r^2',...
                'Ratio of T20 to T30 %%',...
                'Dynamic range: Peak to last 10%% (dB)'};
            t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
            set(t1,'ColumnWidth',{60});
            
            if bpo == 1
                
                
                dat2 = [out.EDTlow(ch),out.EDTmid(ch),out.EDThigh(ch),out.BR_EDT(ch),out.TR_EDT(ch);...
                    out.T20low(ch),out.T20mid(ch),out.T20high(ch),out.BR_T20(ch),out.TR_T20(ch);...
                    out.T30low(ch),out.T30mid(ch),out.T30high(ch),out.BR_T30(ch),out.TR_T30(ch)];
                cnames2 = {'Low Freq','Mid Freq','High Freq','Bass Ratio','Treble Ratio'};
                rnames2 = {'EDT', 'T20', 'T30'};
                t2 =uitable('Data',dat2,'ColumnName',cnames2,'RowName',rnames2);
                set(t2,'ColumnWidth',{90});
                [~,tables] = disptables(f,[t1 t2],{['Ch' num2str(ch) ' - Reverb parameters'],['Ch' num2str(ch) ' - Reverb parameters summary']});
            else
                [~,tables] = disptables(f,t1,{['Ch' num2str(ch) ' - Reverb parameters']});
            end
            %out = [];
            out.tables = [out.tables tables];
        end
        
        %--------------------------------------------------------------------------
        % PLOTTING
        %--------------------------------------------------------------------------
        
        if doplot
            % Define number of rows and columns (for y axis labelling)
            [r,c] = subplotpositions(bands+1,0.4);
            
            
            % preallocate
            levdecayend = zeros(1,chans, bands);
            
            for dim2 = 1:chans
                for dim3 = 1:bands
                    irstart(1,dim2,dim3) = find(levdecay(:,dim2,dim3) <= 0, 1, 'first'); % 0 dB
                    tstart(1,dim2,dim3) = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                    edtend(1,dim2,dim3) = find(levdecay(:,dim2,dim3) <= -10, 1, 'first'); % -10 dB
                    t20end(1,dim2,dim3) = find(levdecay(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
                    t30end(1,dim2,dim3) = find(levdecay(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
                    levdecayend(1,dim2,dim3) = length(levdecay(:,dim2,dim3)); % time at last sample
                    levdecayend(1,dim2,dim3) = levdecayend(1,dim2,dim3)./fs;
                end
            end
            
            
            for ch = 1:chans
                
                figure('Name',['Channel ', num2str(ch), ', Level Decay and Regression Lines'])
                
                for band = 1:bands
                    
                    
                    subplot(r,c,band)
                    
                    
                    hold on
                    
                    % plot the level decay(s) on a single subplot
                    plot(((1:len)-1)./fs, levdecay(:,ch,band),'Color',[0.2 0.2 0.2], ...
                        'LineStyle',':','DisplayName','Level Decay')
                    if (noisecomp == 0) || (noisecomp == 1) || (noisecomp == 2)
                        % linear regression for T30
                        plot(((tstart(1,ch,band):t30end(1,ch,band))./fs), ...
                            (tstart(1,ch,band):t30end(1,ch,band)).* ...
                            q(1,ch,band)+q(2,ch,band), ...
                            'Color',[0 0 0.6],'DisplayName', 'T30')
                        
                        % linear regression for T20
                        plot(((tstart(1,ch,band):t20end(1,ch,band))./fs), ...
                            (tstart(1,ch,band):t20end(1,ch,band)).* ...
                            p(1,ch,band)+p(2,ch,band), ...
                            'Color',[0 0.6 0],'DisplayName','T20')
                        
                        % linear regression for EDT
                        plot(((irstart(1,ch,band):edtend(1,ch,band))./fs), ...
                            (irstart(1,ch,band):edtend(1,ch,band)).* ...
                            o(1,ch,band)+o(2,ch,band), ...
                            'Color',[0.9 0 0],'DisplayName','EDT')
                    end
%                     if (noisecomp == 3) || (noisecomp == 4)
%                         
%                         plot(t(irstart(1,ch,band):dataend(1,ch,band),ch,band)...
%                             ,10*log10((q(1,ch,band).* ...
%                             exp(-abs(q(2,ch,band))...
%                             .*t(irstart(1,ch,band):dataend(1,ch,band),ch,band)) + ...
%                             (q(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
%                             -t(irstart(1,ch,band):dataend(1,ch,band),ch,band))))), ...
%                             'Color',[0.7 0 0],'DisplayName', 'Nonlinear Fit')
%                         
%                         plot(t(irstart(1,ch,band):dataend(1,ch,band),ch,band),...
%                             10*log10(q(1,ch,band).* ...
%                             exp(-abs(q(2,ch,band))...
%                             .*t(irstart(1,ch,band):dataend(1,ch,band),ch,band))),...
%                             'LineStyle',':','Color',[0 0 0.7],'DisplayName', 'Linear fit no B3')
%                     end
%                     if (noisecomp == 5)
%                         
%                         
%                         plot(t(tstart(1,ch,band):t30end(1,ch,band),ch,band)...
%                             ,10*log10((q(1,ch,band).* ...
%                             exp(-abs(q(2,ch,band))...
%                             .*t(tstart(1,ch,band):t30end(1,ch,band),ch,band)) + ...
%                             (q(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
%                             -t(tstart(1,ch,band):t30end(1,ch,band),ch,band))))), ...
%                             'Color',[0 0 0.7],'DisplayName', 'T30')
%                         
%                         plot(t(tstart(1,ch,band):t30end(1,ch,band),ch,band),...
%                             10*log10(q(1,ch,band).* ...
%                             exp(-abs(q(2,ch,band))...
%                             .*t(tstart(1,ch,band):t30end(1,ch,band),ch,band))),...
%                             'LineStyle',':','Color',[0 0 0.7],'DisplayName', 'T30 no B3')
%                         
%                         plot(t(tstart(1,ch,band):t20end(1,ch,band),ch,band),...
%                             10*log10((p(1,ch,band).* ...
%                             exp(-abs(p(2,ch,band))...
%                             .*t(tstart(1,ch,band):t20end(1,ch,band),ch,band)) + ...
%                             (p(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
%                             -t(tstart(1,ch,band):t20end(1,ch,band),ch,band))))), ...
%                             'Color',[0 0.7 0],'DisplayName', 'T20')
%                         
%                         plot(t(tstart(1,ch,band):t20end(1,ch,band),ch,band),...
%                             10*log10(p(1,ch,band).* ...
%                             exp(-abs(p(2,ch,band))...
%                             .*t(tstart(1,ch,band):t20end(1,ch,band),ch,band))),...
%                             'LineStyle',':','Color',[0 0.7 0],'DisplayName', 'T20 no B3')
%                         
%                         plot(t(irstart(1,ch,band):edtend(1,ch,band),ch,band)...
%                             ,10*log10((o(1,ch,band).* ...
%                             exp(-abs(o(2,ch,band))...
%                             .*t(irstart(1,ch,band):edtend(1,ch,band),ch,band)) + ...
%                             (o(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
%                             -t(irstart(1,ch,band):edtend(1,ch,band),ch,band))))), ...
%                             'Color',[0.7 0 0],'DisplayName', 'EDT')
%                         
%                         plot(t(irstart(1,ch,band):edtend(1,ch,band),ch,band),...
%                             10*log10(o(1,ch,band).* ...
%                             exp(-abs(o(2,ch,band))...
%                             .*t(irstart(1,ch,band):edtend(1,ch,band),ch,band))),...
%                             'LineStyle',':','Color',[0.7 0 0],'DisplayName', 'EDT no B3')
%                         
%                         
%                     end
                    
                    if exist('crosspoint','var')
                        plot(crosspoint(1,ch,band)./fs,...
                            levdecay(crosspoint(1,ch,band),ch,band),...
                            'mo')
                    end
                    
                    % x axis label (only on the bottom row of subplots)
                    if band > (c*r - c)
                        xlabel('Time (s)')
                    end
                    
                    % y axis label (only on the left column of subplots)
                    if mod(band-1, c) == 0
                        ylabel('Level (dB)')
                    end
                    
                    xlim([0 levdecayend(1,ch,band)])
                    ylim([-65 0])
                    
%                     if (noisecomp == 0) || (noisecomp == 1) || (noisecomp == 2)|| (noisecomp == 5)
                        % text on subplots
                        text(levdecayend(1,ch,band)*0.45,-5,...
                            ['EDT ',num2str(0.01*round(100*EDT(1,ch,band)))],'Color',[0.9 0 0])
                        
                        text(levdecayend(1,ch,band)*0.45,-10,...
                            ['T20 ',num2str(0.01*round(100*T20(1,ch,band)))],'Color',[0 0.6 0])
                        
                        text(levdecayend(1,ch,band)*0.45,-15,...
                            ['T30 ',num2str(0.01*round(100*T30(1,ch,band)))],'Color',[0 0 0.6])
%                     elseif (noisecomp == 3) || (noisecomp == 4)
%                         text(levdecayend(1,ch,band)*0.45,-5,...
%                             ['T ',num2str(0.01*round(100*EDT(1,ch,band)))],'Color',[0 0 0.6])
%                     end
                    
                    
                    % subplot title
                    title([num2str(bandfc(band)),' Hz'])
                    
                end % for band
                
                % DIY legend
                subplot(r,c,band+1)
                
                plot([0.1,0.4], [0.8,0.8],'Color',[0.2 0.2 0.2], ...
                    'LineStyle',':','DisplayName','Level Decay')
                xlim([0,1])
                ylim([0,1])
                hold on
                text(0.5,0.8,'Decay');
%                 if (noisecomp == 0) ||(noisecomp == 1) ||(noisecomp == 2) ||(noisecomp == 5)
                    plot([0.1,0.4], [0.6,0.6],'Color',[0.9 0 0],'DisplayName','EDT')
                    text(0.5,0.6,'EDT');
                    plot([0.1,0.4], [0.4,0.4],'Color',[0 0.6 0],'DisplayName','T20')
                    text(0.5,0.4,'T20');
                    plot([0.1,0.4], [0.2,0.2],'Color',[0 0 0.6],'DisplayName', 'T30')
                    text(0.5,0.2,'T30');
%                 elseif(noisecomp == 3) ||(noisecomp == 4)
%                     plot([0.1,0.4], [0.6,0.6],'Color',[0.7 0 0],'DisplayName','Nlin')
%                     text(0.5,0.6,'Non-lin');
%                     plot([0.1,0.4], [0.4,0.4],'LineStyle',':','Color',[0 0 0.7],'DisplayName','Lin')
%                     text(0.5,0.4,'Linear');
%                 end
                set(gca,'YTickLabel','',...
                    'YTick',zeros(1,0),...
                    'XTickLabel','',...
                    'XTick',zeros(1,0))
            end
            hold off
            
            if nonlin > 0
                for ch = 1:chans
                    
                    figure('Name',['Channel ', num2str(ch), ', Level Decay Nonlinear Regression'])
                    
                    for band = 1:bands
                        
                        
                        subplot(r,c,band)
                        
                        
                        hold on
                        
                        % plot the level decay(s) on a single subplot
                        plot(((1:len)-1)./fs, levdecay(:,ch,band),'Color',[0.2 0.2 0.2], ...
                            'LineStyle',':','DisplayName','Level Decay')
                        
                        
                            
                        plot(t(irstart(1,ch,band):dataend(1,ch,band),ch,band)...
                            ,10*log10((xx(1,ch,band).* ...
                            exp(-abs(xx(2,ch,band))...
                            .*t(irstart(1,ch,band):dataend(1,ch,band),ch,band)) + ...
                            (xx(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
                            -t(irstart(1,ch,band):dataend(1,ch,band),ch,band))))), ...
                            'Color',[0.7 0 0],'DisplayName', 'Nonlinear Fit')
                        
                        plot(t(irstart(1,ch,band):dataend(1,ch,band),ch,band),...
                            10*log10(xx(1,ch,band).* ...
                            exp(-abs(xx(2,ch,band))...
                            .*t(irstart(1,ch,band):dataend(1,ch,band),ch,band))),...
                            'Color',[0 0 0.7],'DisplayName', 'Linear fit no B3')
                        
                        
                        
                        if exist('crosspoint','var')
                            plot(crosspoint(1,ch,band)./fs,...
                                levdecay(crosspoint(1,ch,band),ch,band),...
                                'mo')
                        end
                        
                        % x axis label (only on the bottom row of subplots)
                        if band > (c*r - c)
                            xlabel('Time (s)')
                        end
                        
                        % y axis label (only on the left column of subplots)
                        if mod(band-1, c) == 0
                            ylabel('Level (dB)')
                        end
                        
                        xlim([0 levdecayend(1,ch,band)])
                        ylim([-65 0])
                        

                        text(levdecayend(1,ch,band)*0.45,-5,...
                            ['T ',num2str(0.01*round(100*Tnonlin(1,ch,band)))],'Color',[0 0 0.6])

                        
                        
                        % subplot title
                        title([num2str(bandfc(band)),' Hz'])
                        
                    end % for band
                    
                    % DIY legend
                    subplot(r,c,band+1)
                    
                    plot([0.1,0.4], [0.8,0.8],'Color',[0.2 0.2 0.2], ...
                        'LineStyle',':','DisplayName','Level Decay')
                    xlim([0,1])
                    ylim([0,1])
                    hold on
                    text(0.5,0.8,'Decay');
                    plot([0.1,0.4], [0.6,0.6],'Color',[0.7 0 0],'DisplayName','Nlin')
                    text(0.5,0.6,'Non-lin');
                    plot([0.1,0.4], [0.4,0.4],'LineStyle',':','Color',[0 0 0.7],'DisplayName','Lin')
                    text(0.5,0.4,'Linear');   
                    
                    
                    set(gca,'YTickLabel','',...
                        'YTick',zeros(1,0),...
                        'XTickLabel','',...
                        'XTick',zeros(1,0))
                end
                hold off
            end
            
            
        end
        if isstruct(data)
            doresultleaf(levdecay,'Level [dB]',{'Time'},...
                'Time',      ((1:len)-1)./fs,  's',           true,...
                'channels',  chanID,           'categorical', [],...
                'Frequency', num2cell(bandfc), 'Hz',          false,...
                'name','Reverb_time');
        end
    end
    out.funcallback.name = 'ReverberationTime_IR1.m';
    out.funcallback.inarg = {fs,startthresh,bpo,doplot,filterstrength,phasemode,noisecomp,autotrunc,f_low,f_hi,nonlin,customstart,customend};
else
    out = [];
end
% eof

