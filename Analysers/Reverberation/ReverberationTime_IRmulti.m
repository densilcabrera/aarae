function out = ReverberationTime_IRmulti(data,fs,startthresh,bpo,doplot,filterstrength,phasemode,noisecomp,autotrunc,f_low,f_hi,averagingmethod,SNR,additionalaudio)
% This function calculates reverberation time parameters from multiple
% impulse responses, using the an average of their squared decay functions.
% This may yield a more representative reverberation time than averaging
% values calculated from individual impulse responses.
%
% The intended use for this is spatial averaging in room acoustics
% measurements. It may be instructive to compare the results from this
% function with reverberation times of individual impulse responses
% (calculared from AARAE's ReverberationTime_IR1).
%
% All impulse responses in the audio input, and in additional inputs via
% dialog box (if used), are squared and combined to yield the result, apart
% from the cases mentioned below. This includes all channels, and higher
% dimensions (except for dimension 3, which is reserved for bands). If
% multiband audio is input, the bands are mixed (and the filterbank is
% applied to the mixed bands). Silent cycles are not included in the
% averaging (by detecting that properties.relgain == -inf). Also a
% signal-to-noise criterion is applied, and, after applying the
% filterbank,those channels/bands that do not meet the criterion are not
% included in the average.
%
% To use this function from the AARAE GUI, you should select ONE audio leaf
% (if you select more, then the function will run multiple times). If you
% wish to load additional audio, this is done via the dialog box of this
% function. Multiple additional audio inputs can be chosen (but if the
% primary audio input is chosen again, it will not be used twice if the
% analysis: audio with the same .name field string as the primary input is
% ignored).
%
% Impulse responses should normally be reasonably well prepared (by
% appropriate truncation) prior to calling this function. Bear in mind that
% the length of the averaged impulse response will be taken from the
% minimum length of the audio inputs. Also all of the individual IRs are
% shifted (independently) to start at the beginning of the analysis,
% approximately in synchrony. The IRs can be gain adjusted within this
% function (no adjustment, normalized, or equal energy) - see
% 'averagingmethod' below.
%
% This function is adapted from ReverberationTime_IR1 (by Grant Cuthbert &
% Densil Cabrera).
% Version 1.00 (25 June 2015)
%
%
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
% 3: Nonlinear fitting of entire reverse-integrated decay, following Xiang
% 4: Nonlinear fitting of entire reverse-integrated decay in dB, similar to
% Xiang
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
% averagingmethod = how the multiple impulse responses are gain adjusted
% for averaging:
% 0: no gain adjustment prior to averaging
% 1: normalize all IRs prior to averaging
% 2: adjust IRs for equal energy prior to averaging (default)
%
% SNR = the minimum signal-to-noise ratio, in dB, required for an impulse
% response in any band to be included in the average. The SNR value is
% taken as the level difference between the IR peak and the power of the
% last 10% of the IR (assumed to be noise) of each channel and band.
%
% additionalaudio = whether or not additional audio is input (by calling
% AARAE's utility function choose_audio
% 0: only primary audio input
% 1: additional audio

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2013-2015, Grant Cuthbert and Densil Cabrera
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
if nargin < 14, SNR = 10; end % minimum raw SNR in each chan & band (if this requirement is not met, the chan/band IR will not be used)
if nargin < 13, averagingmethod = 2; end
if nargin < 12, additionalaudio = 0; end
if nargin < 10, f_low = 125; f_hi = 8000; end
if nargin < 9, autotrunc = 1; end
if nargin < 8, noisecomp = 1; end
if nargin < 7, phasemode = -1; end
if nargin < 6, filterstrength = 2; end
if nargin < 5, doplot = 1; end
if nargin < 4, bpo = 1; end
if nargin < 3
    % startthresh = -20;
    % Dialog box for settings
    prompt = {'Input additional impulse responses (0 | 1)', ...
        'Ensemble averaging method: raw (0); normalized (1) ; or constant energy (2)', ...
        'Minimum signal-to-noise ratio for inclusion in the average (dB)',...
        'Threshold for IR start detection', ...
        'Bands per octave (0 | 1 | 3)', ...
        'Lowest centre frequency (Hz)', ...
        'Highest centre frequency (Hz)', ...
        'Filter strength', ...
        'Zero phase (0), Maximum phase (-1) or Minimum phase (1) filters',...
        'Noise compensation: None (0), Subtract final 10% Chu (1), Extrapolate Late Decay (2), Nonlinear curve fit of whole decay Xiang (3), Nonlinear curve fit of whole decay in dB (4), Nonlinear curve fit of evaluation ranges in dB (5)', ...
        'Automatic end truncation: None (0), Lundeby (1)',...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = [1 60];
    def = {'0','2','10','-20',num2str(bpo),num2str(f_low),num2str(f_hi),...
        '2','-1','1','1','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if length(answer) < 12, answer = []; end
    if ~isempty(answer)
        additionalaudio = str2num(answer{1,1});
        averagingmethod = str2num(answer{2,1});
        SNR = str2num(answer{3,1});
        startthresh = str2num(answer{4,1});
        bpo = str2num(answer{5,1});
        f_low = str2num(answer{6,1});
        f_hi = str2num(answer{7,1});
        filterstrength = str2num(answer{8,1});
        phasemode = str2num(answer{9,1});
        noisecomp = str2num(answer{10,1});
        autotrunc = str2num(answer{11,1});
        doplot = str2num(answer{12,1});
    else
        out = [];
        return
    end
end

if ~exist('bpo','var')
    bpo=1;
    f_low = 125;
    f_hi = 8000;
end


if isstruct(data)
    [mindim1,dim2,dim3,dim4,dim5,dim6] = size(data.audio);
    fs = data.fs;
    if dim3>1
        data.audio = sum(data.audio,3); % mixdown bands
    end
    if dim4 > 1
        if isfield(data,'properties')
            if isfield(data.properties,'relgain')
                if isinf(data.properties.relgain(1))
                    % check for silent cycle
                    dim4 = dim4-1;
                    if dim4 < 1
                        h=warndlg('All of the audio seems to be from a silent cycle (properties.relgain = -inf) - there is nothing suitable for analysis.','AARAE info','modal');
                        uiwait(h)
                        out = []; % consider including function callback in the output
                        return
                    else
                        % discard silent cycle
                        data.audio = data.audio(:,:,1,2:end,:,:);
                    end
                end
            end
        end
    end
    data.audio = reshape(data.audio,mindim1,dim2*dim4*dim5*dim6);
    if isfield(data,'name') % Get the AARAE name if it exists
        name1 = data.name;
    else
        name1 = [];
    end
else
    [mindim1,dim2,dim3,dim4,dim5,dim6] = size(data);
    %fs = fs;
    if dim3>1
        data = sum(data,3);
    end
    data = reshape(data,mindim1,dim2*dim4*dim5*dim6);
    name1 = [];
end

% Get additional audio and match fs and length to primary audio input
if additionalaudio
    additionalaudiodata = choose_audio('multiple');
    if ~isempty(additionalaudiodata)
    if iscell(additionalaudiodata)
        for i = 1:size(additionalaudiodata,1)
            % check the size and fs of each input.
            % if fs does not match then resample to that of primary input
            % if audio is identical to primary input, then discard
            fs2 = additionalaudiodata{i,1}.fs;
            if fs2 ~= fs
                gcd_fs = gcd(fs,fs2); % greatest common denominator
                additionalaudiodata{i,1}.audio =...
                    resample(additionalaudiodata{i,1}.audio,fs/gcd_fs,fs2/gcd_fs);
            end
            if size(additionalaudiodata{i,1},4) > 1
                if isfield(additionalaudiodata{i,1},'properties')
                    if isfield(additionalaudiodata{i,1}.properties,'relgain')
                        if isinf(additionalaudiodata{i,1}.properties.relgain(1))
                            % discard silent cycle
                            if size(additionalaudiodata{i,1},4)>1
                                additionalaudiodata{i,1}.audio = additionalaudiodata{i,1}.audio(:,:,:,2:end,:,:);
                            else
                                additionalaudiodata{i,1}.audio = additionalaudiodata{i,1}.audio*0;
                            end
                        end
                    end
                end
            end
            
            [len2,dim2b,dim3,dim4b,dim5b,dim6b] = size(additionalaudiodata{i,1}.audio);
            mindim1 = min([mindim1,len2]); % update maximum length
            dim2 = dim2+dim2b;
            dim4 = dim4+dim4b-1;
            dim5 = dim5+dim5b-1;
            dim6 = dim6+dim6b-1;
            % mixdown dimension 3
            if dim3 > 1
                additionalaudiodata{i,1}.audio = sum(additionalaudiodata{i,1}.audio,3);
            end
            additionalaudiodata{i,1}.audio = ...
                reshape(additionalaudiodata{i,1}.audio,len2,dim2b*dim4b*dim5b*dim6b);
            if isfield(additionalaudiodata{i,1},'name') && ~isempty(name1)
                if strcmp(additionalaudiodata{i,1}.name,name1)
                    additionalaudiodata{i,1}.audio = 0*additionalaudiodata{i,1}.audio; % zero redundant input
                end
            end
        end
    else
        fs2 = additionalaudiodata.fs;
        if fs2 ~= fs
            gcd_fs = gcd(fs,fs2); % greatest common denominator
            additionalaudiodata.audio =...
                resample(additionalaudiodata.audio,fs/gcd_fs,fs2/gcd_fs);
        end
        
        if size(additionalaudiodata,4) > 1
            if isfield(additionalaudiodata,'properties')
                if isfield(additionalaudiodata.properties,'relgain')
                    if isinf(additionalaudiodata.properties.relgain(1))
                        % discard silent cycle
                        if size(additionalaudiodata,4)>1
                            additionalaudiodata.audio = additionalaudiodata.audio(:,:,:,2:end,:,:);
                        else
                            additionalaudiodata.audio = additionalaudiodata.audio*0;
                        end
                    end
                end
            end
        end
        [len2,dim2b,dim3,dim4b,dim5b,dim6b] = size(additionalaudiodata.audio);
        mindim1 = min([mindim1,len2]); % update maximum length
        dim2 = dim2+dim2b;
        dim4 = dim4+dim4b-1;
        dim5 = dim5+dim5b-1;
        dim6 = dim6+dim6b-1;
        % mixdown dimension 3
        if dim3 > 1
            additionalaudiodata.audio = sum(additionalaudiodata.audio,3);
        end
        additionalaudiodata.audio = ...
            reshape(additionalaudiodata.audio,len2,dim2b*dim4b*dim5b*dim6b);
        if isfield(additionalaudiodata,'name') && ~isempty(name1)
            if strcmp(additionalaudiodata.name,name1)
                additionalaudiodata.audio = 0*additionalaudiodata.audio; % zero redundant input
            end
        end
    
    end

    % Create one big matrix, reshaped to 2 dimensions
    ir = zeros(mindim1,dim2*dim4*dim5*dim6);
    
    % load primary audio into matrix
    if isstruct(data)
        ir(:,1:size(data.audio,2)) = data.audio(1:mindim1,:);
        chancount = size(data.audio,2);
    else
        ir(:,1:size(data,2)) = data(1:mindim1,:);
        chancount = size(data,2);
    end
    
    % load additional audio into matrix
    if iscell(additionalaudiodata)
        for i = 1:size(additionalaudiodata,1)
            ir(:,chancount+1:chancount+size(additionalaudiodata{i,1}.audio,2)) = additionalaudiodata{i,1}.audio(1:mindim1,:);
            chancount = chancount + size(additionalaudiodata{i,1}.audio,2);
        end
    else
        ir(:,chancount+1:end) = additionalaudiodata.audio(1:mindim1,:);
        %chancount = chancount + size(additionalaudiodata.audio,2);
    end
    
    % need to improve chanID so it identifies the input sources
    chanID = cellstr([repmat('Chan',size(ir,2),1) num2str((1:size(ir,2))')]);
    clear additionalaudiodata
    else
        additionalaudio = 0;
    end
end
   
if ~additionalaudio
    if isstruct(data)
        ir(:,1:size(data.audio,2)) = data.audio(1:mindim1,:);
    else
        ir(:,1:size(data,2)) = data(1:mindim1,:);
    end
    % need to improve chanID so it identifies the input sources
    chanID = cellstr([repmat('Chan',size(ir,2),1) num2str((1:size(ir,2))')]);
end






if ~isempty(ir) && ~isempty(fs) && ~isempty(startthresh) && ~isempty(bpo) && ~isempty(doplot) && ~isempty(filterstrength) && ~isempty(phasemode) && ~isempty(noisecomp) && ~isempty(autotrunc) && ~isempty(f_low) && ~isempty(f_hi)
    %--------------------------------------------------------------------------
    % START TRUNCATION & ALIGNMENT
    %--------------------------------------------------------------------------
    
    [len,chans,bands] = size(ir);
    
    if (bands>1)
        multibandIR = 1; % not possible
    else
        multibandIR = 0;
    end
    
    % Get last 10 % for Chu noise compensation if set (before the circular
    % shift is done as part of startpoint detection/allignment)
    
    ir_end10 = ir(round(0.9*len):end,:);
    
    
    if noisecomp == 2;
        autotrunc = 1;
    end
    
    % Preallocate
    m = zeros(1,chans); % maximum value of the IR
    startpoint = zeros(1,chans); % the auto-detected start time of the IR
    
    for dim2 = 1:chans
        m(dim2) = max(ir(:,dim2).^2); % maximum value of the IR
        startpoint(1,dim2) = find(ir(:,dim2).^2 >= m(1,dim2)./ ...
            (10^(abs(startthresh)/10)),1,'first'); % Define start point
        
        %for cases where startpoint was not found
        startpoint(isempty(startpoint))=1;
        
        %startpoint = min(startpoint,[],3);
        if startpoint(dim2) >1
            
            % zero the data before the startpoint
            ir(1:startpoint(dim2)-1,dim2) = 0;
            
            % rotate the zeros to the end (to keep a constant data length)
            ir(:,dim2) = circshift(ir(:,dim2),-(startpoint(dim2)-1));
            
        end % if startpoint
    end % for dim2
    
    %     early50 = ir(1:1+floor(fs*0.05),:); % Truncate Early80
    %     early80 = ir(1:1+floor(fs*0.08),:); % Truncate Early80
    %     late50 = ir(ceil(fs*0.05):end,:); % Truncate Late50
    %     late80 = ir(ceil(fs*0.08):end,:); % Truncate Late80
    
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
            %             early50oct = octbandfilter_viaFFT(early50,fs,bandfc,order,0,1000,0,phasemode);
            %             early80oct = octbandfilter_viaFFT(early80,fs,bandfc,order,0,1000,0,phasemode);
            %             late50oct = octbandfilter_viaFFT(late50,fs,bandfc,order,0,1000,0,phasemode);
            %             late80oct = octbandfilter_viaFFT(late80,fs,bandfc,order,0,1000,0,phasemode);
            
            ir_end10oct = octbandfilter_viaFFT(ir_end10,fs,bandfc,order,0,1000,0,phasemode);
            
            
        else
            order = [36,24] * filterstrength;
            iroct = thirdoctbandfilter_viaFFT(ir,fs,bandfc,order,0,1000,0,phasemode);
            %             early50oct = thirdoctbandfilter_viaFFT(early50,fs,bandfc,order,0,1000,0,phasemode);
            %             early80oct = thirdoctbandfilter_viaFFT(early80,fs,bandfc,order,0,1000,0,phasemode);
            %             late50oct = thirdoctbandfilter_viaFFT(late50,fs,bandfc,order,0,1000,0,phasemode);
            %             late80oct = thirdoctbandfilter_viaFFT(late80,fs,bandfc,order,0,1000,0,phasemode);
            
            ir_end10oct = thirdoctbandfilter_viaFFT(ir_end10,fs,bandfc,order,0,1000,0,phasemode);
            
        end
        
        % check the number of bands again, in case some were not ok to filter
        bands = size(iroct,3);
        fc = fc(1:bands);
        bandfc = bandfc(1:bands);
        
        
        % SQUARE IRs
        iroct = iroct.^2;
        ir_end10oct = ir_end10oct.^2;
        
        
        
        
        switch averagingmethod
            case 1
                % normalize
                iroct = iroct ./ repmat(max(iroct),[len,1,1]);
                ir_end10oct = ir_end10oct ./ repmat(max(iroct),[size(ir_end10oct,1),1,1]);
            case 2
                % equal power (same as equal energy because len is
                % constant)
                iroct = iroct ./ repmat(mean(iroct),[len,1,1]);
                ir_end10oct = ir_end10oct ./ repmat(mean(iroct),[size(ir_end10oct,1),1,1]);
        end

        
%         iroct(isinf(iroct)) = 0;
%         iroct(isnan(iroct)) = 0;
%         ir_end10oct(isinf(ir_end10oct)) = 0;
%         ir_end10oct(isnan(ir_end10oct)) = 0;
        
        % ZERO THE IRs THAT HAVE INADEQUATE SNR
        % this is to avoid unduly contaminating the average with noise.
        % however it carries a risk of eliminating too much data
        for ch = 1:chans
            for b = 1:bands
                if mean(iroct(:,ch,b)) == 0 || sum(isnan(iroct(:,ch,b))) > 0 || sum(isinf(iroct(:,ch,b))) > 0
                    iroct(:,ch,b) = 0;
                    ir_end10oct(:,ch,b) = 0;
                elseif pow2db(max(iroct(:,ch,b))) < pow2db(mean(ir_end10oct(:,ch,b)))+SNR
                    iroct(:,ch,b) = 0;
                    ir_end10oct(:,ch,b) = 0;
                end
            end
        end
        if sum(sum(sum(iroct))) == 0
            h=warndlg('Inadequate signal to noise ratio - there is nothing suitable for analysis.','AARAE info','modal');
            uiwait(h)
            out = []; % consider including function callback in the output
            return
        end
        
        
%         PNR = nan(1,ch,b);
%         for ch = 1:chans
%             for b = 1:bands
%                 if max(iroct(:,ch,b)) > 0
%                     PNR(1,ch,b) = 10*log10(max(iroct(:,ch,b))./mean(ir_end10oct(:,ch,b)));
%                 end
%             end
%         end
        % apply noise subtraction (Chu) to individual channels
        if noisecomp == 1
            iroct = iroct - repmat(mean(ir_end10oct),[len,1,1]);
            iroct(iroct<0)=0;
        end

        
        % AVERAGE SQUARED IRs
        iroct = mean(iroct,2);
        % peak-to-noise ratio (using last 10% as 'noise'
        PNR = 10*log10(max(iroct)./mean(mean(ir_end10oct,1),2));
        chans = 1;
        
        
        
        if autotrunc == 1
            % Lundeby crosspoint
            [crosspoint, Tlate, okcrosspoint] =...
                lundebycrosspoint(iroct(1:(end-max(max(max(startpoint)))),:,:), fs,fc);
            
            % autotruncation
            for ch = 1:chans
                for b = 1:bands
                    if crosspoint(1,ch,b) < len
                        iroct(crosspoint(1,ch,b):end,ch,b) = 0;
                        %                         if crosspoint(1,ch,b) > fs*0.05
                        %                             late50oct(round(crosspoint(1,ch,b)-fs*0.05+1):end,ch,b)=0;
                        %                         end
                        %                         if crosspoint(1,ch,b) > fs*0.08
                        %                             late80oct(round(crosspoint(1,ch,b)-fs*0.08+1):end,ch,b)=0;
                        %                         end
                    end
                end
            end
        else
            % zero the data that should be zero from startpoint detection (but
            % isn't due to filter ringing)
            % this is especially important for non-linear curve fitting
            for ch = 1:chans
                iroct(end-startpoint(1,ch,1):end,ch,:) = 0;
            end
        end
        
        %----------------------------------------------------------------------
        % CALCULATE ENERGY PARAMETERS
        %----------------------------------------------------------------------
        
        %         early50oct = mean(sum(early50oct.^2),2);
        %         early80oct = mean(sum(early80oct.^2),2);
        %         late50oct = mean(sum(late50oct.^2),2);
        %         late80oct = mean(sum(late80oct.^2),2);
        %         alloct = sum(iroct); % this has already been squared and channel-averaged
        %
        %
        %         C50 = 10*log10(early50oct ./ late50oct); % C50
        %         C80 = 10*log10(early80oct ./ late80oct); % C80
        %         D50 = (early50oct ./ alloct); % D50
        %         D80 = (early80oct ./ alloct); % D80
        
        % time values of IR in seconds
        %         tstimes = (0:(length(iroct)-1))' ./ fs;
        %
        %         Ts = (sum(iroct .* ...
        %             repmat(tstimes,[1,chans,bands])))./alloct; % Ts
        %     end
        
        
        
        
    end % if multibandIR == 0
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    
    
    
    %--------------------------------------------------------------------------
    % REVERBERATION TIME (Reverse Integration, Linear Regression)
    %--------------------------------------------------------------------------
    
    %***************************************
    % Derive the reverse-integrated decay curve(s)
    
    
    iroct(iroct<0) = 0;
    
    if noisecomp == 2
        % extrapolate decay beyond crosspoint
        tau = 0.8*Tlate./log(1e6); % 0.8 is more conservative than 1
        for ch = 1:chans
            for b = 1:bands
                if okcrosspoint(1,ch,b)
                    endmean = mean(iroct(crosspoint(1,ch,b)-50:crosspoint(1,ch,b)-1,ch,b));
                    iroct((crosspoint(1,ch,b)):end,ch,b) = ...
                        endmean .* exp(((((crosspoint(1,ch,b)):len)...
                        -(crosspoint(1,ch,b)+25))./fs)./-tau(1,ch,b));
                end
            end
        end
    end
    
    % Reverse integrate squared IR, and express the result in decibels
    levdecay = 10*log10(flip(cumsum(flip(iroct,1)),1));
    
    for dim2 = 1:chans
        for dim3 = 1:bands
            % Adjust so that the IR starts at 0 dB
            levdecay(:,dim2,dim3) = levdecay(:,dim2,dim3) - levdecay(1,dim2,dim3);
        end
    end
    
    if (noisecomp == 0) || (noisecomp == 1) || (noisecomp == 2)
        %***************************************
        % Derive decay times, via linear regressions over the appropriate
        % evaluation ranges
        
        % preallocate
        [o,p,q] = deal(zeros(2, chans, bands));
        [EDT,T20,T30,EDTr2,T20r2,T30r2] = deal(zeros(1, chans, bands));
        
        for dim2 = 1:chans
            for dim3 = 1:bands
                
                % find the indices for the relevant start and end samples
                irstart = find(levdecay(:,dim2,dim3) <= 0, 1, 'first'); % 0 dB
                tstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                edtend = find(levdecay(:,dim2,dim3) <= -10, 1, 'first'); % -10 dB
                t20end = find(levdecay(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
                t30end = find(levdecay(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
                
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
                
                
                
            end % dim3
        end % dim2
        
    elseif (noisecomp == 3) || (noisecomp == 4) || (noisecomp == 5)
        % nonlinear curve fit following Xiang 1995
        
        
        % nonlinear curve fitting seed coefficients
        b1 = 1*rand;
        b2 = -20*rand;
        b3 = 1e-7*rand;
        b = [b1;b2;b3];
        maxloopcount = 10; % maximum number of fitting attempts for each parameter
        jthreshold = 0.9; % value less than 1, the higher the value, the stricter the threshold
        
        % preallocate
        t = zeros(len, chans, bands);
        o = zeros(3, chans, bands);
        p = zeros(3, chans, bands);
        q = zeros(3, chans, bands);
        levdecaymodelEDT = zeros(len, chans, bands);
        levdecaymodelT20 = zeros(len, dim2, dim3);
        levdecaymodelT30 = zeros(len, dim2, dim3);
        edtmodlen = zeros(1, chans, bands);
        t20modlen = zeros(1, chans, bands);
        t30modlen = zeros(1, chans, bands);
        EDT = zeros(1, chans, bands);
        T20 = zeros(1, chans, bands);
        T30 = zeros(1, chans, bands);
        EDTr2 = zeros(1, chans, bands);
        T20r2 = zeros(1, chans, bands);
        T30r2 = zeros(1, chans, bands);
        dataend1 = zeros(1, chans, bands);
        
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
                irstart = find(levdecay(:,dim2,dim3) <= 0, 1, 'first'); % 0 dB
                edtend = find(levdecay(:,dim2,dim3) <= -10, 1, 'first'); % -10 dB
                t(1:dataend(:,dim2,dim3),dim2,dim3) = ((1:dataend(:,dim2,dim3))-1)./fs;
                
                irlen = (t(dataend(:,dim2,dim3),dim2,dim3)); % duration of ir
                
                
                % decay curve model
                if noisecomp == 3
                    % use entire data length for method 3
                    % (modelling does not work well if just EDT, T20 and T30
                    % lengths are used using this method)
                    y = 10.^(levdecay((irstart:dataend(1,dim2,dim3)),dim2,dim3)./10);
                    x = t((irstart:dataend(1,dim2,dim3)),dim2,dim3); % discrete time
                    modelfun = @(B,x)(abs(B(1)).*exp(-abs(B(2)).*x)...
                        + abs(B(3)).*(irlen-x));
                elseif noisecomp == 4
                    % use entire length for method 4, but curve-fit in decibels
                    y = levdecay((irstart:dataend(1,dim2,dim3)),dim2,dim3);
                    x = t((irstart:dataend(1,dim2,dim3)),dim2,dim3); % discrete time
                    modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
                        + abs(B(3)).*(irlen-x)));
                elseif noisecomp == 5
                    % use edt length for method 5
                    y = levdecay((irstart:edtend),dim2,dim3);
                    x = t((irstart:edtend),dim2,dim3); % discrete time
                    modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
                        + abs(B(3)).*(irlen-x)));
                end
                
                % derive coefficients
                fitted = 0;
                loopcount = 0;
                while ~fitted
                    loopcount = loopcount+1;
                    [o(:,dim2,dim3),~,j] = nlinfit(x,y,modelfun,b);
                    jflag = min(sum(j.^2));
                    if jflag > jthreshold || loopcount == maxloopcount,
                        fitted = 1;
                        b = o(:,dim2,dim3)';
                    else
                        b = [100*rand;-100*rand;0.1*rand];
                        disp('*')
                    end
                    
                end
                o(1,dim2,dim3) = abs(o(1,dim2,dim3));
                o(2,dim2,dim3) = -abs(o(2,dim2,dim3));
                o(3,dim2,dim3) = abs(o(3,dim2,dim3));
                disp(['EDT nonlinear parameters: ',num2str(o(1,dim2,dim3)),...
                    ' ',num2str(o(2,dim2,dim3)),' ',num2str(o(3,dim2,dim3))])
                % create acoustic parameter-specific model using coefficients
                levdecaymodelEDT(:,dim2,dim3) = ...
                    10*log10(o(1,dim2,dim3).*exp(-abs(o(2,dim2,dim3)).*t(:,dim2,dim3)));
                levdecaymodelEDT(:,dim2,dim3) = ...
                    levdecaymodelEDT(:,dim2,dim3)-max(levdecaymodelEDT(:,dim2,dim3));
                
                startmodelEDT = ...
                    find(levdecaymodelEDT(:,dim2,dim3) <= 0, 1, 'first'); % 0 dB
                endmodelEDT = ...
                    find(levdecaymodelEDT(:,dim2,dim3) <= -10, 1, 'first'); % -10 dB
                
                % length of specific model
                edtmodlen(:,dim2,dim3) = length(t(startmodelEDT:endmodelEDT,dim2,dim3));
                
                EDT(1,dim2,dim3) = 6.*(edtmodlen(1,dim2,dim3))./fs; % EDT in seconds
                
                EDTr2(1,dim2,dim3) = corr(levdecay(irstart:edtend,dim2,dim3), ...
                    (irstart:edtend)' * o(1,dim2,dim3) + ...
                    o(2,dim2,dim3)).^2; % correlation coefficient, EDT
                
                
            end % dim3
        end % dim2
        
        % no need for further fitting etc for noisecomp 3 and 4
        if (noisecomp == 3) || (noisecomp ==4)
            [p,q] = deal(o);
            %[T20,T30] = deal(EDT);
            disp('Noise compensation methods 3 and 4 use a non-linear regression over the entire reverse-integrated decay.')
            disp('Hence the evaluation ranges of EDT, T20 and T30 are not used in their full sense,')
            disp('and their values are likely to be very similar to each other.')
        end
        
        % T20 nonlinear model
        for dim2 = 1:chans
            for dim3 = 1:bands
                tstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                t20end = find(levdecay(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
                if (isempty(t20end)) || (isempty(tstart))|| (t20end > dataend(:,dim2,dim3))
                    T20(1,dim2,dim3) = NaN;
                    T20r2(1,dim2,dim3) = NaN;
                    p(:,dim2,dim3) = [0;0;0];
                    %disp('insufficient dynamic range for T20')
                else
                    irlen = (t(dataend(:,dim2,dim3),dim2,dim3)); % finite length of IR (ULI)
                    
                    % decay curve model
                    
                    if noisecomp == 5
                        % vector of response (dependent variable) values
                        x = t((tstart:t20end),dim2,dim3); % discrete time
                        y = levdecay((tstart:t20end),dim2,dim3);
                        modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
                            + abs(B(3)).*(irlen-x)));
                        
                        % derive coefficients
                        fitted = 0;
                        loopcount = 0;
                        while ~fitted
                            loopcount = loopcount+1;
                            [p(:,dim2,dim3),~,j] = nlinfit(x,y,modelfun,b);
                            jflag = min(sum(j.^2));
                            if jflag > jthreshold || loopcount == maxloopcount,
                                fitted = 1;
                                b = p(:,dim2,dim3)';
                            else
                                b = [100*rand;100*rand;rand];
                                disp('*')
                            end
                            
                        end
                        p(1,dim2,dim3) = abs(p(1,dim2,dim3));
                        p(2,dim2,dim3) = -abs(p(2,dim2,dim3));
                        p(3,dim2,dim3) = abs(p(3,dim2,dim3));
                        
                        disp(['T20 nonlinear parameters: ',num2str(p(1,dim2,dim3)),...
                            ' ',num2str(p(2,dim2,dim3)),' ',num2str(p(3,dim2,dim3))])
                    end
                    % create acoustic parameter-specific model using coefficients
                    levdecaymodelT20(:,dim2,dim3) = ...
                        10*log10(p(1,dim2,dim3).*exp(-abs(p(2,dim2,dim3)).*t(:,dim2,dim3)));
                    levdecaymodelT20(:,dim2,dim3) = ...
                        levdecaymodelT20(:,dim2,dim3)-max(levdecaymodelT20(:,dim2,dim3));
                    
                    startmodelT20 = ...
                        find(levdecaymodelT20(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                    endmodelT20 = ...
                        find(levdecaymodelT20(:,dim2,dim3) <= -25, 1, 'first'); % -25 dB
                    t20modlen(:,dim2,dim3) = length(t(startmodelT20:endmodelT20,dim2,dim3));
                    
                    T20(1,dim2,dim3) = 3.*(t20modlen(1,dim2,dim3))./fs; % T20 in seconds
                    
                    T20r2(1,dim2,dim3) = corr(levdecay(tstart:t20end,dim2,dim3), ...
                        (tstart:t20end)'*p(1,dim2,dim3) ...
                        + p(2,dim2,dim3)).^2; % correlation coefficient, T20
                end
            end % dim 3
        end % dim 2
        
        
        % T30 nonlinear model
        for dim2 = 1:chans
            for dim3 = 1:bands
                %tstart = find(levdecay(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                t30end = find(levdecay(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
                if (isempty(t30end)) || (isempty(tstart)) || (t30end > dataend(:,dim2,dim3))
                    T30(1,dim2,dim3) = NaN;
                    T30r2(1,dim2,dim3) = NaN;
                    q(:,dim2,dim3) = [0;0;0];
                    %disp('insufficient dynamic range for T30')
                else
                    
                    irlen = (t(dataend(:,dim2,dim3),dim2,dim3)); % finite length of IR (ULI)
                    
                    
                    % decay curve model
                    
                    if noisecomp == 5
                        % vector of response (dependent variable) values
                        x = t((tstart:t30end),dim2,dim3); % discrete time
                        y = levdecay((tstart:t30end),dim2,dim3);
                        modelfun = @(B,x)(10*log10(abs(B(1)).*exp(-abs(B(2)).*x)...
                            + abs(B(3)).*(irlen-x)));
                        
                        
                        % derive coefficients
                        fitted = 0;
                        loopcount = 0;
                        while ~fitted
                            loopcount = loopcount+1;
                            [q(:,dim2,dim3),~,j] = nlinfit(x,y,modelfun,b);
                            jflag = min(sum(j.^2));
                            if jflag > jthreshold || loopcount == maxloopcount,
                                fitted = 1;
                                b = q(:,dim2,dim3)';
                            else
                                b = [100*rand;100*rand;rand];
                                disp('*')
                            end
                            
                        end
                        q(1,dim2,dim3) = abs(q(1,dim2,dim3));
                        q(2,dim2,dim3) = -abs(q(2,dim2,dim3));
                        q(3,dim2,dim3) = abs(q(3,dim2,dim3));
                        disp(['T30 nonlinear parameters: ',num2str(q(1,dim2,dim3)),...
                            ' ',num2str(q(2,dim2,dim3)),' ',num2str(q(3,dim2,dim3))])
                        
                    end
                    % create acoustic parameter-specific model using coefficients
                    
                    levdecaymodelT30(:,dim2,dim3) = ...
                        10*log10(q(1,dim2,dim3).*exp(-abs(q(2,dim2,dim3)).*t(:,dim2,dim3)));
                    
                    levdecaymodelT30(:,dim2,dim3) = ...
                        levdecaymodelT30(:,dim2,dim3)-max(levdecaymodelT30(:,dim2,dim3));
                    
                    startmodelT30 = ...
                        find(levdecaymodelT30(:,dim2,dim3) <= -5, 1, 'first'); % -5 dB
                    endmodelT30 = ...
                        find(levdecaymodelT30(:,dim2,dim3) <= -35, 1, 'first'); % -35 dB
                    t30modlen(:,dim2,dim3) = length(t(startmodelT30:endmodelT30,dim2,dim3));
                    
                    T30(1,dim2,dim3) = 2.*(t30modlen(1,dim2,dim3))./fs; % T30 in seconds
                    
                    T30r2(1,dim2,dim3) = corr(levdecay(tstart:t30end,dim2,dim3), ...
                        (tstart:t30end)'*q(1,dim2,dim3) ...
                        + q(2,dim2,dim3)).^2; % correlation coefficient, T30
                end
            end % dim 3
        end % dim 2
    end
    
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
    %     out.C50 = permute(C50,[2,3,1]);
    %     out.C80 = permute(C80,[2,3,1]);
    %     out.D50 = permute(D50,[2,3,1]);
    %     out.D80 = permute(D80,[2,3,1]);
    %     out.Ts = permute(Ts,[2,3,1]);
    out.EDTr2 = permute(EDTr2,[2,3,1]);
    out.T20r2 = permute(T20r2,[2,3,1]);
    out.T30r2 = permute(T30r2,[2,3,1]);
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
            %             dat1 = [out.EDT(ch,:);out.T20(ch,:);out.T30(ch,:);out.C50(ch,:);...
            %                 out.C80(ch,:);out.D50(ch,:); ...
            %                 out.D80(ch,:);out.Ts(ch,:); ...
            %                 out.EDTr2(ch,:);out.T20r2(ch,:);out.T30r2(ch,:);...
            %                 out.T20T30ratio(ch,:)];
            %             cnames1 = num2cell(bandfc);
            %             rnames1 = {'Early decay time (s)',...
            %                 'Reverberation time T20 (s)',...
            %                 'Reverberation time T30 (s)',...
            %                 'Clarity index C50 (dB)',...
            %                 'Clarity index C80 (dB)',...
            %                 'Definition D50',...
            %                 'Definition D80',...
            %                 'Centre time Ts (s)',...
            %                 'Correlation coefficient EDT r^2',...
            %                 'Correlation coefficient T20 r^2',...
            %                 'Correlation coefficient T30 r^2',...
            %                 'Ratio of T20 to T30 %%'};
            dat1 = [out.EDT(ch,:);out.T20(ch,:);out.T30(ch,:);...
                out.EDTr2(ch,:);out.T20r2(ch,:);out.T30r2(ch,:);...
                out.T20T30ratio(ch,:);out.PNR(ch,:)];
            cnames1 = num2cell(bandfc);
            rnames1 = {'Early decay time (s)',...
                'Reverberation time T20 (s)',...
                'Reverberation time T30 (s)',...
                'Correlation coefficient EDT r^2',...
                'Correlation coefficient T20 r^2',...
                'Correlation coefficient T30 r^2',...
                'Ratio of T20 to T30 %%',...
                'Dynamic range: peak to final 10% (dB)'};
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
                
                figure('Name','Power-averaged Level Decay and Regression Lines')
                
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
                    if (noisecomp == 3) || (noisecomp == 4)
                        
                        plot(t(irstart(1,ch,band):dataend(1,ch,band),ch,band)...
                            ,10*log10((q(1,ch,band).* ...
                            exp(-abs(q(2,ch,band))...
                            .*t(irstart(1,ch,band):dataend(1,ch,band),ch,band)) + ...
                            (q(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
                            -t(irstart(1,ch,band):dataend(1,ch,band),ch,band))))), ...
                            'Color',[0.7 0 0],'DisplayName', 'Nonlinear Fit')
                        
                        plot(t(irstart(1,ch,band):dataend(1,ch,band),ch,band),...
                            10*log10(q(1,ch,band).* ...
                            exp(-abs(q(2,ch,band))...
                            .*t(irstart(1,ch,band):dataend(1,ch,band),ch,band))),...
                            'LineStyle',':','Color',[0 0 0.7],'DisplayName', 'Linear fit no B3')
                    end
                    if (noisecomp == 5)
                        
                        
                        plot(t(tstart(1,ch,band):t30end(1,ch,band),ch,band)...
                            ,10*log10((q(1,ch,band).* ...
                            exp(-abs(q(2,ch,band))...
                            .*t(tstart(1,ch,band):t30end(1,ch,band),ch,band)) + ...
                            (q(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
                            -t(tstart(1,ch,band):t30end(1,ch,band),ch,band))))), ...
                            'Color',[0 0 0.7],'DisplayName', 'T30')
                        
                        plot(t(tstart(1,ch,band):t30end(1,ch,band),ch,band),...
                            10*log10(q(1,ch,band).* ...
                            exp(-abs(q(2,ch,band))...
                            .*t(tstart(1,ch,band):t30end(1,ch,band),ch,band))),...
                            'LineStyle',':','Color',[0 0 0.7],'DisplayName', 'T30 no B3')
                        
                        plot(t(tstart(1,ch,band):t20end(1,ch,band),ch,band),...
                            10*log10((p(1,ch,band).* ...
                            exp(-abs(p(2,ch,band))...
                            .*t(tstart(1,ch,band):t20end(1,ch,band),ch,band)) + ...
                            (p(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
                            -t(tstart(1,ch,band):t20end(1,ch,band),ch,band))))), ...
                            'Color',[0 0.7 0],'DisplayName', 'T20')
                        
                        plot(t(tstart(1,ch,band):t20end(1,ch,band),ch,band),...
                            10*log10(p(1,ch,band).* ...
                            exp(-abs(p(2,ch,band))...
                            .*t(tstart(1,ch,band):t20end(1,ch,band),ch,band))),...
                            'LineStyle',':','Color',[0 0.7 0],'DisplayName', 'T20 no B3')
                        
                        plot(t(irstart(1,ch,band):edtend(1,ch,band),ch,band)...
                            ,10*log10((o(1,ch,band).* ...
                            exp(-abs(o(2,ch,band))...
                            .*t(irstart(1,ch,band):edtend(1,ch,band),ch,band)) + ...
                            (o(3,ch,band).*(t(dataend(1,ch,band),ch,band)...
                            -t(irstart(1,ch,band):edtend(1,ch,band),ch,band))))), ...
                            'Color',[0.7 0 0],'DisplayName', 'EDT')
                        
                        plot(t(irstart(1,ch,band):edtend(1,ch,band),ch,band),...
                            10*log10(o(1,ch,band).* ...
                            exp(-abs(o(2,ch,band))...
                            .*t(irstart(1,ch,band):edtend(1,ch,band),ch,band))),...
                            'LineStyle',':','Color',[0.7 0 0],'DisplayName', 'EDT no B3')
                        
                        
                    end
                    
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
                    
                    if (noisecomp == 0) || (noisecomp == 1) || (noisecomp == 2)|| (noisecomp == 5)
                        % text on subplots
                        text(levdecayend(1,ch,band)*0.45,-5,...
                            ['EDT ',num2str(0.01*round(100*EDT(1,ch,band)))],'Color',[0.9 0 0])
                        
                        text(levdecayend(1,ch,band)*0.45,-10,...
                            ['T20 ',num2str(0.01*round(100*T20(1,ch,band)))],'Color',[0 0.6 0])
                        
                        text(levdecayend(1,ch,band)*0.45,-15,...
                            ['T30 ',num2str(0.01*round(100*T30(1,ch,band)))],'Color',[0 0 0.6])
                    elseif (noisecomp == 3) || (noisecomp == 4)
                        text(levdecayend(1,ch,band)*0.45,-5,...
                            ['T ',num2str(0.01*round(100*EDT(1,ch,band)))],'Color',[0 0 0.6])
                    end
                    
                    
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
                if (noisecomp == 0) ||(noisecomp == 1) ||(noisecomp == 2) ||(noisecomp == 5)
                    plot([0.1,0.4], [0.6,0.6],'Color',[0.9 0 0],'DisplayName','EDT')
                    text(0.5,0.6,'EDT');
                    plot([0.1,0.4], [0.4,0.4],'Color',[0 0.6 0],'DisplayName','T20')
                    text(0.5,0.4,'T20');
                    plot([0.1,0.4], [0.2,0.2],'Color',[0 0 0.6],'DisplayName', 'T30')
                    text(0.5,0.2,'T30');
                elseif(noisecomp == 3) ||(noisecomp == 4)
                    plot([0.1,0.4], [0.6,0.6],'Color',[0.7 0 0],'DisplayName','Nlin')
                    text(0.5,0.6,'Non-lin');
                    plot([0.1,0.4], [0.4,0.4],'LineStyle',':','Color',[0 0 0.7],'DisplayName','Lin')
                    text(0.5,0.4,'Linear');
                end
                set(gca,'YTickLabel','',...
                    'YTick',zeros(1,0),...
                    'XTickLabel','',...
                    'XTick',zeros(1,0))
            end
            hold off
        end
        if isstruct(data)
            doresultleaf(levdecay,'Level [dB]',{'Time'},...
                'Time',      ((1:len)-1)./fs,  's',           true,...
                'channels',  chanID,           'categorical', [],...
                'Frequency', num2cell(bandfc), 'Hz',          false,...
                'name','Reverb_time');
        end
    end
    out.funcallback.name = 'ReverberationTime_IRmulti.m';
    out.funcallback.inarg = {fs,startthresh,bpo,doplot,filterstrength,phasemode,noisecomp,autotrunc,f_low,f_hi,averagingmethod,SNR,additionalaudio};
else
    out = [];
end
% eof

