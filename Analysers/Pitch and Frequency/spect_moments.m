function [OUT, varargout] = spect_moments(IN,spectrogramwinlen,meandims,fscale,yscale,fs)
% SPECT_MOMENTS calculates the statistical spectral moments from a signal's
% spectrum.
% m = spect_moments(x, fs) calculates the first 4 spectral moments:
% centroid / bandwidth (variance) / skewness / kurtosis
% FFT is used as the method of calculation
% Inputs:
% x - one/two channels signal
% fs - sampling frequency
% Outputs:
% 1st moment - centroid | derived from the frequency-weighted mean of the
% critical band distribution, associated with ?brightness?.
% 2nd moment - bandwidth | derived from the dispersion centred at the
% spectral centroid, how wide or narrow the spectrum is.
% 3rd moment - skewness | a measure of asymmetry in the distribution.
% 4th moment - kurtosis | a measure of the width of the peaks in the
% distribution (?peakedness?).

% *************************************************************************
% GET INPUT
if isstruct(IN)
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    
    
    if isfield(IN,'chanID') % Get the channel ID if it exists
        chanID = IN.chanID;
    end
else
    audio = IN;
end




% *************************************************************************
% DIALOG BOX
if nargin ==1
    % defaults
    spectrogramwinlen = 0;
    meandims = 2;
    fscale = 1;
    yscale = 1;
    param = inputdlg({'Time-spectrum window length in samples (use 0 to not do a spectrogram)';...
        'Dimensions to average';...
        'Frequency scale: linear [1],logarithmic [2],critical band [3]';...
        'Magnitude power: e.g. magnitude [1], squared magnitude [2], Stevens power law [0.6]'},...
        'Spectrum Moments Settings',... % This is the dialog window title.
        [1 60],...
        {num2str(spectrogramwinlen);...
        num2str(meandims);...
        num2str(fscale);...
        num2str(yscale)}); % default settings
    if length(param) < 4, param = []; end
    if ~isempty(param)
        spectrogramwinlen = str2num(char(param{1}));
        meandims = str2num(char(param{2}));
        fscale = str2num(char(param{3}));
        yscale = str2num(char(param{4}));
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end



% *************************************************************************
% PROCESS DATA
if ~isempty(audio) && ~isempty(fs) && ~isempty(spectrogramwinlen) && ~isempty(fscale) && ~isempty(yscale)
    
    % the following is unnecessary in aarae
    % if length(audio) ~= size(audio,1)
    %     audio = audio';
    % end
    
    [len,chans,bands,dim4,dim5,dim6] = size(audio);
    
    NFFT = 2 ^ nextpow2(len);
    
    
    MXtmp = abs(fft(audio,NFFT));
    NumUniquePts = ceil((NFFT+1)/2);
    MXtmp = MXtmp(1:end/2,:,:,:,:,:);
    
    % dimension mean
    for dim = 2:6
        if ~isempty(find(meandims==dim, 1))
            MXtmp = (mean(MXtmp.^2,dim)).^0.5; % rms
        end
    end
    [~,chans,bands,dim4,dim5,dim6] = size(MXtmp);
    
    MXtmp = MXtmp / length(audio);
    MXtmp = MXtmp .^ yscale; 
    MXtmp = MXtmp * 2;  % Account for throwing out second half of FFTX above
    MXtmp(1) = MXtmp(1) / 2;                            % Account for DC uniqueness
    if ~rem(NFFT,2)
        MXtmp(length(MXtmp)) = MXtmp(length(MXtmp)) / 2;  % Account for Nyquist uniqueness
    end
    
    f_n = (1:NumUniquePts-1)' * fs / NFFT; % not sure that this is right (what about 0 Hz?)

    
    [centroid,bandwidth,skewness,kurtosis] = calculatemoments(MXtmp,f_n);
    
    
    
    % *************************************************************************
    % NEED TO CHANGE THIS TO BE ABLE TO DEAL WITH MULTIDIMENSIONAL DATA.
    % This can be done by using RowName to identify the data coordinate
    fig1 = figure('Name','My results table');
    table1 = uitable('Data',[centroid bandwidth skewness kurtosis],...
        'ColumnName',{'Centroid','Bandwidth','Skewness','Kurtosis'},...
        'RowName',{'Results'});
    [~,tables] = disptables(fig1,table1);
    
    
    % *************************************************************************
    % TIME-SPECTRUM ANALYSIS
    if spectrogramwinlen > 0
        %THIS STILL NEEDS TO BE WRITTEN
        % use spectrogram to get the data
        % find spectral moments for each frame
        % generate AARAE results leaf (big data presumably)
        % generate chart of parameters over time
    end
    
    
    % *************************************************************************
    if isstruct(IN)
        OUT.tables = tables;
        OUT.funcallback.name = 'spect_moments.m';
        OUT.funcallback.inarg = {spectrogramwinlen,meandims,fscale,yscale}; %fs is not needed in callback
    else
        OUT = [];
        varargout{1} = centroid;
        varargout{2} = bandwidth;
        varargout{3} = skewness;
        varargout{4} = kurtosis;
    end
else
    OUT = []; % return empty output if any required input apart from meandims is empty
end
end %eof

function [centroid,bandwidth,skewness,kurtosis] = calculatemoments(y,f)
% this local function could be made non-local if that was useful
centroid = sum(f .* y) ./ sum(y);
bandwidth = sum((f - centroid).^2 .* y) ./ sum(y);
bandwidth = sqrt(bandwidth);
skewness = sum((f - centroid).^3 .* y) ./ sum(y);
skewness = skewness / bandwidth.^3;
kurtosis = sum((f - centroid).^4 .* y) ./ sum(y);
kurtosis = kurtosis /  bandwidth.^4;
end

%**************************************************************************
% Copyright (c) 2014, Ella Manor ella.manor@sydney.edu.au
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