function OUT = narrowbandfilterbank(IN,fs,filtdur,flo,fhi,componentsperband,orderout,cutoffattenuation,phasemode,pad)
% This filterbank is based on an fft and ifft of the input signal. The
% filters are specified in relation to the spectrum components of the fft.
%
% Note that it is possible for this filterbank to return a very large
% number of bands, so this function should be used with caution!
%
% Although filtering is done in the frequency domain, this function avoids
% temporal aliasing (from spectral wrap-around) by zeropadding in the time
% domain. This is inefficient from a computational perspective, but it
% provides dramatic improvements in the time response of the minimum and
% maximum phase filters.

if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
else
    audio = IN;
end

[len,chans,bands,dim4,dim5,dim6] = size(audio);

% if the input is already multiband, then mixdown first
if bands > 1
    audio = sum(audio,3);
    bands = 1;
    disp('Multi-band audio has been mixed-down prior to octave band filtering')
end

if nargin < 3
    param = inputdlg({'FIR filter duration (s), or use 0 to make it equal to the input wave duration';...
        'Lowest centre frequency (Hz)';...
        'Highest centre frequency (Hz)';...
        'Number of spectrum components within each filter bandwidth (odd number >= 3, includes cutoff components)';...
        'Filter out-of-band order';...
        'Filter attenuation at cutoff component frequencies (e.g., 0, 3 or 6 dB)';...
        'Linear phase [0], Minimum phase [1] or Maximum phase [-1]';...
        'Remove zero-padding [0] or retain zero-padding [1]'},...
        'Filterbank settings',...
        [1 60],...
        {'1';'20';'200';'3';'12';'3';'0';'0'});
    if isempty(param)
        OUT = [];
        return
    end
    param = str2num(char(param));
    filtdur = param(1);
    flo = param(2);
    fhi = param(3);
    componentsperband = param(4);
    orderout = param(5);
    cutoffattenuation = param(6);
    phasemode = param(7);
    pad = param(8);
else
    if ~exist('fs','var'),fs = 48000; end
    if ~exist('filtdur','var'),filtdur = 0; end
    if ~exist('flo','var'),flo = 20; end
    if ~exist('fhi','var'),fhi = 200; end
    if ~exist('componentsperband','var'),componentsperband = 3; end
    if ~exist('orderout','var'),orderout = 12; end
    if ~exist('cutoffattenuation','var'),cutoffattenuation = 3; end
    if ~exist('phasemode','var'),phasemode = 0; end
end

if filtdur == 0, filtdur = len/fs; end
% if phasemode ~= 0 || phasemode ~=10
%     filtdur = 2*filtdur;
% end
filterlen = 2*ceil(filtdur*fs/2); % even length fft & filter
f = fs*((1:filterlen)'-1)./filterlen; % list of frequencies
indlo = find(abs(f(1:end/2)-flo) == min(abs(f(1:end/2)-flo)),1,'first');
indhi = find(abs(f(1:end/2)-fhi) == min(abs(f(1:end/2)-fhi)),1,'first');

if componentsperband < 3, componentsperband = 3; end
if rem(componentsperband,2) == 0, componentsperband=componentsperband+1; end

len = 2*ceil((len+filterlen)/2); % even total length
% FFT with zero padding to avoid circular convolution
switch phasemode
    case {1, 11}
        % zero pad after for min phase
        audio = fft(audio,len);
    case {-1,-11}
        % zero pad before for max phase
        audio = fft([zeros(filterlen,chans,bands,dim4,dim5,dim6);audio],len);
    case {0,10}
        % symetric zero pad for linear phase
        audio = fft([zeros(filterlen/2,chans,bands,dim4,dim5,dim6);audio],len);
end

findex = indlo:componentsperband-2:indhi;
bands = length(findex);
bandfiltered = zeros(len,chans,bands,dim4,dim6,dim6);
cutoffmag = db2mag(-cutoffattenuation);
for b = 1:bands
    % magnitude envelope
    mag = zeros(filterlen,1);
    lowcutoff = findex(b) - floor(componentsperband/2);
    hicutoff = findex(b) + floor(componentsperband/2);
    mag(lowcutoff+1:hicutoff-1) = 1; % 0 dB
    [mag(lowcutoff), mag(hicutoff)] = deal(cutoffmag);
    mag(2:lowcutoff-1) = ...
        (f(2:lowcutoff-1)./ f(lowcutoff) ).^(orderout) .* cutoffmag;
    mag(hicutoff+1:filterlen/2+1) = ...
        (f(hicutoff+1:filterlen/2+1) ./ f(hicutoff)).^(-orderout) .* cutoffmag;
    mag(filterlen/2+2:end) = flipud(mag(2:filterlen/2));
    
    switch phasemode
        case 0
            bandfiltered(:,:,b,:,:,:) = real(ifft(repmat(...
                fft([zeros(round((len-filterlen)/2),1);ifft(mag)],len),...
               [1,chans,1,dim4,dim5,dim6]) .* audio));
        case 10
            mag(filterlen/2:end) = 0;
            bandfiltered(:,:,b,:,:,:) = ifft(repmat(...
                fft([zeros(round((len-filterlen)/2),1);ifft(mag)],len),...
               [1,chans,1,dim4,dim5,dim6]) .* audio);
        case 1
            mag = minphasefreqdomain(mag,120,1);
            bandfiltered(:,:,b,:,:,:) = real(ifft(repmat(...
                fft([ifft(mag);zeros(len-filterlen,1)],len),...
                [1,chans,1,dim4,dim5,dim6]) .* audio));
        case 11
            mag = minphasefreqdomain(mag,120,1);
            mag(filterlen/2:end) = 0;
            bandfiltered(:,:,b,:,:,:) = ifft(repmat(...
                fft([zeros(len-filterlen,1);ifft(mag)],len),...
                [1,chans,1,dim4,dim5,dim6]) .* audio);
        case -1
            mag = conj(minphasefreqdomain(mag,120,1));
            bandfiltered(:,:,b,:,:,:) = real(ifft(repmat(...
                fft([zeros(len-filterlen,1);ifft(mag)],len),...
                [1,chans,1,dim4,dim5,dim6]) .* audio));
        case -11
            mag = conj(minphasefreqdomain(mag,120,1));
            mag(filterlen/2:end) = 0;
            bandfiltered(:,:,b,:,:,:) = ifft(repmat(...
                fft([zeros(len-filterlen,1);ifft(mag)],len),...
                [1,chans,1,dim4,dim5,dim6]) .* audio);
    end
end

% crop the ends (remove zeropadding)
if pad ~= 1
    switch phasemode
        case {1, 11}
            bandfiltered = bandfiltered(1:end-filterlen,:,:,:,:,:);
        case {-1,-11}
            bandfiltered = bandfiltered(filterlen+1:end,:,:,:,:,:);
        case {0,10}
            bandfiltered = bandfiltered(filterlen/2+1:end-filterlen/2,:,:,:,:,:);
    end
end

clear audio

if isstruct(IN)
    OUT = IN;
    OUT.audio = bandfiltered;
    OUT.bandID = f(findex);
    OUT.funcallback.name = 'narrowbandfilterbank.m';
    OUT.funcallback.inarg = {fs,filtdur,flo,fhi,componentsperband,orderout,cutoffattenuation,phasemode,pad};
else
    OUT = bandfiltered;
end

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