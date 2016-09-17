function OUT = bandpass(IN, flo, fhi, order, fs, phase)
% This function applies a zero phase bandpass filter in the frequency
% domain to the input audio.
%
% The high and low frequency cutoff frequencies can be specified, which
% should be in the range between 0 Hz and the Nyquist frequency, and the
% high cutoff frequency should be greater than the low cutoff frequency.
%
% The filter response at cutoff frequency is -3 dB.
%
% The filter pseudo-order is used to shape the magnitude response - such
% that the filter skirt slopes are 6n per octave (n being the filter
% pseudo-order). The default filter order of 24 yields a slope of 144
% dB/octave. However, if the band limits are close to the 0 Hz or Nyquist
% frequencies (and the audio input is brief), then this slope might not be
% realized in practice.
%
% The 0 Hz component is always zeroed by this function.
%
% Code by Densil Cabrera
% version 1.00 (19 August 2014)

if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
elseif  nargin > 1
    audio = IN;
end

if nargin ==1
    param = inputdlg({'High frequency limit (Hz)';...
        'Low frequency limit (Hz)';...
        'Filter pseudo-order in-band';...
        'Filter pseudo-order out-of-band';...
        'Maximum phase [-1], Zero phase [0] or Mimimum phase [1]'},...
        'Band Limit Settings',...
        [1 60],...
        {num2str(fs/2);'0';'24';'24';'0'});
    
    param = str2num(char(param));
    
    if length(param) < 5, param = []; end
    
    if ~isempty(param)
        fhi = param(1);
        flo = param(2);
        order = [param(3) param(4)];
        phase = param(5);
    else
        OUT = [];
        return
    end
end

if ~exist('phase','var')
    phase = 0;
end

if fhi > flo ...
        && fhi <= fs/2 ...
        && flo >= 0
    ok = 1;
else
    ok = 0;
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(order) && ok ==1
    [len,chans,bands,dim4,dim5,dim6] = size(audio);
    
    if length(order) == 1, order = [order order]; end
%     if mod(len,2) == 1
%         % use an even length signal
%         audio = [audio;ones(1,chans,bands,dim4,dim5,dim6)];       
%     end
    fftlen = 2*len; % note that fftlen must be even
    
    % list of fft component frequencies
    f = ((1:fftlen)'-1) * fs / fftlen;
    
    % index of low cut-off
    indlo = find(abs(f(1:end/2)-flo) == min(abs(f(1:end/2)-flo)),1,'first');
    
    % index of high cut-off
    indhi = find(abs(f(1:end/2)-fhi) == min(abs(f(1:end/2)-fhi)),1,'first');
    
    % centre frequency index
    fc = (flo * fhi).^0.5; % geometric mean frequency
    indfc = find(abs(f(1:end/2)-fc) ...
        == min(abs(f(1:end/2)-fc)),1,'first');
    
    
    % magnitude envelope
    mag = zeros(fftlen,1); % preallocate and set DC to 0
    
    % below centre frequency
    % out-of-band, using exact 6n dB/oct skirts
    mag(2:indlo-1) = ...
        (f(2:indlo-1)./ flo ).^(order(2)) ./2.^0.5;
    
    % in-band
    mag(indlo:indfc-1) = ...
        (1 ./ (1 + (f(indlo:indfc-1)./ flo ).^(-2*order(1)))).^0.5;
    
    % from centre frequency to Nyquist frequency
    % in-band
    mag(indfc:indhi) = ...
        (1 ./ (1 + (f(indfc:indhi)./ fhi ).^(2*order(1)))).^0.5;
    
    % out-of-band, using exact 6n dB/oct skirts
    mag(indhi+1:fftlen/2+1) = ...
        (f(indhi+1:fftlen/2+1) ./ fhi).^(-order(2)) ./ 2.^0.5;
    
    % normalize gain to 0 dB at fc
    % note that if flo and fhi are too close together (or the in-band
    % order is too low) then this will change the gain at the cutoff
    % frequencies
    mag = mag ./ mag(indfc);
    
    % above Nyquist frequency
    mag(fftlen/2+2:end) = flipud(mag(2:fftlen/2));
    
    % phase
    switch phase
        case 1
            mag = minphasefreqdomain(mag);
        case -1
            mag = conj(minphasefreqdomain(mag));
    end
    
    
    % apply the filter
    audio = ifft(fft(audio,fftlen).* repmat(mag,[1,chans,bands,dim4,dim5,dim6]));
    audio = audio(1:len,:,:,:,:,:);
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
        OUT.funcallback.name = 'bandpass.m';
        OUT.funcallback.inarg = {flo,fhi,order,fs,phase};
    else
        OUT = audio;
    end
    
else
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2014, Densil Cabrera
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