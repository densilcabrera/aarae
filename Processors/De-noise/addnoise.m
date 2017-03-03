function OUT = addnoise(IN, method, signoiseratio,SNRtype,fexponent,fhigh,flow,windowfun,selection)
% This function adds noise to the input audio



    
if nargin == 1
    str = {'1. Add noise that has the same power spectrum as the signal';
        '2. Add steady pink, white or similar noise, using unweighted SNR';
        '3. Add steady pink, white or similar noise, using A-weighted SNR';
        '4. Add steady pink, white or similar noise, using SNR weighted by any available filter'};
    [method,ok] = listdlg('PromptString','Select the calculation method',...
        'SelectionMode','single',...
        'ListString',str,...
        'ListSize', [400,400]);
    if ~ok
        OUT = [];
        return
    end
end

if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
    [len,chans,bands,dim4,dim5,dim6] = size(audio);
end


switch method
    case 1
        if nargin < 3
            param = inputdlg({'Signal to noise ratio [dB]';...
                'SNR is RMS [0] or peak-to-noise [1] ratio'},...
                'Noise input parameters',1, ...
                {'10';'0'});
            param = str2num(char(param));
            if length(param) ~= 2
                OUT = [];
                return
            end
            if ~isempty(param)
                signoiseratio = param(1);
                SNRtype = param(2);
            end
        end
        X = fft(audio);
        if rem(len,2) == 0 % even number of spectrum components
            X = abs(X(2:len/2,:,:,:,:,:)).* exp(1i*rand(len/2-1,chans,bands,dim4,dim5,dim6).*2*pi);
            X = ifft([zeros(1,chans,bands,dim4,dim5,dim6);...
                X; zeros(1,chans,bands,dim4,dim5,dim6); flip(conj(X))]);
        else % odd number of spectrum components
            X = abs(X(2:ceil(len/2),:,:,:,:,:)).* exp(1i*rand(floor(len/2),chans,bands,dim4,dim5,dim6).*2*pi);
            X = ifft([zeros(1,chans,bands,dim4,dim5,dim6);...
                X; flip(conj(X))]);
        end
        if SNRtype == 0
            audio = audio + X ./ db2mag(signoiseratio);
        else
            SNR0 = max(abs(audio))./rms(X);
            audio = audio + X.* repmat(SNR0./db2mag(signoiseratio) , [len,1,1,1,1,1]);
        end
                    
        
    case {2,3,4}
        if nargin < 6
            param = inputdlg({'Signal to noise ratio [dB]';...
                'SNR is RMS [0] or peak-to-noise [1] ratio';...
                'Spectral slope [dB/octave]';...
                'High cutoff frequency [Hz]';...
                'Low cutoff frequency [Hz]';...
                'Tukey window ratio [0:1]'},...
                'Noise input parameters',1, ...
                {'10';'0';'0';num2str(0.98*fs/2);'20';'0'});
            param = str2num(char(param));
            if length(param) ~= 6
                OUT = [];
                return
            end
            if ~isempty(param)
                signoiseratio = param(1);
                SNRtype = param(2);
                fexponent = (param(3)-3)/3;
                fhigh = param(4);
                flow = param(5);
                windowfun = param(6);
            end
        end
        
        display = 0;
        X = audio;
        for d3 = 1:bands
            for d4 = 1:dim4
                for d5 = 1:dim5
                    for d6 = 1:dim6
                        Xnoise = noise(fexponent, len/fs, fs,...
                            flow, fhigh, chans, windowfun, display);
                        if size(Xnoise.audio,1) == len
                            X(:,:,d3,d4,d5,d6) = Xnoise.audio;
                        elseif size(Xnoise.audio,1) < len
                            X(1:size(Xnoise.audio,1),:,d3,d4,d5,d6) =Xnoise.audio;
                        else
                            X(:,:,d3,d4,d5,d6) =Xnoise.audio(1:len,:);
                        end
                    end
                end
            end
        end
        if method == 2
            if SNRtype == 0
                SNR0 = rms(audio)./rms(X);
            else
                SNR0 = max(abs(audio))./rms(X);
            end
                        
        elseif method == 3
            if SNRtype == 0
                SNR0 = rms(Aweight(audio,fs))./rms(Aweight(X,fs));
            else
                SNR0 = max(abs(audio))./rms(X);
            end
                        
        else

            weights = what([cd '/Processors/Filters']);
            if ~exist('selection','var')
                if ~isempty(weights.m)
                    [selection,ok] = listdlg('ListString',cellstr(weights.m),'SelectionMode','single');
                else
                    warndlg('No weighting filters found!','AARAE info')
                    return
                end
                if ok ~= 1, return; end
            end
            [~,funname] = fileparts(weights.m{selection,1});
            S = feval(funname,IN);
            N.audio = X;
            N.fs = fs;
            N = feval(funname,N);
            if SNRtype == 0
                SNR0 = rms(S.audio)./rms(N.audio);
            else
                SNR0 = max(abs(S.audio))./rms(N.audio);
            end
        end
        audio = audio + X.* repmat(SNR0./db2mag(signoiseratio) , [len,1,1,1,1,1]);
end




    
if isstruct(IN)
    OUT = IN; 
    OUT.audio = audio; 
    OUT.funcallback.name = 'addnoise.m'; 
    switch method
        case 1
            OUT.funcallback.inarg = {method, signoiseratio,SNRtype};
        case {2, 3}
            OUT.funcallback.inarg = {method, signoiseratio,SNRtype,fexponent,fhigh,flow,windowfun};
        otherwise
            OUT.funcallback.inarg = {method, signoiseratio,SNRtype,fexponent,fhigh,flow,windowfun,selection};
    end
else
    OUT = audio;
end


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