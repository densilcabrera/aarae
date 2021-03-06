function [OUT, varargout] = Golay_process(IN, fs, audio2, dostack)
% This function is used to analyse signals that were recorded using the
% Golay generator within AARAE. The output is an implulse response.
%
% code by Densil Cabrera
% version 1.01 (1 August 2014)

if nargin < 4, dostack = true; end
% dostack=false is used to bypass the stacking of multicycle measurements
% and complementary signal subtraction when this function is used by
% AARAE's '*' button (convolve audio with audio2) because the stacking is
% done within that function in various ways.

try
    if ~isstruct(IN) && nargin >=3  
        ok = true;
    elseif isfield(IN.properties,'GolayN')
        ok = true;
    elseif isfield(IN,'funcallback') && strcmp(IN.funcallback.name,'Golay.m')
        ok = true;
    else
        ok = false;
    end
catch
    ok = false;
end

if ok
    if isstruct(IN)
        audio = IN.audio;
        fs = IN.fs;
        if isfield(IN,'audio2')
            audio2 = IN.audio2;
        else
            disp('The original Golay test signal is required to be in audio2 - not found')
            OUT = [];
            return
        end
    else
        audio = IN;
    end
     
    if ~isempty(audio) && ~isempty(fs)
        % find the relevent indices from audio2 - the original test signal
        lasta = find(audio2 == 0,1,'first')-1;
        firstb = find(audio2(lasta+1:end) ~= 0,1,'first')+lasta;
        lastb = length(audio2);

        [len,chans,bands,dim4,dim5,dim6] = size(audio);
        
        if len < lastb
            disp('Recorded audio is too short for Golay processing')
            disp('It must be at least the same duration as the test signal')
            OUT = [];
            return
        end
        
        if ~isempty(lasta) && ~isempty(firstb)
            a = audio2(1:lasta);
            b = audio2(firstb:end);
        else
            disp('Audio 2 is not in the required format')
            disp('This function should be used with test signals that were generated by the ''Golay'' generator in AARAE')
            OUT = [];
            return
        end
        
        if length(a) ~= length(b)
            disp('Audio 2 is not in the required format')
            disp('This function should be used with test signals that were generated by the ''Golay'' generator in AARAE')
            OUT = [];
            return
        end
        
        if dostack
            % average phase-complementary pairs of signals (if they exist)
            complementarysignals = false;
            if isfield(IN,'properties')
                if isfield(IN.properties,'complementarysignals')
                    if IN.properties.complementarysignals && isfield(IN.properties,'complementarysignalsoffset')
                        complementarysignals = true;
                        if isfield(IN.properties,'startflag')
                            startflag = IN.properties.startflag;
                        else
                            startflag = 1;
                        end
                        for d=1:length(startflag)
                            startindex1 = startflag(d);
                            endindex1 = startindex1 + IN.properties.complementarysignalsoffset-1;
                            startindex2 = startflag(d) + IN.properties.complementarysignalsoffset;
                            endindex2 = startindex2 + IN.properties.complementarysignalsoffset-1;
                            audio(startindex1:endindex1,:,:,:,:,:) = ...
                                mean(cat(7,audio(startindex1:endindex1,:,:,:,:,:),...
                                -audio(startindex2:endindex2,:,:,:,:,:)),7);
                        end
                    end
                end
            end
            
            % Stack IRs in dimension 4 if AARAE's multi-cycle mode was used
            if isfield(IN,'properties')
                if isfield(IN.properties,'startflag') && dim4==1
                    startflag = IN.properties.startflag;
                    dim4 = length(startflag);
                    if complementarysignals
                        len2 = IN.properties.complementarysignalsoffset;
                    else
                        len2 = startflag(2)-startflag(1);
                    end
                    audiotemp = zeros(len2,chans,bands,dim4);
                    for d=1:dim4
                        audiotemp(:,:,:,d,:,:) = ...
                            audio(startflag(d):startflag(d)+len2-1,:,:,1,:,:);
                    end
                end
            end
            if exist('audiotemp','var')
                audio = audiotemp;
            end
        end
        
        [~,chans,bands,dim4,dim5,dim6] = size(audio);
        
        aa = audio(1:lasta,:,:,:,:,:);
        bb = audio(firstb:lastb,:,:,:,:,:);
        
        
        % cross-spectrum, sum and scale
        ir = ifft(repmat(conj(fft(a)),[1,chans,bands,dim4,dim5,dim6]) .* fft(aa) ...
            + repmat(conj(fft(b)),[1,chans,bands,dim4,dim5,dim6]) .* fft(bb)) ...
            ./ (2*lasta);
        
        
        if isstruct(IN)
            OUT = IN;
            OUT.audio = ir;
            OUT.funcallback.name = 'Golay_process.m';
            OUT.funcallback.inarg = {};
        else
            OUT = ir;
        end
        varargout{1} = fs;
    else
        OUT = [];
    end
else
    OUT = [];
    warndlg('The audio to be processed was not recorded from a signal generated using the Golay generator within AARAE.','AARAE info')
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