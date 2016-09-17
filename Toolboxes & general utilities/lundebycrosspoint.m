function [crosspoint, Tlate, ok, prelimcrosspoint] = lundebycrosspoint(IR2, fs, fc)
% This function implements the crosspoint finding algorithm
% described in:
% A. Lundeby, T.E. Vigran, H. Bietz and M. Vorlaender, "Uncertainties of
% Measurements in Room Acoustics", Acustica 81, 344-355 (1995)
%
% INPUTS:
% IR2 is a squared IR - can be multi channel and multiband
%
% fs is sampling rate in Hz
%
% fc (optional) is the band centre frequencies
%
% OUTPUT:
% crosspoint is an index corresponding to the IR crosspoint. Each band &
% channel has its own crosspoint (chan in dim 2, band in dim 3).
%
% Tlate is the late reverberation time.
%
% ok indicates whether the algorithm was fully executed (for each
% band/channel etc) - in some cases it may not be possible to reasonably
% execute it (e.g. due to poor signal-to-noise ratio).
%
% prelimcrosspoint is the initial preliminary crosspoint index (rather than
% the final crosspoint after refinement).
%
% The main difference in this implementation is that the windows used for
% smoothing the squared decay are implemented with a fftfilt (and so do not
% involve any downsampling). Presumably this should slightly improve the
% curve fitting.
%
%
% Code by Densil Cabrera
% version 0.04 (Dec 2015)


%**************************************************************************
% Copyright (c) 2013, Densil Cabrera
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


[len,chans,bands,dim4,dim5,dim6] = size(IR2);

% clear fc if it does not match bands or is just a count of bands
if exist('fc','var')
    if length(fc)>1
        if length(fc) ~= bands
            clear fc;
        elseif fc(2)-fc(1)==1
            clear fc;
        end
    elseif fc == 1
        clear fc;
    end
end


Tlate = zeros(1,chans,bands,dim4,dim5,dim6);
ok = true(1,chans,bands,dim4,dim5,dim6);
[crosspoint,prelimcrosspoint] = deal(round(0.9*len)*ones(1,chans,bands,dim4,dim5,dim6));

% The following nested for-loops are a quick and dirty way of making this
% function compatible with up to 6 dimensional audio input. (Originally the
% function was written just for 3 dimensional audio input.)
for d4 = 1:dim4
    for d5 = 1:dim5
        for d6 = 1:dim6
            IR2a = IR2(:,:,:,d4,d5,d6);
            
            
            % window length for step 1
            if ~exist('fc','var') && bands > 1
                winlen = round(0.001*fs*linspace(50,10,bands));
            elseif exist('fc','var') && bands > 1
                winlen = zeros(1,bands);
                fc_low = fc<=125;
                winlen(fc_low) = round(0.001*fs*50);
                fc_hi = fc>=8000;
                winlen(fc_hi) = round(0.001*fs*10);
                fc_mid = ~(fc_low | fc_hi);
                if sum(fc_mid) >1
                    winlen(fc_mid) = round(0.001*fs*linspace(50,10,sum(fc_mid)));
                elseif sum(fc_mid) ==1
                    winlen(fc_mid) = round(0.001*fs*25);
                end
            else
                winlen(1:bands) = round(0.001*fs*25);
            end
            winlen = repmat(permute(winlen,[1,3,2]),[1,chans,1]);
            
            
            % 1. AVERAGE SQUARED IR IN LOCAL TIME INTERVALS
            IR2smooth = IR2a; %just for preallocation
            for b = 1:bands
                IR2smooth(:,:,b) = fftfilt(ones(winlen(1,1,b),1)./winlen(1,1,b),IR2a(:,:,b));
            end
            IR2smooth(IR2smooth<=0)=1e-300;
            IR2smoothdB =  10*log10(IR2smooth);
            maxIR2smoothdB = max(IR2smoothdB);
            maxind = ones(1,chans,bands);
            for ch = 1:chans
                for b = 1:bands
                    maxind(1,ch,b) = find(IR2smoothdB(:,ch,b) == maxIR2smoothdB(1,ch,b),1,'first');
                    IR2smoothdB(:,ch,b) = IR2smoothdB(:,ch,b) - maxIR2smoothdB(1,ch,b);
                end
            end
            
            
            
            % 2. ESTIMATE BACKGROUND NOISE LEVEL USING THE TAIL
            IR2tail = mean(IR2a(round(0.9*end):end,:,:));
            IR2taildB = 10*log10(IR2tail) - maxIR2smoothdB;
            
            
            
            % 3. ESTIMATE SLOPE OF DECAY FROM 0 dB TO NOISE LEVEL
            o = zeros(2,chans, bands);
            [crosspoint1,tend] = deal(len*ones(1,chans,bands));
            for ch = 1:chans
                for b = 1:bands
                    
                    index = find(IR2smoothdB(maxind(1,ch,b):end,ch,b) ...
                        <= IR2taildB(1,ch,b), 1, 'first')+maxind(1,ch,b)-1;
                    if ~isempty(index)
                        tend(1,ch,b) = index;
                    else
                        ok(1,ch,b,d4,d5,d6) = false;
                    end
                    
                    o(:,ch,b) = polyfit((maxind(1,ch,b):tend(1,ch,b))', ...
                        IR2smoothdB(maxind(1,ch,b):tend(1,ch,b),ch,b),1)';
                    
                    % 4. FIND PRELIMINARY crosspoint1
                    crosspoint1(1,ch,b) = ...
                        -round((o(2,ch,b)-IR2taildB(1,ch,b))/o(1,ch,b));
                    
                end
            end
            
            for ch=1:chans
                for b = 1:bands
                    if crosspoint1(1,ch,b) < round(maxind(1,ch,b)+0.05*fs)
                        crosspoint1(1,ch,b) = round(maxind(1,ch,b)+0.05*fs);
                    end
                end
            end
            crosspoint1(crosspoint1>round(0.9*len)) = round(0.9*len);
            prelimcrosspoint(1,:,:,d4,d5,d6) = crosspoint1;
            
            % 5. FIND NEW LOCAL TIME INTERVAL LENGTH
            slopelength = crosspoint1 - maxind;
            
            winlen = zeros(1,chans,bands);
            for ch = 1:chans
                if ~exist('fc','var') && bands > 1
                    winlen(1,ch,:) = round(slopelength(1,ch,:) ./ (-IR2taildB(1,ch,:) ./permute(linspace(3,10,bands),[1,3,2])));
                elseif exist('fc','var') && bands > 1
                    winlen(1,ch,fc_low) = round(slopelength(1,ch,fc_low) ./ (-IR2taildB(1,ch,fc_low) ./3));
                    winlen(1,ch,fc_hi) = round(slopelength(1,ch,fc_hi) ./ (-IR2taildB(1,ch,fc_hi) ./10));
                    if sum(fc_mid) >1
                        winlen(1,ch,fc_mid) = round(slopelength(1,ch,fc_mid) ./ (-IR2taildB(1,ch,fc_mid) ./permute(linspace(3,10,sum(fc_mid)),[1,3,2])));
                    elseif sum(fc_mid) ==1
                        winlen(1,ch,fc_mid) = round(slopelength(1,ch,fc_mid) ./ (-IR2taildB(1,ch,fc_mid) ./6));
                    end
                else
                    winlen(1,ch,:) = round(slopelength(1,ch,:) ./ (-IR2taildB(1,ch,:) ./6));
                end
            end
            
            
            
            % Fix out-of-range winlen
            winlen(isinf(winlen)) = round(0.001*fs*25); %25 ms
            winlen(isnan(winlen)) = round(0.001*fs*25);
            winlen(winlen>0.2*len) = round(0.2*len);
            winlen(winlen<10) = 10;
            
            
            
            % 6. AVERAGE SQUARED IR IN NEW LOCAL TIME INTERVALS
            for ch = 1:chans
                for b = 1:bands
                    IR2smooth(:,ch,b) = fftfilt(ones(winlen(1,ch,b),1)./winlen(1,ch,b),IR2a(:,ch,b));
                end
            end
            IR2smooth(IR2smooth<=0)=1e-300;
            IR2smoothdB =  10*log10(IR2smooth);
            maxIR2smoothdB = max(IR2smoothdB);
            maxind = ones(1,chans,bands);
            for ch = 1:chans
                for b = 1:bands
                    index = find(IR2smoothdB(:,ch,b) == maxIR2smoothdB(1,ch,b),1,'first');
                    if ~isempty(index)
                        maxind(1,ch,b) = index;
                    else
                        ok(1,ch,b,d4,d5,d6) = false;
                    end
                    IR2smoothdB(:,ch,b) = IR2smoothdB(:,ch,b) - maxIR2smoothdB(1,ch,b);
                end
            end
            
            
            
            
            % ITERATE STEPS 7-9
            for iter = 1:5
                
                
                % 7. ESTIMATE THE BACKGROUND NOISE LEVEL
                noisefloorindex = zeros(1,chans,bands);
                for ch = 1:chans
                    for b = 1:bands
                        if ok(1,ch,b,d4,d5,d6)
                            try
                            noisefloorindex(1,ch,b) = ...
                                round((o(2,ch,b)+IR2taildB(1,ch,b)-10)/o(1,ch,b));
                            catch
                                noisefloorindex(1,ch,b) = round(0.9*len);
                            end
                        end
                    end
                end
                noisefloorindex(noisefloorindex > round(0.9*len)) = round(0.9*len);
                %noisefloorindex(noisefloorindex < round(0.05*fs)) = round(0.05*fs);
                
                for ch = 1:chans
                    for b = 1:bands
                        if ok(1,ch,b,d4,d5,d6)
                            try
                                IR2tail(1,ch,b) = mean(IR2a(noisefloorindex(1,ch,b):end,ch,b));
                                IR2taildB(1,ch,b) = 10*log10(IR2tail(1,ch,b))- maxIR2smoothdB(1,ch,b);
                            catch
                                ok(1,ch,b,d4,d5,d6) = false;
                            end
                        end
                    end
                end
                
                
                
                % 8. ESTIMATE THE LATE DECAY SLOPE
                for ch = 1:chans
                    for b = 1:bands
                        if ok(1,ch,b,d4,d5,d6)
                            if IR2taildB(1,ch,b) < -35
                                LateSlopeEnddB = IR2taildB(1,ch,b) + 10;
                                LateSlopeStartdB = IR2taildB(1,ch,b) + 30;
                            elseif IR2taildB(1,ch,b) < -30
                                LateSlopeEnddB = IR2taildB(1,ch,b) + 8.5;
                                LateSlopeStartdB = IR2taildB(1,ch,b) + 25;
                            elseif IR2taildB(1,ch,b) < -25
                                LateSlopeEnddB = IR2taildB(1,ch,b) + 7;
                                LateSlopeStartdB = IR2taildB(1,ch,b) + 20;
                            elseif IR2taildB(1,ch,b) < -20
                                LateSlopeEnddB = IR2taildB(1,ch,b) + 6;
                                LateSlopeStartdB = IR2taildB(1,ch,b) + 16;
                            elseif IR2taildB(1,ch,b) < -15
                                LateSlopeEnddB = IR2taildB(1,ch,b) + 5;
                                LateSlopeStartdB = IR2taildB(1,ch,b) + 15;
                            else
                                ok(1,ch,b,d4,d5,d6) = false; % give up - insufficient SNR
                            end
                        end
                        if ok(1,ch,b,d4,d5,d6)
                            tend = find(IR2smoothdB(maxind(1,ch,b):end,ch,b) <= LateSlopeEnddB, 1, 'first')+maxind(1,ch,b)-1; %
                            tstart = find(IR2smoothdB(maxind(1,ch,b):end,ch,b) <= LateSlopeStartdB, 1, 'first')+maxind(1,ch,b)-1; %
                            if ~(isempty(tend) || isempty(tstart))
                                
                                o(:,ch,b) = polyfit((tstart:tend)', ...
                                    IR2smoothdB(tstart:tend,ch,b),1)';
                                Tlate(1,ch,b,d4,d5,d6) = 2*((o(2,ch,b)-60)/o(1,ch,b) ...
                                    -(o(2,ch,b))/ ...
                                    o(1,ch,b))/fs; % Late reverberation time
                                
                                % 9. FIND NEW crosspoint1
                                crosspoint1(1,ch,b) = -round((o(2,ch,b)-IR2taildB(1,ch,b))/o(1,ch,b));
                            else
                                ok(1,ch,b,d4,d5,d6) = false;
                            end
                        end
                    end
                end
            end
            
            % DOUBLE-CHECK THAT crosspoint1 IS REALLY OK
            for ch = 1:chans
                for b = 1:bands
                    if crosspoint1(1,ch,b) <= maxIR2smoothdB(1,ch,b)+0.05*fs || ... % at least 50 ms decay
                            crosspoint1(1,ch,b) > len
                        crosspoint1(1,ch,b) = len; % give up - no crosspoint1 found
                        ok(1,ch,b,d4,d5,d6) = false;
                    end
                end
            end
            crosspoint(:,:,:,d4,d5,d6) = crosspoint1;
        end
    end
end

%crosspoint(~ok) = round(0.9*len);
crosspoint(~ok) = len;