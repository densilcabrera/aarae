function X = RoomNoiseSTIAnalysis(X,IRduration,Noisetype,SubtractNoise)
% This function analyses a multiple test signal designed for assessing
% room acoustics, background noise and STI (direct and indirect) in rooms.
% It was originally written for measurements in lecture theatres.
%
% Use the RoomNoiseSTITestSignal generator to generate a test signal for
% this analyser.
% The test signal consists of STIPA (30 s), followed by silence (35 s, followed by an
% exponential sweep (60 s).
% Default STIPA band values are for the NTI Talkbox
%
% This analyser derives:
% * STIPA (from the STIPA signal recording)
% * octave band speech level (from the STIPA recording)
% * octave band noise level (from the silence recording), including RC and
%    RNC values
% * Reverberation time, clarity index and related parameters (from the
%    sweep recording)
% * Indirect STI (from the speech level, noise level and sweep)
% The analyser returns an impulse response, and includes effective
% background noise of the impulse response in the second channel
%
% Due to the large amount of output data, it may be easiest to use the log
% file to read the results.
%
% Currently this generator and its analyser are only tested with single
% channel, single cycle measurements.

if ~isfield(X, 'STIPAsignal')
    disp('This is not the right test signal for this analyser')
    X = [];
    return
end

if nargin == 1
    param = inputdlg({'Duration of impulse response [s]';...
        'Background noise type: Leq (0), L90 fast (1) or L90 slow (2)';...
        'Subtract noise level from speech level [0 | 1]'},...
        'Analysis parameters',1,{'3';'0';'1'});
    param = str2num(char(param));
    if length(param) < 3, X = []; return; end
    if ~isempty(param)
        IRduration = param(1);
        Noisetype = param(2);
        SubtractNoise = param(3);
    else
        X = [];
        return
    end
end

cal = X.cal;
[len,chans,~,cycles,outchans,dim6] = size(X.audio);

% Find and compensate for latency of the system, write lags to log file
Cindex = zeros(chans,cycles,outchans,dim6);
startindex = X.properties.startindex;
STIPArecording = zeros(startindex(2),chans,cycles,outchans,dim6);
Silencerecording = zeros(startindex(3)-startindex(2),chans,cycles,outchans,dim6);
Sweeprecording = zeros(len-startindex(3),chans,cycles,outchans,dim6);
for ch = 1:chans
    for d4 = 1:cycles
        for d5 = 1:outchans
            for d6 = 1:dim6
                [C, LAGS] = xcorr(X.audio(1:10*X.fs,ch,1,d4,d5,d6), X.STIPAsignal(1:10*X.fs));
                Cindex(ch,d4,d5,d6) = find(C == max(C),1,'first');
                logtext(['LAG, ' num2str(Cindex(ch,d4,d5,d6)) ', samples, Chan ' num2str(ch) ' Cycle ' num2str(d4) ' Outchan ' num2str(d5) ' D6 ' num2str(d6)]);
                STIPArecording(:,ch,1,d4,d5,d6) = X.audio(LAGS(Cindex(ch,d4,d5,d6)):LAGS(Cindex(ch,d4,d5,d6))+X.properties.startindex(2)-1,ch,1,d4,d5,d6);
                Silencerecording(:,ch,1,d4,d5,d6) = X.audio(LAGS(Cindex(ch,d4,d5,d6))+X.properties.startindex(2):LAGS(Cindex(ch,d4,d5,d6))+X.properties.startindex(3)-1,ch,1,d4,d5,d6);
                Sweeprecording(1:(len-(LAGS(Cindex(ch,d4,d5,d6))+X.properties.startindex(3))+1),ch,1,d4,d5,d6) = X.audio(LAGS(Cindex(ch,d4,d5,d6))+X.properties.startindex(3):end,ch,1,d4,d5,d6);
            end
        end
    end
end
X.properties.lags = Cindex; % include lags as a new field in output structure

% Process sweep to make an IR, which will be returned by this function
startindex = 1444000-24000; % this can be refined
endindex = startindex+round(IRduration*X.fs);
X.audio = convolveaudiolocal(Sweeprecording,X.audio2);
X.audio = X.audio(startindex:endindex,:,:,:,:,:);

% Analyse silence as if it is a sweep, to allow for effective SNR analysis
IRsilence = convolveaudiolocal(Silencerecording,X.audio2);
IRsilence = IRsilence(startindex:endindex,:,:,:,:,:);

for ch = 1:chans
    for d4 = 1:cycles
        for d5 = 1:outchans
            for d6 = 1:dim6
                % Analyse STIPA signal to determine speech level
                STIPAresult = STIPA_direct(STIPArecording(:,ch,1,d4,d5,d6),48000,cal(ch),[],1,1,0);
                speechlevel = permute(STIPAresult.Level,[2,3,1]);
                
                % Analyse octave band sound pressure level of silence
                %X = octave_band_level_barplot(Silencerecording,48000,0,1,16,16000,0,0,0);
                NOISEresult = RoomNoiseEvaluation_NC_RNC(Silencerecording(:,ch,1,d4,d5,d6),48000,cal(ch),1);
                switch Noisetype
                    case 1
                        noiselevel = NOISEresult.L90fast(4:end);
                    case 2
                        noiselevel = NOISEresult.L90slow(4:end);
                    otherwise
                        noiselevel = NOISEresult.Leq(4:end);
                end
                
                % Subtract noise from speech - could cause problems if noise is not
                % constant! In benign cases this will have little effect.
                if SubtractNoise
                    subtractok = speechlevel>noiselevel;
                    speechlevel(subtractok) = pow2db(db2pow(speechlevel(subtractok)) - db2pow(noiselevel(subtractok)));
                    speechlevel(~subtractok) = -10; % -10 dB effectively no speech sound
                end
                
                X1 = X;
                X1.audio = X1.audio(:,ch,1,d4,d5,d6);
                % Reverberation parameters and indirect STI
                RTresult = ReverberationTimeIR2(X1,-20,1,125,8000,1,-1,1,2,-5,-15,1);
                
                STI_IRresult = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);
                
                % STI with various NC curve background noise levels
                noiselevel = [44 37 31 27 24 22 22]; % NC25
                 STI_IRresultNC25 = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);
                noiselevel = [48 41 35 32 29 28 27];% NC30
                 STI_IRresultNC30 = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);
                noiselevel = [52 45 40 36 34 33 32];% NC35
                 STI_IRresultNC35 = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);
                noiselevel = [56 50 44 41 39 38 37];% NC40
                 STI_IRresultNC40 = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);
                noiselevel = [60 54 49 46 44 43 42];% NC45
                 STI_IRresultNC45 = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);
                noiselevel = [64 58 54 51 49 48 47];% NC50
                 STI_IRresultNC50 = STI_IR(X1,48000,speechlevel,noiselevel,2011,1,1,2,1,0);

                
                % harmonic distortion - potentially could be included in
                % future
                % X1.audio = cat(4,IRsilence(:,ch,1,d4,d5,d6),X1.audio);
                % X1.properties.relgain = [-inf 0];
                % THDresult = ;
                

                % Concatenate tables
                if ~isfield(X,'tables'), X.tables = []; end
                X.tables = [X.tables STIPAresult.tables NOISEresult.tables ...
                    RTresult.tables STI_IRresult.tables ...
                    STI_IRresultNC25.tables STI_IRresultNC30.tables STI_IRresultNC35.tables ...
                    STI_IRresultNC40.tables STI_IRresultNC45.tables STI_IRresultNC50.tables];
            end
        end
    end
end


X.audio = cat(4,X.audio, IRsilence);

X.funcallback.name = 'RoomNoiseSTIAnalysis.m';
X.funcallback.inarg = {IRduration,Noisetype,SubtractNoise};
X = rmfield(X,{'STIPAsignal','audio2'}); % remove extra audio fields to reduce data size
end % eof


function IR = convolveaudiolocal(audio1,audio2)
fftlen = length(audio1)+length(audio2)-1;
IR = ifft(fft(audio1,fftlen) .* fft(audio2,fftlen));
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