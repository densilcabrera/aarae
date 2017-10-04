function OUT = RoomNoiseSTITestSignal(duration, fs, octavebandlevel)
% This function generates a multiple test signal designed for assessing
% room acoustics, background noise and STI (direct and indirect) in rooms.
% It was originally written for measurements in lecture theatres.
%
% The test signal consists of STIPA (30 s), followed by silence (35 s, followed by an
% exponential sweep (60 s).
% Default STIPA band values are for the NTI Talkbox
%
% Analysis should be done with the RoomNoiseSTIAnalysis analyser, which
% derives:
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
% Currently this generator and its analyser are only tested with single
% channel, single cycle measurements.


if nargin == 0
    param = inputdlg({'Duration of the wave [s]';...
        '125 Hz octave band level (dB)';...
        '250 Hz octave band level (dB)';...
        '500 Hz octave band level (dB)';...
        '1000 Hz octave band level (dB)';...
        '2000 Hz octave band level (dB)';...
        '4000 Hz octave band level (dB)';...
        '8000 Hz octave band level (dB)';...
        'Broadband gain (dB)'; ...
        'Sampling frequency [samples/s]'}, ...
        'Noise input parameters',1,{'30'; ...
        '-0.7';'0';'-3';'-9.4';'-16.5';'-23.4';'-31.6'; ...
        '0';'48000'});
    param = str2num(char(param));
    if length(param) < 10
        OUT = [];
        return
    end
    duration = param(1);
    octavebandlevel(1) = param(2);
    octavebandlevel(2) = param(3);
    octavebandlevel(3) = param(4);
    octavebandlevel(4) = param(5);
    octavebandlevel(5) = param(6);
    octavebandlevel(6) = param(7);
    octavebandlevel(7) = param(8);
    Lgain = param(9);
    fs = param(10);
    octavebandlevel = octavebandlevel + Lgain;
end
if exist('fs','var') && fs<44100, fs = 44100; end


% Generate sine sweep
OUT = exponential_sweep(30,80,12500,fs,0,150);

% Generate STIPA
OUT2 = STIPA_signal(duration,fs,octavebandlevel);

% Generate silence
silence = zeros(fs*35,1);


OUT.audio = [OUT2.audio; silence; OUT.audio];
OUT.STIPAsignal = OUT2.audio; % write the STIPA signal as an extra field
OUT.properties.startindex =[1, length(OUT2.audio)+1,length(OUT2.audio)+length(silence)+1];
OUT.tag = 'RoomNoiseSTI';
OUT.funcallback.name = 'RoomNoiseSTITestSignal.m'; 
OUT.funcallback.inarg = {}; 





end % End of function

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