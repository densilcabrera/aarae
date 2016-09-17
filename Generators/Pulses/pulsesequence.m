function OUT = pulsesequence(dur,fs,reverse,threshold,choice,parameter)
% This function generates a set of indices which are used to make pulses in
% an otherwise silent waveform.
%
% MAIN INPUTS
% dur is the duration of the wave in seconds.
%
% fs is the audio sampling rate in Hz.
%
% reverse: whether or not to reverse the waveform in time.
%
% threshold is used in generating an approximate inverse filter. Spectrum
% componenents in the inverted spectrum that are greater than the minimum
% spectrum component's level plus threshold are zeroed. However, the DC
% component is ignored in this processing.
%
% TYPES OF PULSE SERIES (choice)
%
% PERIODIC
% A periodic sequence of pulses with the following parameters:
%   1. Period between pulses (in samples)
%   2. Jitter standard deviation (in samples), created using Matlab's randn
%
% CUMULATIVE SUM OF A PERIODIC SEQUENCE
% A periodic sequence is generated, jitter applied (if > 0), and then the
% sequence is cumulatively summed at least once. This sequence of pulses
% has the following parameters:
%   1. Initial step (in samples)
%   2. Cumulative sum order (number of times the sequence is cumulatively
%   summed)
%   3. Jitter standard deviation (in samples), created using Matlab's
%   randn, which is applied prior to the cummulative summation.
%
%
% UNIFORM RANDOM
% Unique indices from a uniform random distribution are used to create
% pulses. This has the following parameter:
%   1. Density - a value greater than 0 and less than 1

if nargin < 3
    param = inputdlg({'Duration of the wave (s)';...
        'Sampling rate (Hz)';...
        'Reverse [0 | 1]';...
        'Inverse spectrum threshold (dB)'},...
        'Settings',... % This is the dialog window title.
        [1 60],...
        {'1';'48000';'0';'80'}); % default answers
    param = str2num(char(param));
    if length(param) < 4, param = []; end
    if ~isempty(param)
        dur = param(1);
        fs = param(2);
        reverse = param(3);
        threshold = param(4);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end

len = round(dur*fs);
audio = zeros(len,1);

if ~exist('choice', 'var')
    % Choice of series type
    S = {'Periodic';...
        'Cumsum of periodic';...
        'Exponential';...
        'Uniform Random';...
        'Primes';...
        'Fibonacci'};
    
    [choice,ok] = listdlg('ListString',S,...
        'SelectionMode','Single',...
        'Name','HELLO');
    
    if ~ok
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end

switch choice
    case 1
        % Periodic
        if ~exist('parameter','var')
            param = inputdlg({'Period between pulses (samples)';...
                'Jitter standard deviation (samples)'},...
                'Arithmetic pulse series settings',... % This is the dialog window title.
                [1 60],...
                {num2str(round(fs/20));'0'}); % default answers
            param = str2num(char(param));
            if length(param) < 2, param = []; end
            if ~isempty(param)
                period = param(1);
                jitter = param(2);
            else
                % get out of here if the user presses 'cancel'
                OUT = [];
                return
            end
        else
            period = parameter(1);
            if length(parameter)>1
                jitter = parameter(2);
            else
                jitter = 0;
            end
        end
        if period < 1, period = 1; end
        indices = 1:period:len;
        if jitter > 0
            indices = jitter*randn(1,length(indices)) + indices;
        end
        indices = round(indices);
        audio(unique(indices(indices>0 & indices<=len))) = 1;
        parameter = [period jitter];
        tag = 'PeriodicPulses';
    case 2
        % Cumulative sum (integral) of periodic series of indices
        if ~exist('parameter','var')
            param = inputdlg({'Cumulative sum order';...
                'Initial step (samples)';...
                'Jitter standard deviation (samples)'},...
                'Cumulative sum pulse series settings',... % This is the dialog window title.
                [1 60],...
                {'1';'1';'0'}); % default answers
            param = str2num(char(param));
            if length(param) < 3, param = []; end
            if ~isempty(param)
                order = param(1);
                step = param(2);
                jitter = param(3);
            else
                % get out of here if the user presses 'cancel'
                OUT = [];
                return
            end
        else
            order = parameter(1);
            if length(parameter)>1
                step = parameter(2);
                if length(parameter)>2
                    jitter = parameter(3);
                else
                    jitter = 0;
                end
            else
                step = 1;
                jitter = 0;
            end
        end
        indices = 1:step:len;
        if jitter > 0
            indices = jitter*randn(1,length(indices)) + indices;
        end
        for o = 1:order
            indices = cumsum(indices);
        end
        indices = round(indices);
        audio(unique(indices(indices>0 & indices<=len))) = 1;
        parameter = [step order jitter];
        tag = 'CumsumPeriodicPulses';
    case 3
        % exponential
        disp('not implemented')
        OUT = [];
        return
    case 4
        % uniform random
        if ~exist('parameter','var')
            param = inputdlg({'Pulse density (between 0 and 1)'},...
                'Arithmetic pulse series settings',... % This is the dialog window title.
                [1 60],...
                {'0.5'}); % default answers
            param = str2num(char(param));
            if length(param) < 1, param = []; end
            if ~isempty(param)
                density = param(1);
            else
                % get out of here if the user presses 'cancel'
                OUT = [];
                return
            end
        else
            density = parameter(1);
        end
        indices = round(len*rand(round(len*density),1));
        audio(unique(indices(indices>0 & indices<=len))) = 1;
        parameter = density;
        tag = 'RandomPulses';
    case 5
        % primes
        audio(isprime(2:len+1))=1;
        parameter = [];
        tag = 'PrimesPulses';
    case 6
        % fibonacci
        n = 1:100; % unlikely to need more than this!
        phi = (1+sqrt(5))/2;
        indices = round((phi.^(n+1) - (1-phi).^(n+1))/(2*phi-1));
        audio(indices(indices<=len)) = 1;
        parameter = [];
        tag = 'FibonacciPulses';
        
end

if reverse == 1
    audio = flipud(audio);
end

% Generate a kind of inverse filter (probably not particularly useful)
audio2 = 1 ./ fft(audio);
% exclude very large spectrum components (ignore DC in applyng the threshold)
audio2(abs(audio2) > abs(db2mag(threshold))*min(abs(audio2(2:end)))) = 0;
audio2 = ifft(audio2);

OUT.audio = audio;
OUT.audio2 = audio2;
OUT.fs = fs;
OUT.tag = tag;
OUT.properties.parameter = parameter;
OUT.funcallback.name = 'pulsesequence.m';
OUT.funcallback.inarg = {dur,fs,reverse,threshold,choice,parameter};


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
