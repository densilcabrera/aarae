function OUT = OctSpacedSweep(durfinal,fs,reverse,tukeyratio,hannexp)
% This function generates an octave-spaced sweep, essentially following the
% concept of the Shepard-Risset glissando
%
% The main aim of this is to demonstrate the chroma vs pitch-height
% distinction, or just for fun

if nargin == 0 
    
    param = inputdlg({'Waveform duration (s)';... 
        'Time per octave (s)';...
        'Sampling rate (Hz)';...
        'Ascending [0] or descending [1] sweep';...
        'Overall fade-in, fade-out Tukey window ratio';...
        'Hann window exponent for spectrum envelope'},...
        'Octave Spaced Sweep',... % This is the dialog window title.
        [1 60],... 
        {'60';'10';'48000';'0';'0.1';'1'}); % defaults
    
    param = str2num(char(param));  
    if length(param) < 6, param = []; end 
    if ~isempty(param) 
        durfinal = param(1);
        octdur = param(2);
        fs = param(3);
        reverse = param(4);
        tukeyratio = param(5);
        hannexp = param(6);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end




if ~isempty(param) || nargin ~= 0
    


    noct = 10; % number of octaves
    dur = octdur*noct;
    octlen = round(fs*octdur); 
    end_freq = fs/2;
    start_freq = end_freq / 2^noct;
    SI = 1/fs;
    ampl = 0.1;
    
    w1 = 2*pi*start_freq; w2 = 2*pi*end_freq;
    K = (dur*w1)/(log(w2/w1));
    L = log(w2/w1)/dur;
    t = (0:round(dur/SI)-1)*SI;
    phi = K*(exp(t*L) - 1);
    S = (ampl*sin(phi))';
    if reverse == 1, S = flip(S); end
    
    S = S .* hann(length(S)).^hannexp;
    Smix = S;
    for n = 1:noct-1
        Smix = Smix + circshift(S,n*octlen);
    end
    
    numcat = ceil(durfinal./(length(Smix)/fs));
    if numcat > 1
        Smix = repmat(Smix,[numcat 1]);
    end
    Smix = Smix(1:floor(durfinal*fs));
    Smix = Smix .* tukeywin(length(Smix),tukeyratio);
   
    OUT.audio = Smix; % You NEED to provide the audio you generated.
    %OUT.audio2 = ?;     You may provide additional audio derived from your function.
    OUT.fs = fs;       % You NEED to provide the sampling frequency of your audio.
    OUT.tag = 'OctSpacedSweep';      %You may assign it a name to be identified in AARAE.
    %OUT.properties.prop1 = prop1;
    %OUT.properties.prop2 = prop2; You may provide additional info
    %OUT.properties.prop3 = prop3; about your generated signal in
    %                              a structure type field called
    %                              .properties
    OUT.funcallback.name = 'OctSpacedSweep.m'; % Provide AARAE
    % with the name of your function
    OUT.funcallback.inarg = {durfinal,fs,reverse,tukeyratio,hannexp}; % assign all of the
    % input parameters that could be used to call the function
    % without dialog box to the output field param (as a cell
    % array).
    
    
    
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
    OUT = [];
end

end % End of function

%**************************************************************************
% Copyright (c) 2016, Densil Cabrera
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