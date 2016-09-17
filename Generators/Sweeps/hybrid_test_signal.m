function OUT = hybrid_test_signal
% This function generates a sequence of test signals from some of AARAE's
% generators, which can be analysed to yield an impulse response using
% AARAE's Hybrid_process processor.
%
% You can use variations of the one signal type and/or a variety of signals
% to generate the one impulse response.
%
% By combining a variety of test signals into a single measurement your
% can: 
% * explore the variation between the IR results from each;
% * manage signal-to-noise ratio by combining the IRs (e.g. using average
% or median)
%
% Code by Densil Cabrera
% Version 0 (beta) (11 May 2014)


if nargin == 0 
    
    param = inputdlg({'Number of samples between signals';...
        'Number of OATSP signals';... 
        'Number of exponential sweep signals';...
        'Number of ''sweep from signal'' signals';...
        'Number of MLS signals';...
        'Number of complementary Golay signal pairs'},...
        'Choose signals to generate',... 
        [1 60],... 
        {'48000';'1';'0';'0';'0';'0'}); 
    
    param = str2num(char(param)); 
    
    if length(param) < 5, param = []; end 
 
    if ~isempty(param) 
        silencelen = param(1);
        doOATSP = param(2);
        doExpSweep = param(3);
        doSweepfromSignal = param(4);
        doMLS = param(5);
        doGolay = param(6);
    end
else
    param = [];
end




    
    
    if ~isempty(param) || nargin ~= 0
        OUT.audio = [];
        OUT.audio2 = [];
        OUT.properties.hybridindex = 1;
        OUT.properties.hybridindex2 = 1;
        OUT.properties.signals = [doOATSP,doExpSweep,doSweepfromSignal,doMLS,doGolay];
        OUT.fs = [];
        if doOATSP > 0
            for n = 1: doOATSP
                OUTtemp = OATSP;
                OUT = addsignal(OUTtemp,OUT,silencelen);
            end
        end
        
        if doExpSweep > 0
            for n = 1: doExpSweep
                OUTtemp = exponential_sweep;
                OUT = addsignal(OUTtemp,OUT,silencelen);
            end
        end  

        if doSweepfromSignal > 0
            for n = 1: doSweepfromSignal
                OUTtemp = sweep_from_signal_call;
                OUT = addsignal(OUTtemp,OUT,silencelen);
            end
        end  
        
        if doMLS > 0
            for n = 1: doMLS
                OUTtemp = mls;
                OUT = addsignal(OUTtemp,OUT,silencelen);
                OUT.properties.MLSn(n) = OUTtemp.properties.n;
                OUT.properties.MLScycles(n) = OUTtemp.properties.cycles;
            end
        end  
        
       if doGolay > 0
            for n = 1: doGolay
                OUTtemp = Golay;
                OUT = addsignal(OUTtemp,OUT,silencelen);
            end
        end  
        

            OUT.tag = 'Hybrid';     % You may assign it a name to be identified in AARAE.
            OUT.funcallback.name = 'hybrid_test_signal.m';
            %OUT.funcallback.inarg = {silencelen};
            OUT.funcallback.inarg = {};


    else
        % AARAE requires that in case that the user doesn't input enough
        % arguments to generate audio to output an empty variable.
        OUT = [];
    end
    
end % End of function

function OUT = addsignal(OUTtemp,OUT,silencelen)
OUT.audio = [OUT.audio; OUTtemp.audio;zeros(silencelen,1)];
OUT.audio2 = [OUT.audio2; OUTtemp.audio2];
OUT.properties.hybridindex = ...
    [OUT.properties.hybridindex, ...
    OUT.properties.hybridindex(end) + length(OUTtemp.audio) + silencelen];
OUT.properties.hybridindex2 = ...
    [OUT.properties.hybridindex2, ...
    OUT.properties.hybridindex2(end) + length(OUTtemp.audio2)];
if ~isempty(OUT.fs)
    if OUT.fs ~= OUTtemp.fs
        disp('Sampling rate differs between signals!')
    end
else
    OUT.fs = OUTtemp.fs;
end

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