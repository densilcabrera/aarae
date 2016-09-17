function OUT = RLB2004(IN,fs)
% RLB weighting filter, based on:
% G. Souldore, "Evaluation of Objective Loudness Meters",
% Presented at the 116th Audio Engineering SocietyConvention,
% 2004 May 8-11 Berlin, Germany.
%
% The RLB filter is a simple high-pass filter, which was 
% found to be an effective weighting filter for estimating the 
% loudness of various program materials prepared for broadcast audio.
%
% This processor was ported from PsySound3 by Densil Cabrera
% version 1.00 (28 March 2014).

processed = [];
if nargin < 2
    if isstruct(IN)
        audio = IN.audio;
        fs = IN.fs;
    else
        audio = IN;
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
else
    audio = IN;
end
if ~isempty(audio) && ~isempty(fs)
    %if isdir([cd '/Processors/Filters/' num2str(fs) 'Hz'])
    if false % bypass ths code
        content = load([cd '/Processors/Filters/' num2str(fs) 'RLB2004Filter.mat']);
        filterbank = content.filterbank;
        processed = filter(filterbank,1,audio);
        
    else
        % This core part of the function was written by Farhan Rizwi for
        % the PsySound3 project.
        % Filter coefficients from [1], pg 12
        % These are defined for 48k
        b = [1 -2 1];
        a = [1 -1.99004745483398 0.99007225036621];
        if fs ~= 48000
            poles = roots(a);
            
            % Make polynomial after fixing up the roots
            %
            % z = exp(s*T) --> s = ln(z)/T
            %
            % s = ln(z1)/T1 = ln(z2)/T2  -->  z2 = exp(ln(z1)*T2/T1)
            %
            a = poly(exp(log(poles)*48e3/fs));
            
            % Note that the two zeros at 1 remain there.
            % Note also, that the negligible high frequency gain adjustment
            % is ignored.
        end
        processed = filter(b,a,audio);
    end
    if isstruct(IN)
        OUT = IN;
        OUT.audio = processed;
        OUT.funcallback.name = 'RLB2004.m';
        OUT.funcallback.inarg = {};
    else
        OUT = processed;
    end
end