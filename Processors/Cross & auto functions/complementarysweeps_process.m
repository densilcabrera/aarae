function OUT = complementarysweeps_process(IN)
% This function stacks a pair of phase-complementary sweeps into a single
% sweep, which can subsequently be converted into an impulse response using
% the usual process (with the '*' button).
%
% Use the complementarysweeps generator to generate a pair of exponential
% sweeps, where the second one is delayed and phase-inverted.
%
% Note that this function will return no output if the required fields are
% not present in the input recording (notably
% IN.properties.complementarysweepsgap). It will also return no output if
% the input audio is too short.


if ~isfield(IN,'properties') 
    OUT = [];
    return
end

% check for the distinctive field from complementarysweeps generator
if ~isfield(IN.properties,'complementarysweepsgap') 
    OUT = [];
    return
else
    gap = IN.properties.complementarysweepsgap;
end

if ~isfield(IN,'audio2')
    OUT = [];
    return
else
    sweeplen = length(IN.audio2) + round(IN.fs*gap);
end

% stack multicycle measurements in dimension 4
% if isfield(IN,'properties')
%     if isfield(IN.properties,'startflag')
%         IN = stackaudio_aarae(IN,1);
%         IN = rmfield(IN.properties , 'startflag');
%         %IN = rmfield(IN , 'relgain');
%     end
% end

if size(IN.audio,1) >= 2*sweeplen
    % combine the two sweeps by subtraction
    audio = IN.audio(1:sweeplen,:,:,:,:,:) ...
        - IN.audio(sweeplen+1:2*sweeplen,:,:,:,:,:);
else
    OUT = [];
    return
end

OUT = IN;
OUT.audio = audio;
OUT.funcallback.name = 'complementarysweeps_process.m';
OUT.funcallback.inarg = {};