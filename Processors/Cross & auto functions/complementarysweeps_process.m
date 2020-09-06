function OUT = complementarysweeps_process(IN)
% This function stacks a pair of phase-complementary sweeps into a single
% sweep, which can subsequently be converted into an impulse response using
% the usual process (with the '*' button).

% check for the distinctive field from complementarysweeps generator
if ~isfield(IN,'properties')
    OUT = [];
    return
end

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