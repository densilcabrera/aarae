function OUT = pseudorand1(duration,seed,fs)
% Generates a pseudo-random function over a particular time range.
% For a given 'seed', the same 'white' noise wave will be generated every
% time this function is called.
%
% Based on:
% S. K. Park and K. W. Miller, "Random number generators: Good ones are
% hard to find", Communications of the ACM, 31 (1988), pp. 1192-1201.

if nargin == 0
    param = inputdlg({'Duration of the wave [s]';...
        'seed';...
        'Sampling frequency [samples/s]'}, ...
        'Impulse input parameters',1,{'1';num2str(2147483646);'48000'});
    param = str2num(char(param));
    if length(param) < 3, param = []; end
    if ~isempty(param)
        duration = param(1);
        seed = param(2);
        fs = param(3);
    end
else
    param = [];
end
if ~isempty(param) || nargin ~= 0
    % Parameters for pseudo-random number generator, as recommended in
% S. K. Park and K. W. Miller, "Random number generators: Good ones are
% hard to find", Communications of the ACM, 31 (1988), pp. 1192-1201.
a = 16807; % 7^5
m = 2147483647; % 2^31-1;
y = seed * ones(ceil(duration * fs),1);

% generate the pseudo-random number sequence
for k = 1:(ceil(duration * fs)-1)
    y(k+1) = mod(a*y(k),m);
end
y = 2*(y-1) ./(max(y)-1) -1; % normalize to -1 / 1

    tag = ['Pseudorand1' num2str(duration)];
    
    OUT.audio = y;
    OUT.audio2 = flipud(y);
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.properties.dur = duration;
    OUT.properties.seed = seed;
    OUT.funcallback.name = 'pseudorand1.m';
    OUT.funcallback.inarg = {duration,seed,fs};
else
    OUT = [];
end

end % End of function