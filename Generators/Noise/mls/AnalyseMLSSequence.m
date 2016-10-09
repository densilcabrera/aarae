function ir = AnalyseMLSSequence(signal,offset,reps,N,DCCoupling, impalign)

% AnalyseFullSequence
%
% Detects alignment impulse and analyses MLS bursts. USE IN CONJUNTION WITH
% GenerateMLSSequence; MLS SIGNAL HAS 2^N-1 DIFFERENT POSSIBILITIES FOR ANY
% GIVEN N SO IT'S ESSENTIAL YOU USE THE RIGHT ONE!
%
%   ir = AnalyseMLSSequence(signal,[offset],[reps],[N],[DCCoupling],[impalign])
%
%   Inputs
%       signal:     NxM recorded sequence from system under test, where N
%                   is the length of the recording and M is the number of
%                   channels.
%       offset: 	(Optional) The impulse detection finds the point where 
%                   the signal exceeds a certain rate of change. Pulse 
%                   spreading may cause the detection to be erroneous by a 
%                   few samples (rarely more than 10) and this may be 
%                   corrected with the offset parameter. Positive offset 
%                   causes negative time shift (default = 0).
%       reps:       (Optional) The number of reps of each amplitude 
%                   (default = 4). Note change: Initial burst and tail are 
%                   now summed to form one 'repetition' - no bursts are 
%                   discarded.
%       N:          (Optional) Order of the MLS, where P=2^N-1. 
%                   (default = 14).
%       DCCoupling: (Optional) Set to true for DC recovery method. Set to 
%                   false for loudspeaker measurements (default = true).
%       impalign:   1 to include alignment impulse, 0 to exclude
%   Outputs:    
%       ir:         An NxM vector of impulse responses, where N is the 
%                   order of the MLS and M is the number of channels.
%	References:
%					[1] J. Borish and J. B. Angell, "An Efficient Algorithm for Measuring the Impulse Response Using Pseudorandom Noise," J. Audio Eng. Soc., vol. 31, pp. 478–489, 1983.
%					[2] D. D. Rife and J. Vanderkooy, "Transfer-Function Measurement with Maximum-Length Sequences," J. Audio Eng. Soc., vol. 37, pp. 419–444, 1989.%
%					[3] M. Cohn and A. Lempel, "On Fast M-Sequence Transforms," IEEE Trans. Inf. Theory, IT-23, pp. 135–137, 1977.

%**************************************************************************
% Author:           M. R. P. Thomas 
% Date:             21 Feb 2007
% Last modified:    09 May 2007
%**************************************************************************

narginchk(1,6); % Change to allow some defaults

if (nargin < 6)
    impalign = 1; 
end
if (nargin < 5)
    DCCoupling = 1;
end
if (nargin < 4)
    N = 14;
end
if (nargin < 3)
    reps = 4;
end
if (nargin < 2)
    offset = 0;
end

P=2^N-1;
mls = GenerateMLS(N,1);
tagS = GeneratetagS(mls,P,N);
tagL = GeneratetagL(mls,P,N);

if (impalign)
    impulset = zeros(1,length(signal(1,:)));
    for i=1:length(signal(1,:)) % Go through all channels
        k=findpeaks_mls(signal(:,i),'',2^N);
        impulset(i) = k(1);
    end
    impulseindex = min(impulset)
    absstartindex = impulseindex + P - offset;
else
    absstartindex = 1;
end


sz = size(signal);
signal = [signal;zeros(P,sz(2))];   % Append P zeros to end in case of overrun

% Go through all channels
for i=1:sz(2)
    acc = 0;
    startindex = absstartindex;
    for j=1:1:reps+1;
        acc = acc + signal(startindex:startindex+P-1,i);
        startindex = startindex + P;
    end
    
    mean(:,i) = acc/(reps); %-1 because we are ignoring the first burst
    ir(:,i) = AnalyseMLS(mean(:,i),tagS,tagL,N,DCCoupling);
end

%% Analyses an MLS sequence. 
function impulseresp = AnalyseMLS(signal,tagS,tagL,N,DCCoupling)
% signal: Must be the same length as the MLS sequence and be taken from the
% beginning of a SECOND MLS sequence to approximate circular convolution.
%
% mls: a P-length MLS sequence (Where P=2^N-1)
%
% fs: sampling frequency in Hz
%
% DCCoupling: True if device under test is DC coupled, false otherwise.

P = 2^N-1;
perm = PermuteSignal(signal, tagS, P, DCCoupling);
had = FastHadamard(perm, P+1, N);
resp = PermuteResponse(had, tagL, P);
impulseresp = resp;

%% GeneratetagL
function tagL = GeneratetagL(mls, P, N)
% Generates array for the rearrangement of samples after a Hadamard Transform
%
% tagL = GeneratetagL(mls, P, N);
%
% mls:		The MLS signal for which tagL is valid.
% P:		Length of the MLS signal
% N:		Order of the MLS signal
% tagL:		1x(P+1) vector of indices.


% Convert {-1,1} to binary
binmls = (mls-1)./-2; 

S = GeneratetagS(mls,P,N);

% Find which values of the tagS vector are powers of 2
for i=1:1:P
    for j=1:1:N
        if (S(i) == 2^(j-1))
            index(j) = i;
        end
    end
end

powerindices = 0:1:N-1;
powers = 2.^powerindices;

for i=1:1:N
  L(i,1:mod(index(i),P)) = binmls(mod(index(i),P):-1:1);
  L(i,mod(index(i),P)+1:P) = binmls(P:-1:mod(index(i),P)+1);
end

tagL = powers*L;

%% GeneratetagS
function tagS = GeneratetagS(mls, P, N)
% Generates array for the rearrangement of samples before a Hadamard Transform
%
% tagS = GeneratetagS(mls, P, N);
%
% mls:		The MLS signal for which tagS is valid.
% P:		Length of the MLS signal
% N:		Order of the MLS signal
% tagS:		1x(P+1) vector of indices.


% Convert {-1,1} to binary
binmls = (mls-1)./-2; 

% Make S matrix by making first line mls and shifting every subsequent row
% RIGHT up to N then multiply each row by the correct power of 2.
powerindices = N-1:-1:0;
powers = 2.^powerindices;

for i=1:1:N
   S(i, 1:i-1) = binmls(P-i+2:P);
   S(i, i:P) = binmls(1:P-i+1);
end

tagS = powers*S;

%% PermuteSignal
function perm = PermuteSignal(signal, tagS, P, dcCoupled)
% Rearranges input signal according to tagS.
%
% perm = PermuteSignal(signal, tagS, P, dcCoupled);
%
% signal:       Signal to be arranged
% tagS:         1xP vector of indces from GeneratetagS
% P:            Length of MLS, where P=2N-1
% DCCoupled:	Set to true for DC recovery method. Set to false for 
%               loudspeaker	measurements.
% perm:     	The rearranged signal


%DC coupling:
if (dcCoupled == 1)
    dc = 0;
    for i=1:1:P
        dc = dc + signal(i);
    end
    perm(1) = -dc;
else
    perm(1) = 0;    % Not sure if this is the right thing to do yet...
end

%for i=1:1:P
%    perm(tagS(i)+1) = signal(i);
%end
perm(tagS+1) = signal;

%% PermuteResponse
function resp = PermuteResponse(perm, tagL, P)
% Rearranges output of Hadamard Transform according to tagL.
%
% resp = PermuteResponse(perm, tagL, P)
%
% perm:     Output of Hadamard Transform
% tagL:     1xP vector of indices from GeneratetagL.
% P:        Length of MLS, where P=2N-1
% Resp:     Rearranged signal


fact = 1/(P+1);

%Ignore first element
perm = perm(2:end);

%for i=1:1:P
%    resp(i) = perm(tagL(i))*fact;
%end

resp = perm(tagL).*fact;

resp(P+1) = 0;

%% FastHadamard
function y = FastHadamard(x, P1, N)
% Applies a Fast Hadamard transform to a 1-D signal
%
% y = FastHadamard(x, P, N)
%
% x:	Signal to be transformed
% P:	Length of MLS, where P=2N-1
% N:	Order of MLS
% y: 	Transformed signal

k1 = P1;
for k=1:1:N
    k2 = k1/2;
    for j=1:1:k2
        for i=j:k1:P1
            i1 = i + k2;
            temp = x(i) + x(i1);
            x(i1) = x(i) - x(i1);
            x(i) = temp;
        end
    end
    k1 = k1/2;
end

y = x;