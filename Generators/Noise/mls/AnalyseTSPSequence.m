function ir = AnalyseTSPSequence(signal,offset,reps,N, impalign, lincirc)

% AnalyseTSPSequence
%
%   ir = AnalyseTSPSequence(signal,[offset],[reps],[N],[impalign])
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
%                   (default = 4).
%       N:          (Optional) Order of the TSP, where P=2^N. 
%                   (default = 14).
%       impalign:   1 to include alignment impulse, 0 to exclude
%       lincirc:    (Optional) 'lin' assumes 2^N gap between bursts and
%                   applies linear deconvolution. 'circ' assumes no gap
%                   between bursts and applies circular deconvolution.
%   Outputs:    
%       ir:         An NxM vector of impulse responses, where N is the 
%                   order of the TSP and M is the number of channels.
%   References:
%       Y. Suzuki, F. Asani, H.-Y. Kim and T. Sone, "An optimum
%       computer-generated pulse signal suitable for the measurement of
%       very long impulse responses," J. Acoust. Soc. Am. Vol. 97(2), pp.
%       1119-1123, 1995.
%
%**************************************************************************
% Author:           Wrapper by M. R. P. Thomas, functional code F. Asano
% Date:             17 Aug 2009
% Last modified:    17 Aug 2009
%**************************************************************************

error(nargchk(1,6,nargin)); % Change to allow some defaults

if (nargin < 6)
    lincirc = 'lin';
end
if (nargin < 5)
    impalign = 1;
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

if ~(strcmp(lincirc,'lin') || strcmp(lincirc,'circ'))
    error('Argument strcmp must be either ''lin'' or ''circ''');
end

normalize = 1;
P=2^N;

if (impalign)
    impulset = zeros(1,length(signal(1,:)));
    for i=1:length(signal(1,:)) % Go through all channels
        k=findpeaks_mls(signal(:,i),'',2^N);
        impulset(i) = k(1);
    end
    impulseindex = min(impulset)
    absstartindex = impulseindex + P - offset -1;
else
    absstartindex = 1;
end


sz = size(signal);
signal = [signal;zeros(P,sz(2))];   % Append P zeros to end in case of overrun

[tsp invtsp] = GenerateTSP(N,normalize);

if(strcmp(lincirc,'lin'))
    for i=1:sz(2)
        acc = 0;
        startindex = absstartindex;
        for j=1:1:reps;
            acc = acc + signal(startindex:startindex+2*P-1,i);
            startindex = startindex + 2*P;
        end

        mymean(:,i) = acc/(reps);
        
        tmpir = fftfilt(invtsp,mymean(:,i));
        ir(:,i) = tmpir(P+1:2*P);
    end
elseif(strcmp(lincirc,'circ'))
    %error('Circular deconvolution not yet implemented');

    % For circular convolution
    % Go through all channels
    for i=1:sz(2)
        acc = 0;
        startindex = absstartindex;
        for j=1:1:reps+1;
            acc = acc + signal(startindex:startindex+P-1,i);
            startindex = startindex + P;
        end

        mymean(:,i) = acc/(reps);

        ir(:,i) = real(ifft(fft(mymean(:,i)).*fft(invtsp)));
    end
end
