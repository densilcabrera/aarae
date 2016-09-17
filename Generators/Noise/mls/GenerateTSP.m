function [tsp invtsp] = GenerateTSP(N, normalize)

% GenerateTSP
%
% Generates a TSP sequence of order N
%
%   tsp = GenerateTSP(N)
%
%   Inputs:
%       N:              TSP order  
%       normalize:      (Optional) If nonzero, normalize forward TSP 
%                       sequence to minimise quantization noise and 
%                       maximise analogue SNR (Default = 0); 
%   Outputs:
%       tsp:            TSP burst of effective length 2^N
%       invtsp:         Inverse TSP burst such that conv(tsp,invtsp) is an
%                       impulse at 2^N+1.
%
% References:
%       Y. Suzuki, F. Asani, H.-Y. Kim and T. Sone, "An optimum
%       computer-generated pulse signal suitable for the measurement of
%       very long impulse responses," J. Acoust. Soc. Am. Vol. 97(2), pp.
%       1119-1123, 1995.
%
%**************************************************************************
% Author:           Wrapper by M. R. P. Thomas, functional code by F. Asano
% Date:             17 Aug 2009
% Last modified:    17 Aug 2009
%**************************************************************************

nargchk(nargin,2,1);
if(nargin<2)
    normalize=0;
end

%-------------------------------------
%  Constants & Defualts
%-------------------------------------
default_tsp_length = [512 1024 2048 4096 8192 16384];
average_error_design = [10 12 14 14 15 15];
max_error_design = [7 10 12 13 14 15];
default_stretch = max_error_design;

MAX_MAG = 10000;

%-------------------------------------
%  Parameter Input
%-------------------------------------
%N = input('Tsp Length: ');
N = 2^N;
if mod(log2(N),1) ~= 0
	fprintf( 'TSP Length must be power of 2\n' );
	return;
end

default_stretch_id = 0;
for k=1:length(default_tsp_length)
	if N == default_tsp_length(k)
		default_stretch_id = k;
		break;
	end
end 

if default_stretch_id>0
	tsp_stretch = default_stretch(default_stretch_id);
elseif N>max(default_tsp_length)
	tsp_stretch = 15;
else
	tsp_stretch = input( 'TSP Stretch (integer from 1 to 15):' );
end

%-----------------------------------------------------------------
%  Report Paramters
%-----------------------------------------------------------------
M = (tsp_stretch/32)*N;
M_ratio = (2*M)/N*100;
%
% fprintf( '******** Parameters **********\n' );
% fprintf( 'Length of TSP : %d\n', N );
% fprintf( 'Stretch: %d\n', tsp_stretch );
% fprintf( 'M : %d \n', M );
% fprintf( 'Effective Length of TSP : %d (%4.1f [percent])\n', 2*M, M_ratio );

%-----------------------------------------------------------------
%  Design of TSP & Inverse TSP
%-----------------------------------------------------------------
H = zeros(1,N);

H(1) = 1;
H(N/2+1) = exp(j*M*pi);
%
for k=1:N/2-1
	H(k+1) = exp(j*4*M*pi*k^2/N^2);
	H(N-k+1) = conj(H(k+1));
end

G = 1./H;

%----------------------------------------------------------------
%  IFFT
%----------------------------------------------------------------
h = real(ifft(H,N));
% MRPT: Optional normalization of forward sequence.
if(normalize)
    hlev = max(abs(h));
    h = h/hlev;
    g = hlev*real(ifft(G,N));
else
    g = real(ifft(G,N));
end

%----------------------------------------------------------------
%  Rotation
%----------------------------------------------------------------
buf=h;
for k=1:N
	kk = k+(N/2-M);
	if kk>N
		kk=kk-N;
	end
	h(k)=buf(kk);
end

buf=g;
for k=1:N
	kk = k-(N/2-M);
	if kk<=0 
		kk=kk+N;
	end
	g(k)=buf(kk);
end

tsp = h';
invtsp = g';