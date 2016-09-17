function outSig = fftfiltmat(FltMat,inpSig,nmbFft)
%FFTFILTMAT Overlap-add method of FIR filtering with a matrix of filters.
%
%  Y = FFTFILT(F,X) filters the input signals in matrix X with the matrix
%  of FIR filters F using the overlap/add method, using internal parameters
%  (FFT size and block length) which guarantee efficient execution.
%
%  Y = FFTFILT(F,X,N) allows you to have some control over the internal
%  parameters, by using an FFT of at least N points.
% 
%  F must be a 3D array of dimension nmbTap x nmbOut x nmbInp, where
%  nmbInp = size(X,2) and nmbOut = size(Y,2).
%
%  See also FFTFILT 

%  N.Epain, 2010.

% Dimensions
nmbSmp = size(inpSig,1) ;
nmbInp = size(inpSig,2) ;
nmbOut = size(FltMat,2) ;

% Tests
if nmbInp ~= size(FltMat,3)
    error('The number of cols in inpSig must be equal to size(FltMat,3)') ;
end

% Initialise output signals
outSig = zeros(nmbSmp,nmbOut) ;

% Process the signals
if nargin < 3
    for I = 1 : nmbInp
        outSig = outSig + fftfilt(FltMat(:,:,I),inpSig(:,I)) ;
    end   
else
    for I = 1 : nmbInp
        outSig = outSig + fftfilt(FltMat(:,:,I),inpSig(:,I),nmbFft) ;
    end
end
