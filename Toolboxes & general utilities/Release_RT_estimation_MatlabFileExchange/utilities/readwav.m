%--------------------------------------------------------------------------
% Read WAV-file and resample to given sampling frequency
%--------------------------------------------------------------------------
% function [x,bits] = readwav(wav_file,fs,do_norm,scaling)
%
% Input
%    wav_file:   Name and path of file
%    fs:         desired sampling frequency in Hz
%    do_norm     '1': normalize input speech to 'scaling' (default)
%    scaling     scaling for normalization (default: 0.99)
%
% Output
%    x:          Audio Matrix (Channel x Data)
%    bits:       number of bits
%--------------------------------------------------------------------------
%
% Copyright (c) 2012, Marco Jeub
% Institute of Communication Systems and Data Processing
% RWTH Aachen University, Germany
% Contact information: jeub@ind.rwth-aachen.de
%
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the RWTH Aachen University nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function [x,bits] = readwav(wav_file,fs,do_norm,scaling)
%--------------------------------------------------------------------------
if nargin < 4
    scaling = 0.99;
end
if nargin < 3
    do_norm = 1;
end;

[x, fs_in, bits] = wavread(wav_file);
% Resample if necessary
if fs ~= fs_in
    x=resample(x,fs,fs_in);
    if fs >= fs_in
        disp('readwav: Given sampling frequency is higher than the fs of the file, Proceed with upsampling');
    end;      
end;
x = x'; % transpose
% ensure even number of samples
even_size = 2*floor(length(x)/2);
x = x(:,1:even_size);  

% Normalization
if do_norm == 1
    x_max = max(abs(x(:)));
    x=x./x_max.*scaling;
end;
   
%--------------------------------------------------------------------------

