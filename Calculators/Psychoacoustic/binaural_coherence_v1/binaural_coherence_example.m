%--------------------------------------------------------------------------
% A Semi-Analytical Model for the Binaural Coherence of Noise Fields
%--------------------------------------------------------------------------
% -> Example simulation script
%--------------------------------------------------------------------------
% Related paper:
% M. Jeub, M. Dörbecker, P. Vary: "A Semi-Analytical Model for the Binaural
% Coherence of Noise Fields", IEEE Signal Processing Letters, Vol.18, No.3,
% March 2011 (DOI 10.1109/LSP.2011.2108284)
%--------------------------------------------------------------------------
% Copyright (c) 2011, Marco Jeub
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
% Version 1.0
%--------------------------------------------------------------------------
clear all;close all;clc

%--------------------------------------------------------------------------
% Simulation settings
%--------------------------------------------------------------------------
fs = 16e3;      % Sampling frequency (Hz)
c = 340;        % Sound velocity (m/s)
N = 256;        % Number of frequency bins
frq_vec = linspace(1,fs/2,N);   % Frequency vector
beta = 2*pi*frq_vec/c; % Vector containing the values 2*pi*f/c for which
% the coherence has to be calculated

d_mic  = 0.17;    % Distance between the two microphones (m)
d_head = 0.17;    % Distance between discs, i.e. head diameter (m)
r_head = 0.075;   % Radius of the obstacle between the microphones (m)
%                        for r_head=0 the obstacle vanishes and the
%                        calculated coherence is equal to the free-field
%                        coherence without head shadowing
%--------------------------------------------------------------------------
% Calculate coherence with head shadowing
%--------------------------------------------------------------------------
[bin_coh_3d,bin_coh_2d] = binaural_coherence(d_mic,d_head,r_head,beta,fs);

%--------------------------------------------------------------------------
% Calculate coherence w/o head shadowing (free-field)
%--------------------------------------------------------------------------
% 3D noise field (Eq.2)
coh_3d =  sinc(2*frq_vec*d_mic/c);
% 2D noise field (Eq.3)
nu = 0;
coh_2d = besselj(nu,2*pi*frq_vec*d_mic/c);

%--------------------------------------------------------------------------
% Plot results (magnitude squared coherence)
%--------------------------------------------------------------------------
figure,hold on;
plot(frq_vec,bin_coh_3d.^2,'-r');
plot(frq_vec,bin_coh_2d.^2,'-k');
plot(frq_vec,coh_3d.^2,'--r');
plot(frq_vec,coh_2d.^2,'--k');
legend('3D model with head','2D model with head',...
    '3D model (free-field)','2D model (free-field)');
xlabel('Frequency (Hz)');
ylabel('Coherence');
axis([0 4000 0 1])
grid on,box on
%--------------------------------------------------------------------------