function [bin_coh_3d,bin_coh_2d] = binaural_coherence(d_mic, d_head,...
    r_head, beta, fs,c)
%--------------------------------------------------------------------------
% A Semi-Analytical Model for the Binaural Coherence of Noise Fields
%--------------------------------------------------------------------------
% Syntax:
%   [bin_coh_3d,bin_coh_2d] = binaural_coherence(d_mic, d_head,
%                                                r_head, beta, fs)
%
% Input arguments:
%   d_mic:  Distance between the two microphones (m)
%   d_head: Distance between discs, i.e. head diameter (m)
%   r_head: Radius of disc/head (m)
%                        for r_head=0 the obstacle vanishes and the
%                        calculated coherence is equal to the free-field
%                        coherence without head shadowing
%   beta:  vector containing the values 2*pi*f/c for which the coherence
%          has to be calculated
%   fs:    sampling frequency (Hz); required to determine the number of
%          summations
%
% Output arguments:
%   bin_coh_3d : Binaural coherence (conventional definition, 3 dim.)
%   bin_coh_2d : Binaural coherence (2 dim. definition)
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

%--------------------------------------------------------------------------
% Check input arguments
%--------------------------------------------------------------------------
if nargin < 5
    error('Not enough input arguments');
end;

if ~(d_head <= d_mic)
    error('<d_mic> must be equal to or greater than <d_head>');
end;

%--------------------------------------------------------------------------
% Inits
%--------------------------------------------------------------------------
if nargin < 6
    c   = 340;    % sound velocity, set fixed to 340m/s
end
r_q = 0;

% Number of intervals:
%    The intervals dr and d_phi_r have to be chosen such that the
%    corresponding surface elements are small compared to the sound
%    wavelength. Here, it is proposed that the maximum length of every
%    surface element should be one-tenth of the wavelength.

%   N_r = r_head/delta(r)=r_head/(lambda/10) = 10*r_head*fs/c
N_r = ceil(10*r_head*fs/c);     % -> Eq.(18)

%   N_phi = 2pi/delta(tilde(phi))=2*pi*r_head/(lambda/10))
%          = 10*2*pi*r_head*fs/c
N_phi = ceil(10*2*pi*r_head*fs/c);  % -> Eq.(18)
N_theta = 360/10;   % used for Eq.(22)

% Summation intervals:
dr        = r_head/N_r;     % -> Eq.(18)
d_phi_r   = 2*pi/N_phi;     % -> Eq.(18)
d_theta   = pi/2/N_theta;   % ->Eq.(22)

%--------------------------------------------------------------------------
% Calculate acoustic transfer functions -> Eqs.(18,19,20)
%--------------------------------------------------------------------------
H1 = zeros(N_theta, length(beta));
H2 = zeros(N_theta, length(beta));

for cnt = 1:N_theta,
    theta=d_theta*(cnt-.5);
    
    cos_uq = cos(theta);    % -> Eq.(17.1)
    rq_rm  = r_q + d_mic/2 * cos(theta);    % -> Eq.(14)
    H_1    = exp(-1i*beta*rq_rm);
    H_2    = exp( 1i*beta*rq_rm);
    
    for r = dr:dr:r_head    % -> Eq.(18)
        l_m    = sqrt(r*r+.25*(d_mic-d_head)*(d_mic-d_head));  % -> Eq.(16)
        cos_um = (d_mic-d_head)/(2*l_m);    % -> Eq.(17.2)
        for phi_r = d_phi_r:d_phi_r:2*pi,
            % -> Eq.(15)
            l_q  = r_q + d_head/2*cos(theta) - r*sin(theta)*cos(phi_r);
            H_1  = H_1 -1i*beta/(4*pi).*(exp(-1i*beta*(l_q+l_m))/l_m).* ...
                (cos_uq + (1-1i./(beta*l_m))*cos_um) *r*d_phi_r*dr ;
        end
    end
    H2(cnt,:) = H_2;
    H1(cnt,:) = H_1;
end

%--------------------------------------------------------------------------
% Calculate coherence -> Eq.(22)
%--------------------------------------------------------------------------
numerator = zeros(size(beta));
denumerator0 = zeros(size(beta));

% 2D noise field (without sin(theta))
for cnt=1:N_theta,
    numerator = numerator + 2*real(H2(cnt,:).*conj(H1(cnt,:)));
    denumerator0 = denumerator0 + (H1(cnt,:).*conj(H1(cnt,:)) + ...
        H2(cnt,:).*conj(H2(cnt,:)));
end
bin_coh_2d = ((real(numerator)) ./ abs(denumerator0));

numerator = zeros(size(beta));
denumerator0 = zeros(size(beta));

% 3D noise field (with sin(theta))
for cnt=1:N_theta,
    theta=d_theta*(cnt-.5);
    numerator = numerator + 2*real(H2(cnt,:).*conj(H1(cnt,:))) *sin(theta);
    denumerator0 = denumerator0 + (H1(cnt,:).*conj(H1(cnt,:)) + ...
        H2(cnt,:).*conj(H2(cnt,:))) * sin(theta);
end
bin_coh_3d = ((real(numerator)) ./ abs(denumerator0));
%--------------------------------------------------------------------------