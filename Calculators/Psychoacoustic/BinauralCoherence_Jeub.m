function OUT = BinauralCoherence_Jeub(d_head,r_head,d_mic,c)
% This function ports Marco Jeub's Binaural Coherence of Noise Fields
% calculator to AARAE. The following is a quote from the Matlab file
% exchange site:
%
% A Matlab reference implementation of a novel semi-analytical signal
% processing model for the binaural coherence of homogeneous isotropic
% noise fields is presented.
%
% The model is derived from a simplified geometrical model of the human head,
% where the shadowing between the left and right ear is modeled by two
% non-reflecting circular plates. Based on Kirchhoff?s diffraction theory,
% it is shown how the corresponding coherence is calculated. This model can
% be used as part of various binaural signal processing algorithms, such as
% speech enhancement for digital hearing aids or binaural speech
% transmission systems. In experiments using an artificial head in a highly
% reverberant environment, it is confirmed that the proposed theoretical
% model shows a good match with the coherence obtained from measurements.
%
% Within this implmentation, arbitrary dimensions for head and microphone
% distances can be employed and no acoustic measurements are required.
%
% Related Paper:
% M. Jeub, M. Dörbecker, P. Vary: "A Semi-Analytical Model for the Binaural
% Coherence of Noise Fields", IEEE Signal Processing Letters, Vol.18, No.3,
% March 2011 (DOI 10.1109/LSP.2011.2108284)

if nargin == 0 % If the function is called within the AARAE environment it
    % won't have any input arguments, this is when the inputdlg
    % function becomes useful.
    
    param = inputdlg({'Head width diameter (m)';...
        'Head median plane diameter (m)';...
        'Distance between microphones (m)';...
        'Speed of sound (m/s)'},...
        'Settings',... % This is the dialog window title.
        [1 60],...
        {'0.17';'0.17';'0.17';'340'});
    
    param = str2num(char(param));
    
    if length(param) < 4, param = []; end
    if ~isempty(param)
        d_head = param(1);
        r_head = param(2)*0.5;
        d_mic  = param(3);
        c = param(4);
    end
else
    param = [];
end

if ~isempty(param) || nargin ~= 0
    fs = 16e3;      % Sampling frequency (Hz)
    
    N = 256;        % Number of frequency bins
    frq_vec = linspace(1,fs/2,N);   % Frequency vector
    beta = 2*pi*frq_vec/c; % Vector containing the values 2*pi*f/c for which
    % the coherence has to be calculated
    
    
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
    figure('Name','Magnitude Squared Coherence')
    hold on;
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
    
    % Create output structure
%     if nargin == 0
%         OUT.f = frq_vec;
%         OUT.bin_coh_3d = bin_coh_3d.^2;
%         OUT.bin_coh_2d = bin_coh_2d.^2;
%         OUT.coh_3d = coh_3d.^2;
%         OUT.coh_2d = coh_2d.^2;
%     end
    

doresultleaf([bin_coh_3d.^2',bin_coh_2d.^2',coh_3d.^2',coh_2d.^2'],'BinauralCoherence',{'Frequency','CoherenceType'},...
                 'Frequency',      frq_vec,                'Hz',           true,...
                 'CoherenceType',  [{'bin_coh_3d'},{'bin_coh_2d'},{'coh_3d'},{'coh_2d'}],           'categorical', [],...
                 'name','BinCohereJeub');
    %
    % Input arguments:
    % #1: Your data variable. It can be multidimensional, make sure you
    %     specify what each dimension is.
    % #2: What is your data variable representing? is it level? is it
    %     reverb time? make sure you label it appropriately and assign
    %     units to it, this second argument is a single string.
    % #3: This is a cell array where each cell contains the name of each
    %     dimension.
    %
    % #4: Name of your first dimension. (String)
    % #5: Matrix that matches your first dimension, in this case time.
    % #6: Units of your first dimension. (String)
    % #7: Can this dimension be used as a category? (true, false, [])
    %
    % Replicate arguments 4 to 7 for as many dimensions as your data
    % variable has.
    %
    % The second last input argument is the string 'name', this helps the
    % function identify that the last input argument is the name that will
    % be displayed in AARAEs categorical tree under the results leaf.



    OUT.funcallback.name = 'BinauralCoherence_Jeub.m';
    OUT.funcallback.inarg = {d_head,r_head,d_mic,c};
else
    OUT = [];
end

%**************************************************************************
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

% Adapted by Densil Cabrera for AARAE