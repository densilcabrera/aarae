%--------------------------------------------------------------------------
% Script to demonstrate the use of blind reverberation time estimator
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% The used algorithm is presented in
% Heinrich W. Löllmann, Emre Yilmaz, Marco Jeub and Peter Vary:
% "An Improved Algorithm for Blind Reverberation Time Estimation"
% International Workshop on Acoustic Echo and Noise Control (IWAENC),
% Tel Aviv, Israel, August 2010.
% (reference avaiable at www.ind.rwth-aachen.de/~bib/loellmann10a)
%
% Note:
% The approach for a fast tracking of changing RTs by means of a second
% histogram is not implemented in this version to further reduce the
% complexity of the algorithm.The parameter setting for this program is not
% identical to the parameter setting used for the simulations examples in
% the reference paper.
%
% The algorithm allows to estimate the RT within a range of 0.2s to 1.2s
% and assumes that source and receiver are not within the critical
% distance. A denoising is not performed by this function and has to be
% done in advance.
%
%--------------------------------------------------------------------------
% Copyright (c) 2012, Heinrich Loellmann and Marco Jeub
% Institute of Communication Systems and Data Processing
% RWTH Aachen University, Germany
% Contact information: loellmann@ind.rwth-aachen.de
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

% clear all
% close all

addpath speech_file
addpath utilities
addpath functions

addpath AIR

% room impulse responses taken from the AIR database
% The complete database is available at www.ind.rwth-aachen.de/~air

example = 5;  % examples 1 to 5 for different rooms (RTs)


% Parameters for AIR Database
% (see help load_air.m for further information)

simpar.rir_type = 1;
simpar.channel = 1;
simpar.azimuth = 0;
simpar.head = 0;
simpar.rir_no = 3;
simpar.room = 4;

switch example
    
    case 5  % Stairway (with dummy head), Source-Mic. distance: 3m
        simpar.rir_type = 1;
        simpar.room = 5;
        simpar.channel = 1;
        simpar.azimuth = 90;
        simpar.head = 1;
        simpar.rir_no = 3;
        simpar.phone_pos = 1;
        
    case 4  % Lecture room, Source-Mic. distance: 5.5m        
        simpar.room = 4;
        
    case 3  % Meeting room, Source-Mic. distance: 1.9m 
        simpar.room = 3;
        
    case 2  % Office, Source-Mic. distance: 3m      
        simpar.room = 2;
        
    case 1  % Low reverberant booth, Source-Mic. distance: 1.5m 
        simpar.room = 1;
        
end

speech_file = 'TSP1min48.wav'; % clean speech file
simpar.fs = 16e3;              % sampling frequency

fig = 1;                       % figure handle for plots

%--------------------------------------------------------------------------
% Generation of reverberant speech
%--------------------------------------------------------------------------
s = readwav(speech_file, simpar.fs);
[ h, h_info ] = load_air(simpar);
s_rev = fftfilt(h,s);

% h_info  % information about the room impulse response

%--------------------------------------------------------------------------
% Parameters for overall processing schemes
% They are not identical to the frame sizes and frame shift used for the
% actual RT estimation to demonstrate the integration into a given block
% processing scheme.
%--------------------------------------------------------------------------
simpar.block_size = round(20e-3 * simpar.fs);  % block length
simpar.overlap = round(simpar.block_size/2);   % overlap

%--------------------------------------------------------------------------
% RT Calculation based on RIR by means of the Schroeder method
%--------------------------------------------------------------------------
rt_schroeder = RT_schroeder(h,simpar.fs,[-5,-35],1,1,1,fig);

%--------------------------------------------------------------------------
% Blind RT Estimation based on reverberant speech signal
%--------------------------------------------------------------------------
tic
[rt_est,rt_est_mean,rt_est_dbg] = ML_RT_estimation(s_rev,simpar);
toc

fprintf('Schroeder T60: %3.2fs\n',rt_schroeder);
fprintf('Estimated T60: %3.2fs\n',rt_est_mean);
%--------------------------------------------------------------------------
% Plot estimated RT and 'true' RT obtained by Schroeder method
%--------------------------------------------------------------------------
fr2sec_idx = linspace(1,length(s_rev)/simpar.fs,length(rt_est));
figure(fig+1)
clf
hold on
plot(fr2sec_idx,rt_est,'-r')
line([0 fr2sec_idx(end)],[rt_schroeder rt_schroeder])
grid on,box on
xlabel('Time [s]'),ylabel('RT [s]');
legend('Estimated T60','Schroeder T60','location','southeast');

%--------------------------------------------------------------------------
