function [RT,EDC_log,t_EDC,noise_floor] = RT_schroeder(h,fs,region,method,plot_on,delay_comp, fig)
%--------------------------------------------------------------------------
% Measuring the Reverberation Time using the Schroeder Method
%--------------------------------------------------------------------------
%
% Input:      h:  inpulse response
%             fs: sampling frequency in [Hz]
%             Optional:
%                region(2): Region in the EDC where the RT is computed [dB]
%                           default [-5 -35]
%                method:    detection method
%                            1: least square fitting (default)
%                            2: line between region(1) and (2)
%                plot_on:   plots the detection mechanism
%                            0: no plot (default), 1: plot
%                delay_comp: 0: no delay compensation (default)
%                            1: compensate sound propagation delay
%             fig : figure handle for plot (default = 1)
%
% Output:     RT: reveberation time in [s]
%             EDC_log: energy decay curve (log normalized)
%             t_EDC: time vector for EDC curve
%             noise_floor: noise floor
%             (at time instance the fitting curve reaches -60dB)
%--------------------------------------------------------------------------
%
% Copyright (c) 2012, Marco Jeub and Heinrich Loellmann
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

if nargin <7, fig = 1; end
if nargin < 6, delay_comp = 0; end
if nargin < 5, plot_on = 0; end
if nargin < 4, method = 1; end
if nargin < 3, region = [-5,-35]; end

if nargin < 2; error('at least the RIR and its sampling frequency must be given');end;
%--------------------------------------------------------------------------

h_length=length(h);
%--------------------------------------------------------------------------
% Compensate sound propagation
% (exclude parts of the RIR before the direct path)
%--------------------------------------------------------------------------
% 
if delay_comp == 1
    [~,prop_delay]=max(h);
    h(1:h_length-prop_delay+1)=h(prop_delay:h_length);
    h_length=length(h);
end

%--------------------------------------------------------------------------
% Energy decay curve
%--------------------------------------------------------------------------
h_sq=h.^2;
h_sq = fliplr(h_sq);
EDC = cumsum(h_sq);
EDC = fliplr(EDC);

% normalize to '1'
EDC_norm=EDC./max(EDC);

%--------------------------------------------------------------------------
% Estimate the reverberation time
%--------------------------------------------------------------------------
if method == 1  % least square fitting
    % first value of the EDC decaying 60dB (10^-6)
    EDC_log = 10*log10(EDC_norm);
    %EDC_60 = find (EDC_log <= -60, 1, 'first');
    EDC_reg1  = find (EDC_log <= region(1), 1, 'first');
    EDC_reg2 = find (EDC_log <= region(2), 1, 'first');
    
    EDC_reg12 = EDC_log(EDC_reg1:EDC_reg2);
    x=1:length(EDC_reg12);
    p = polyfit(x,EDC_reg12,1); % linear least square fitting
    
    x=1:length(EDC_reg12);
    y=p(1)*x+p(2);
    
    y0=y(1)-p(1)*EDC_reg1;
    
    % intersection of polyfit line with -60dB
    x_rt = (-60-y0)/p(1);
    
    RT=x_rt/fs;   % Reverberation time in [s]
    
    % fitting line from 0 to -60dB
    x=1:x_rt;
    y = p(1)*x+y0;
end;
%--------------------------------------------------------------------------
if method == 2  % simple line between the 2 thresholds in region
    % first value of the EDC decaying 60dB (10^-6)
    EDC_log = 10*log10(EDC_norm);
    %EDC_60 = find (EDC_log <= -60, 1, 'first');
    EDC_reg1  = find (EDC_log <= region(1), 1, 'first');
    EDC_reg2 = find (EDC_log <= region(2), 1, 'first');
    
    % Line Slope between the Points given by region(1) and region(2)
    m=( EDC_log(EDC_reg2) - EDC_log(EDC_reg1)) / ( EDC_reg2 - EDC_reg1);
    % Line
    x=1:h_length;
    y=EDC_log(EDC_reg1)+m*(x-EDC_reg1);
    % point of intersection at -60dB
    x2=(-60-y(1)) / m;
    
    RT=x2/fs;   % Reverberation time in [s]
   
end

%--------------------------------------------------------------------------
% Noise floor
%--------------------------------------------------------------------------
% time instance the fitting curve 'y' reaches -60dB
% noise_floor = EDC_log(length(y));
noise_floor = [];
%--------------------------------------------------------------------------
% Optional Plots
%--------------------------------------------------------------------------
t_EDC=linspace(0,h_length/fs,h_length);
if plot_on == 1
    linewidth = 2;
    fontsize = 14;
    y_length=length(y);
    t_y=linspace(0,y_length/fs,y_length);
    
    figure(fig)
    clf
    plot(t_EDC,EDC_log,'LineWidth',linewidth);hold on;
    plot(t_y,y,'-r','LineWidth',linewidth);
    line( [0 y_length/fs],[region(1)  region(1)],'Color','black','LineStyle','--','LineWidth',linewidth);
    line( [0 y_length/fs],[region(2)  region(2)],'Color','black','LineStyle','--','LineWidth',linewidth);
    line( [0 y_length/fs],[-60 -60],'Color','black','LineStyle','--','LineWidth',linewidth);
    axis([0 y_length/fs -65 0]);
    set(gca,'fontsize',fontsize)
    xlabel('Time in (s)','FontSize',fontsize);ylabel('Energy Decay in (dB)','FontSize',fontsize);title('Reverberation time estimation - Schroeder method','FontSize',fontsize);
    text(.6*RT ,-10,['T_6_0 (s) = ',num2str(RT)],'FontSize',fontsize,'color','k');
   
end
%--------------------------------------------------------------------------