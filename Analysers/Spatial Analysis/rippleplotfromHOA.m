function OUT = rippleplotfromHOA(IN,fs,planechoice,start_time,end_time,max_order,hif,lof,smoothlen,valtype,plottype,domain)
% This function creates a ripple plot from HOA encoded data, usually for
% impulse response analysis.
% 
% The ripple plot is a representation of the impulse response (or HOA
% waveform) in a plane: either the horizontal plane, median plane or
% transverse plane.
%
% If the plot is time domain, then the radius from the centre represents
% time (between the start time and the end time, which are user inputs).
% The start time is in the centre, and end time at the edge.
%
% If the plot is frequency domain, then the radius from the centre
% represents frequency (between the low and high cutoff frequencies, which
% are user inputs). The low cutoff frequency is in the centre, and the high
% cutoff frequency at the edge.
%
% In the time domain, values can be represented as the waveform itself, the
% waveform's envelope, or the waveform's envelope expressed in decibels.
% These can be displayed directly, or their distribution can be displayed
% using 100-bin histograms (rotated around the circle).
%
% In the frequency domain, values can be represented as the spectrum
% magnitude or the spectrum magnitude expressed in decibels (in frequency
% domain, there is no difference between the VALUES input of 0 and 1,
% except that the former cannot be smoothed). The choice of window
% prior to FFT is: rectangular, half Hann (i.e with maximum at the start)
% and full Hann (i.e. with zero at the start).
%
% A smoothing filter can be applied to the plotted data, which can increase
% readability.
%
% When the data size is large, and user parameters are set appropriately,
% the data is downsampled prior to plotting.
%
% This function uses the HOAToolbox, by Nicolas Epain.
%
% Code by Densil Cabrera and Luis Miranda
% Version 1.02 (5 September 2014)

% TO DO:
% work out a good way of displaying radial axis values

if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    hoaSignals = IN.audio;
    fs = IN.fs;
    if isfield(IN,'cal')
        hoaSignals = cal_reset_aarae(hoaSignals,0,IN.cal);
    end
else
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2double(char(fs));
    end
    hoaSignals = IN;
end

[len,chans,bands,dim4,dim5,dim6] = size(hoaSignals);

if abs(chans^0.5 - round(chans^0.5)) >1e-20
    h=warndlg('This audio does not appear to be in HOA format. Unable to analyse with rippleplotfromHOA.','AARAE info','modal');
    uiwait(h)
    OUT = [];
    return
end

maxtime = len/fs;
if nargin < 12, domain = 0; end
if nargin < 11, plottype = 1; end
if nargin < 10, valtype = 40; end
if nargin < 9, smoothlen = round(4*fs/1000); end
if nargin < 8, lof = 0; end
if nargin < 7, hif = fs/2; end
if nargin < 6, max_order=round(chans.^0.5-1); end
if nargin < 5, 
    end_time = 0.1;
    if end_time > maxtime, end_time = maxtime; end
end
if nargin < 4, start_time = 0; end
if nargin < 3, planechoice = 0; end
if isstruct(IN)
    if nargin < 3
        param = inputdlg({'Plane: Horizontal [0]; Median [1]; Transverse [2]',...
            'Start time [s]'...
            'End time [s]'...
            'Maximum order'...
            'High frequency cutoff (Hz)'...
            'Low frequency cutoff (Hz)'...
            'Smoothing filter length (samples)'...
            'VALUES: Raw amplitude [0], Envelope in amplitude [1], Envelope in dB with a range of [2:100] dB, Schroeder reverse integration [-1]'...
            'PLOT TYPE: Mesh [1]; Surfn [2]; Surfl [3]'...
            'Domain or transformation: Time [0]; Frequency domain rectangular window [1]; Frequency domain half Hann window [2]; Frequency domain full Hann window [3]; Histogram of time domain [4];'},...
            'Input parameters',1,...
            {'0',...
            num2str(start_time),...
            num2str(end_time),...
            num2str(max_order),...
            num2str(fs/2),...
            '0',...
            num2str(round(4*fs/1000)),...
            '40',...
            '3',...
            num2str(domain)});
        if isempty(param) || isempty(param{1,1}) || isempty(param{2,1}) || isempty(param{3,1}) || isempty(param{4,1}) || isempty(param{5,1}) || isempty(param{6,1})
            OUT = [];
            return;
        else
            planechoice = str2double(param{1,1});
            start_time = str2double(param{2,1});
            
            end_time = str2double(param{3,1});
            
            max_order = str2double(param{4,1});
            
            hif = str2double(param{5,1});
            lof = str2double(param{6,1});
            smoothlen = round(str2double(param{7,1}));
            valtype = str2double(param{8,1});
            
            plottype = str2double(param{9,1});
            domain= str2double(param{10,1});
            if isnan(start_time) || isnan(end_time) || isnan(max_order) || isnan(hif) || isnan(lof), OUT = []; return; end
        end
    end
end

start_sample = round(start_time*fs);
if start_sample == 0, start_sample = 1; end

end_sample = round(end_time*fs);
if end_sample >= length(hoaSignals), end_sample = length(hoaSignals); end

nchans = (max_order+1)^2;
if chans>nchans % delete unused channels
    hoaSignals = hoaSignals(:,1:nchans);
elseif chans<nchans % limit max_order to available channels
    max_order = round(chans.^0.5-1);
    if nargin < 2
        warndlg(['Maximum available order for this audio input is ' num2str(max_order) '.'],'AARAE info','modal');
    end
end



if (hif < fs/2 && hif > lof) || (lof > 0 && lof < hif)
    % bandlimit the spectrum
    filterorder = 48;
    hoaSignals = bandpass(hoaSignals, lof, hif, filterorder, fs);
end






hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order) ;


% 1 deg resolution
step = 2*pi/360;
switch planechoice
    case 1 % median plane
        elev_for_directplot1 = -pi/2:step:pi/2;
        elev_for_directplot2 = pi/2-step:-step:-pi/2;
        azim_for_directplot = [zeros(size(elev_for_directplot1)),...
            pi*ones(size(elev_for_directplot2))];
        elev_for_directplot = [elev_for_directplot1,elev_for_directplot2];
    case 2 % transverse plane
        elev_for_directplot1 = -pi/2:step:pi/2;
        elev_for_directplot2 = pi/2-step:-step:-pi/2;
        azim_for_directplot = [pi/2*ones(size(elev_for_directplot1)),...
            3*pi/2*ones(size(elev_for_directplot2))];
        elev_for_directplot = [elev_for_directplot1,elev_for_directplot2];
    otherwise % horizontal plane
        azim_for_directplot = 0:step:2*pi; % azimuth
        elev_for_directplot = zeros(size(azim_for_directplot)); % elevation
end
numberofdirections = numel(azim_for_directplot);

Y = SphericalHarmonicMatrix(hoaFmt,azim_for_directplot,elev_for_directplot);

direct_sound_HOA = hoaSignals(start_sample:end_sample,:,:);
len = length(direct_sound_HOA);
beamsignals = zeros(length(direct_sound_HOA),numberofdirections,bands);

for i = 1:chans;
    for j = 1:numberofdirections;
        for b = 1:bands
            beamsignals(:,j,b) = beamsignals(:,j,b)+(direct_sound_HOA(:,i,b).*Y(i,j));
        end
    end
end

if domain == 1 || domain == 2 || domain == 3
    
    if mod(len,2)==1
        len = len+1;
        beamsignals = [beamsignals;zeros(1,numberofdirections,bands,dim4,dim5,dim6)];
    end
    
    if domain == 2
        window = hann(len*2);
        window = window((end/2+1):end);
    elseif domain == 3
        window = hann(len);
    else
        window = ones(len,1);
    end
    % normalise window
    window = window ./ sum(window);
    
    beamsignals = abs(fft(beamsignals .* repmat(window,[1,numberofdirections,bands,dim4,dim5,dim6]),len));
        % list of fft component frequencies
    f = ((1:len)'-1) * fs / len;
    
    % index of low cut-off
    indlo = find(abs(f(1:end/2)-lof) == min(abs(f(1:end/2)-lof)),1,'first');
    
    % index of high cut-off
    indhi = find(abs(f(1:end/2)-hif) == min(abs(f(1:end/2)-hif)),1,'first');
    beamsignals = beamsignals(indlo:indhi,:,:);
end


switch valtype
    case 0
        % resample the wave is hif is low enough - to avoid plotting
        % excessively large data (but still oversampling by a factor of 2)
        if round(fs/(hif*4)) > 1 && length(beamsignals)>5000 && domain==0
            beamsignals = resample(beamsignals,1,round(fs/(hif*4)));
        end
        % no smoothing here for freq domain
    case 1
        beamsignals = abs(beamsignals);
        if smoothlen > 0 && smoothlen < length(beamsignals)/4
            bcoef = hann(smoothlen) ./ sum(hann(smoothlen));
            beamsignals = filtfilt(bcoef,1,beamsignals);
            if smoothlen/fs >=0.002  && length(beamsignals)>5000 && domain ==0
                beamsignals = resample(beamsignals,1,round(1000*smoothlen/fs));
            end
        end
    case -1 % Schroeder reverse integration
        beamsignals = 10*log10(flipdim(cumsum(flipdim(beamsignals,1).^2),1));
        dBrange = 60;
        for b = 1:bands
            beamsignals(beamsignals(:,:,b) < max(max(beamsignals(:,:,b)))-dBrange)...
                = max(max(beamsignals(:,:,b)))-dBrange;
        end
        case -2 % differenced Schroeder reverse integration
        beamsignals = 10*log10(flipdim(cumsum(flipdim(beamsignals,1).^2),1));
        dBrange = 60;
        for b = 1:bands
            beamsignals(beamsignals(:,:,b) < max(max(beamsignals(:,:,b)))-dBrange)...
                = max(max(beamsignals(:,:,b)))-dBrange;
        end
        beamsignals =diff(beamsignals);
    otherwise
        % values are in dB, with valtype specifying the range of the data
        beamsignals = beamsignals.^2;
        if smoothlen > 0  && smoothlen < length(beamsignals)/4
            bcoef = hann(smoothlen) ./ sum(hann(smoothlen));
            beamsignals = filtfilt(bcoef,1,beamsignals);
            if smoothlen/fs >=0.002 && length(beamsignals)>5000 && domain ==0
                beamsignals = resample(beamsignals,1,round(1000*smoothlen/fs));
            end
            beamsignals(beamsignals<=0) = 1e-99;
        end
        beamsignals = 10*log10(beamsignals);
        % apply level range
        dBrange = valtype; % (just for clarity)
        for b = 1:bands
            beamsignals(beamsignals(:,:,b) < max(max(beamsignals(:,:,b)))-dBrange)...
                = max(max(beamsignals(:,:,b)))-dBrange;
        end
end

if domain == 4
    nbins = 100;
    values = zeros(nbins,numberofdirections,bands);
    for b = 1:bands
        values(:,:,b) = hist(beamsignals(:,:,b),nbins);
    end
    beamsignals = values;
    clear values
end


for b = 1:bands
    figure('color','white');
    switch plottype
        % case 1 is taken care of by 'otherwise'
        case 2
            polarplot3d(beamsignals(:,:,b),'plottype','surfn','AxisLocation','bottom');
            hold on
            if valtype == 0 || valtype == 1 || domain == 4
                    valrange = max(max(beamsignals(:,:,b)))-min(min(beamsignals(:,:,b)));
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-0.2*valrange,max(max(beamsignals(:,:,b)))+0.2*valrange],'LineWidth',2,'Color',[0.8,0.8,0]);
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-0.2*valrange,max(max(beamsignals(:,:,b)))+0.2*valrange],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            else
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-10,max(max(beamsignals(:,:,b)))+10],'LineWidth',2,'Color',[0.8,0.8,0]);
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-10,max(max(beamsignals(:,:,b)))+10],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            end
        case 3
            [x,y,z] = polarplot3d(beamsignals(:,:,b),'plottype','off');
            hold on
            surfl(x,y,z,'light');
            shading interp
            colormap(jet)
            grid on
            hold on
            if valtype == 0 || valtype == 1 || domain == 4
                    valrange = max(max(z))-min(min(z));
                    plot3([0,0],[0,0],[min(min(z))-0.2*valrange,max(max(z))+0.2*valrange],'LineWidth',2,'Color',[0.8,0.8,0]);
                    plot3([0,0],[0,0],[min(min(z))-0.2*valrange,max(max(z))+0.2*valrange],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            else
                    plot3([0,0],[0,0],[min(min(z))-10,max(max(z))+10],'LineWidth',2,'Color',[0.8,0.8,0]);
                    plot3([0,0],[0,0],[min(min(z))-10,max(max(z))+10],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            end
        case 4
            [x,y,z] = polarplot3d(beamsignals(:,:,b),'plottype','off');
            hold on
            contourf(x,y,z);
        otherwise
            % mesh
            polarplot3d(beamsignals(:,:,b),'plottype','mesh','AxisLocation','bottom');
            hold on
            if valtype == 0 || valtype == 1 || domain == 4
                    valrange = max(max(beamsignals(:,:,b)))-min(min(beamsignals(:,:,b)));
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-0.2*valrange,max(max(beamsignals(:,:,b)))+0.2*valrange],'LineWidth',2,'Color',[0.8,0.8,0]);
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-0.2*valrange,max(max(beamsignals(:,:,b)))+0.2*valrange],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            else
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-10,max(max(beamsignals(:,:,b)))+10],'LineWidth',2,'Color',[0.8,0.8,0]);
                    plot3([0,0],[0,0],[min(min(beamsignals(:,:,b)))-10,max(max(beamsignals(:,:,b)))+10],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            end
    end
    
    titlestring = [num2str(start_time), '-', num2str(end_time),' s, '];
    
    switch planechoice
        case 1
            xlabel('x')
            ylabel('z')
            titlestring = [titlestring, 'Median Plane'];
        case 2
            xlabel('y')
            ylabel('z')
            titlestring = [titlestring, 'Transverse Plane'];
        otherwise
            xlabel('x')
            ylabel('y')
            titlestring = [titlestring, 'Horizontal Plane'];
    end
    
    switch domain
        case {1,2,3}
            titlestring = [titlestring, ', Frequency domain: ',...
            num2str(lof), ' Hz - ', num2str(hif), ' Hz (low f at centre) '];
        case 4
            titlestring = [titlestring,...
                ', Rotated histogram of time series values '];
        otherwise
            titlestring = [titlestring, ', Time domain (start time at centre) '];
    end
    
    % Band title
    if isstruct(IN)
        if isfield(IN,'bandID')
            title([num2str(IN.bandID(b)), ' Hz, ', titlestring])
        else
            title(['Band ', num2str(b),', ', titlestring])
        end
    else
        title(['Band ', num2str(b),', ', titlestring])
    end
    
    if domain ~= 4
    switch valtype
        case {0,1}
            zlabel('Amplitude')
        case -1
            zlabel('Reverse-integrated decay (dB)')
        otherwise
            zlabel('Level (dB)')
    end
    else
        zlabel('Histogram count')
    end
    
    colorbar;
    
end

if isstruct(IN)
    OUT.funcallback.name = 'rippleplotfromHOA.m';
    OUT.funcallback.inarg = {fs,planechoice,start_time,end_time,max_order,hif,lof,smoothlen,valtype,plottype,domain};
else
    OUT = hoaSignals;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2014, Densil Cabrera and Luis Miranda
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%