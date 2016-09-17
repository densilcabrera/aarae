function OUT = sphereplotfromHOA(IN,fs,sphere_cover,start_time,end_time,max_order,hif,lof,plottype)
% This function creates a plot of spatial energy distribution from spherical
% harmonic encoded (HOA) multichannel audio input data.
%
% If audio is multiband, then a separate plot is created for each band.
%
% Diffusivity of the soundfield is calculated using the Gover and HOA
% covariance methods (but only when an approximately even distribution of
% directions is used, i.e., not when a rectangular mesh is used).
%
% Currently three plotting methods are supported:
% * Robinson projection (of the sphere onto a 2D surface), which plots
%   using colormap based on the decibel scale
% * Surface plot, which plots a surface using magnitude units and Jet
%   colormap for value (interpolated, with a light source)
% * Lit surface plot, which plots a surface using magnitude units and a
%   Copper colormap for z (more angles are used in this plot than in the 
%   surface plot)
%
% This function uses the HOAToolbox, by Nicolas Epain.
%
% Code by Daniel Jimenez, Luis Miranda and Densil Cabrera
% Version 1.02 (30 August 2014)



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

if abs(size(hoaSignals,2)^0.5 - round(size(hoaSignals,2)^0.5)) >1e-20
    h=warndlg('This audio does not appear to be in HOA format. Unable to analyse with sphereplotfromHOA.','AARAE info','modal');
    uiwait(h)
    OUT = [];
    return
end

if nargin < 9, plottype = 0; end
if nargin < 8, lof = 0; end
if nargin < 7, hif = fs/2; end
if nargin < 6, max_order=round(size(hoaSignals,2).^0.5-1); end
if nargin < 5, end_time = min([0.1, length(hoaSignals)/fs]); end
if nargin < 4, start_time = 0; end
if nargin < 3, sphere_cover = 1; end
if isstruct(IN)
    if nargin < 3
        param = inputdlg({'Spatial density of points (0.25:1)',...
            'Start time [s]'...
            'End time [s]'...
            'Maximum order'...
            'High frequency cutoff (Hz)'...
            'Low frequency cutoff (Hz)'...
            'PLOT TYPE: Trisurf plot of amplitude [0]; Robinson Projection in dB [1]'},... %; 3D colormap surface in amplitude [2]; 3D lit surface in amplitude (slow) [3]'},...
            'Input parameters',1,...
            {num2str(sphere_cover),...
            num2str(start_time),...
            num2str(end_time),...
            num2str(max_order),...
            num2str(fs/2),...
            '0',...
            '1'});
        if isempty(param) || isempty(param{1,1}) || isempty(param{2,1}) || isempty(param{3,1}) || isempty(param{4,1}) || isempty(param{5,1}) || isempty(param{6,1})
            OUT = [];
            return;
        else
            sphere_cover = round(130*str2double(param{1,1}));
            if sphere_cover>130,sphere_cover = 130;end
            if sphere_cover<4,sphere_cover = 4;end
            start_time = str2double(param{2,1});
            
            end_time = str2double(param{3,1});
            
            max_order = str2double(param{4,1});
            
            hif = str2double(param{5,1});
            lof = str2double(param{6,1});
            plottype = str2double(param{7,1});
            if isnan(sphere_cover) || isnan(start_time) || isnan(end_time) || isnan(max_order) || isnan(hif) || isnan(lof), OUT = []; return; end
        end
    end
end

start_sample = round(start_time*fs);
if start_sample == 0, start_sample = 1; end
end_sample = round(end_time*fs);
if end_sample >= length(hoaSignals), end_sample = length(hoaSignals); end
nchans = (max_order+1)^2;
if size(hoaSignals,2)>nchans % delete unused channels
    hoaSignals = hoaSignals(:,1:nchans);
elseif size(hoaSignals,2)<nchans % limit max_order to available channels
    max_order = round(size(hoaSignals,2).^0.5-1);
    if nargin < 2
        warndlg(['Maximum available order for this audio input is ' num2str(max_order) '.'],'AARAE info','modal');
    end
end

            
if (hif < fs/2 && hif > lof) || (lof > 0 && lof < hif)
    % bandlimit the spectrum
    filterorder = 48;
    hoaSignals = bandpass(hoaSignals, lof, hif, filterorder, fs);
end

switch plottype
    case {2, 3}
        % rectangular mesh
        if plottype == 2
            step = pi/(6 * max_order * sphere_cover/130);
        elseif plottype == 3
            step = pi/(12 * max_order * sphere_cover/130);
        end
        elev_for_directplot = -pi/2:step:pi/2; % elevation
        azim_for_directplot = 0:2*step:2*pi; % azimuth
        [azim_for_directplot, elev_for_directplot] =...
            meshgrid(azim_for_directplot, elev_for_directplot);
        [len1, len2] = size(azim_for_directplot);
        
        azim_for_directplot = reshape(azim_for_directplot,numel(azim_for_directplot),1);
        elev_for_directplot = reshape(elev_for_directplot,numel(elev_for_directplot),1);
        numberofdirections = numel(azim_for_directplot);
        
        hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order) ;
        
        Y = SphericalHarmonicMatrix(hoaFmt,azim_for_directplot,elev_for_directplot);
        
        hoa2SpkCfg_for_directplot.filters.gainMatrix = Y;
        direct_sound_HOA = hoaSignals(start_sample:end_sample,:,:);
        bands = size(hoaSignals,3);
        beamsignals_for_directPlot = zeros(length(direct_sound_HOA),numberofdirections,bands);
        
        
        
        for i = 1:size(hoaSignals,2);
            for j = 1:numberofdirections;
                for b = 1:bands
                    beamsignals_for_directPlot(:,j,b) = beamsignals_for_directPlot(:,j,b)+(direct_sound_HOA(:,i,b).*hoa2SpkCfg_for_directplot.filters.gainMatrix(i,j));
                end
            end
        end

        

        
    otherwise
        % case 1
        % approximately even sphere covering
        if plottype == 0
            sphere_cover = 3002;
        end
        [azim_for_directplot,elev_for_directplot] = SphereCovering(sphere_cover);
        numberofdirections = sphere_cover;
        beams_for_directplot = GenerateSpkFmt('sphCoord',[azim_for_directplot elev_for_directplot ones(numberofdirections,1)]);
        
        hoa2SpkOpt_for_directplot = GenerateHoa2SpkOpt('decodType','basic','sampFreq',fs,'filterLength',128,'transFreq',6000,...
            'transWidth',200,'spkDistCor',false);
        
        hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order) ;
        
        hoa2SpkCfg_for_directplot = Hoa2SpkDecodingFilters(hoaFmt,beams_for_directplot,hoa2SpkOpt_for_directplot);
        
        %[~, firstarrivalall_hoaSignals] = max(hoaSignals);
        
        %firstarrival_hoaSignals = min(firstarrivalall_hoaSignals);
        direct_sound_HOA = hoaSignals(start_sample:end_sample,:,:);
        bands = size(hoaSignals,3);
        beamsignals_for_directPlot = zeros(length(direct_sound_HOA),numberofdirections,bands);
        
        
        
        for i = 1:size(hoaSignals,2);
            for j = 1:numberofdirections;
                for b = 1:bands
                    beamsignals_for_directPlot(:,j,b) = beamsignals_for_directPlot(:,j,b)+(direct_sound_HOA(:,i,b).*hoa2SpkCfg_for_directplot.filters.gainMatrix(j,i));
                end
            end
        end
        
end




[Goverdif,HOAdif] = deal(zeros(bands,1));
for b = 1:bands
    switch plottype
        case {2, 3}
            figure('Name','Plot of Magnitude')
            azim_for_directplot = reshape(azim_for_directplot,len1,len2);
            elev_for_directplot = reshape(elev_for_directplot,len1,len2);
            values = reshape(rms(abs(beamsignals_for_directPlot(:,:,b))),len1,len2);
            [x,y,z] = sph2cart(azim_for_directplot,elev_for_directplot,values);
            hold on
            % plot axis poles
            polelength = 1.1*max(max(abs([x;y;z])));
            plot3([0,0],[0,0],[-polelength, polelength],'LineWidth',2,'Color',[0.8,0.8,0]);
            plot3([0,0],[0,0],[-polelength, polelength],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            
            plot3([0,0],[-polelength, polelength],[0,0],'LineWidth',2,'Color',[0.8,0.8,0]);
            plot3([0,0],[-polelength, polelength],[0,0],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','>');
            
            plot3([-polelength, polelength],[0,0],[0,0],'LineWidth',2,'Color',[0.8,0.8,0]);
            plot3([-polelength, polelength],[0,0],[0,0],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','>');
            
            
            if plottype == 2
                surf(x,y,z,values,'FaceColor','interp','FaceLighting','gouraud',...
                'EdgeLighting','gouraud','DiffuseStrength',1,'EdgeAlpha',0.1,...
                'AmbientStrength',0.6);
                camlight right
            elseif plottype ==3
                surfl(x,y,z,'light');
                shading interp
                colormap(copper)
            end
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            set(gca, 'DataAspectRatio', [1 1 1])
            set(gca, 'PlotBoxAspectRatio', [1 1 1])
            % Band title
            if isstruct(IN)
                if isfield(IN,'bandID')
                    title([num2str(IN.bandID(b)), ' Hz'])
                else
                    title(['Band ', num2str(b)])
                end
            else
                title(['Band ', num2str(b)])
            end
            
        case 0

            
            [x1,y1,z1] = sph2cart(azim_for_directplot,...
                elev_for_directplot,...
                ones(size(azim_for_directplot))); % mesh with radius of 1
            K = convhulln([x1,y1,z1]);
            values = rms(abs(beamsignals_for_directPlot(:,:,b)))';
            
            [x,y,z] = sph2cart(azim_for_directplot,...
                elev_for_directplot,...
                values);
            
            
            figure
            
            trisurf(K,x,y,z,values,'FaceColor','interp','FaceLighting','gouraud',...
                'EdgeLighting','gouraud','DiffuseStrength',1,'EdgeAlpha',0.1,...
                'AmbientStrength',0.6);
            camlight right
            %colormap(copper)
            colormap(autumn)
            hold on
            % plot axis poles
            polelength = 1.1*max(max(abs([x;y;z])));
            plot3([0,0],[0,0],[-polelength, polelength],'LineWidth',2,'Color',[0.8,0.8,0]);
            plot3([0,0],[0,0],[-polelength, polelength],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','^');
            
            plot3([0,0],[-polelength, polelength],[0,0],'LineWidth',2,'Color',[0.8,0.8,0]);
            plot3([0,0],[-polelength, polelength],[0,0],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','>');
            
            plot3([-polelength, polelength],[0,0],[0,0],'LineWidth',2,'Color',[0.8,0.8,0]);
            plot3([-polelength, polelength],[0,0],[0,0],'LineWidth',2,'Color',[0,0,0],'LineStyle','--','Marker','>');
            
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            set(gca, 'DataAspectRatio', [1 1 1])
            set(gca, 'PlotBoxAspectRatio', [1 1 1])
            
             Goverdif(b) = GoverDiffuseness(direct_sound_HOA(:,:,b),hoaFmt);
            HOAdif(b) = HoaDiffuseness(direct_sound_HOA(:,:,b),hoaFmt);
            
            % Band title
            if isstruct(IN)
                if isfield(IN,'bandID')
                    title([num2str(IN.bandID(b)), ' Hz, Gover dif. ',...
                        num2str(Goverdif(b)),', HOA dif. ', num2str(HOAdif(b))])
                else
                    title(['Gover dif. ',...
                        num2str(Goverdif(b)),', HOA dif. ', num2str(HOAdif(b))])
                end
            else
                title(['Gover dif. ',...
                    num2str(Goverdif(b)),', HOA dif. ', num2str(HOAdif(b))])
            end
            
        otherwise
            PlotRobinsonProject([azim_for_directplot,elev_for_directplot],mag2db(rms(abs(beamsignals_for_directPlot(:,:,b)))'));
            % Diffuseness calculations and display (only done for
            % approximately even spherical distribution)
            Goverdif(b) = GoverDiffuseness(direct_sound_HOA(:,:,b),hoaFmt);
            HOAdif(b) = HoaDiffuseness(direct_sound_HOA(:,:,b),hoaFmt);
            
            % Band title
            if isstruct(IN)
                if isfield(IN,'bandID')
                    title([num2str(IN.bandID(b)), ' Hz, Gover dif. ',...
                        num2str(Goverdif(b)),', HOA dif. ', num2str(HOAdif(b))])
                else
                    title(['Gover dif. ',...
                        num2str(Goverdif(b)),', HOA dif. ', num2str(HOAdif(b))])
                end
            else
                title(['Gover dif. ',...
                    num2str(Goverdif(b)),', HOA dif. ', num2str(HOAdif(b))])
            end
    end
    
end

if isstruct(IN)
    OUT.funcallback.name = 'sphereplotfromHOA.m';
    OUT.funcallback.inarg = {fs,sphere_cover,start_time,end_time,max_order,hif,lof,plottype};
else
    OUT = hoaSignals;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2014, Daniel Jimenez, Luis Miranda and Densil Cabrera
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