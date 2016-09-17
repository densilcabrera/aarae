function OUT = sphereplotfromMics(IN,percentilethreshold,valuescale,cartgrid,sphgrid)
% This function creates a spherical polar plot based on measurements from
% micrphones surrounding a sound source. Plotting is done by meshing the
% microphone positions (usually with some interpolation - see below). The
% resulting 'mango plots' allow the directional radiation patterns of sound
% sources to be examined.
%
% The input must have chanID data that identifies the position of each
% microphone. These chanIDs can be created using AARAE's makechanID
% function. However, more likely, the chanIDs are created as part of the
% processor MicArraysSpecialProcess (e.g., when making sequential
% measurements using a turntable). Note that you will need to define a
% microphone array (as a structure within the MicArraysSpecial folder
% within Processors/Beamforming) to use MicArraysSpecialProcess.
%
% Input can be multiband, and it is assumed that the intput is limited to 3
% dimensions (time, channels, bands).
%
% While it is ideal for the microphones to be evenly distributed over a
% sphere, in practice that is often not possible. Therefore this function
% provides some extra points in the mesh by interpolation (so that large
% separations between points are not joined by a straight line, but instead
% by a bent line following an interpolated radius).
%
% Code by Densil Cabrera
% version 1 (26 September 2014)

if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    audio = IN.audio;
    fs = IN.fs;
    if isfield(IN,'chanID')
        [coords,format] = readchanID(IN.chanID);
        if format ~= 2 && format ~= 3 && format ~= 4
            warndlg('ChanID is not in the correct format for this function')
            OUT = [];
            return
        end
    end
    % reset cal to 0 dB if it exists (this should already have been done by
    % the processor)
    if isfield(IN,'cal')
        audio = cal_reset_aarae(audio,94,IN.cal);
    end
else
    OUT = [];
    return
end






if nargin == 1
    param = inputdlg({'Percentile mesh distance threshold for interpolation (100 means no interpolation)';...
        'Scale: Amplitude (1), squared amplitude (2) or decibels (3)';...
        'Cartesian grid (0 | 1)';...
        'Spherical grid: None (0), outer sphere (1), circles (2) outer sphere and circles (3)';},...
        'Mango plot settings',... % This is the dialog window title.
        [1 60],...
        {'90';'0';'1';'2'}); % preset answers
    
    param = str2num(char(param));
    
    if length(param) < 2, param = []; end
    if ~isempty(param) % If they have, you can then assign the dialog's
        % inputs to your function's input parameters.
        percentilethreshold = param(1);
        valuescale = round(param(2));
        cartgrid = round(param(3));
        sphgrid = round(param(4));
    else
        OUT = [];
        return
    end
end
percentilethreshold = abs(percentilethreshold);
valuescale = round(valuescale);
if valuescale > 3, valuescale = 3; end
if valuescale < 1, valuescale = 1; end





% values to plot
values = permute(rms(audio),[2,3,1]);
if valuescale == 2
    values = values.^2;
elseif valuescale == 3
    values = 20*log10(values);
    values = values - max(max(values)) + 40;
    values(values<0) = 0;
end





% spherical coordinates
[az,el] = cart2sph(coords(:,1),coords(:,2),coords(:,3));

% coordinates on a unit sphere
[x,y,z] = sph2cart(az,el,ones(length(az),1));

% indices for meshing
K = convhulln([x,y,z]);






% Find points that are far enough apart to require interpolation of values,
% and interpolate, using a maximum of two iterations
if percentilethreshold < 100
    for it = 1:2
        originalxyzlen = length(x);
        
        % find 'distances' between adjacent points
        d1 = ((x(K(:,1)) - x(K(:,2))).^2 ...
            + (y(K(:,1)) - y(K(:,2))).^2 ...
            + (z(K(:,1)) - z(K(:,2))).^2).^0.5;
        
        d2 = ((x(K(:,3)) - x(K(:,2))).^2 ...
            + (y(K(:,3)) - y(K(:,2))).^2 ...
            + (z(K(:,3)) - z(K(:,2))).^2).^0.5;
        
        d3 = ((x(K(:,1)) - x(K(:,3))).^2 ...
            + (y(K(:,1)) - y(K(:,3))).^2 ...
            + (z(K(:,1)) - z(K(:,3))).^2).^0.5;
        
        if it == 1
            % percentile threshold (only calculate on first iteration)
            dthresh = prctile([d1;d2;d3],percentilethreshold);
        end
        d1ind = find(d1 > dthresh);
        d2ind = find(d2 > dthresh);
        d3ind = find(d3 > dthresh);
        
        % interpolate between points
        for n = 1:length(d1ind)
            x = [x;mean([x(K(d1ind(n),1)),x(K(d1ind(n),2))])];
            y = [y;mean([y(K(d1ind(n),1)),y(K(d1ind(n),2))])];
            z = [z;mean([z(K(d1ind(n),1)),z(K(d1ind(n),2))])];
        end
        for n = 1:length(d2ind)
            x = [x;mean([x(K(d2ind(n),3)),x(K(d2ind(n),2))])];
            y = [y;mean([y(K(d2ind(n),3)),y(K(d2ind(n),2))])];
            z = [z;mean([z(K(d2ind(n),3)),z(K(d2ind(n),2))])];
        end
        for n = 1:length(d3ind)
            x = [x;mean([x(K(d3ind(n),1)),x(K(d3ind(n),3))])];
            y = [y;mean([y(K(d3ind(n),1)),y(K(d3ind(n),3))])];
            z = [z;mean([z(K(d3ind(n),1)),z(K(d3ind(n),3))])];
        end
        
        % recreate the spherical mesh with the interpolated points
        [az,el] = cart2sph(x,y,z);
        [x,y,z] = sph2cart(az,el,ones(length(az),1));
        K = convhulln([x,y,z]);
        
        
        % extend values with the interpolated values
        if originalxyzlen < length(x)
            values = [values;...
                zeros((length(x)-length(values)),size(values,2))];
            for n = (originalxyzlen+1):length(x)
                % for each new coordinate, find all of its neighbours
                [r,c] = find(K == n);
                % find the distances between it and known values
                neighborvals = NaN(length(r),2,size(values,2));
                neighbordist = NaN(length(r),2);
                for m = 1:length(r)
                    if K(r(m),1+mod(c(m),3)) <= originalxyzlen+n-1
                        neighborvals(m,1,:) = values(K(r(m),1+mod(c(m),3)),:);
                        neighbordist(m,1) = ((x(n)-x(K(r(m),1+mod(c(m),3)))).^2 ...
                            + (y(n)-y(K(r(m),1+mod(c(m),3)))).^2 ...
                            + (z(n)-z(K(r(m),1+mod(c(m),3)))).^2).^0.5;
                    end
                    if K(r(m),1+mod(c(m)+1,3)) <= originalxyzlen+n-1
                        neighborvals(m,2,:) = values(K(r(m),1+mod(c(m)+1,3)),:);
                        neighbordist(m,2) = ((x(n)-x(K(r(m),1+mod(c(m)+1,3)))).^2 ...
                            + (y(n)-y(K(r(m),1+mod(c(m)+1,3)))).^2 ...
                            + (z(n)-z(K(r(m),1+mod(c(m)+1,3)))).^2).^0.5;
                    end
                end
                neighborvals = permute([neighborvals(:,1,:);neighborvals(:,2,:)],[1,3,2]);
                neighbordist = [neighbordist(:,1);neighbordist(:,2)];
                idx = ~isnan(neighbordist);
                neighbordist = neighbordist(idx);
                neighborvals = neighborvals(idx,:);
                idx = neighborvals(:,1) > 0;
                neighbordist = neighbordist(idx);
                neighborvals = neighborvals(idx,:);
                neighborinvdist = 1./neighbordist;
                %  values(n,:) = mean(neighborvals);
                switch valuescale
                    case 1
                        values(n,:) = (sum(neighborvals.^2 .* repmat(neighborinvdist,[1,size(values,2)]))...
                            ./sum(repmat(neighborinvdist,[1,size(values,2)]))).^0.5;
                    case 2
                        values(n,:) = (sum(neighborvals .* repmat(neighborinvdist,[1,size(values,2)]))...
                            ./sum(repmat(neighborinvdist,[1,size(values,2)])));
                    case 3
                        values(n,:) = 10*log10(sum(10.^(neighborvals./10) .* repmat(neighborinvdist,[1,size(values,2)]))...
                            ./sum(repmat(neighborinvdist,[1,size(values,2)])));
                        
                end
            end
        end
    end
end




% plotting
for b = 1:size(values,2)
    % Cartesian coordinates with values
    [x,y,z] = sph2cart(az,el,values(:,b));
    

    switch valuescale
        case 1
            figure('Name','Polar plot of magnitude')
        case 2
            figure('Name','Polar plot of energy')
        case 3
            figure('Name','Polar plot of level (40 dB scale)')
    end
    
    trisurf(K,x,y,z,values(:,b),'FaceColor','interp','FaceLighting','gouraud',...
        'EdgeLighting','gouraud','DiffuseStrength',1,'EdgeAlpha',0.1,...
        'AmbientStrength',0.6);
    camlight right
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
    
    if sphgrid > 0 && sphgrid <= 4
        spheregrid('radius',max(max(abs([x;y;z]))),...
            'degstep',30,...
            'gridtype',sphgrid-1,...
            'radiusdiv',4,...
            'LineWidth', 0.5,...
            'Color',[0,0,1],...
            'LineStyle',':'); % spheregrid is an aarae utility function
    end
    if cartgrid == 1
        grid on
    else
        grid off
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    set(gca, 'DataAspectRatio', [1 1 1])
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
    
    if isfield(IN,'bandID')
        title([num2str(IN.bandID(b)), ' Hz'])
    else
        title(['Band ', num2str(b)])
    end
    
end



OUT.funcallback.name = 'sphereplotfromMics.m';
OUT.funcallback.inarg = {percentilethreshold,valuescale,cartgrid,sphgrid};




%**************************************************************************
% Copyright (c) 2014, Densil Cabrera
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
%**************************************************************************
