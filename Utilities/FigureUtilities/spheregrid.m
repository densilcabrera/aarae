function spheregrid(varargin)
% This function plots a spherical grid to help with 3D plots.
%
% This function can be used to add to an existing figure (that is how it
% must be used when called from an AARAE function).
%
% INPUT PARAMETERS:
% The following input parameters are parsed from the varargin (see example
% that follows).
%
% radius is the outer radius of the grid lines
%
% degstep is the step in degrees between grid lines
% degstep is not relevant if gridtype == 1
%
% radiusdiv is the number of grid divisions along the radius
% radiusdiv is not relevant if gridtype == 0
%
% gridtype of 0 just plots a sphere (with lines of lattitude and longitude
% gridtype of 1 just plots circles on each of the Cartesian planes
% gridtype of 2 plots a combination of 0 and 1
%
% LineWidth, Color and LineStyle are Matlab line properties (as defined by
% Matlab)
%
%
% EXAMPLE
% spheregrid('radius',10,'degstep',45,'gridtype',0)
%
% code by Densil Cabrera
% 16 October 2014

p = inputParser;
addParameter(p,'radius',1,@isnumeric);
addParameter(p,'degstep',30,@isnumeric);
addParameter(p,'gridtype',2,@isnumeric);
addParameter(p,'radiusdiv',4,@isnumeric);
addParameter(p,'LineWidth',0.5,@isnumeric);
addParameter(p,'LineStyle',':');
addParameter(p,'Color',[0.5 0.5 0.5]);
parse(p,varargin{:});


radius = p.Results.radius;
degstep = p.Results.degstep;
gridtype = p.Results.gridtype;
radiusdiv = p.Results.radiusdiv;
LineWidth = p.Results.LineWidth;
Color = p.Results.Color;
LineStyle = p.Results.LineStyle;

if gridtype == 0 || gridtype == 2
    % ***** Plot a sphere at the specified radius   *****
    % plot lines of lattitude    
    numberoflines = floor(180/degstep);
    radstep = pi*degstep/180;
    az = (0:2*pi/360:2*pi)';
    for n = 1:numberoflines
        el = ((n-1)*radstep-pi/2);
        [x,y,z] = sph2cart(az,el*ones(size(az)),radius*ones(size(az)));
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
    
    
    % plot lines of longitude
    numberoflines = floor(360/degstep);
    el = -pi/2:pi/180:pi/2;
    for n = 1:numberoflines
        az = ((n-1)*radstep);
        [x,y,z] = sph2cart(az*ones(size(el)),el,radius*ones(size(el)));
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
end




if gridtype == 1 || gridtype == 2
    % ***** Plot planar circles   *****
    
    % horizontal
    numberoflines = radiusdiv-1;
    az = (0:2*pi/360:2*pi)';
    el = zeros(size(az));
    for n = 1:numberoflines
        radii = ones(size(az)) * n*radius/(numberoflines+1);
        [x,y,z] = sph2cart(az,el,radii);
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
    if gridtype == 1
        radii = ones(size(az)) * radius;
        [x,y,z] = sph2cart(az,el,radii);
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
    
    % median
    el = [-pi/2:pi/180:pi/2, (pi/2-pi/180):-pi/180:-pi/2];
    az = [zeros(1,181),pi*ones(1,180)];
    for n = 1:numberoflines
        radii = ones(size(az)) * n*radius/(numberoflines+1);
        [x,y,z] = sph2cart(az,el,radii);
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
    if gridtype == 1
        radii = ones(size(az)) * radius;
        [x,y,z] = sph2cart(az,el,radii);
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
    
    % transverse
    az = [pi/2*ones(1,181),3*pi/2*ones(1,180)];
    for n = 1:numberoflines
        radii = ones(size(az)) * n*radius/(numberoflines+1);
        [x,y,z] = sph2cart(az,el,radii);
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
    if gridtype == 1
        radii = ones(size(az)) * radius;
        [x,y,z] = sph2cart(az,el,radii);
        plot3(x,y,z,'Color',Color,'LineWidth',LineWidth,'LineStyle',LineStyle);
        hold on
    end
end





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
