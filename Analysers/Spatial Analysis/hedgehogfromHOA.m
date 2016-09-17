function OUT = hedgehogfromHOA(IN,fs,octfiltering, start_time,end_time,max_order,hif,lof,db_range,direct_limit,early_limit,mag_or_dB,bandfc)
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
% clean up - replace sphere cover with
if nargin < 10, db_range = 25; end
if nargin < 8, lof = 0; end
if nargin < 7, hif = fs/2; end
if nargin < 6, max_order=round(size(hoaSignals,2).^0.5-1); end
if nargin < 5, end_time = min([0.1, length(hoaSignals)/fs]); end
if nargin < 4, start_time = 0; end
if nargin < 3, octfiltering = 1; end
if isstruct(IN)
    if nargin < 3
        param = inputdlg({'Broadband [0], or Octave band [1]',... %'Spatial density of points (0.25:1)'
            'Start time [s]'...
            'End time [s]'...
            'Direct-early sound time limit [ms]'...
            'Early-late sound time limit [ms]'...
            'Maximum order'...
            'High frequency cutoff (Hz)'...
            'Low frequency cutoff (Hz)'...
            'Magnitude [1], or dB [2]'...
            'Peak detection range (dB):'},... %; 3D colormap surface in amplitude [2]; 3D lit surface in amplitude (slow) [3]'},...
            'Input parameters',1,...
            {num2str(octfiltering),...
            num2str(start_time),...
            num2str(end_time),...
            '10',...
            '50',...
            num2str(max_order),...
            num2str(fs/2),...
            '0',...
            '1',...
            '25'});
        if isempty(param) || isempty(param{1,1}) || isempty(param{2,1}) || isempty(param{3,1}) || isempty(param{4,1}) || isempty(param{5,1}) || isempty(param{6,1}) || isempty(param{7,1}) || isempty(param{8,1}) || isempty(param{9,1})|| isempty(param{10,1})
            OUT = [];
            return;
        else
            %sphere_cover = round(130*str2double(param{1,1}));
            % possibly move sphere cover to its own section
            
            octfiltering = str2double(param{1,1});
            start_time = str2double(param{2,1});
            
            end_time = str2double(param{3,1});
            direct_limit = str2double(param{4,1});
            early_limit = str2double(param{5,1});
            max_order = str2double(param{6,1});
            
            hif = str2double(param{7,1});
            lof = str2double(param{8,1});
            mag_or_dB = str2double(param{9,1});
            db_range = str2double(param{10,1});
            if isnan(octfiltering) || isnan(start_time) || isnan(end_time) || isnan(max_order) || isnan(hif) || isnan(lof), OUT = []; return; end
        end
    end
end

start_sample = round(start_time*fs);
if start_sample == 0, start_sample = 1; end


hoaSignals = autocropstart_aarae(hoaSignals,-20,3);

% zero padding
[~,chans,bands] = size(hoaSignals);
hoaSignals = [zeros(2,chans,bands);hoaSignals;zeros(2,chans,bands)];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OCTAVE BAND PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if octfiltering == 0 && ...
        ((hif < fs/2 && hif > lof) || (lof > 0 && lof < hif))
    % bandlimit the spectrum
    filterorder = 48;
    hoaSignals = bandpass(hoaSignals, lof, hif, filterorder, fs);
end


if ~exist('bandfc','var')
    bandfc = [];
end
if size(hoaSignals,3) == 1 && octfiltering == 1
    if isempty(bandfc)
        [hoaSignals bandfc] = octbandfilter_viaFFT(hoaSignals,fs);
    else
        hoaSignals = octbandfilter_viaFFT(hoaSignals,fs,bandfc);
    end
end
bands = size(hoaSignals,3);

plottype = 2;
sphere_cover = round(130*1);
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



number_of_directions = 362; %horzcat error

[azimuths,elevations] = SphereCovering(number_of_directions);

spkFmt = GenerateSpkFmt('sphCoord',[azimuths elevations ones(number_of_directions,1)]);

decodType    = 'basic' ;
sampFreq     = fs ;
filterLength = 128 ;
transFreq    = 'kr' ;
transWidth   = 200 ;
spkDistCor   = false ;
hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order) ;

hoa2SpkOpt = GenerateHoa2SpkOpt('decodType',decodType,'sampFreq',sampFreq,'filterLength',filterLength,'transFreq',transFreq,...
    'transWidth',transWidth,'spkDistCor',spkDistCor);
hoa2SpkCfg = Hoa2SpkDecodingFilters(hoaFmt,spkFmt,hoa2SpkOpt);

no_direct_sound_sigs = hoaSignals;
omni_signal = no_direct_sound_sigs(:,1,:);

beamsignals = zeros(size(no_direct_sound_sigs,1),number_of_directions,bands);

for i = 1:size(hoaSignals,2);
    for j = 1:number_of_directions;
        for b = 1:bands;
            beamsignals(:,j,b) = beamsignals(:,j,b)+(no_direct_sound_sigs(:,i,b).*hoa2SpkCfg.filters.gainMatrix(j,i));
        end
    end
end
for b = 1:bands;
    beamsignals(:,:,b) = beamsignals(:,:,b)./max(max(abs(beamsignals(:,:,b))));
end

%% Hedgehog plot
preroll = round((20./1000).*fs);

direct_sound_trim = preroll + ((20./1000).*fs);

%no_direct_sound_sigs = hoaSignals(direct_sound_trim:end,:,:);


IRreverse_integration = flip(cumsum(flip(omni_signal.^2,1)),1);

no_decay_IR = omni_signal ./ (IRreverse_integration.^0.5);

% no_decay_IR_squared = no_decay_IR.^2;

if end_sample > size(no_decay_IR,1)
    end_sample = size(no_decay_IR,1);
end



no_decay_IR_squared = no_decay_IR(1:end_sample,1,:).^2;
for b = 1:bands
    
    no_decay_IR_squared(:,1,b) = no_decay_IR_squared(:,1,b)./ max(abs(no_decay_IR_squared(:,1,b)));
    
    
    [pks,locs] = findpeaks(no_decay_IR_squared(:,1,b),'MinPeakHeight',db2pow(-abs(db_range))); %User input?
    %pkslocs is a cell array - pks in column 1, locs in column 2.
    pkslocs{b,1} = pks;
    pkslocs{b,2} = locs;
    % figure
    % plot(no_decay_IR_squared), hold on
    % plot(locs,pks,'k^','markerfacecolor',[1 0 0]), hold off
    
    
    for I = 1:length(pks);
        mag_peak = beamsignals(locs(I)-2:locs(I)+2,:,b);
        mag_square_sum = sqrt(sum(mag_peak.^2,1));
        [max_mag, ID] = max(mag_square_sum);
        directions_refl(I,:) = [azimuths(ID) elevations(ID)];
        magnitude_refl(I,:) = max_mag;
    end
    
    magnitude_refl_norm = magnitude_refl./max(magnitude_refl);
    magnitude_refl_db = mag2db(magnitude_refl);
    magnitude_refl_db = magnitude_refl_db + abs(min(magnitude_refl_db)); %Scaling?
    
    locs_time = locs/fs;
    
    switch mag_or_dB
        case 1
            [xplotrefs, yplotrefs, zplotrefs] = sph2cart(directions_refl(:,1), directions_refl(:,2), magnitude_refl_norm);%or that
            [xsphere, ysphere, zsphere] = sphere(20);
        case 2
            [xplotrefs, yplotrefs, zplotrefs] = sph2cart(directions_refl(:,1), directions_refl(:,2), magnitude_refl_db);%either this
            [xsphere, ysphere, zsphere] = sphere(20);
            xsphere = xsphere .* max(magnitude_refl_db);% dB
            ysphere = ysphere .* max(magnitude_refl_db);% db
            zsphere = zsphere .* max(magnitude_refl_db);% db
    end
    
    
    
    % Early and late time periods
    for I = 1:length(locs_time);
        if locs_time(I) <= direct_limit/1000; % I Think this is the time parameters?
            color_refls(I,:) = [1 0 0];
        elseif locs_time(I) >direct_limit/1000 && locs_time(I) <= early_limit/1000;
            color_refls(I,:) = [0 0 1];
        else
            color_refls(I,:) = [0 1 0];
            %alpha(0.2);
        end
    end
    
    
    figure1 = figure;
    
    %create axes
        axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],...
            'DataAspectRatio',[1 1 1],...
            'CameraViewAngle',10.339584907202,...
            'Position',[0.05,0.3,0.96,0.65]);
        view(axes1,[-45 35]);
        grid(axes1,'on');
        hold(axes1,'all');
    if ~isempty(bandfc)
        title([num2str(bandfc(b)), ' Hz'])
    elseif isfield(IN,'bandID')
        title([num2str(IN.bandID(b)), ' Hz'])
    else
        title(['Band ', num2str(b)])
    end
    
        % Create surf
        surf(xsphere,ysphere,zsphere,'Parent',axes1,'EdgeLighting','flat',...
            'FaceLighting','none',...
            'LineStyle',':',...
            'FaceColor','none');
    
    hold on
    for i = 1:length(xplotrefs);
        plot3([0,xplotrefs(i)],[0,yplotrefs(i)],[0,zplotrefs(i)],'LineWidth',1.5,'Color',color_refls(i,:))
        hold on
    end
    
    hold off
    
    
    axes2 = axes('Parent',figure1,...
    'Position',[0.0687393951356129 0.0789431157078216 0.903339198452805 0.166478129713424]);
    
    directendindex = round(1+fs*direct_limit/1000);
    earlyendindex = round(1+fs*early_limit/1000);
    normgain = 1/max(omni_signal(1:directendindex,1,b));
%     plot((0:length(no_decay_IR_squared(:,1,1))-1)./fs,pow2db(no_decay_IR_squared(:,1,b)),'Parent',axes2)
    plot((pkslocs{b,2}-1)./fs,pow2db(pkslocs{b,1}),'k.','markerfacecolor',[1 1 0],'Parent',axes2) 
    hold on   
    plot((0:directendindex-1)./fs,mag2db(normgain*omni_signal(1:directendindex,1,b)),'Parent',axes2,'color',[1,0.5,0.5])
    hold on
    plot((directendindex:earlyendindex-1)./fs,mag2db(normgain*omni_signal(directendindex+1:earlyendindex,1,b)),'Parent',axes2,'color',[0.5,0.5,1])
    plot((earlyendindex:size(no_decay_IR_squared,1)-1)./fs,mag2db(normgain*omni_signal(earlyendindex+1:size(no_decay_IR_squared,1),1,b)),'Parent',axes2,'color',[0.3,0.8,0.3])
    xlabel('Time (s)')
    ylim([-abs(db_range) 0])
    hold off
end


% end


if isstruct(IN)
    OUT.funcallback.name = 'hedgehogfromHOA.m';
    OUT.funcallback.inarg = {fs,octfiltering,start_time,end_time,max_order,hif,lof,db_range,direct_limit,early_limit,mag_or_dB,bandfc};
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