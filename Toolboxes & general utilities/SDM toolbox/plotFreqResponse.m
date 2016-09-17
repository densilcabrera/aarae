function plotFreqResponse(H, f, color, lw, dBdynamics, maxdB, ...
    name, plotStyle, flagGrid, dBSpacing)
% plotSpatio(H, color, lw, res, plane, mindB, maxdB, name)
% Plots polar response of given H vector
%
% H     : Absolute squared spectrum [nfft 1]
% f     : frequencies corresponding to H [nfft 1]
% color : color in matlab RGB colors, for example [0 0 1] or 'k'
% lw    : the width of line, for example 2, value
% dbDynamics : dynamic range in dB, value
% maxdB : maximum dB value in H, value
% name  : title of the plot
% flagGrid  : true if grid is visible, false is off
% plotStyle : 'Line' or 'Fill' -style in the visualization 
% dBSpacing : grid of dB values in the visuzalization
% DOAspacing : grid of frequency values

% SDM toolbox : plotFreqResponse
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

f(1) = 1;
% --- Plot the frequency response ----
switch lower(plotStyle)
    case 'line'
        if lw <= 0 || isempty(lw)
            warning('plotPolarResponse : linewidth 0 or not specified, using linewidth 1')
            plot(log10(f),10*log10(H)-maxdB,'-','Linewidth',1,'Color', color)
        else
            plot(log10(f),10*log10(H)-maxdB,'-','Linewidth',lw,'Color', color)
        end
    case 'fill'
        % use gray color for the edge if edge is shown
        if lw > 0
            fill(log10([f f(end:-1:1)]),[10*log10(H)-maxdB; -dBdynamics*ones(size(H))]',...
               color, 'EdgeColor',[0.7 0.7 0.7],'Linewidth', lw)
        elseif isempty(lw) || lw == 0
             fill(log10([f f(end:-1:1)]),[10*log10(H)-maxdB; -dBdynamics*ones(size(H))]', ...
                color, 'EdgeColor','none')
        end
        %set(gca,'Xscale','log');
end

% --- EOF Plot the frequency response ----

% This part adds orientation and bearing description and grid lines to the plot
if flagGrid % is true
    
    if numel(dBSpacing) > 1
        lim = -dBSpacing(end:-1:1);
    else
        lim = 0:dBSpacing:dBdynamics; % dB grid points (at every 10 dB)
        lim = lim(end:-1:1);
    end
    set(gca,'Ytick',sort(-lim))
    ylim([-dBdynamics 3])
    set(gca,'Xtick',log10([30:10:100 200:100:1000 2000:1000:10000 ...
        20000]))
    set(gca,'XtickLabel',{'30','','','','','','','100','','','','','','','','',...
        '1000','','','','','','','','','10000','20000'}')
    xlim(log10([30 20000]))
    grid on
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]')
end
box off

% Name of the plot
title(name);
hold on;