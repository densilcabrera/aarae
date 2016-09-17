function parameterVisualization(P, p)
% parameterVisualization(P, p)
% Visualizes the applied parameters in p for P.
%
% USAGE:
% 
% P     : Impulse response (Pressure) in the center of the microphone array
%          Size [N 1]
% p     : struct from createVisualizationStruct.m or can be manually given
% 
% EXAMPLES:
% 
% For further examples for setting the parameters in p, see 
% demo*.m
%
% References
% [1] J. Pätynen, S. Tervo, T. Lokki, "Analysis of concert hall acoustics via
% visualizations of time-frequency and spatiotemporal responses",
% In J. Acoustical Society of America, vol. 133, no. 2, pp. 842-857, 2013.

% SDM toolbox : parameterVisualization
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

if nargin < 2
    error(['SDM toolbox : spatioTemporalVisualization : At least impulse response'
        ' and microphone array must be defined > SDMPar(IR, micArray)']);
end

disp(['Started visualization of parameters.' ]);tic

ir_threshold = 0.2; % treshold level for the beginning of direct sound

HS = zeros(360/p.res + 1, length(p.t)); 

numOfSources = length(P);


T = cell(numOfSources,1);
for s = 1:numOfSources
    % meaningful time indices for a room
    t = round(p.t/1000*p.fs);
    
    % Find the direct sound
    ind = find(abs(P{s})/max(abs(P{s})) > ir_threshold,1,'first');
    pre_threshold = round(0.001*p.fs); % go back 1 ms
    t_ds = ind - pre_threshold; % direct sound begins from this sample
    
    % make sure that the time index of direct sound is greater than or equal to 1
    t_ds = max(t_ds, 1);
    
    % make sure the integration does not exceed the vector length
    t_end = length(P{s});

    % Iterate through different time windows
    
    for k = 1:length(t)
        switch lower(p.DOI)
              case {'backward'} % USEFUL FOR SMALL ROOMS
                % --- BACKWARD INTEGRATION from t to t(end) ----
                t2 = min(t_ds+t(k),t_end);
                T{s}(k,:) = [t2, t_end];
                
            case {'forward'} % USEFUL FOR LARGE ROOMS
                % ----- FORWARD INTEGRATION from direct sound to t ------
                t2 = min(t_ds+t(k),t_end);
                T{s}(k,:) = [t_ds t2];
        end
    end
end

% --- Plot the polar responses for each time window ----
% find the maximum for plotting
figure;

for s = 1:numOfSources
    subplot(numOfSources,1,s)
    hold on
    
    maxval = max(abs(P{s}));
    
    heightk = linspace(maxval, 2*maxval,length(t));
    
    switch lower(p.DOI)
        case {'backward'}
            for k = 1:length(t)
                plotParameters([T{s}(length(t)-k+1,1) T{s}(length(t)-k+1,2)]/p.fs, p.plotStyle, ...
                    p.colors(length(t)-k+1,:), ...
                    p.linewidth(max(1,k)),...
                    heightk(length(t)-k+1))
                
            end
        case {'forward'}

            for k = 1:length(t)
                plotParameters([T{s}(length(t)-k+1,1) T{s}(length(t)-k+1,2)]/p.fs, p.plotStyle, ...
                    p.colors(length(t)-k+1,:), ...
                    p.linewidth(max(1,length(t)-k+1)),...
                    heightk(length(t)-k+1))
            end
    end
    plot((1:length(P{s}))/p.fs, P{s},'-','Color',mean(p.colors));
    hold off
    ylabel('Pressure [.]')
    if s == 1
        title('Applied time windowing and impulse responses')
    end
end

% Create names for the legend
legName = cell(1,length(t));

switch lower(p.DOI)
    case {'backward'}
        for k = 1:length(t)
            legName{length(t)-k+1} = [num2str(p.t(k)) ' - ' num2str(round(t_end/p.fs*1000)) ' ms'];
        end
    case {'forward'}
        for k = 1:length(t)
            legName{length(t)-k+1} = ['0 - ' num2str(p.t(k)) ' ms'];
        end
end
legName{length(t)+1} = 'IR';
legend(legName);
hold off
xlabel('Time [s]')


disp(['Ended visualization of parameters in ' num2str(toc) ' seconds.' ])

end

function plotParameters(T, plotStyle,color,lw, maxval)
% plotParameters(T, plotStyle,color,lw, maxvalue)
% Plots time parameters
%
% H     : Directional energy histogram, i.e., energy response in linear
%        scale (not in dB)
% plane : 'lateral', 'median' or 'transverse'
% color : color in matlab RGB colors, for example [0 0 1] or 'k'
% lw    : the width of line, for example 2
% res   : resolution
% dbDynamis : dynamic range in dB
% name  : title of the plot
% X     : origin in the visualization, e.g., [0 0]
% flagGrid  : true if grid is visible, false is off
% plotStyle : 'Line' or 'Fill' -style in the visualization 
% dBSpacing : grid of dB values in the visuzalization
% DOAspacing : grid of DOA values in the visuzalization

% SDM toolbox :  plotPolarResponse
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Copyleft

% Plot the time windows
switch lower(plotStyle)
    case 'line'
         if lw <= 0 || isempty(lw)
            warning('plotPolarResponse : linewidth 0 or not specified, using linewidth 1')
            plot([T(1) T(1) T(2) T(2)],[0 maxval maxval 0],...
                '-','Linewidth',1,'Color', color);
        else
             plot([T(1) T(1) T(2) T(2)],[0 maxval maxval 0],...
                 '-','Linewidth',lw,'Color', color)
        end
    case 'fill'
        % use gray color for the edge if edge is shown
        if lw > 0
            fill([T(1) T(2) T(2) T(1)],[0 0 maxval maxval], color, 'EdgeColor',[0.7 0.7 0.7],'Linewidth', lw)
        elseif isempty(lw) || lw == 0
            fill([T(1) T(2) T(2) T(1)],[0 0 maxval maxval], color, 'EdgeColor','none')
        end
end

end

