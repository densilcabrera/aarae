function timeFrequencyVisualization(P, p)
% spatioTemporaVisualization(P, DOA)
% Shows and calculates the spatio temporal diagram for given plane and 
% time indices.
%
% USAGE:
%  
% plane : 'lateral', 'median' or 'transverse'
% DOA   : 3-D locations (Cartesian) of each image-source, from SDMpar
% p     : struct from createTimeFreqStruct.m or can be manually given
% 
% EXAMPLES:
% 
% For further examples for setting the parameters in p, see 
% createTimeFreqStruct.m or demo*.m
%
% References
% [1] J. Pätynen, S. Tervo, T. Lokki, "Analysis of concert hall acoustics via
% visualizations of time-frequency and spatiotemporal responses",
% In J. Acoustical Society of America, vol. 133, no. 2, pp. 842-857, 2013.

% SDM toolbox : timeFrequencyVisualization
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

if nargin < 2
    error(['SDM toolbox : timeFrequencyVisualization : At least impulse response'
        ' and microphone array must be defined > SDMPar(IR, micArray)']);
end

disp(['Started time-frequency visualization.' ]);tic

ir_threshold = 0.2; % treshold level for the beginning of direct sound
nfft = p.fs*0.2;
f = (0:nfft-1)/nfft*p.fs;
HS = zeros(nfft, length(p.t)); 

numOfSources = length(P);

% Find the direct sound
for s = 1:numOfSources
    % meaningful time indices for a room
    t = round(p.t/1000*p.fs);
    
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
                H = abs(fft(P{s}(t2:t_end),nfft)).^2;
                
            case {'forward'} % USEFUL FOR LARGE ROOMS
                % ----- FORWARD INTEGRATION from direct sound to t ------
                t2 = min(t_ds+t(k),t_end);
                H = abs(fft(P{s}(t_ds:t2),nfft)).^2;
                
        end
        % Smoothing of the energy distribution with octave band smoothing
        HS(:,k) = HS(:,k) + reshape(octSmooth(H, p.fs, 3),[nfft 1]);
    end
end

% --- Plot the frequency responses for each time window ----
% find the maximum for plotting
maxdB = 10*log10(max(HS(:)));

flagGrid = false;
figure;
hold on
switch lower(p.DOI)
    case {'backward'}
        for k = 1:length(t)
            if k == length(t) && p.showGrid
                flagGrid = true;
            end
            plotFreqResponse(HS(:,k),f, p.colors(k,:), ...
                p.linewidth(max(1,k)), ...
                p.dBDynamics, maxdB, p.name, p.plotStyle, flagGrid, ...
                p.dBSpacing);
        end
    case {'forward'}
        HS = HS(:,end:-1:1);
        for k = 1:length(t)
            if k == length(t) && p.showGrid
                flagGrid = true;
            end
            plotFreqResponse(HS(:,k), f, p.colors(length(t)-k+1,:), ...
                p.linewidth(max(1,length(t)-k+1)), ...
                p.dBDynamics, maxdB, p.name, p.plotStyle, flagGrid, ...
                p.dBSpacing);
                    
        end
end
% Create names for the legend
legName = cell(1,length(t));

switch lower(p.DOI)
    case {'backward'}
        for k = 1:length(t)
            legName{k} = [num2str(p.t(k)) ' - ' num2str(round(t_end/p.fs*1000)) ' ms'];
        end
    case {'forward'}
        for k = 1:length(t)
            legName{length(t)-k+1} = ['0 - ' num2str(p.t(k)) ' ms'];
        end
end
legend(legName);

hold off

%set(gca,'Xscale','log');

disp(['Ended time-frequency visualization in ' num2str(toc) ' seconds.' ])

end

function [y] = octSmooth(x,fs,n)
% Smooth the frequency response over 1/n fractional octave bands
% 
% x : Matrix of frequency domain signal s, x = fft(s,N);
%     size [N numberOfsignals]
% y : smoothed spectrum, same size as x
% fs : sampling frequency in Hz
% n : 1/n octave-band for smoothing

% SDM toolbox : octSmooth
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi
% Copyleft

N = size(x,1);
f = (0:N-1)/N*fs;
f1 = round(f/(2^(1/(2*n)))/fs*N)+1;
f2 = round(f*(2^(1/(2*n)))/fs*N)+1;

f1(f1 > N-1) = N;
f2(f2 > N-1) = N;

y = zeros(N,size(x,2));

for i = 1:N
    y(i,:) = mean(abs(x(f1(i):f2(i),:)));
end

if nargout < 1
    figure(1)
    clf
    semilogx(f, 20*log10(abs(x)));
    hold all
    semilogx(f, 20*log10(abs(y)));
    hold off    
end

end



