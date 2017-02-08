function OUT = dissipation(frequencies, temperatures, relhumidities, pressures, doplot)
% This function generates values for the atmospheric attenuation 
% coefficient (alpha) in decibels per metre. 
%
% When it is used within AARAE (or with no input arguments),
% it also generates an FIR filter (as an audio waveform)
% based on the atmospheric attenuation for a specified
% distance.
%
% Code by Densil Cabrera
% version 1.02 (9 December 2015)

% Inputs are frequencies (Hz), temperatures (deg C), relative humidities
% (%),  and atmospheric presures (kPa).
%
% Inputs can be single values, or vectors for any or all of the input
% arguments. The output is an up to 4D matrix with dimensions determined by
% the length of each vector, i.e., f x t x h x p.
%
% doplot enables or disables 2D and 3D plots (when the data are 1D and 2D)
% Plots are not generated for 0D, 3D or 4D data.

%
if nargin < 5, doplot = true; end
if nargin < 4, pressures = 90:5:110; end
if nargin < 3, relhumidities = 0:5:100; end
if nargin < 2, temperatures = 0:5:40; end
%if nargin < 1, frequencies = 10.^((13:43)./10); end

if nargin == 0
    % Dialog box mode
    % in this mode, only single values for pressure, relhumidity and
    % temperature are input, and the function generates a filter impulse
    % response that approximates the calculated dissipation.
    
    % Default values
    fs = 48000;
    distance = 20;
    temperatures = 20;
    relhumidities = 50;
    pressures = 100;
    filterlen = 100;
    minphase = 1;
    
    
    param = inputdlg({'Audio sampling rate (Hz)';...
        'Distance (m)'; ...
        'Temperature (deg C)'; ...
        'Relative humidity (%)'; ...
        'Atmospheric pressure (kPa)';...
        'Filter length (samples)';...
        'Minimum phase (0 | 1)'}, ...
        'Calculation parameters',1, ...
        {num2str(fs);num2str(distance); num2str(temperatures); ...
        num2str(relhumidities); num2str(pressures); num2str(filterlen);...
        num2str(minphase)});
    param = str2num(char(param));
    if length(param) < 7, param = []; end
    if ~isempty(param)
        fs = param(1);
        distance = param(2);
        temperatures = param(3);
        relhumidities = param(4);
        pressures = param(5);
        filterlen = param(6);
        minphase = param(7);
    else
        OUT = [];
        return
    end
    
    frequencies = linspace(0, fs/2,1000);
end

freq_l = length(frequencies);
temp_l = length(temperatures);
humid_l = length(relhumidities);
pres_l = length(pressures);

frequencies = repmat(reshape(frequencies, freq_l,1), [1 temp_l humid_l pres_l]);
temperatures = repmat(reshape(temperatures, 1, temp_l), [freq_l 1 humid_l pres_l]);
relhumidities = repmat(reshape(relhumidities, [1 1 humid_l]),[freq_l temp_l 1 pres_l]);
pressures = repmat(reshape(pressures, [1 1 1 pres_l]), [freq_l temp_l humid_l 1]);


% Dissiplation calculation
temperatures = temperatures + 273.15; % temparatures in Kelvin
T_on_T0 = temperatures ./ 293.15; % temperatures with respect to 20deg
psat_on_pr = 10.^(-6.8346.*(273.16./temperatures).^1.261+4.6151);
h = relhumidities.*psat_on_pr./(pressures./101.325); % molar concentration of water vapour (%)
pa_on_pr = pressures ./101.325;
FrO = pa_on_pr .*(24+40400.*h.*(0.02+h)./(0.391+h)); % relaxation freq of oxygen (Hz)
FrN = pa_on_pr.*T_on_T0.^-0.5.*(9+280.*h.*exp(-4.17.*(T_on_T0.^(-1/3)-1))); % relaxation freq of nitrogen (Hz)
alpha = 8.686.*frequencies.^2.*((0.0000000000184.*pa_on_pr.^-1.*T_on_T0.^0.5)...
    +T_on_T0.^(-5/2).*(0.01275.*exp(-2239.1./temperatures).*(FrO+frequencies.^2./FrO).^-1 ...
    +0.1068.*exp(-3352./temperatures).*(FrN+frequencies.^2./FrN).^-1));
%out.alpha = alpha;
% 
% speed of sound
h=h./100;
Mw = 29 - 11*h;
yw = (7+h)./(5+h);
c = 331.3 * (1 + (temperatures-273.15) ./ 273.15).^0.5...
    .* 4.5513 .* (yw./Mw).^0.5;
%disp(['Speed of sound: ' num2str(c(1)) ' m/s'])

% characteristic acoustic impedance
% rho = pressures ./ (287.05*temperatures) +
% z = rho .* c;

clear T_on_T0 psat_on_pr h pa_on_pr FrO FrN
frequencies = squeeze(frequencies(:,1,1,1));
temperatures = squeeze(temperatures(1,:,1,1)) - 273.15;
relhumidities = squeeze(relhumidities(1,1,:,1));
pressures = squeeze(pressures(1,1,1,:));



% Generate a filter for dissipation
if nargin == 0
    magnitude = 10.^((-alpha * distance) / 20);
    h = fir2(filterlen,frequencies./(fs/2),magnitude)';
    if minphase == 1
        [~, h] = rceps(h);
    end
    OUT.audio = h;
    OUT.fs = fs;
    OUT.properties.temperature = temperatures;
    OUT.properties.relhumidity = relhumidities;
    OUT.properties.pressure = pressures;
    OUT.properties.distance = distance;
    if minphase == 1
        OUT.properties.phasetype = 'minimum';
    else
        OUT.properties.phasetype = 'linear';
    end
    
else
    OUT = alpha;
end
if isstruct(OUT)
    OUT.funcallback.name = 'dissipation.m';
    OUT.funcallback.inarg = {frequencies, temperatures, relhumidities, pressures, doplot};
end

% PLOTTING
if doplot
    plotlogic = false(4,1);
    if freq_l > 1, plotlogic(1)=true; end
    if temp_l > 1, plotlogic(2)=true; end
    if humid_l > 1, plotlogic(3)=true; end
    if pres_l > 1, plotlogic(4)=true; end
    
    if plotlogic(1) && ~plotlogic(2) && ~plotlogic(3) && ~plotlogic(4)
        figure('Name','Dissipation')
        plot(frequencies, squeeze(alpha))
        xlabel('Frequency (Hz)')
        ylabel('alpha (dB/m)')
        title(['Atmospheric attenuation, temp ', num2str(temperatures), ' deg, rh ', ...
            num2str(relhumidities), ' %, pres ', num2str(pressures), ' kPa'])
        text(1000,max(squeeze(alpha))*0.97,['c = ' num2str(c(1)) ' m/s'])
    end
    if ~plotlogic(1) && plotlogic(2) && ~plotlogic(3) && ~plotlogic(4)
        figure('Name','Dissipation')
        plot(temperatures, squeeze(alpha))
        xlabel('Temperature (C)')
        ylabel('alpha (dB/m)')
        title(['Atmospheric attenuation, freq ', num2str(frequencies), ' Hz, rh ', ...
            num2str(relhumidities), ' %, pres ', num2str(pressures), ' kPa'])
    end
    if ~plotlogic(1) && ~plotlogic(2) && plotlogic(3) && ~plotlogic(4)
        figure('Name','Dissipation')
        plot(relhumidities, squeeze(alpha))
        xlabel('Relative Humidity (%)')
        ylabel('alpha (dB/m)')
        title(['Atmospheric attenuation, freq ', num2str(frequencies), ' Hz, temp ', ...
            num2str(temperatures), ' deg, pres ', num2str(pressures), ' kPa'])
    end
    if ~plotlogic(1) && ~plotlogic(2) && ~plotlogic(3) && plotlogic(4)
        figure('Name','Dissipation')
        plot(pressures, squeeze(alpha))
        xlabel('Pressure (kPa)')
        ylabel('alpha (dB/m)')
        title(['Atmospheric attenuation, freq ', num2str(frequencies), ' Hz, temp ', ...
            num2str(temperatures), ' deg, rh ', num2str(relhumidities), ' %'])
    end
    
    
    if plotlogic(1) && plotlogic(2) && ~plotlogic(3) && ~plotlogic(4)
        figure('Name','Dissipation')
        surf(temperatures, frequencies, squeeze(alpha))
        ylabel('Frequency (Hz)')
        xlabel('Temperature (C)')
        zlabel('alpha (dB/m)')
        title(['Atmospheric attenuation,  rh ', ...
            num2str(relhumidities), ' %, pres ', num2str(pressures), ' kPa'])
    end
    if plotlogic(1) && ~plotlogic(2) && plotlogic(3) && ~plotlogic(4)
        figure('Name','Dissipation')
        surf(relhumidities, frequencies, squeeze(alpha))
        ylabel('Frequency (Hz)')
        xlabel('Relative Humidity (%)')
        zlabel('alpha (dB/m)')
        title(['Atmospheric attenuation,  temp ', ...
            num2str(temperatures), ' deg, pres ', num2str(pressures), ' kPa'])
    end
    if plotlogic(1) && ~plotlogic(2) && ~plotlogic(3) && plotlogic(4)
        figure('Name','Dissipation')
        surf(pressures, frequencies, squeeze(alpha))
        ylabel('Frequency (Hz)')
        xlabel('Pressure (kPa)')
        zlabel('alpha (dB/m)')
        title(['Atmospheric attenuation,  temp ', ...
            num2str(temperatures), ' deg, rh ', num2str(relhumidities), ' %'])
    end
    if ~plotlogic(1) && plotlogic(2) && plotlogic(3) && ~plotlogic(4)
        figure('Name','Dissipation')
        surf(relhumidities, temperatures, squeeze(alpha))
        ylabel('Temperature (deg)')
        xlabel('Relative Humidity (%)')
        zlabel('alpha (dB/m)')
        title(['Atmospheric attenuation,  freq ', ...
            num2str(frequencies), ' Hz, pres ', num2str(pressures), ' kPa'])
    end
    if ~plotlogic(1) && plotlogic(2) && ~plotlogic(3) && plotlogic(4)
        figure('Name','Dissipation')
        surf(pressures, temperatures, squeeze(alpha))
        ylabel('Temperature (deg)')
        xlabel('Pressure (kPa)')
        zlabel('alpha (dB/m)')
        title(['Atmospheric attenuation,  freq ', ...
            num2str(frequencies), ' Hz, rh ', num2str(relhumidities), ' %'])
    end
    if ~plotlogic(1) && ~plotlogic(2) && plotlogic(3) && plotlogic(4)
        figure('Name','Dissipation')
        surf(pressures, relhumidities, squeeze(alpha))
        ylabel('Relative Humidity (%)')
        xlabel('Pressure (kPa)')
        zlabel('alpha (dB/m)')
        title(['Atmospheric attenuation,  freq ', ...
            num2str(frequencies), ' Hz, temp ', num2str(temperatures), ' deg'])
    end
end

