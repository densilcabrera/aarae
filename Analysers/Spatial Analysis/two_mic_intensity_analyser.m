function [OUT, varargout] = two_mic_intensity_analyser(IN, applygain,nfft)
% ASM (2016) release - this is a work in progress.
% Simple intensity analyser. IN should be a 2 channel signal. Channel 1 is
% pressure and channel 2 is particle velocity. See the two mic intensity
% beamforming processor.
% The outputs are time domain and frequency domain. Pressure, velocity
% (time only), Intensity, direction and impedance.
% Currently a work in progress.

if nargin ==1 

    param = inputdlg({'apply gain cal 1 = yes';... % These are the input box titles in the
        'fft power e.g. 2^n. 0 = length of signal (needs better implimentation)'},...% inputdlg window.
        'Intensity Analyser Settings',... % This is the dialog window title.
        [1 30],... % You can define the number of rows per
        ...        % input box and the number of character
        ...        % spaces that each box can display at once
        ...        % per row.
        {'1';'0'}); % And the preset answers for your dialog.
    
    param = str2num(char(param));
    
    if length(param) < 2, param = []; end % You should check that the user
    % has input all the required
    % fields.
    if ~isempty(param) % If they have, you can then assign the dialog's
        % inputs to your function's input parameters.
        applygain = param(1);
        nfft = 2^param(2);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end


% *************************************************************************
if isstruct(IN) % You should check that the function is being called within
    % the AARAE environment, if so, you can extract the
    % information you need to run your processor. Audio in
    % AARAE is in a Matlab structure (hence the isstruct test).
    
    % AARAE audio data can have up to 6 dimension, but most analysers cannot
    % handle all of those dimensions. The following utility function call
    % reduces the number of dimensions in the audio field of the input
    % structure if they are greater than the second input argument (3 in
    % the example below). For further information see the comments in
    % choose_from_higher_dimensions in AARAE's utilities directory. Unless
    % there is a reason not to, it is usually a good idea to support
    % analysis of the first three dimensions.
    IN = choose_from_higher_dimensions(IN,3,1);
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    
    % The following field values might not be needed for your anlyser, but
    % in some cases they are (delete if not needed). Bear in mind that
    % these fields might not be present in the input structure, and so a
    % way of dealing with missing field values might be needed. Options
    % include: using default values, asking for user input via a dialog
    % box, analysing the sound to derive values (probably impractical), and
    % exiting from the function
    if isfield(IN,'cal') % Get the calibration offset if it exists
        cal = IN.cal;
        % Note that the aarae basic processor cal_reset_aarae could be
        % called here. However, in this template it is called later.
    elseif ~isfield(IN, 'cal') && applygain ==1
        h=warndlg('Calibration data missing - no calibration will be applied','AARAE info','modal');
        uiwait(h)
        cal = 0;
        % Here is an example of how to exit the function with a warning
        % message
        
    else
        h=warndlg('Calibration data missing - please calibrate prior to calling this function.','AARAE info','modal');
        uiwait(h)
        OUT = []; % you need to return an empty output
        return % get out of here!
    end
    
    % chanID is a cell array of strings describing each channel
    if isfield(IN,'chanID') % Get the channel ID if it exists
        chanID = IN.chanID;
    else
        % or make a chanID using AARAE's utility function
        chanID = makechanID(size(audio,2),0);
    end
    
    % The name of the audio input could be useful in generating figures
    % (for example, in the title of a figure). This is a string.
    if isfield(IN,'name') % Get the AARAE name if it exists
        name = IN.name;
    else
        name = [];
    end
    
    
elseif ~isempty(param) || nargin > 1
    % If for example you want to enable your function to
    % run as a standalone MATLAB function, you can use
    % the IN input argument as an array of type double
    % and assign it to the audio your function is going
    % to process.
    audio = IN;
    % for the sake of this example, input_3 is being used for fs and
    % input_4 for cal. These values would be provided by an AARAE
    % structure, but need to be provided as additional inputs if the audio
    % is not input as a structure.
    fs = input_3;
    cal = input_4;
    chanID = makechanID(size(audio,2),0);
    bandID = 1:size(audio,3);
    name = '';
end
% *************************************************************************





% Check that the required data exists for analysis to run
if ~isempty(audio) && ~isempty(fs) && ~isempty(cal)
    
    % the audio's length, number of channels, and number of bands
    [len,chans,bands] = size(audio);
    if param(2)==0
        nfft = len;
    end
    if bands > 1
        error('Error:bands not allowed (at the moment)')
    end
    if applygain ==1 %apply gain cal. See processor 'two_mic_intensity'
        pressure = audio(:,1).*10^(cal(:,1)/20).*2e-5;
        velocity = audio(:,2).*10^(cal(:,1)/20).*2e-5;
    else
        pressure = audio(:,1);
        velocity = audio(:,2);
    end
    impedance = pressure./velocity;
    %impedance = intensity./velocity.^2;
    
    %calculate time domain vectors
    intensity = pressure.*velocity;
    %impedance = pressure.^2./intensity;
    %phasediff = angle(pressure + 1i*velocity);
    pressure_h = hilbert(pressure);
    velocity_h = hilbert(velocity);
    % phasediff = angle(pressure_h); p2 = angle(velocity_h); % Instantaneous phase (wrapped)
    phasediff = wrapToPi(angle(pressure_h) - angle(velocity_h));
    % direction
    % forward = 0 (or 2pi but not using wrapToPi)
    % reverse = pi
    % standing wave = pi/2
    
    pmax = max(abs(pressure));
    vmax = max(abs(velocity));
    imax = max(abs(intensity));
    
    % %     if snap == 1 %    0-pi/2-pi-(pi+pi/2)-pi/2
    % %         for i = 1:length(phasediff)
    % %             if abs(phasediff(i)) <=(pi/2)/2
    % %                 phasediff(i)=0;
    % %             elseif abs(phasediff(i))<=pi-(pi/2)/2
    % %                 phasediff(i)=pi/2;
    % %             elseif abs(phasediff(i)) <= pi +pi/4
    % %                 phasediff(i) = pi;
    % %             elseif abs(phasediff(i)) <= 2*pi-pi/4
    % %                 phasediff(i) = pi/2;
    % %             else
    % %                 phasediff(i) = 0;
    % %             end
    % %         end
    % %     end
    
    % % %     avint = cumsum(intensity);
    % % %     for i = 1:len
    % % %         averageintensity(i,1) = avint(i,1)/i;
    % % %     end
    % % %     cumsumpress = cumsum(pressure.^2);
    % % %     for i = 1:len
    % % %         runningpressure(i,1) = sqrt(cumsumpress(i,1)/i);
    % % %     end
    % % %     P = fft(pressure);
    % % %     U = fft(velocity);
    % % %     for i = 1:len;
    % % %         if i <=1000
    % % %             P(i) = P(i);
    % % %             U(i) = U(i);
    % % %         else
    % % %             P(i) = 0;
    % % %             U(i) = 0;
    % % %         end
    % % %     end
    % I = real(ifft(P.*U)); %fft multiplication - no USB.
    
    % average Intensity level and pressure level.
    % consider adding max levels?
    Ieq = 10*log10(abs(mean(intensity)/1e-12));
    Leq = 10*log10(mean(pressure.^2)/(2e-5)^2);
    %     disp(Leq)
    %     disp((Ieq))
    
    %%
    
    
    % *** WRITING TO THE LOG FILE ***
    % There is usually no reason to write to the log file, because a lot of
    % information is automatically written to it (depending on the outputs,
    % including tables). However, if you wish to write to the log file,
    % this is easily done using AARAE's logtext utility:
    % logtext('%% This is a test of writing to the log file.\n');
    
    % *** MAKING PLOTS ***
    % You may also include figures to display your results as plots.
    
    
    %%% plot time domain
    duration = len/fs;
    t = linspace(0,duration,len);
    figure1 = figure('Name', ['two mic intensity: ' IN.name]);
    set(figure1, 'Position', [100, 100, 1049, 1100]);
    filename = (['filename: ' IN.name]);
    results = (['Average Intensity Level:  ' num2str(Ieq) ' dB  |  Average Sound Pressure Level (Leq): ' num2str(Leq) ' dB']);
    annotation(figure1,'textbox',...
        [0.25 0.95 0.5 0.05],...
        'String',{filename, 'Instructions: Set time range with "X axis limits" button (bottom left)', results},...
        'FitBoxToText','off', 'interpreter', 'none', 'FontSize', 11);
    %pressure plot
    subplot(5,1,1)
    title('Pressure');
    %% Uncomment the following line to preserve the X-limits of the axes
    % xlim(axes1,[str2num(xmin) str2num(xmax)]);
    ylim([-pmax pmax]);
    box('on');
    grid('on');
    hold('on');
    % Create plot
    plot(t,pressure,'DisplayName','Pressure','LineWidth',1,...
        'Color',[0 0.447 0.741]);
    % Create ylabel
    ylabel('Pa');
    % Create legend
    legend show
    
    subplot(5,1,2)
    title('Particle Velocity')
    %% Uncomment the following line to preserve the X-limits of the axes
    % xlim(axes2,[str2num(xmin) str2num(xmax)]);
    ylim([-vmax vmax]);
    box on
    grid on
    hold on
    plot(t,velocity,'DisplayName','Velocity','LineWidth',1,...
        'Color',[0.257142857142857 0.9 0]);
    % Create ylabel
    ylabel('m/s');
    xlabel('Time (s)');
    % Create legend
    legend show
    
    subplot(5,1,3)
    title('Intensity')
    %% Uncomment the following line to preserve the X-limits of the axes
    % xlim(axes3,[20 20.01]);
    ylim([-imax imax]);
    box on
    grid on
    hold on
    plot(t,intensity, 'LineWidth',1, 'DisplayName', 'Intensity','Color',[0.9 0 0] );
    % Create ylabel
    ylabel('W/m^2');
    xlabel('Time (s)');
    % Create legend
    legend show
    
    %%direction
    % Create axes
    subplot(5,1,4)
    title('Direction')
    ylim([-3.5 3.5]);
    box on
    grid on
    hold on
    plot(t,phasediff, 'LineWidth',1, 'DisplayName',...
        'Direction','Color',[0.9 0.5 0]);
    ax = gca;
    ax.YTick = [-3.14 -1.57 0 1.57 3.14];
    ax.YTickLabel = {'reverse (\pi)', 'standing wave (\pi/2)',...
        'forward','standing wave (\pi/2)','reverse (\pi)'};
    xlabel('Time (s)');
    subplot(5,1,5)
    ylim([-1000 1000]);
    box on
    grid on
    hold on
    plot(t,impedance, 'LineWidth',1, 'DisplayName', 'Impedance','Color',[0.9 0.2 0.2] );
    title('Impedance')
    % Create ylabel
    ylabel('N*s/m^3');
    % Create legend
    legend show
    % Create xlabel
    xlabel('Time (s)');
    uicontrol('Style', 'pushbutton', 'String', 'X axis limits',...
        'Position', [0 0 65 30],...
        'Callback', 'setXaxislimits');
    legend show;
    hold off
    
    
    %% find frequency range for IRs and automatically truncate x axis?
    
    %% plot intensity spectrum
    fHz = linspace(0,fs,nfft);
    nyq = ceil(length(fHz)/2);
    Pfft = fft(pressure,nfft);
    Ufft = fft(velocity,nfft);
    Ifft = Pfft.*Ufft; %intensity fft?
    directionfft = angle(Pfft)-angle(Ufft);
    IMPfft = Pfft./Ufft;
    
    %     Ires = real(abs(Ifft);
    %     Iangle = angle(Ifft);
    figure2 = figure('Name', ['two mic intensity spectrum: ' IN.name]);
    set(figure2, 'Position', [100, 100, 1049, 1100]);
    subplot(4,1,1)
    semilogx(fHz(1:nyq), 10*log10(abs(Pfft(1:nyq).^2./len./(2e-5)^2)), 'DisplayName','Pressure','LineWidth',1,...
        'Color',[0 0.447 0.741]);
    xlabel('Frequency [Hz]'); ylabel('dB [ref: 2e-5]');
    title('Pressure')
    xlim([0 fs/2]);
    legend show;
    subplot(4,1,2)
    semilogx(fHz(1:nyq), 10*log10(abs(Ifft(1:nyq)./len./1e-12)),...
        'DisplayName','Intensity','LineWidth',1,'Color',...
        [0.9 0 0])
    xlabel('Frequency [Hz]'); ylabel('dB [ref: 1e-12]');
    title('Intensity')
    xlim([0 fs/2]);
    legend show;
    subplot(4,1,3)
    semilogx(fHz(1:nyq), directionfft(1:nyq),'LineWidth',1, 'DisplayName',...
        'Direction','Color',[0.9 0.5 0])
    grid on
    title('Direction')
    ax = gca;
    ax.YTick = [-6.28 -4.71 -3.14 -1.57 0 1.57 3.14 4.71 6.28];
    ax.YTickLabel = {'forward (2\pi)', 'standing wave (\pi/2)',...
        'reverse (\pi)', 'standing wave (\pi/2)','forward',...
        'standing wave (\pi/2)','reverse (\pi)','standing wave (\pi/2)',...
        'forward (2\pi)'};
    xlim([0 fs/2])
     xlabel('Frequency [Hz]');
    legend show;
    subplot(4,1,4)
    semilogx(fHz(1:nyq), real(IMPfft(1:nyq)),'LineWidth',1, 'DisplayName', 'real(Z)','Color',[0.9 0.2 0.2])
    hold on
    semilogx(fHz(1:nyq), imag(IMPfft(1:nyq)),'LineWidth',1, 'DisplayName', 'imag(Z)','Color',[0.5 0.5 0.5])
    xlabel('Frequency [Hz]'); ylabel('N*s/m^3');
    ylim([-1000 1000])
    xlim([0 fs/2])
    title('Impedance')
    legend show;
    
    uicontrol('Style', 'pushbutton', 'String', 'X axis limits',...
        'Position', [0 0 65 30],...
        'Callback', 'setXaxislimits');
    
    
    
    
    % *** CREATING A RESULTS LEAF (FOR BIG NON-AUDIO DATA) ***
    % You may output data to be plotted in a variety of charts, including
    % lines, mesh, surf, imagesc, loglog, and others depending on the
    % number of dimensions of your data using the doresultleaf.m function:
    % E.g.:
    %
    %     doresultleaf(figure1,'Type [units]',{'Time','channels','Frequency'},...
    %         'Time',      t,                's',           true,...
    %         'channels',  chanID,           'categorical', [],...
    %         'Frequency', num2cell(bandfc), 'Hz',          false,...
    %         'name','two_mic_intensity');
    %     %
    % Input arguments:
    % #1: Your data variable. It can be multidimensional, make sure you
    %     specify what each dimension is.
    % #2: What is your data variable representing? is it level? is it
    %     reverb time? make sure you label it appropriately and assign
    %     units to it, this second argument is a single string.
    % #3: This is a cell array where each cell contains the name of each
    %     dimension.
    %
    % #4: Name of your first dimension. (String)
    % #5: Matrix that matches your first dimension, in this case time.
    % #6: Units of your first dimension. (String)
    % #7: Can this dimension be used as a category? (true, false, [])
    %
    % Replicate arguments 4 to 7 for as many dimensions as your data
    % variable has.
    %
    % The second last input argument is the string 'name', this helps the
    % function identify that the last input argument is the name that will
    % be displayed in AARAEs categorical tree under the results leaf.
    
    
    
    % And once you have your result, you should set it up in an output form
    % that AARAE can understand.
    if isstruct(IN)
        % The output of this function, running within AARAE, is a
        % structure. As a minimum, the structure requires the function
        % callback. However, in the unusual event that you wish to
        % replicate the input structre before adding to it, it is easiest
        % to do this first:
        
        % *** OUTPUTTING AUDIO ***
        % Most analysers do not output audio. If you wish to output audio,
        % first consider whether your analyser should be a processor
        % instead. If it should be an analyser, then audio can be output as
        % follows:
        % OUT = IN; % You can replicate the input structure for your output
        % OUT.audio = audio; % And modify the fields you processed
        % However, for an analyser, you might not wish to output audio (in
        % which case the two lines above might not be wanted.
        
        % *** FUNCTION CALLBACKS ***
        % The following outputs are required so that AARAE can repeat the
        % analysis without user interaction (e.g. for batch processing),
        % and for writing the log file (and hence for generating
        % workflows).
        OUT.funcallback.name = 'two_mic_intensity_analyser.m'; % Provide AARAE
        % with the name of your function
        OUT.funcallback.inarg = {applygain, nfft};
        % assign all of the input parameters that could be used to call the
        % function without dialog box to the output field param (as a cell
        % array) in order to allow batch analysing. Do not include the
        % first (audio) input here. If there are no input arguments (apart
        % from the audio input), then use:
        % OUT.funcallback.inarg = {};
        
        
        % If you generated tables using AARAE's disptables (as described
        % above, you can attach the tables to the output structure
        %       OUT.tables = tables;
        
        
        % *** OUTPUTTING NEW FIELDS ***
        % Or simply output the fields you consider necessary after
        % processing the input audio data, AARAE will figure out what has
        % changed and complete the structure. But remember, it HAS TO BE a
        % structure if you're returning more than one field:
        
        %         OUT.properties.duration = duration;
        %         OUT.properties.maximum = maximum;
        %         OUT.properties.minimum = minimum;
        % The advantages of providing outputs as subfields of properties
        % are that AARAE has a button that opens a window to display
        % properties, and that properties are also written to the log file
        % of the AARAE session. Outputing fields without making them
        % subfields of properties is possible, but this makes the output
        % data harder to access.
        % Note that the above outputs might be considered to be redundant
        % if OUT.tables is used, as described above. Generally the
        % properties fields output is suitable for small data only (e.g.
        % single values). Tables is best for small and medium data, while a
        % results leaf is is best for big data.
        
        
        % AARAE will only display the output in the main GUI if it includes
        % tables, audio or other fields. It will not display if only
        % function callbacks are returned. Result leaves are not created by
        % the functions main output, but instead are created by calling
        % AARAE's doresultleaf as described above, and these will be
        % displayed in AARAE's main GUI regardless of the function's main
        % outputs.
    else
        % You may increase the functionality of your code by allowing the
        % output to be used as standalone and returning individual
        % arguments instead of a structure.
        OUT = audio;
    end
    % some additional outputs (not used by AARAE)
    %     varargout{1} = duration;
    %     varargout{2} = maximum;
    %     varargout{3} = minimum;
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2016, Jonothan Holmes
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
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
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