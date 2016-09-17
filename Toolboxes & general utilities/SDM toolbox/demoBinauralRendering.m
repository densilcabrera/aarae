%% Example usage of SDM toolbox for analysis, synthesis, and binaural reproduction. 
% The data are room impulse responses measured in a
% control room type space in home environment.

% SDM toolbox : demoBinauralRendering
% Sakari Tervo & Jukka Pätynen, Aalto University, 2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

%% Load the impulse response and the source signal
% 1s long impulse response measured at 192 kHz
% Loudspeakers are custom-made passive 2-way loudspeakers at approximately
% 1 meter distance from the microphone array.
% IRs contain left and right channel of the audio system.

% Download a spatial room impulse response
ir_filename = 'IR_home_control_room';
if ~exist([ir_filename '.mat'],'file')
    disp(['Downloading an example IR ' ir_filename ' from the database.'])
    url_ir = ['https://mediatech.aalto.fi/~tervos/' ir_filename '.mat'];
    websave([ir_filename '.mat'],url_ir);
end
% Download a stereofile from free music archive
audio_filename = 'paper_navy_swan_song';
if ~exist([audio_filename ,'.mp3'],'file')
    disp('Downloading an example music file from free music archive.')
    url_of_the_song = 'https://mediatech.aalto.fi/~tervos/demoJAES/samples/Song1_CR1.mp3';
    outfilename = websave([audio_filename '.mp3'],url_of_the_song);
end
% If websave not supported, you have to download IRs and source signals
% manually from the urls given below
% 'https://mediatech.aalto.fi/~tervos/IR_home_control_room.mat'
% 'https://mediatech.aalto.fi/~tervos/demoJAES/samples/Song1_CR1.mp3'

%% Download CIPIC HRTFS
% You can also download CIPIC HRTFs manually by inserting the url in your
% browser. Save it to your working folder and unzip. It should be located
% in : currentfolder\CIPIC_hrtf_database
cipic_dirname = 'CIPIC_hrtf_database';
if ~exist(cipic_dirname,'dir')
    answer = input('Download and unzip CIPIC HRTFS (170 MB)? Yes/No : ','s');
    if strcmpi(answer,'yes');
        disp('Downloading and extracting example CIPIC HRTF database, ~170 MB')
        disp('This may take a while')
        % this is the original url of CIPIC HRTF database
        % url_cipic = 'http://interface.cipic.ucdavis.edu/data/CIPIC_hrtf_database.zip';
        url_cipic = 'https://mediatech.aalto.fi/~tervos/CIPIC_hrtf_database.zip';
        tic
        hrtfnames = unzip(url_cipic,'.'); % <--- Unzip to this (.) folder
        tc = toc;
        disp(hrtfnames')
        disp(['Ended downloading and unzipping in ' num2str(tc) ' seconds'])
    else
        disp('Binaural rendering requires HRTFS, you can specify them in createSynthesisStruct.m')
    end
end

%% Read the data
% Read impulse response
load([ir_filename '.mat'])
% Read stereo signal
S = audioread([audio_filename '.mp3']);
% Choose 10 seconds and resample
Sr = resample(S(1:44.e3*10,:),480,441);

%% Create SDM struct for analysis with a set of parameters
% Parameters required for the calculation
% Load default array and define some parameters with custom values
fs = 192e3;
a = createSDMStruct('DefaultArray','GRASVI25','fs',fs);

%% Calculate the SDM coefficients
% Solve the DOA of each time window assuming wide band reflections, white
% noise in the sensors and far-field (plane wave propagation model inside the array)
DOA{1} = SDMPar(ir_left, a);

% Here we are using the top-most microphone as the estimate for the
% pressure in the center of the array
P{1} = ir_left(:,5);

% Same for right channel
DOA{2} = SDMPar(ir_right, a);
P{2} = ir_right(:,5);


%% Create a struct for visualization with a set of parameters
% Load default setup for very small room and change some of the variables
v = createVisualizationStruct('DefaultRoom','VerySmall',...
    'name','Home Control Room','fs',fs);
% For visualization purposes, set the text interpreter to latex
set(0,'DefaultTextInterpreter','latex')

%% Draw analysis parameters and impulse responses
parameterVisualization(P, v);
%% Draw  time frequency visualization
timeFrequencyVisualization(P, v)
%% Draw the spatio temporal visualization for each section plane
v.plane = 'lateral';
spatioTemporalVisualization(P, DOA, v)
v.plane = 'transverse';
spatioTemporalVisualization(P, DOA, v)
v.plane = 'median';
spatioTemporalVisualization(P, DOA, v)

%% Create a struct for synthesis with a set of parameters
% Compare 2.1, 5.1, 22.2 and Aalto 24 loudspeaker setup via binaural synthesis
lspSetupNames = {'2.1','5.1','22.2','AALTO_24'};
for lspSetups = 1:length(lspSetupNames)
    % load default 5.1 setup and define some parameters with custom values
    s = createSynthesisStruct('Binaural',true,...
        'DefaultArray',lspSetupNames{lspSetups},...
        'HRTFset',28,... % <---- CIPIC HRTF subject number
        'snfft',length(P{1}),...
        'fs',192e3,...
        'c',343);
    % You always need to define 'snfft'
    
    %% Synthesize the spatial impulse response with NLS as binaural
    
    Hbin = cell(1,2);
    H = cell(1,2);
    for channel = 1:2
        [~, Hbin{channel}] = synthesizeSDMCoeffs(P{channel},DOA{channel}, s);
    end
    
    %% Convolution with the source signal
    
    % Resample H to 48e3 [Hz] sampling frequency for auralization
    Hbin{1} = resample(Hbin{1},1,4);
    Hbin{2} = resample(Hbin{2},1,4);
    
    Y = zeros(size(Sr,1),2);
    for channel = 1:2 % measured left and right channel
        for lsp = 1:2 % left and right ear
            % Convolution with Matlab's overlap-add
            Y(:,lsp) = Y(:,lsp) +  fftfilt(Hbin{channel}(:,lsp),Sr(:,channel));
        end
    end
    % Y contains the auralization of the spatial IRs with S
    
    %% Saving the auralization to a file
    % Save the file to the default folder with a custom filename.
    % Save the result as wav, as wav can handle upto 256 channels.
    disp('Started Auralization');tic
    savename = [ir_filename '_' audio_filename '_Binaural_' lspSetupNames{lspSetups} '.wav'];
    if max(abs(Y(:))) > 1
        Y = Y/max(abs(Y(:)))*.9;
        disp('Sound normalized, since otherwise would have clipped')
    end
    disp(['Ended Auralization in ' num2str(toc) ' seconds.'])
    disp('Started writing the auralization wav file')
    disp([savename  ' on the disk.']);tic
    audiowrite(savename,Y/10,s.fs/4) % <---- save the result as wav
    info = audioinfo(savename);
    disp('Wrote ... ');
    disp(info)
    disp(['... in ' num2str(toc) ' seconds'])
end
%% Playback using Matlab or other applications

% <--- EOF demoBinauralRendering.m




