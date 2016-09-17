%% Example usage of SDM toolbox for analysis, synthesis, and multi-channel reproduction.
% The data are room impulse responses measured in a living room.
% Sakari Tervo & Jukka Pätynen, Aalto University, 2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi
% Copyleft

% SDM toolbox : demoDefaultSetups
% Sakari Tervo & Jukka Pätynen, Aalto University, 2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

%% Load the impulse response and the source signal
% 1s long impulse response measured at 192 kHz
% IRs contains left and right channels
% loudspeakers are large size active loudspeakers ( Genelec <something> )
% Distance to loudspeakers about > 2 meters

% Download a spatial room impulse response
ir_filename = 'IR_living_room';
if ~exist([ir_filename '.mat'],'file')
     disp(['Downloading an example IR ' ir_filename ' from the database.'])
    url_ir = ['https://mediatech.aalto.fi/~tervos/' ir_filename '.mat'];
    websave([ir_filename '.mat'],url_ir);
end
%% Download a stereofile (originally from free music archive)
audio_filename = 'paper_navy_swan_song';
if ~exist([audio_filename ,'.mp3'],'file')
    disp('Downloading an example music file from free music archive.')
    url_of_the_song = 'https://mediatech.aalto.fi/~tervos/demoJAES/samples/Song1_CR1.mp3';
    outfilename = websave([audio_filename '.mp3'],url_of_the_song);
end

% If websave not supported, you have to download IRs and source signals
% manually from the urls given below
% 'https://mediatech.aalto.fi/~tervos/IR_living_room.mat'
% 'https://mediatech.aalto.fi/~tervos/demoJAES/samples/Song1_CR1.mp3'

%% Read the data
% Read impulse response
load([ir_filename '.mat'])
% Read stereo signal
S = audioread([audio_filename '.mp3']);


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
    'name','My Living Room','fs',fs);
% For visualization purposes, set the text interpreter to latex
set(0,'DefaultTextInterpreter','latex')

%% Draw analysis parameters and impulse responses
parameterVisualization(P, v);

%% Draw time frequency visualization
% Drawing only the lateral plane
timeFrequencyVisualization(P, v)

%% Draw the spatio temporal visualization
spatioTemporalVisualization(P, DOA, v)

%% Create synthesis struct with the given parameters
% Load default 5.1 setup and define some parameters with custom values
s = createSynthesisStruct('defaultArray','5.1',...
    'snfft',length(P{1}),...
    'fs',192e3,...
    'c',343);
% You always need to define 'snfft'

%% Synthesize the spatial impulse response with NLS
H = cell(1,2);
for channel = 1:2
    H{channel} = synthesizeSDMCoeffs(P{channel},DOA{channel}, s);
end

%% Convolution with the source signal

% Choose 10 seconds and resample
Sr = resample(S(1:44.e3*10,:),480,441);
numOfLsp = size(s.lspLocs,1);
Y = zeros(size(Sr,1),numOfLsp);

% Resample H to 48e3 [Hz] sampling frequency for auralization
H{1} = resample(H{1},1,4);
H{2} = resample(H{2},1,4);

for channel = 1:2;
    for lsp = 1:numOfLsp
        % Convolution with Matlab's overlap-add
        Y(:,lsp) = Y(:,lsp) +  fftfilt(H{channel}(:,lsp),Sr(:,channel));
    end
end
% Y contains the auralization of the spatial IRs with S

%% Saving the auralization to a file
% Save the file to the default folder with a custom filename.
% Save the result as wav, as wav can handle upto 256 channels.
disp('Started Auralization');tic
savename = [ir_filename '_' audio_filename '.wav'];
if max(abs(Y(:))) > 1
    Y = Y/max(abs(Y(:)))*.9;
    disp('Sound normalized, since otherwise would have clipped')
end
disp(['Ended Auralization in ' num2str(toc) ' seconds.'])
disp('Started writing the auralization wav file')
disp([savename  ' on the disk.']);tic
audiowrite(savename,Y/10,s.fs/4)
info = audioinfo(savename);
disp('Wrote ... ');
disp(info)
disp(['... in ' num2str(toc) ' seconds'])
%% Playback using Matlab or other applications

% <--- EOF demoBinauralRendering.m


