%%%Create spectrum movie from *.wav files
%%%It will take a while to make the video.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wav2mov is a simple srcipt to turn a audio *.wav file into a video of
%evolving spectrum
%
%
%
% INPUTS:
%       The script will ask for a *.wav file. For files of about 5 minutes
%       it takes about 15 minutes in my computer.
%
%
%
% OUTPUTS:
%       Once it is finished in the current folder you will find an *.avi file with the generated video already compresed to upload for ecample in youtube.
%
%
%
%This function was written by :
%                             Héctor Corte
%                             B.Sc. in physics 2010
%                             M.Sc. in physics of complex systems 2012
%                             Battery Research Laboratory
%                             University of Oviedo
%                             Department of Electrical, Computer, and Systems Engineering
%                             Campus de Viesques, Edificio Departamental 3 - Office 3.2.05
%                             PC: 33204, Gijón (Spain)
%                             Phone: (+34) 985 182 559
%                             Fax: (+34) 985 182 138
%                             Email: cortehector@uniovi.es
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First load the wav
[name, path]=uigetfile('*.wav');
[signal fs]= wavread([path name]);
[a,b]=size(signal);
duration = a / fs; %Lenght in seconds.
nSec=duration;
% Initialization of movie object.
%Video
frameRate=10;
mmObject=signalblks.MultimediaFileWriter([name,'.avi'],'AudioInputPort',1,'VideoInputPort',1,...
    'SampleRate',fs,'FrameRate',frameRate);
mmObject.AudioCompressor='GSM 6.10';
mmObject.VideoCompressor='MJPEG Compressor';
%Audio
nFrames=nSec*frameRate;
audioFrames=buffer(signal(:,1),round(length(signal(:,1))/nFrames));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now the iterations process to generate the movie.

%We precharge a number of elements to acelearte process.
numelelements=ceil(a/nFrames);
mygraph =figure('visible','off','position',[100,100,800,400]);
ini=1;
fin=(ini+ceil(a/nFrames)*50)*(ini+ceil(a/nFrames)*50<a)+a*(ini+ceil(a/nFrames)*50>=a);
%Each channel has it's own spectrum
[~,F,T,P] =  spectrogram(signal(ini:fin,1),200,80,200,10);
surf(T,F(1:40,1),10*log10(P(1:40,:)),'edgecolor','none'); axis tight;
[S,F,T,P] =  spectrogram(signal(ini:fin,2),200,80,200,10);
surf(T,F(1:40,1),10*log10(P(1:40,:)),'edgecolor','none'); axis tight;
%xlabel('Time (Seconds)'); ylabel('Hz');
view(145,65)
axis off
colormap hot
%To have better resolution, we print the current figure into a tiff
%image
print('-dtiff','-r90','temporal')
[X,map] = imread('temporal.tif');
F=im2frame(X,map);


for ind=1:nFrames
    
    %This is the same as before, calculate de spectrum of the two channels
    %and plot it.
    [~,F,T,P] =  spectrogram(signal(ini:fin,1),200,80,200,10);
    surf(T,F(1:40,1),10*log10(P(1:40,:)),'edgecolor','none'); axis tight;
    hold on
    [S,F,T,P] =  spectrogram(signal(ini:fin,2),200,80,200,10);
    surf(T,-F(1:40,1),10*log10(P(1:40,:)),'edgecolor','none'); axis tight;
    %Make some graphical changes
    hold off
    axis off
    view(125,55)
    %Save a temporal imagen and
    print('-dtiff','-r70','temporal')
    [X,map] = imread('temporal.tif');
    F=im2frame(X,map);
    
    %Add frames and audio to current video object.
    step(mmObject,F.cdata,audioFrames(:,ind));
    
    %Prepare for next step.
    ini=(ini+ceil(a/nFrames))*(ini+ceil(a/nFrames)<a)+a*(ini+ceil(a/nFrames)>=a);
    fin=(ini+ceil(a/nFrames)*30)*(ini+ceil(a/nFrames)*30<a)+a*(ini+ceil(a/nFrames)*30>=a);
    
end
set(gcf,'visible','on');
close (gcf)
close (mmObject)