function MultiWavePlot_aarae(in, fs, answer)
% This function plots a multichannel and/or multiband audio wave as lines
% stacked vertically on a single plot.
%
% It calls Christopher Hummersone's multiwaveplot function (2010).

if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1);
    audio = in.audio;
    fs = in.fs;
    if isfield(in,'bandID')
        bandID = in.bandID;
    end
    if isfield(in,'chanID')
        chanID = in.chanID;
    end
else
    audio = in;
end

if nargin < 3, answer = 'Channels'; end

[len, chans, bands] = size(audio);

t = (0:(len-1))./fs;

if ~exist('bandID','var')
    bandID = 1:bands;
end

if ~exist('chanID','var')
    chanID = 1:chans;
end



if chans > 1 && bands > 1
    if nargin < 3
        answer = questdlg('Channels or Bands in each figure?', ...
            'Settings', ...
            'Channels', ...
            'Bands', ...
            'Channels');
    end
    switch answer
        case 'Channels'
            for b = 1:bands
                figure('Name', ['Multiwave Plot, Band ', num2str(bandID(b))])
                multiwaveplot(t,1:bands,squeeze(audio(:,:,b))',1);
            end

        otherwise
            for c = 1:chans
                figure('Name', ['Multiwave Plot, Channel ', num2str(chanID(c))])
                multiwaveplot(t,1:chans,squeeze(audio(:,c,:))',1);
            end

    end

    
elseif chans == 1 && bands > 1
    figure('Name', ['Multiwave Plot, Channel ', num2str(chanID(1))])
    audio = permute(audio,[3,1,2]);
    multiwaveplot(t,1:bands,audio,1);
else
    figure('Name', ['Multiwave Plot, Band ', num2str(bandID(1))])
    multiwaveplot(t,1:chans,audio',1);
end



