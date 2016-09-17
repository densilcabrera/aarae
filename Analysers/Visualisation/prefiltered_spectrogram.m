function out = prefiltered_spectrogram(in, fs, zscaling)
% Uses pre-filtered data to generate a spectrogram.
%
% Make sure that the input data is multiband.
if nargin < 3
    %dialog box for settings
    prompt = {'Amp (0), Amp^2 (1), |Hilb| (2), Level (3)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        
        zscaling = str2num(answer{1,1});
    end
end
if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1);
    audio = in.audio;
    fs = in.fs;
    if isfield(in, 'cal')
       cal = in.cal; % get the calibration offset from the input structure
    end
    if isfield(in, 'bandID')
        bandfreq = in.bandID;
    end
else
    audio = in;
end

S = size(audio); % size of the IR
ndim = length(S); % number of dimensions
switch ndim
    case 1
        len = S(1); % number of samples in IR
        chans = 1; % number of channels
        bands = 1; % number of bands
    case 2
        len = S(1); % number of samples in IR
        chans = S(2); % number of channels
        bands = 1; % number of bands
    case 3
        len = S(1); % number of samples in IR
        chans = S(2); % number of channels
        bands = S(3); % number of bands
end

% apply calibration factor if it has been provided
if exist('cal','var')
    audio = audio .* 10.^(repmat(cal,[len,1,bands])/20);
end

if zscaling == 1
    audio = audio .^2;
    titletext = 'Squared amplitude';
elseif zscaling == 2;
    for b = 1:bands
        audio(:,:,b) = abs(hilbert(audio(:,:,b)));
    end
    titletext = 'Hilbert envelope';
elseif zscaling == 3;
    audio = 10*log10(audio.^2);
    titletext = 'Level (dB)';
else
    titletext = 'Amplitude';
end

t = ((1:len)'-1)/fs;

if zscaling == 0
    colors = selectcolormap('bipolar1');
else
    colors = 'jet';
end

if ~exist('bandfreq','var')
    bandfreq = 1:bands;
end

figure('Name','Prefiltered Spectrogram');

for ch = 1:chans
    subplot(chans,1,ch)
    imagesc(t,1:bands,squeeze(audio(:,ch,:))')
    colormap(colors)
    set(gca,'YDir','normal');
    set(gca,'YTickLabel',num2cell(bandfreq),'YTick',1:bands);
    if bands > 1, ylim([1 bands]); end
    ylabel('Band');
    xlabel('Time (s)');
    if zscaling == 3
        Lmax =  max(max(audio));
        Lmin = Lmax - 60;
        set(gca, 'Clim', [Lmin Lmax]);
    end
end

out.funcallback.name = 'prefiltered_spectrogram.m';
out.funcallback.inarg = {fs,zscaling};