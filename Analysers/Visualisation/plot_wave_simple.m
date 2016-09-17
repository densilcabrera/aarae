function out = plot_wave_simple(in, fs, yscaling, transposesubplots)
% Creates a simple plot of a waveform, or a matrix of waveforms, as a
% function of time.
%
% By default:
%  If the waveform is multichannel, then channels are in rows.
%  If the waveform is multiband, then bands are in columns.
%  If 'Transpose subplots' is set in the dialog box, then these dimensions
%  are transposed (channels in columns, bands in rows).
%
% The waveform can be plotted 'as is' or it can be transformed via:
%  * squaring
%  * absolute value of the Hilbert transform
%  * squared amplitude expressed in decibels
if nargin < 4, transposesubplots = 0; end
if nargin < 3, yscaling = 0; end
if nargin < 3
    %dialog box for settings
    prompt = {'Transpose subplots (0 | 1)', 'Amp (0), Amp^2 (1), |Hilb| (2), Level (3)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        transposesubplots = str2num(answer{1,1});
        yscaling = str2num(answer{2,1});
    end
end
if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1);
    audio = in.audio;
    fs = in.fs;
    if isfield(in, 'cal')
       cal = in.cal; % get the calibration offset from the input structure
    end
else
    audio = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end
if ~isempty(audio) && ~isempty(fs) && ~isempty(yscaling) && ~isempty(transposesubplots)
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

    if yscaling == 1
        audio = audio .^2;
        titletext = 'Squared amplitude';
    elseif yscaling == 2;
        for b = 1:bands
            audio(:,:,b) = abs(hilbert(audio(:,:,b)));
        end
        titletext = 'Hilbert envelope';
    elseif yscaling == 3;
        audio = 10*log10(audio.^2);
        titletext = 'Level (dB)';
    else
        titletext = 'Amplitude';
    end

    t = ((1:len)'-1)/fs;

    figure('Name', titletext)

    % get some colors
    colors = selectcolormap('rainbow256');
    if bands == 1
        colorsb = colors(128,:);
    else
    colorsb = colors(round(1:255/(bands-1):256),:);
    end

    k = 1; % subplot counter

    if ~transposesubplots
        for ch = 1:chans
            for b = 1:bands
                subplot(chans,bands,k)
                plot(t,audio(:,ch,b),'Color',colorsb(b,:))
                if ch == chans
                    xlabel('Time (s)')
                end
                if b == 1
                    ylabel(['Chan ',num2str(ch)])
                end
                if ch == 1
                    title(['Band ',num2str(b)])
                end
                k = k+1;
            end
        end
    else
        for b = bands:-1:1
            for ch = 1:chans
                subplot(bands,chans,k)
                plot(t,audio(:,ch,b),'Color',colorsb(b,:))
                if b == 1
                    xlabel('Time (s)')
                end
                if b == bands
                    title(['Chan ',num2str(ch)])
                end
                if ch == 1
                    ylabel(['Band ',num2str(b)])
                end
                k = k+1;
            end
        end
    end

    out.transposesubplots = transposesubplots;
    out.yscaling = yscaling;
    out.funcallback.name = 'plot_wave_simple.m';
    out.funcallback.inarg = {fs,yscaling,transposesubplots};
else
    out = [];
end


