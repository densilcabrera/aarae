function out = histogram_of_wave_simple(in,transposesubplots,yscaling,nbars,normdist)
% Creates a histogram plot of a waveform, or a matrix of waveforms, as a
% function of value.
%
% By default:
%  If the waveform is multichannel, then channels are in rows.
%  If the waveform is multiband, then bands are in columns.
%  If 'Transpose subplots' is set in the dialog box, then these dimensions
%  are transposed (channels in columns, bands in rows).
%
% The histogram can be generated from the raw amplitude values, or it can
%  be transformed via squaring or by the absolute value of the Hilbert
%  transform.
%
% Optionally a normal distribution curve can be plotted.

if nargin < 2
    %dialog box for settings
    prompt = {'Transpose subplots (0 | 1)', ...
        'Amp (0), Amp^2 (1), |Hilb| (2)', ...
        'Number of histogram bars', ...
        'Plot normal distribution'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0','0','100','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        transposesubplots = str2num(answer{1,1});
        yscaling = str2num(answer{2,1});
        nbars = str2num(answer{3,1});
        normdist = str2num(answer{4,1});
    end
end
if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1);
    audio = in.audio;
    if isfield(in, 'cal')
       cal = in.cal; % get the calibration offset from the input structure
    end
else
    audio = in;
end

[len,chans,bands] = size(audio); 


% apply calibration factor if it has been provided
if exist('cal','var')
    audio = cal_reset_aarae(audio,0,cal);
end

if normdist
    dist = 'normal'; 
else
    dist = [];
end

if yscaling == 1
    audio = audio .^2;
    titletext = 'Squared amplitude';
elseif yscaling == 2;
    for b = 1:bands
        audio(:,:,b) = abs(hilbert(audio(:,:,b)));
    end
    titletext = 'Hilbert envelope';
% elseif yscaling == 3;
%     audio = 10*log10(audio.^2);
%     titletext = 'Level (dB)';
else
    titletext = 'Amplitude';
end


figure('Name', [titletext, ' histogram'])

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
            if normdist
                H = histfit(audio(:,ch,b),nbars,dist);
                set(H(2),'Color',[0.5 0.5 0.5])
            else
                hist(audio(:,ch,b),nbars)
            end
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',colorsb(b,:),'EdgeColor',colorsb(b,:))
            if ch == chans
                xlabel(titletext)
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
            if normdist
                H = histfit(audio(:,ch,b),nbars,dist);
                set(H(2),'Color',[0.5 0.5 0.5])
            else
                hist(audio(:,ch,b),nbars,'Color',colorsb(b,:))
            end
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',colorsb(b,:),'EdgeColor',colorsb(b,:))
            if b == 1
                xlabel(titletext)
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

% out.transposesubplots = transposesubplots;
% out.yscaling = yscaling;
% out.nbars = nbars;
out.funcallback.name = 'histogram_of_wave_simple.m';
out.funcallback.inarg = {transposesubplots,yscaling,nbars,normdist};


