function OUT = spectrogram_2chancolor(in, fs, transposesubplots, winlen, NOVERLAP, dBrange, maxfreq, downsamp, wfchoice, hue)
% Generates a 2-channel spectrogram using complementary hues for the two
% channels
%
% Code by Densil Cabrera
% Version 1.0 (19 October 2013)
if isstruct(in)
    in = choose_from_higher_dimensions(in,2,1);
    audio = in.audio;
    fs = in.fs;
else
    audio = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end
Nyquist = fs/2;
if nargin < 3
    %dialog box for settings
    prompt = {'Transpose subplots (0 | 1)', ...
        'Window length (samples)', ...
        'Number of overlapping samples', ...
        'Decibel range', ...
        'Highest frequency to display (Hz)', ...
        'Downsample factor', ...
        'Rectangular, Hann, Blackman-Harris window (r|h|b)',...
        'Hues (0|1|2|3|4|5)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0','2048','1024','90',num2str(Nyquist), '1','h','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        transposesubplots = str2num(answer{1,1});
        winlen = str2num(answer{2,1});
        NOVERLAP = str2num(answer{3,1});
        dBrange = str2num(answer{4,1});
        maxfreq = str2num(answer{5,1});
        downsamp = str2num(answer{6,1});
        wfchoice = answer{7,1};
        hue = str2num(answer{8,1});
    else
        OUT = [];
        return
    end
end

% sum bands, if multiband
audio = mean(audio,3);

[~, chans] = size(audio);
if chans == 1
    % replicate single channel audio
    audio = repmat(audio,[1, 2]);
else
    % first two channels only
    audio = audio(:,1:2);
end

if wfchoice == 'h'
    wf = window(@hann,winlen);
    w = 'Hann window ';
elseif wfchoice == 'b'
    wf = window(@blackmanharris,winlen);
    w = 'Blackman-Harris window ';
else
    wf = ones(winlen,1);
    w = 'Rectangular window ';
end

if ~(downsamp == 1)
    if downsamp < 1, downsamp = 1/downsamp; end
    audio = resample(audio, 1,downsamp);
    fs = fs / downsamp;
    Nyquist = fs/2;
end

if maxfreq > Nyquist, maxfreq = Nyquist; end


[~,~,~,P] = spectrogram(audio(:,1),wf,NOVERLAP,[],fs);
y1 = 10*log10(P);
[~,F,T,P] = spectrogram(audio(:,2),wf,NOVERLAP,[],fs);
y2 = 10*log10(P);


titletext = '2-Channel Spectrogram';

figure('Name', titletext)



m= max([max(max(y1)) max(max(y2))]);




rangestring = [num2str(dBrange), ' dB range'];
z1 = ((y1 - m)+dBrange)./dBrange;
z1 = z1.*(sign(z1) .* (sign(z1)+1))./2;
z2 = ((y2 - m)+dBrange)./dBrange;
z2 = z2.*(sign(z2) .* (sign(z2)+1))./2;

z1 = repmat(z1, [1 1 3]);
z2 = repmat(z2, [1 1 3]);
switch hue
    case 0
        x1 = 0.1;
        x2 = 0.5;
        x3 = 0.9;
        chanstring = ' (ch1 blue, ch2 orange)';
    case 1
        x1 = 0.9;
        x2 = 0.5;
        x3 = 0.1;
        chanstring = ' (ch1 orange, ch2 blue)';
    case 2
        x1 = 0.5;
        x2 = 0.9;
        x3 = 0.1;
        chanstring = ' (ch1 green, ch2 violet)';
    case 3
        x1 = 0.5;
        x2 = 0.1;
        x3 = 0.9;
        chanstring = ' (ch1 violet, ch2 green)';
    case 4
        x1 = 0.1;
        x2 = 0.9;
        x3 = 0.5;
        chanstring = ' (ch1 green, ch2 magenta)';
    otherwise
        x1 = 0.9;
        x2 = 0.1;
        x3 = 0.5;
        chanstring = ' (ch1 magenta, ch2 green)';
end
z1(:,:,1) = z1(:,:,1) .* x1;
z1(:,:,2) = z1(:,:,2) .* x2;
z1(:,:,3) = z1(:,:,3) .* x3;
z2(:,:,1) = z2(:,:,1) .* (1-x1);
z2(:,:,2) = z2(:,:,2) .* (1-x2);
z2(:,:,3) = z2(:,:,3) .* (1-x3);

image(T,F,(z1+z2));
set(gca,'YDir','normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)')
set(gca,'TickDir', 'out');
title([w,num2str(winlen),'-pt window (fs=',num2str(fs),' Hz), ',...
    num2str(NOVERLAP/winlen*100),'% overlap, ', rangestring, chanstring]);

OUT.funcallback.name = 'spectrogram_2chancolor.m';
OUT.funcallback.inarg = {fs,transposesubplots,winlen,NOVERLAP,dBrange,maxfreq,downsamp,wfchoice,hue};