function OUT = apply_tukeywindow(in, fs, tukeyratio, halfwin, doplay)
% Applies a Tukey window (or tapered cosine window) to the audio data.
% This can be thought of as providing a fade-in and fade-out, with the 
% middle part of the audio unaffected.
% 
% tukeyratio is the ratio of fadein+fadeout to the total window duration
% (the remainder of the window has unity gain).
%
% Half a window can be applied (i.e., just fade-in or just fade-out) by
% setting halfwin to -1 (fade-in only) or 1 (fade-out only).

if nargin < 5, doplay = 0; end
if nargin < 4, halfwin = 0; end
if nargin < 3
    %dialog box for settings
    prompt = {'Tukey window ratio (0:1)', ...
        'Half window (-1 | 0 | 1)', ...
        'Play and/or display (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0.1', '0', '0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        tukeyratio = str2num(answer{1,1});
        halfwin = str2num(answer{2,1});
        doplay = str2num(answer{3,1});
    else
        OUT = [];
        return
    end
end
if isstruct(in)
    data = in.audio;
    fs = in.fs;
else
    data = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(data) && ~isempty(fs) && ~isempty(tukeyratio) && ~isempty(halfwin) && ~isempty(doplay)
    [len,chans,bands,dim4,dim5,dim6] = size(data); % size of the audio
    

    if abs(halfwin) == 1
        wlen = 2*len;
    else
        wlen = len;
    end

    wf = tukeywin(wlen,tukeyratio);
    wtag = ['Tukey ',num2str(tukeyratio),' '];

    switch halfwin
        case -1
            wf = wf(1:len);
        case 1
            wf = flipud(wf(1:len));
    end


    y = data .* repmat(wf,[1,chans,bands,dim4,dim5,dim6]);
    if isstruct(in)
        OUT.audio = y;
        OUT.funcallback.name = 'apply_tukeywindow.m';
        OUT.funcallback.inarg = {fs,tukeyratio,halfwin,doplay};
    else
        OUT = y;
    end

    if doplay
        figure('Name', 'Window function applied to wave')
        plot(((1:len)'-1)./fs,mean(mean(mean(mean(mean(data,3),2),4),5),6), 'Color', [0.7 0.7 0.7])
        hold on
        plot(((1:len)'-1)./fs,mean(mean(mean(mean(mean(y,3),2),4),5),6), 'Color', [0.5 0.5 0.5])

        plot(((1:len)'-1)./fs,wf, 'r')
        title([wtag 'window function'])
        xlabel('Time (s)')
        ylabel('Value')
        hold off

        % Loop for replaying, saving and finishing
        choice = 0;

        % loop until the user presses the 'Done' button
        while choice < 5
            choice = menu('What next?', ...
                'Play audio', 'Plot all channels and bands', ...
                'Plot spectrogram', 'Save wav file', 'Discard', 'Done');
            switch choice
                case 1
                    ysumbands = mean(mean(mean(sum(y,3),4),5),6);
                    sound(ysumbands./max(max(abs(ysumbands))),fs)
                case 2
                    plot_wave_simple(y, fs)

                case 3
                    spectrogram_simple(y,fs);
                case 4
                    [filename, pathname] = uiputfile({'*.wav'},'Save as');
                    ysumbands = mean(mean(mean(sum(y,3),4),5),6);
                    if max(max(abs(ysumbands)))>1
                        ysumbands = ysumbands./max(max(abs(ysumbands)));
                        disp('Wav data has been normalized to avoid clipping')
                    end
                    if ischar(filename)
                        audiowrite([pathname,filename], ysumbands, fs);
                    end
                case 5
                    y = [];
            end
        end
    end
else
    OUT = [];
end