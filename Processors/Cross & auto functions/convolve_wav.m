function OUT = convolve_wav(in,fs,wave2,fs2,more_options)
% This function convolves the selected audio data with a user-selected 
% wav file, using the frequency domain multiplication method.
%
% Code by Densil Cabrera
% version 1.01 (15 December 2013)
if nargin < 5, more_options = 0; end
% Use a dialog box to select a wav file.
if nargin < 4
    selection = choose_audio;
    if ~isempty(selection)
        wave2 = selection.audio;
        fs2 = selection.fs;
    else
        OUT = [];
        return
    end
end
if isstruct(in)
    wave1 = squeeze(sum(in.audio(:,:,:),3)); % sum the 3rd dimension if it exists
    fs = in.fs;
else
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
    wave1 = in;
end

if ~isempty(wave1) && ~isempty(fs) && ~isempty(wave2) && ~isempty(fs2)
    if ~(fs2 == fs)
        gcd_fs = gcd(fs,fs2); % greatest common denominator
        wave2 = resample(wave2,fs/gcd_fs,fs2/gcd_fs);
    end
    [len1, chans1] = size(wave1);
    [len2, chans2] = size(wave2);

    outputlength = len1 + len2 - 1;

    if chans1 == chans2
        y = ifft(fft(wave1, outputlength) .* fft(wave2, outputlength));
    elseif chans1 ==1
        y = ifft(fft(repmat(wave1,[1, chans2]), outputlength) ...
            .* fft(wave2, outputlength));
    elseif chans2 ==1
        y = ifft(fft(wave1, outputlength) ...
            .* fft(repmat(wave2,[1,chans1]), outputlength));
    else
        % in this case only the first channel of the wav file's audio is used
        y = ifft(fft(wave1, outputlength) ...
            .* fft(repmat(wave2(:,1),[1,chans1]), outputlength));
    end
    if more_options ~= 0
        % Loop for replaying, saving and finishing - currently bypassed
        choice = 0;
        % loop until the user presses the 'Done' button
        while choice < 4
            choice = menu('What next?', ...
                'Play', ...
                'Save wav file', 'Adjust gain', 'Discard', 'Done');
            switch choice

                case 1
                    sound(y ./ max(max(abs(y))),fs)

                case 2
                    [filename, pathname] = uiputfile({'*.wav'},'Save as');
                    if ischar(filename)
                        if max(max(abs(y))) <= 1
                            audiowrite([pathname,filename], y, fs);
                        else
                            audiowrite([pathname,filename], y./ max(max(abs(y))), fs);
                            disp('The saved audio has been normalized to prevent clipping.')
                        end
                    end
                case 3
                    Lmax = 20*log10(max(max(abs(y))));
                    prompt = {['Gain (dB) or ''n'' (max is ',num2str(Lmax), ' dBFS)']};
                    dlg_title = 'Gain';
                    num_lines = 1;
                    def = {num2str(-Lmax)};
                    answer = inputdlg(prompt,dlg_title,num_lines,def);
                    gain = answer{1,1};
                    if ischar(gain)
                        % normalize each channel individually
                        maxval = max(abs(y));
                        y = y ./ repmat(maxval,[outputlength,1]);
                    else
                        gain = str2num(gain);
                        y = y * 10.^(gain/20);
                    end

                case 4
                    y = [];
            end
        end
    end
    if isstruct(in)
        OUT = in; % replicate input structure, to preserve fields
    end
    if size(y,2) ~= size(wave1,2)
        OUT.chanID = makechanID(size(y,2),0);
    end
    OUT.audio = y;
    OUT.funcallback.name = 'convolve_wav.m';
    OUT.funcallback.inarg = {fs,wave2,fs2,more_options};
else
    OUT = [];
    warndlg('Not enough input arguments','AARAE info')
end