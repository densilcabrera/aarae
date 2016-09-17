function out = reversetimefilter(in)
% This function time-reverses the audio waveform prior to applying octave
% band or 1/3-octave band filtering

prompt = {'Octave or 1/3-octave band? (1 | 3)', ...
        'Zero-pad length at start (samples)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        out = [];
        return
    else
        bpo = str2double(answer{1,1});
        zeropadlen = abs(str2double(answer{2,1}));
        
    end
     
    audio = flipdim(in.audio,1);
    fs = in.fs;
    
    [~, chans, bands] = size(audio);
    if bands > 1
        % mix down bands if more than one present
        audio = mean(audio,3);
        disp('The multi-band signal has been mixed-down prior to filtering')
    end
    
    if zeropadlen > 0
        audio = [audio; zeros(zeropadlen,chans)];
    end
    
    if bpo == 3
        [audio, bandID] = thirdoctbandfilter(audio,fs);
    else
        [audio, bandID] = octbandfilter(audio,fs);
    end
    
    out.audio = flipdim(audio,1);
    out.bandID = bandID;