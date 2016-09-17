function OUT = tilt_spectrum_logf(in, fs, dBperOct, lowf,hif, pivotf, doplay)
% This function applies a tilt to the magnitude spectrum of the input. The
% tilt is linear as a function of log(frequency), and is specified in
% decibels per octave.
%
% The operation is done in the frequency domain, without changing phase.
%
% This can, for example, be used to transform between white and pink noise,
% or assist in the derivation of impulse responses from a logarithmic swept
% sine that is missing its inverse filter, or for spectrally weighting MLS
% and unweighting it in the analysis.
%
% As well as specifying the spectral slope, you can also specify low and
% high frequency limits to the spectrum slope (beyond which the gain
% remains constant). If the slope is not limited in frequency, it can
% create problems with very large gain or attenuation (especially in the
% low frequency range).
%
% The pivot frequency is the frequency at which the filter has 0 dB gain
% (1000 Hz by default).
%
%
% Code by Densil Cabrera version 1.01 (14 August 2014)



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


if nargin < 4, doplay = 0; end
if nargin < 3
    %dialog box for settings
    prompt = {'Spectrum tilt (dB/octave)', ...
        'Low frequency limit (Hz)',...
        'High frequency limit (Hz)',...
        'Pivot frequency (Hz) at which there is 0 dB gain',...
        'Play and/or display (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'3','20',num2str(fs/2),'1000','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        dBperOct = str2num(answer{1,1});
        lowf = str2num(answer{2,1});
        hif = str2num(answer{3,1});
        pivotf = str2num(answer{4,1});
        doplay = str2num(answer{5,1});
    else
        OUT = [];
        return
    end
end

if ~isempty(data) && ~isempty(fs) && ~isempty(dBperOct) && ~isempty(doplay)
    [len,chans,bands,dim4,dim5,dim6] = size(data); % size of the audio
    
    if mod(len,2) == 0
        evenlen = 1;
    else
        evenlen = 0;
    end
    nsamples = 2*ceil(len/2); % use an even length fft

    % exponent used to generate magnitude slope
    fexponent = dBperOct/3;

    % magnitude slope of half spectrum not including DC and Nyquist components
    magslope = ((1:nsamples/2-1)./(nsamples/4)).^(fexponent*0.5)';
    
    % find and apply low and high frequency indices
    if lowf > 0
        flowindx = ceil(nsamples * lowf/fs);
        magslope(1:flowindx) = magslope(flowindx);
    end
    
    if hif > lowf && hif < fs/2
        fhiindx = floor(nsamples * hif/fs);
        magslope(fhiindx:end) = magslope(fhiindx);
    end
    
    % apply pivot frequency
    pivotindx = round(nsamples * pivotf/fs);
    if pivotindx < 1, pivotindx = 1; end
    if pivotindx > length(magslope), pivotindx = length(magslope); end
    magslope = magslope ./ magslope(pivotindx);
    

    % magnitude slope across the whole spectrum
    magslope = [magslope(1); magslope; magslope(end); flipud(magslope)];

    % apply the spectrum slope filter
    y = ifft(fft(data,nsamples) .* repmat(magslope,[1,chans,bands,dim4,dim5,dim6]));
    if evenlen == 0
        y = y(1:end-1,:,:,:,:,:);
    end
    
    if isstruct(in)
        OUT.audio = y;
        OUT.funcallback.name = 'tilt_spectrum_logf.m';
        OUT.funcallback.inarg = {fs,dBperOct,lowf,hif,pivotf,doplay};
    else
        OUT = y;
    end

    if doplay
        % Loop for replaying, saving and finishing
        choice = 0;

        % loop until the user presses the 'Done' button
        while choice < 4
            choice = menu('What next?', ...
                'Play audio', ...
                'Plot spectrogram', 'Save wav file', 'Discard', 'Done');
            switch choice
                case 1
                    ysumbands = sum(y,3);
                    sound(ysumbands./max(max(abs(ysumbands))),fs)
                case 2
                    spectrogram_simple(y,fs);
                case 3
                    [filename, pathname] = uiputfile({'*.wav'},'Save as');
                    ysumbands = sum(y,3);
                    if max(max(abs(ysumbands)))>1
                        ysumbands = ysumbands./max(max(abs(ysumbands)));
                        disp('Wav data has been normalized to avoid clipping')
                    end
                    if ischar(filename)
                        audiowrite([pathname,filename], ysumbands, fs);
                    end
                case 4
                    y = [];
            end
        end
    end
else
    OUT = [];
end