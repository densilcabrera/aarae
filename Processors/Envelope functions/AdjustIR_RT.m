function [OUT, varargout] = AdjustIR_RT(IN,T,freq,iterations,autocrop,openclosed,Tevalrange,filterstrength,AbsRelT,GuideChan,fs)
% This function is designed to change the reverberation time of an impulse
% response (probably a room impulse response) to desired target band values
% using octave or 1/3-octave band processing (or other bandwidths if the
% audio is pre-filtered).
%
% Some background to this function is given in: D. Cabrera, D. Lee, M.
% Yadav and W.L. Martens, "Decay envelope manipulation of room impulse
% responses: Techniques for auralization and sonification," Proceedings of
% Acoustics 2011 (Australian Acoustical Society Conference), Gold Coast,
% Australia, 2011.
%
% IMPORTANT: This function requires Matlab's Curve Fitting Toolbox.
%
%
% INPUTS: 
% IN is an impulse response (either an AARAE audio structure, or a vector
% or matrix). If dimension 3 of the audio data is non-singleton, it is
% assumed that the audio has been pre-filtered (with bands in dimension 3),
% and so this function does not filter the audio in that case. The function
% is fully compatible with AARAE's audio dimension conventions (allowing up
% to 6-dimensional audio). In the case of dimensions 2, 4, 5, and/or 6
% being non-singleton, reverberation time is calculated based on the
% power-sum of the simultaneous decays. This ensures that multi-channel
% impulse responses have the same adjustment made to all channels (and the
% target reverberation time represents the reverberation time of the
% combined channels).
%
% Other inputs are entered via a dialog box if only the first input
% argument is present.
%
% T is a vector containing the target reverberation times (T20) in seconds.
% Note that if T values are unreasonably long (for the input impulse
% response) then the processing might not succeed.
%
% freq is a vector containing the band centre frequencies corresponding to
% the target reverberation times. Please only use standard nominal centre
% frequencies for octave and 1/3-octave bands. This input is not used if
% the audio input is multiband.
%
% iterations specifies the number of iterations used to refine the
% reverberation time adjustment. Usually more than 1 iteration is required
% because the changing the reverberation time changes the length of the
% vector used for T evaluation.
%
% autocrop specifies whether AARAE's autocropstart_aarae function is run
% (to remove samples before the start of the impulse response). Usually
% this is a helpful operation. Autocrop is done using 'method' 3, where the
% same shift is applied to all columns, based on the earliest arival.
%
% openclosed specifies whether the lowest and highest filters are bandpass
% filters (closed, 1) or, respectively, lopass and hipass filters (open,
% 0). Using closed-end filters provides more control, but bandpasses the
% frequency range of the adjusted impulse response.
%
% Tevalrange is the evaluation range in decibels over which reverberation
% time is evaluated (at all stages of the processing). It can be dangerous
% to use a large evaluation range (unless signal to noise ratio is
% excellent), and 20 dB is the suggested value for Tevalrange.
%
% filterstrength is a factor controling the frequency selectivity of the
% filters. A value of 1 gives the equivalent of 12th order filters, a value
% of 2 gives the equivalent of 24th order filters. If there are large
% changes in target reverberation time between bands, then greater filter
% strength may be beneficial. Furthermore, 1/3-octave band filters will
% may often be better with a high filterstrength. However, if reverberation
% time is short, then the filters' own impulse response may introduce
% unwanted artefacts if the filterstrength is too great.
%
% AbsRelT allows the target reverberation times to be specified either in
% absolute terms (reverberation time in seconds) or in relative terms (as a
% factor change in reverberation time relative to the original). Use 0 for
% absolute; use 1 for relative.
%
% This version does filtering in the frequency domain (whereas the version
% described in the paper did time domain filtering). The frequency domain
% filters are adapted from AARAE's fft-based filterbanks, which are very 
% well behaved even at very high equivalent orders. Filtering is zero
% phase so that the impulse response can be reconstructed from the bands
% with minimal loss of its fine temporal structure.
%
% The function will only be effective if it is used intelligently. The
% input impulse response should have a high signal-to-noise ratio,
% especially if reverberation time is being increased. While this function
% does include noise-floor extrapolation, it is best not to rely on that
% for important parts of the impulse response. Some experimentation may be
% required to achieve useful results. If possible, use a high bit depth
% (e.g. 64 bit, rather than 16-bit) for the IR (this will always be the
% case if the IR is measured in AARAE, but may not be the case for imported
% IRs) - this should avoid some occasional curve-fitting issues.
%
% The reverberation times displayed by this function are indicative -
% please check the results with a more comprehensive reverberation time
% analyser. However, don't forget that the analysis introduces its own
% artefacts (e.g. bandpass filter selectivity) which may cause apparent
% differences from the expected values.
%
% Code by Densil Cabrera
% version 2, 9 Sept 2019

v = ver('curvefit');
if isempty(v)
    h=warndlg('This function requires Matlab''s Curve Fitting Toolbox.','AARAE info','modal');
        uiwait(h)
        OUT = []; % you need to return an empty output
        return % get out of here!
end

if nargin == 1
    if size(IN.audio,3)==1 % if audio is multiband already, then we don't have to filter
        param = inputdlg({'Bands per octave [1 | 3]';...
            'Highest centre frequency (Hz)';...
            'Lowest centre frequency (Hz)';...
            'Number of iterations';...
            'Auto-crop start? [0 | 1]';...
            'Open [0] or closed [1] end filters at the extremes';...
            'Reverberation time evaluation range (dB)';...
            'Absolute or relative reverberation time adjustment [0 | 1]';...
            'Guide channels [0 for all, or specify channel(s)]';...
            'Filter strength factor'},...
            'Reverberation Time Adjustment Settings',... % dialog window title.
            [1 60],...
            {'1';'8000';'125';'4';'1';'0';'20';'0';'0';'2'}); % preset answers.
        
        param = str2num(char(param));
        
        if length(param) < 10, param = []; end 
        if ~isempty(param) 
            bandwidths = param(1);
            fhigh = param(2);
            flow = param(3);
            iterations = param(4);
            autocrop = param(5);
            openclosed = param(6);
            Tevalrange = abs(param(7));
            AbsRelT = abs(param(8));
            GuideChan = param(9);
            filterstrength = param(10);
            noctaves = log2(fhigh/flow); % number of octaves
            if bandwidths == 1
                % octave bands
                freq = flow .* 2.^(0:round(noctaves));
            else
                % 1/3-octave bands
                freq = flow .* 2.^(0:1/3:round(3*noctaves)/3);
                bandwidths = 3;
            end
            
            if length(freq) < 3
                openclosed = 1;
                % closed-ended filters are used if there are less than 3
                % bands. (Maybe this is not necessary)
            end
        else
            % get out of here if the user presses 'cancel'
            OUT = [];
            return
        end
    else
        % audio has been pre-filtered
        
        param = inputdlg({'Number of iterations';...
            'Auto-crop start? [0 | 1]';...
            'Reverberation time evaluation range (dB)';...
            'Absolute or relative reverberation time adjustment [0 | 1]';...
            'Guide channels [0 for all, or specify channel(s)]';},...
            'Reverberation Time Adjustment Settings',... % dialog window title.
            [1 60],...
            {'4';'1';'20';'0';'0'}); % preset answers.
        
        param = str2num(char(param));
        
        if length(param) < 5, param = []; end
        if ~isempty(param)
            
            iterations = param(1);
            autocrop = param(2);
            Tevalrange = abs(param(3));
            AbsRelT = abs(param(4));
            GuideChan = param(5);
        else
            % get out of here if the user presses 'cancel'
            OUT = [];
            return
        end
        if ~exist('freq','var')
            freq = [];
        end
    end
else
    param = [];
end

if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
    if size(audio,3) > 1
        if isfield(IN,'bandID') % Get the band ID if it exists
            freqnom = IN.bandID;
        else
            % asssign ordinal band numbers if bandID does not exist
            freqnom = 1:size(audio,3);
        end
    end
else
    audio = IN;
end

if autocrop
    % autocrop method 3 crops based on the shortest latency channel
    % autocrop method 2 (not used) crops based on the mix of channels
    audio = autocropstart_aarae(audio,-20,3);
end

bands = size(audio,3);

if bands ==1
    freqnom = exact2nom_oct(freq); % this is an AARAE utility
    bands = length(freqnom);
end

if isempty(bands)
    h=warndlg('Input error: no frequency bands specified by input settings.','AARAE info','modal');
        uiwait(h)
        OUT = []; % you need to return an empty output
        return % get out of here!
end


% ***********************************
% DIALOG BOX FOR REVERBERATION TIMES
if ~exist('T','var')
    T = ones(1,bands);
    for n = 1:bands
        freqcell{n} = num2str(freqnom(n));
        defaultcell{n} = '1';
    end
    
    maxfield = 15;
    
    if bands >= maxfield
        freqcell1 = freqcell(1:maxfield);
        defaultcell1 = defaultcell(1:maxfield);
        if bands < 2*maxfield
            freqcell2 = freqcell(maxfield+1:end);
            defaultcell2 = defaultcell(maxfield+1:end);
        else
            freqcell2 = freqcell(maxfield+1:2*maxfield);
            defaultcell2 = defaultcell(maxfield+1:2*maxfield);
            freqcell3 = freqcell(2*maxfield+1:end);
            defaultcell3 = defaultcell(2*maxfield+1:end);
        end
    else
        freqcell1 = freqcell;
        defaultcell1 = defaultcell;
    end
    if ~AbsRelT
    param = inputdlg(freqcell1,...
        'Reverberation Times',... % dialog window title.
        [1 60],...
        defaultcell1); % preset answers for dialog.
    else
        param = inputdlg(freqcell1,...
        'Reverb Change Factor',... % dialog window title.
        [1 60],...
        defaultcell1); % preset answers for dialog.
    end
    
    param = str2num(char(param));
    
    if length(param) < min([bands,maxfield]), param = []; end
    if ~isempty(param)
        if bands > maxfield
            T(1:maxfield) = param;
        else
            T = param;
        end
    end
    if bands>=maxfield
        if ~AbsRelT
        param = inputdlg(freqcell2,...
            'Reverberation Times',... % dialog window title.
            [1 60],...
            defaultcell2); % preset answers for dialog.
        else
            param = inputdlg(freqcell2,...
            'Reverb Change Factor',... % dialog window title.
            [1 60],...
            defaultcell2); % preset answers for dialog.
        end
        param = str2num(char(param));
        
        if length(param) < (bands-maxfield), param = []; end
        if ~isempty(param)
            T(maxfield+1:maxfield+length(param)) = param;
        end
    end
    if bands>=2*maxfield
        if ~AbsRelT
        param = inputdlg(freqcell3,...
            'Reverberation Times',... % dialog window title.
            [1 60],...
            defaultcell3); % preset answers for dialog.
        else
            param = inputdlg(freqcell3,...
            'Reverb Change Factor',... % dialog window title.
            [1 60],...
            defaultcell3); % preset answers for dialog.
        end
        param = str2num(char(param));
        
        if length(param) < (bands-2*maxfield), param = []; end
        if ~isempty(param)
            T(2*maxfield+1:end) = param;
        end
    end
end
% ***********************************

[r,c]=size(T);
if r > 1 && c == 1
    T = T';
end




if ~isempty(audio) && ~isempty(fs) && ~isempty(T)
    
    [len,chans,bandsin,dim4,dim5,dim6] = size(audio);
    
    
    
    % If bandsin > 1, then we will process using these instead of filtering.
    % otherwise, filter into bands here
    if bandsin == 1
        if ~exist('bandwidths','var')
            if round(mean(diff(log2(freqnom)))) == 1
                bandwidths = 1;
            else
                bandwidths = 3;
            end
        end
        audio = filterbankwithopenends(audio,fs,bandwidths,freqnom,openclosed,filterstrength);
    end
    
    
    
    
    % ESTIMATE THE POINT AT WHICH NOISE OVERWHELMS THE RIR IN EACH BAND
    % derive a smoothed envelope function for each band
    Nyquist = fs/2;
    lopassfreq = 4; % smoothing filter cutoff frequency in Hz
    halforder = 1; % smoothing filter order
    [num, den] = butter(halforder, lopassfreq/Nyquist, 'low');
    
    for ch = 1:chans
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
                    
                    %envelopes = 10*log10(filtfilt(num, den, audio .^2));
                    envelopes = filtfilt(num, den,...
                        permute(10*log10(audio(:,ch,:,d4,d5,d6) .^2 + 1e-300),[1 3 2]));
                    if bands>2
                        for n = 2:bands-1
                            envelopes(:,1:end-n) = filtfilt(num, den, envelopes(:,1:end-n));
                        end
                    end
                    envelopes = envelopes - repmat(max(envelopes),len,1); % make max = 0 dB
                    maxsample = zeros(1, bands);
                    a = zeros(1, bands);
                    b = zeros(1, bands);
                    for band = 1:bands
                        maxsample(band) = find(envelopes(:,band) == 0, 1, 'last');
                        times = (0:length(envelopes) - maxsample(band))./fs;
                        s = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower',[-1000,0],...
                            'Upper',[0,max(times)],...
                            'Startpoint',[1 1]);
                        f = fittype('10*log10(10^(a*x/10)+b)','options',s);
                        [c,gof] = fit(times',envelopes(maxsample(band):end,band),f);
                        a(band) = c.a;
                        b(band) = c.b;
                        if false
                            figure
                            plot(times',envelopes(maxsample(band):end,band))
                            hold on
                            plot(c,'r')
                            hold off
                        end
                        
                        % remove noise floor
                        noiseplus10dB = floor(((10*log10(b(band))+10)...
                            ./ a(band)) .* fs);
                        
                        snr = (10.^(a(band).*times./10)...
                            ./ (10.^(a(band).*times./10) +b(band))) .^ 0.5;
                        
                        audio(noiseplus10dB+maxsample(band):end,ch,band,d4,d5,d6) =...
                            audio(noiseplus10dB+maxsample(band):end,ch,band,d4,d5,d6)...
                            .* snr(noiseplus10dB+1:end)';
                    end
                end
            end
        end
    end
    
    
    tabledata = zeros(iterations+1,bands);
    
    % adjust reverberation times
    for i = 1:iterations
        
        % CALCULATE Tx
        RT = zeros(1,bands);
        x = Tevalrange;
        %Schroeder reverse integration
        for band = 1:bands
            if GuideChan == 0
            decaycurve = flipud(10*log10(cumsum(flipud(...
                mean(mean(mean(mean(audio(:,:,band,:,:,:),2),4),5),6).^2))+1e-300));
            else
                decaycurve = flipud(10*log10(cumsum(flipud(...
                mean(mean(mean(mean(audio(:,GuideChan,band,:,:,:),2),4),5),6).^2))+1e-300));
            end
            % make IR start time 0 dB
            decaycurve = decaycurve - decaycurve(1);
            
            Tstart = find(decaycurve <= -5, 1, 'first'); % -5 dB
            Tend = find(decaycurve <= -x-5, 1, 'first'); % -x-5 dB
            p = polyfit((Tstart:Tend)', decaycurve(Tstart:Tend),1); %linear regression
            RT(band) = 60/x*((p(2)-x-5)/p(1) - (p(2)-5)/p(1))/fs; % reverberation time
        end
        if AbsRelT && i==1
            T = T.*RT;
        end
        %if i ==1, disp(['T20            (s)  ', num2str(RT)]); end
        tabledata(i,:) = RT;
        tau1 = RT ./ (log(60)); % exponential decay constants of original
        tau2 = T ./ (log(60)); % desired exponential decay constants
        
        gainfunction = permute(exp((repmat((0:len-1)',1, bands)./fs)...
            ./repmat(tau1,len,1)+ (repmat(-(0:len-1)',1, bands)./fs)...
            ./repmat(tau2,len,1)),[1 3 2]);
        audio = audio .* repmat(gainfunction,[1 chans 1 d4 d5 d6]);
    end
    
    
   % RIRoct = permute(mean(mean(mean(mean(audio,2),4),5),6),[1 3 2]);
    % CALCULATE T OF CHANGED WAVE
    Tfinal = zeros(bands,1);
    %Schroeder reverse integration
    for band = 1:bands
        %decaycurve = flipud(10*log10(cumsum(flipud(RIRoct(:,band).^2))+1e-300));
        if GuideChan == 0
            decaycurve = flipud(10*log10(cumsum(flipud(...
                mean(mean(mean(mean(audio(:,:,band,:,:,:),2),4),5),6).^2))+1e-300));
        else
            decaycurve = flipud(10*log10(cumsum(flipud(...
                mean(mean(mean(mean(audio(:,GuideChan,band,:,:,:),2),4),5),6).^2))+1e-300));
        end
        % make IR start time 0 dB
        decaycurve = decaycurve - decaycurve(1);
        Tstart = find(decaycurve <= -5, 1, 'first'); % -5 dB
        Tend = find(decaycurve <= -x-5, 1, 'first'); % -x dB
        p = polyfit((Tstart:Tend)', decaycurve(Tstart:Tend),1); %linear regression
        Tfinal(band) = 60/x*((p(2)-x-5)/p(1) - (p(2)-5)/p(1))/fs; % reverberation time
    end
    %disp(['T',num2str(x),' adjusted   (s)  ', num2str(Tfinal')])
    tabledata(end,:) = Tfinal';
    
    fig1 = figure('Name','Reverberation time evolution from original to final');
    table1 =  uitable('Data',tabledata,...
                      'RowName',cellstr(num2str((1:iterations+1)')),...
                      'ColumnName',cellstr(num2str(freqnom')));
%       disptables(fig1,table1);
    [~,table] = disptables(fig1,table1);
    
    
    if isstruct(IN)
        OUT = IN;
        OUT.audio = sum(audio,3);
        % OUT.tables = table; % this is commented out because in the
        % current version of AARAE it will result in the audio field being
        % empty
        OUT.funcallback.name = 'AdjustIR_RT.m';
        OUT.funcallback.inarg = {T,freq,iterations,autocrop,openclosed,Tevalrange,filterstrength,AbsRelT,GuideChan,fs};
    else
        OUT = sum(audio,3);
    end
    varargout{1} = fs;
else
    OUT = [];
end
end










function filtered = filterbankwithopenends(audio,fs,bpo,frequencies,openclosed,filterstrength)
% adapted from AARAE's octbandfilter_viaFFT.m, but with open ended top and
% bottom bands.


phasemode = 0; % zero phase filter
orderin = 12*filterstrength;  % almost flat passband
orderout = 12*filterstrength; % high filter selectivity
zeropad = 0;
base = 10;


if bpo == 1
    % potential nominal centre frequencies
    nominalfreq = [16,31.5,63,125,250,500,1000, ...
        2000,4000,8000,16000,31500,63000,125000,250000,500000,1000000];
    if base == 10
        exactfreq = 10.^((12:3:60)/10);
    elseif base == 2
        exactfreq = 1000.* 2.^(-6:10);
    else
        exactfreq = nominalfreq; % not a good idea, but included to demonstrate this!
    end
    % possible nominal frequencies
    %         nominalfreq = nominalfreq(exactfreq <= maxfrq);
    %         exactfreq = exactfreq(exactfreq <= maxfrq);
    
else
    % potential nominal centre frequencies
    nominalfreq = [12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,...
        500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,...
        8000,10000,12500,16000,20000,25000,31500,40000,50000,63000,...
        80000,100000,125000,160000,200000,250000,315000,400000,500000,...
        630000,800000,1000000];
    
    exactfreq = 10.^((11:60)/10);
    % possible nominal frequencies
    % nominalfreq = nominalfreq(exactfreq <= maxfrq);
    % exactfreq = exactfreq(exactfreq <= maxfrq);
    
end

S = zeros(size(frequencies));
for i = 1:length(frequencies)
    check = find(nominalfreq == frequencies(i));
    if isempty(check), check = 0; end
    S(i) = check;
end
if all(S)
    exactfreq = exactfreq(S);
    frequencies = sort(frequencies,'ascend');
    exactfreq = sort(exactfreq,'ascend');
    ok = 1;
else
    ok = 0;
    disp('unexpected frequencies')
    % more code needed for this contingency?
end

[~,chans,~,dim4,dim5,dim6] = size(audio);
if zeropad > 0
    audio = [zeros(zeropad,chans,1,dim4,dim5,dim6); audio; zeros(zeropad,chans,1,dim4,dim5,dim6)];
end

len = size(audio,1);
filtered = zeros(len,chans,length(frequencies),dim4,dim5,dim6);

minfftlenfactor = 1000;
minfftlen = 2.^nextpow2((fs/min(frequencies)) * minfftlenfactor);

if len >= minfftlen
    fftlen = 2.^nextpow2(len);
else
    fftlen = minfftlen;
end

audio = fft(audio, fftlen);


for b = 1: length(frequencies)
    
    % list of fft component frequencies
    f = ((1:fftlen)'-1) * fs / fftlen;
    
    
    % index of low cut-off
    if base == 10
        flo = exactfreq(b) / 10.^(0.15/bpo);
    else
        flo = exactfreq(b) / 2.^(0.5/bpo);
    end
    
    
    % index of high cut-off
    if base == 10
        fhi = exactfreq(b) * 10.^(0.15/bpo);
    else
        fhi = exactfreq(b) * 2.^(0.5/bpo);
    end
    
    
    indlo = find(abs(f(1:end/2)-flo) == min(abs(f(1:end/2)-flo)),1,'first');
    indhi = find(abs(f(1:end/2)-fhi) == min(abs(f(1:end/2)-fhi)),1,'first');
    
    
    % centre frequency index
    indfc = find(abs(f(1:end/2)-exactfreq(b)) ...
        == min(abs(f(1:end/2)-exactfreq(b))),1,'first');
    
    
    
    % magnitude envelope
    
    mag = zeros(fftlen,1); % preallocate and set DC to 0
    
    if b > 1 || openclosed
        % below centre frequency
        % out-of-band, using exact 6n dB/oct skirts
        mag(2:indlo-1) = ...
            (f(2:indlo-1)./ flo ).^(orderout) ./2.^0.5;
        % the following alternative uses Butterworth skirts
        %(1 ./ (1 + (f(2:indlo-1)./ flo ).^(-2*orderout))).^0.5;
        
        % in-band
        mag(indlo:indfc-1) = ...
            (1 ./ (1 + (f(indlo:indfc-1)./ flo ).^(-2*orderin))).^0.5;
    else
        mag(2:indfc-1) = 1;
    end
    
    if b < length(frequencies) || openclosed
        % from centre frequency to Nyquist frequency
        % in-band
        mag(indfc:indhi) = ...
            (1 ./ (1 + (f(indfc:indhi)./ fhi ).^(2*orderin))).^0.5;
        
        % out-of-band, using exact 6n dB/oct skirts
        mag(indhi+1:fftlen/2+1) = ...
            (f(indhi+1:fftlen/2+1) ./ fhi).^(-orderout) ./ 2.^0.5;
        % the following alternative uses Butterworth skirts
        %(1 ./ (1 + (f(indhi+1:fftlen/2+1)./ fhi ).^(2*orderout))).^0.5;
    else
        mag(indfc:fftlen/2+1) = 1;
    end
    
    % normalize gain
    mag = mag ./ mag(indfc);
    
    % above Nyquist frequency
    mag(fftlen/2+2:end) = flipud(mag(2:fftlen/2));
    
    if (phasemode == 1) || (phasemode == 11)
        % convert mag to min phase complex coefficients
        mag = minphasefreqdomain(mag);
    elseif (phasemode == -1) || (phasemode == -11)
        % convert mag to max phase complex coefficients
        mag = conj(minphasefreqdomain(mag));
        
    end
    
    if (phasemode == -11) || (phasemode == 10) || (phasemode == 11)
        % zero the upper half of the spectrum for quadrature (complex) filters
        mag(fftlen/2:end) = 0;
    end
    
    if (phasemode == 0) || (phasemode == 10)
        bandfiltered = ifft(repmat(mag,[1,chans,1,dim4,dim5,dim6]) .* audio);
    elseif (phasemode == -1) || (phasemode == 1)
        % real output only for min phase and max phase
        bandfiltered = real(ifft(repmat(mag,[1,chans,1,dim4,dim5,dim6]) .* audio));
    elseif (phasemode == -11) || (phasemode == 11)
        % quadrature min and max phase
        bandfiltered = ifft(repmat(mag,[1,chans,1,dim4,dim5,dim6]) .* audio);
    else
        disp('Phasemode value not recognized');
        
        return
    end
    
    
    
    % truncate waveform and send to filtered waveform matrix
    
    filtered(:,:,b,:,:,:) = bandfiltered(1:len,:,1,:,:,:);
end

end

%**************************************************************************
% Copyright (c) 2015, Densil Cabrera
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%**************************************************************************
