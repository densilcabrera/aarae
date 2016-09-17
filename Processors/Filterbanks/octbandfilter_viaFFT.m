function [OUT,varargout] = octbandfilter_viaFFT(IN,fs,param,order,zeropad,minfftlenfactor,test,phasemode,base)
% This function does zero, linear, minimum and maximum phase filtering, 
% using a single large fft.
%
% Rather than brick-wall filters, this implementation
% has somewhat Butterworth-like magnitude responses. Hence the 'order' input
% argument determines the slope of the filter skirts, and the flatness of 
% the in-band response. If desired, a different order can be used for the 
% in-band and out-of-band frequency ranges. When calling this function
% directly, this can be done by supplying a two value vector instead of a
% single value for 'order' (the first value is for in-band, the second for
% out-of-band).
%
% Filters that meet IEC class 0 criteria can be designed with appropriate
% selection of filter order. The minimum in-band order that meets class 0 
% criteria is 7, and the minimum out-of-band order is 5 (based on testing 
% with a 2 s duration waveform at fs = 48 kHz). Note that the tested filter 
% performance at low frequencies depends on the test signal duration.
%
% param is a list of octave band centre frequencies (nominal freqencies)
%
% phasemode = 0 for zero phase (or linear phase if a delay is used)
% phasemode = 1 for minimum phase
% phasemode = -1 for maximum phase
% phasemode = 10 for zero phase quadrature (complex) filter
% phasemode = 11 for minimum phase quadrature (complex) filter
% phasemode = -11 for maximum phase quadrature (complex) filter
%
% zeropad is the number of samples added both before and after the input
% audio, to capture the filters' (acausal) build-up and decay. Hence a
% non-zero zeropad value will result in a delay, making the filter linear
% rather than zero phase.
%
% If test == 1, then a plot is produced for each band, showing the tested
% filter magnitude response together with IEC 61260 limit curves for class
% 0, class 1 and class 2 octave band filters. These plots are derived from
% an analysis of an impulse centred in time using a audio wave that has the
% same number of samples as the analysed wave.
%
% base is only available at present by calling the function with it as an
% input argument (it is not available in dialog box). This allows base 2 or
% base 10 octave band filters to be used (slightly changing the bandwidth
% and centre frequencies). Enter 2 or 10 (base 10 is default).
%
% Code by Densil Cabrera
% Version 1.06 (24 September 2014)

if ~exist('test','var')
    test = 0;
end

% Spectral resolution can be increased (or reduced) by adjusting the value
% of minfftlenfactor (below). Increasing it increases the minimum FFT 
% length, which may increase computational load (processing time and memory
% utilization). The actual FFT length is always an integer power of 2.

% minimum number of spectral compoments up to the lowest octave band fc:
if ~exist('minfftlenfactor','var')
    minfftlenfactor = 1000; % adjust this to get the required spectral resolution
else
    minfftlenfactor = round(abs(minfftlenfactor));
end

if ~exist('base','var')
    base = 10; % base 10 octave band frequencies by default
end

if ~exist('phasemode','var')
    phasemode = 0;
end


if ~exist('zeropad','var')
    zeropad = 0;
else
    zeropad = round(abs(zeropad));
end

if ~exist('order','var')
    orderin = 12; % default filter in-band pseudo-order
    orderout = 12; % default filter out-of-band pseudo-order
else
    if length(order) == 1
        orderin = abs(order); 
        orderout = abs(order);
    elseif length(order) == 2
        orderin = abs(order(1)); 
        orderout = abs(order(2));
    else
        error('order input argument is incorrectly dimensioned')
    end
end

if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
else
    audio = IN;
    if ~exist('fs','var')
        error('fs must be specified')
    end
end

TOOBIG = 1e7;
if numel(audio) >= TOOBIG;
    if nargin == 1
        warndlg('This audio input is probably too big for octbandfilter_viaFFT. Try AARAE''s octbandfilter processor instead.')
        OUT = [];
        varargout = {};
        return
    else
        % automatically use octbandfilter instead if this function was called by another
        if ~exist('param','var')
            param = [31.5 63 125 250 500 1000 2000 4000 8000 16000];
        end
        if ~exist('phasemode','var')
            phasemode = 1;
        end
        [OUT,centerf] = octbandfilter(audio,fs,param,phasemode);
        varargout = {centerf};
        return
    end
end


maxfrq = fs / 2.^1.51; % maximum possible octave band centre frequency


if exist('param','var')
    % nominal frequencies only allowed
    param = exact2nom_oct(param);
    % minimum and maximum allowable frequencies
    param = param(param >=16  & param <= maxfrq);
    
end

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
nominalfreq = nominalfreq(exactfreq <= maxfrq);
exactfreq = exactfreq(exactfreq <= maxfrq);


if nargin < 3
    param = nominalfreq;
    [S,ok] = listdlg('Name','Octave band filter input parameters',...
        'PromptString','Center frequencies [Hz]',...
        'ListString',[num2str(param') repmat(' Hz',length(param),1)]);
    param = param(S);
    exactfreq = exactfreq(S);
else
    S = zeros(size(param));
    for i = 1:length(param)
        check = find(nominalfreq == param(i));
        if isempty(check), check = 0; end
        S(i) = check;
    end
    if all(S)
        exactfreq = exactfreq(S);
        param = sort(param,'ascend');
        exactfreq = sort(exactfreq,'ascend');
        ok = 1; 
    else
        ok = 0; 
    end
end

if ok == 1 && isstruct(IN) && nargin < 4
    param1 = inputdlg({'Pseudo-Butterworth filter in-band order';... 
        'Filter out-of-band order';... 
        'Number of samples to zero-pad before and after the wave data';...
        'Zero phase [0], Minimum phase [1] or Maximum phase [-1]';...
        'Test filter response [0|1]'},...
        'Filter settings',... 
        [1 60],... 
        {num2str(orderin);num2str(orderout);num2str(zeropad);...
        num2str(phasemode);'0'}); 
    
    param1 = str2num(char(param1)); 
    
    if length(param1) < 5, param1 = []; end 
    
    if ~isempty(param1) 
        orderin = param1(1);
        orderout = param1(2);
        zeropad = param1(3);
        phasemode = param1(4);
        test = param1(5);
    end
end




    
% if the input is already multiband, then mixdown first
if size(audio,3) > 1
    audio = sum(audio,3);
    disp('Multi-band audio has been mixed-down prior to octave band filtering')
end

if ok == 1
    [~,chans,~,dim4,dim5,dim6] = size(audio);
    if zeropad > 0
        audio = [zeros(zeropad,chans,1,dim4,dim5,dim6); audio; zeros(zeropad,chans,1,dim4,dim5,dim6)];
    end
    
    len = size(audio,1);
    filtered = zeros(len,chans,length(param),dim4,dim5,dim6);
        
    if test == 1
        testaudio = zeros(len,1);
        testaudio(round(len/2)) = 1;
        testaudioout = zeros(len,1,length(param));
    end
    
    
    minfftlen = 2.^nextpow2((fs/min(param)) * minfftlenfactor);
    
    if len >= minfftlen
        fftlen = 2.^nextpow2(len);
    else
        fftlen = minfftlen;
    end
    
    audio = fft(audio, fftlen);
    if test == 1
        testspectrum = fft(testaudio, fftlen);
    end
    
    for b = 1: length(param)
        
        % list of fft component frequencies
        f = ((1:fftlen)'-1) * fs / fftlen;
        
        % index of low cut-off
        if base == 10
            flo = exactfreq(b) / 10.^0.15;
        else
            flo = exactfreq(b) / 2.^0.5;
        end
        indlo = find(abs(f(1:end/2)-flo) == min(abs(f(1:end/2)-flo)),1,'first');
        
        % index of high cut-off
        if base == 10
            fhi = exactfreq(b) * 10.^0.15;
        else
            fhi = exactfreq(b) * 2.^0.5;
        end
        indhi = find(abs(f(1:end/2)-fhi) == min(abs(f(1:end/2)-fhi)),1,'first');
        
        % centre frequency index
        indfc = find(abs(f(1:end/2)-exactfreq(b)) ...
            == min(abs(f(1:end/2)-exactfreq(b))),1,'first');
        
       
        
        % magnitude envelope
        
        mag = zeros(fftlen,1); % preallocate and set DC to 0
        
         % below centre frequency
         % out-of-band, using exact 6n dB/oct skirts
         mag(2:indlo-1) = ...
             (f(2:indlo-1)./ flo ).^(orderout) ./2.^0.5;
         % the following alternative uses Butterworth skirts
             %(1 ./ (1 + (f(2:indlo-1)./ flo ).^(-2*orderout))).^0.5;
         
         % in-band
         mag(indlo:indfc-1) = ...
             (1 ./ (1 + (f(indlo:indfc-1)./ flo ).^(-2*orderin))).^0.5;
        
         % from centre frequency to Nyquist frequency
         % in-band
         mag(indfc:indhi) = ...
            (1 ./ (1 + (f(indfc:indhi)./ fhi ).^(2*orderin))).^0.5;
        
        % out-of-band, using exact 6n dB/oct skirts
        mag(indhi+1:fftlen/2+1) = ...
            (f(indhi+1:fftlen/2+1) ./ fhi).^(-orderout) ./ 2.^0.5;
        % the following alternative uses Butterworth skirts
            %(1 ./ (1 + (f(indhi+1:fftlen/2+1)./ fhi ).^(2*orderout))).^0.5;
        
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
            OUT = [];
            return
        end
        
        
        
        % truncate waveform and send to filtered waveform matrix

            filtered(:,:,b,:,:,:) = bandfiltered(1:len,:,1,:,:,:);

        
        
        if test
            testaudiofiltered = ifft(testspectrum .* mag);
            testaudioout(:,1,b) = testaudiofiltered(1:len,:);
        end
        
    end
    clear audio bandfiltered mag testaudiofiltered
    
    if test == 1
        testlevel = 10*log10(abs(fft(testaudioout)).^2);
        frequencies = fs .* ((1:length(testlevel))-1) ./ length(testlevel);
        class0uppergain = [0.15;
                            0.15;
                            0.15;
                            0.15;
                            0.15;
                            -2.3;
                            -18.0;
                            -42.5;
                            -62;
                            -75];
        class0uppergain = [flipud(class0uppergain(2:end));class0uppergain];                
        class0lowergain = [-0.15;
                            -0.2;
                            -0.4;
                            -1.1;
                            -4.5;
                            -4.5;
                            -1e10; % should be -inf, but this plots better
                            -1e10;
                            -1e10;
                            -1e10]; 
        class0lowergain = [flipud(class0lowergain(2:end));class0lowergain];
        
        class1uppergain = [0.3;
                            0.3;
                            0.3;
                            0.3;
                            0.3;
                            -2;
                            -17.5;
                            -42;
                            -61;
                            -70];
        class1uppergain = [flipud(class1uppergain(2:end));class1uppergain];                
        class1lowergain = [-0.3;
                            -0.4;
                            -0.6;
                            -1.3;
                            -5;
                            -5;
                            -1e10; % should be -inf, but this plots better
                            -1e10;
                            -1e10;
                            -1e10]; 
        class1lowergain = [flipud(class1lowergain(2:end));class1lowergain];
        
        class2uppergain = [0.5;
                            0.5;
                            0.5;
                            0.5;
                            0.5;
                            -1.6;
                            -16.5;
                            -41;
                            -55;
                            -60];
        class2uppergain = [flipud(class2uppergain(2:end));class2uppergain];                
        class2lowergain = [-0.5;
                            -0.6;
                            -0.8;
                            -1.6;
                            -5.5;
                            -5.5;
                            -1e10; % should be -inf, but this plots better
                            -1e10;
                            -1e10;
                            -1e10]; 
        class2lowergain = [flipud(class2lowergain(2:end));class2lowergain];
        
        
        for b = 1:length(param)
            testfreq = exactfreq(b) .* 10.^(3/10).^ ...
            [-4;-3;-2;-1;-1/2;-1/2;-3/8;-1/4;-1/8;0;...
                1/8;1/4;3/8;1/2;1/2;1;2;3;4];
            figure('Name',[num2str(param(b)), ' Hz Band'])
            
            semilogx(testfreq,class2uppergain,'Color',[1,0.8,0.8]); 
            hold on
            semilogx(testfreq,class2lowergain,'Color',[1,0.8,0.8]);
            text(testfreq(end-1),-53,'Class 2','Color',[1,0.8,0.8]);
            
            semilogx(testfreq,class1uppergain,'Color',[1,0.5,0.5]);
            semilogx(testfreq,class1lowergain,'Color',[1,0.5,0.5]);
            text(testfreq(end-1),-59,'Class 1','Color',[1,0.5,0.5]);
            
            semilogx(testfreq,class0lowergain,'r');
            semilogx(testfreq,class0uppergain,'r');
            text(testfreq(end-1),-73,'Class 0','Color','r');
            
            semilogx(frequencies(2:round(end/2)), ...
                testlevel(2:round(end/2),1,b),'b');
            semilogx([frequencies(2),1e99],[-3,-3],':k')
            text(testfreq(end-1),-1,'-3 dB','Color','k');
            %legend('show','Location','EastOutside');
            xlim([testfreq(1),testfreq(end)])
            ylim([-90,10])
            xlabel('Frequency (Hz)')
            ylabel('Gain (dB)')
            hold off
        end
    end
    
    if isstruct(IN) && ~isempty(filtered)
        OUT = IN;
        OUT.audio = filtered;
        OUT.bandID = param;
        OUT.funcallback.name = 'octbandfilter_viaFFT.m';
        OUT.funcallback.inarg = {fs,param,[orderin orderout],zeropad,minfftlenfactor,test,phasemode};
    else
        OUT = filtered;
    end
    varargout{1} = param;
else
    OUT = [];
end
end % eof