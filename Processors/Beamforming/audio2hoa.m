function OUT = audio2hoa(IN, fs, mic_coords, cropIR , max_order,filterlength, micType)
% This function reduces raw recordings made from a spherical microphone
% array to higher order Ambisonics (HOA) format, to allow them to be used
% and analysed with processors and analysers that use this generic spatial
% audio format.
%
% To use this processor, a mat file listing the microphone coordinates
% (spherical) for the particular model of microphone used must be in the
% Processors/Beamforming/Microphones directory (see the examples already
% present). At the time of writing, files for the 32-channel Eigenmike and
% the 64-channel Visisonics microphone are available.
%
% This function calls the HOAToolbox by Nicolas Epain.
%
% Code by Daniel Jimenez, Luis Miranda and Densil Cabrera
% Version 1.0 (19 August 2014)

if nargin < 7, micType = 'omni'; end
if nargin < 6, filterlength = 256; end
if nargin < 5, max_order = 4; end
if nargin < 4, cropIR = 0; end
if nargin < 3, mic_coords = importdata([cd '/Processors/Beamforming/Microphones/visisonics_mic_loc.mat']); end
if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
    if nargin < 4
        param = inputdlg({'Maximum order','Apply autocropstart_aarae.m','Filter length (in samples)','Mic Type'},...
                       'Input parameters',1,...
                      {num2str(max_order),num2str(cropIR),num2str(filterlength),micType});
        if isempty(param) || isempty(param{1,1}) || isempty(param{2,1}) || isempty(param{3,1}) || isempty(param{4,1})
            OUT = [];
            return;
        else
            max_order = str2double(param{1,1});
            cropIR = str2double(param{2,1});
            filterlength = str2double(param{3,1});
            switch param{4,1}
                    case {'omni','cardio','super', ...
                            'hyper','eight','measured'}
                        micType = param{4,1};
                    otherwise
                        error('Unknown microphone type') ;
            end       
            if isnan(max_order) || isnan(cropIR) || isnan(filterlength), OUT = []; return; end
        end
    end
    if nargin < 3
        mics = what([cd '/Processors/Beamforming/Microphones']);
        if isempty(mics.mat), warndlg('No microphone coordinates available. To add new mic coordinate files (*.mat) go to /Processors/Beamforming/Microphones','AARAE info','modal'); OUT = []; return; end
        [S,ok] = listdlg('Name','Microphones',...
                         'PromptString','Select microphone coordinates',...
                         'ListString',mics.mat,'SelectionMode','Single');
        if isempty(S), warndlg('No microphone coordinates selected. To add new mic coordinate files (*.mat) go to /Processors/Beamforming/Microphones','AARAE info','modal'); OUT = []; return; end
        mic_coords = importdata([cd '/Processors/Beamforming/Microphones/' mics.mat{S,1}]);
        if size(audio,2) ~= size(mic_coords,1), warndlg('The selected microphone coordinates do not match the number of audio channels.','AARAE info','modal'); OUT = []; return; end
    end
else
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2double(char(fs));
    end
    audio = IN;
end



% Pre-treatment IRs
if cropIR == 1
    trimmedIRs = autocropstart_aarae(audio,-20,2);
else
    trimmedIRs = audio;
end

% reset cal to 0 dB if it exists
if isstruct(IN)
    if isfield(IN,'cal')
        audio = cal_reset_aarae(audio,0,IN.cal);        
    end
end

%[~, firstarrivalall] = max(IRs);

%firstarrival = min(firstarrivalall);

%preroll = round((100./1000).*fs);

%if firstarrival-preroll < 0;
%    preroll = 0;
%end

%lengthIRs_after_arrival = round(analysis_length.*fs);

%if firstarrival+lengthIRs_after_arrival > size(IRs,1);
%    lengthIRs_after_arrival = size(IRs,1) - firstarrival;
%    disp('Analysis length too long - trimmed to end of file');
%end

%trimmedIRs = IRs((firstarrival-preroll):(firstarrival+lengthIRs_after_arrival),:);

%trimmedIRs = trimmedIRs./max(max(trimmedIRs)); % Normalize to direct sound

%clear IRs

% Mic setup

if size(mic_coords,2) ~= 3;
    warndlg('Mic coordinates are three columns with spherical coordinates','AARAE info','modal'); OUT = []; return;
end

micFmt = GenerateMicFmt({'sphCoord',mic_coords,'micType',micType});

% Check that max_order is not impossibly big
if (max_order+1)^2 > size(audio,2)
    % reduce to maximum possible order for the number of input channels
    max_order = floor(size(audio,2).^0.5)-1;
end

% HOA Signals

hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order) ;

mic2HoaOpt.sampFreq      = fs;
mic2HoaOpt.filterType    = 'firMatrix' ;
mic2HoaOpt.filterLength  = filterlength;
mic2HoaOpt.limiterMethod = 'tikh' ;
mic2HoaOpt.limiterLevel  = 6 ;
mic2HoaOpt.higherOrders  = false ;
mic2HoaOpt.subArrayFilt  = false ;
mic2HoaOpt.highFreqEq    = false ;
mic2HoaOpt.lowPassFreq   = 22000 ;

mic2HoaCfg = Mic2HoaEncodingFilters(hoaFmt,micFmt,mic2HoaOpt);

[len,~,bands,dim4,dim5,dim6] = size(trimmedIRs);
hoaSignals = zeros(len,hoaFmt.nbComp,bands,dim4,dim5,dim6);
for d6=1:dim6
    for d5 = 1:dim5
        for d4 = 1:dim4
            for b = 1:bands
                for I = 1:size(trimmedIRs,2);
                    for J = 1:hoaFmt.nbComp;
                        % revisar si se puede hacer matricial
                        hoaSignals(:,J,b,d4,d5,d6) =...
                            hoaSignals(:,J,b,d4,d5,d6) +...
                            fftfilt(mic2HoaCfg.filters.firMatrix(:,J,I),trimmedIRs(:,I,b,d4,d5,d6));
                    end
                end
            end
        end
    end
end


if isstruct(IN)
    OUT=IN;
    OUT.audio = hoaSignals;
    %OUT.fs = fs;
    OUT.cal = zeros(1,size(OUT.audio,2));
    %OUT.chanID = cellstr([repmat('Y ',[size(OUT.audio,2),1]),num2str(hoaFmt.index)]);
    OUT.chanID = makechanID(size(OUT.audio,2),1); %using an aarae utility function
    OUT.funcallback.name = 'audio2hoa.m';
    OUT.funcallback.inarg = {fs, mic_coords, cropIR, max_order,filterlength, micType};
else
    OUT = hoaSignals;
end