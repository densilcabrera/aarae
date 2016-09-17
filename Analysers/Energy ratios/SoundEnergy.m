function out = SoundEnergy(data, fs, startthresh, bpo, doplot)
% This function calculates Clarity Index, Definition and Centre Time from
% an impulse response.
%
% THIS FUNCTION REQUIRES MORE WORK - PLEASE USE ReverberationTime_IR1 to
% calculate Clarity Index, Definition and Centre Time instead of this
% function.

%
%--------------------------------------------------------------------------
% INPUT VARIABLES
%
% IR = .wav impulse response
%
% fs = Sampling frequency
%
% startthresh = Defines beginning of the IR as the first startthresh sample
%               as the maximum value in the IR (e.g if startthresh = -20,
%               the new IR start is the first sample that is >= 20 dB below
%               the maximum
%
% bpo = frequency scale to analyse
%             (1 = octave bands (default); 3 = 1/3 octave bands)
%
% doplot = Output plots (1 = yes; 0 = no)
%
%--------------------------------------------------------------------------
% OUTPUT VARIABLES
%
% C50 = Early (<50 ms) to late energy ratio, in decibels
% C80 = Early (<80 ms) to late energy ratio, in decibels
% D50 = Early (<50 ms) to total energy ratio, as a percentage
% D80 = Early (<80 ms) to total energy ratio, as a percentage
% Ts = Time of the centre of gravity of the squared IR, in seconds
%
%--------------------------------------------------------------------------
%
% ************************  TO DO  ****************************************
% OBVIOUS IMPROVEMENTS
% Make compatible with multichannel input
% consider whether to allow multiband input (maybe with disclaimer)
% replace octave band and 1/3-oct band filters with AARAE's filters
% display table of results
% consider whether to display chart(s)
% include authors (Grant & Densil), BSD license text, & version number
% reformat output structure
%
% *** What could this function do better than ReverberatioTime_IR1? *******
% POSSIBLE ANSWERS: 
% * Detailed plots showing the reasons for the resulting values (e.g using
%       cumulative sums)
% * Interpret values in terms of published theories/models/criteria
% * Introduce non-standard options such as a cross-fade between early and
%       late (could be controlled by a single input parameter)
% * Include some other values such as C10, C35 and C100 (not too many, as
%       there is another AARAE function that does as many as you like)
% * Include end truncation (both automatic and user-controlled?)
% * Determine the end point where C50 etc is within tolerance limit (e.g.
%       0.1 dB, 0.5 dB, 1 dB etc, or percentage error for D50 etc)
% * Allow special filterbanks - e.g. Gammatone (it is better to filter after
%       truncation, and so to call the filterbank from this function)

if nargin < 5, doplot = 1; end
if nargin < 4, bpo = 1; end
if nargin < 3
    startthresh = -20;
    %dialog box for settings
    prompt = {'Threshold for IR start detection', ...
        'Bands per octave (1 | 3)', ...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'-20','1','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        startthresh = str2num(answer{1,1});
        bpo = str2num(answer{2,1});
        doplot = str2num(answer{3,1});
    else
        out = [];
        return
    end
end
if isstruct(data)
    data = choose_from_higher_dimensions(data,3,1); 
    IR = data.audio;
    fs = data.fs;
else
    IR = data;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(IR) && ~isempty(fs) && ~isempty(startthresh) && ~isempty(bpo) && ~isempty(doplot)
    %--------------------------------------------------------------------------
    % TRUNCATION
    %--------------------------------------------------------------------------

    % Check the input data dimensions
    S = size(IR); % size of the IR
    ndim = length(S); % number of dimensions
    switch ndim
        case 1
            len = S(1); % number of samples in IR
            chans = 1; % number of channels
        case 2
            len = S(1); % number of samples in IR
            chans = S(2); % number of channels
        case 3
            warndlg('Sound Energy input IR cannot be multiband','AARAE info','modal');
            %out.error = 'Input IR cannot be multiband';
            out = [];
            return
    end

    % Preallocate
    m = zeros(1, chans); % maximum value of the IR
    startpoint = zeros(1, chans); % the auto-detected start time of the IR

    for dim2 = 1:chans
        m(1,dim2) = max(IR(:,dim2).^2); % maximum value of the IR
        startpoint(1,dim2) = find(IR(:,dim2).^2 >= m(1,dim2)./ ...
            (10^(abs(startthresh)/10)),1,'first'); % Define start point

        if startpoint(1,dim2) > 1
            % zero the data before the startpoint
            IR(1:startpoint(1,dim2)-1,dim2) = 0;

            % rotate the zeros to the end (to keep a constant data length)
            IR(:,dim2) = circshift(IR(:,dim2),-(startpoint(1,dim2)-1));
        end
    end

    Early50 = IR(1:1+floor(fs*0.05),:); % Truncate Early80
    Early80 = IR(1:1+floor(fs*0.08),:); % Truncate Early80
    Late50 = IR(ceil(fs*0.05):end,:); % Truncate Late50
    Late80 = IR(ceil(fs*0.08):end,:); % Truncate Late80


    %--------------------------------------------------------------------------
    % FILTERING
    %--------------------------------------------------------------------------

    if bpo == 3
        bandnumber = 20:37; % filter band numbers (1/3 octaves 100 Hz - 5 kHz)
        bandwidth = 1/3;
        halforder = 2; % half of the filter order
    else
        bandnumber = 21:3:36; % filter band numbers (octave bands 125 Hz - 4 kHz)
        bandwidth = 1;
        halforder = 3; % half of the filter order
    end

    fc = 10.^(bandnumber./10); % filter centre frequencies in Hz
    bands = length(fc);
    f_low = fc./10^(0.15*bandwidth); % low cut-off frequency in Hz
    f_hi = fc.*10^(0.15*bandwidth); % high cut-off frequency in Hz
    Nyquist = fs/2; % Nyquist frequency
    b = zeros(halforder*2+1,length(fc)); % pre-allocate filter coefficients
    a = b; % pre-allocate filter coefficients

    % calculate filter coefficients
    for k = 1:bands
        [b(:,k), a(:,k)]=butter(halforder, [f_low(k)/Nyquist f_hi(k)/Nyquist]);
    end

    % Preallocate
    IRoct = zeros(len,chans,bands);
    Early50oct = zeros(length(Early50),chans,bands);
    Early80oct = zeros(length(Early80),chans,bands);
    Late50oct = zeros(length(Late50),chans,bands);
    Late80oct = zeros(length(Late80),chans,bands);

    % filter IR and Early/Late
    for k = 1:bands
        IRoct(:,:,k) = filter(b(:,k),a(:,k), IR); % IR
        Early50oct(:,:,k) = filter(b(:,k),a(:,k), Early50); % Early50
        Early80oct(:,:,k) = filter(b(:,k),a(:,k), Early80); % Early80
        Late50oct(:,:,k) = filter(b(:,k),a(:,k), Late50); % Late50
        Late80oct(:,:,k) = filter(b(:,k),a(:,k), Late80); % Late80
    end

    %--------------------------------------------------------------------------
    % CALCULATE ENERGY PARAMETERS
    %--------------------------------------------------------------------------

    Early50oct = squeeze(sum(Early50oct.^2));
    Early80oct = squeeze(sum(Early80oct.^2));
    Late50oct = squeeze(sum(Late50oct.^2));
    Late80oct = squeeze(sum(Late80oct.^2));
    all_oct = squeeze(sum(IRoct.^2));

    C50 = 10*log10(Early50oct ./ Late50oct)'; % C50
    C80 = 10*log10(Early80oct ./ Late80oct)'; % C80
    D50 = (Early50oct ./ all_oct)'; % D50
    D80 = (Early80oct ./ all_oct)'; % D80

    TsTimes = (length(IRoct)-1)' ./ fs; % length of IR bands in seconds
    Ts = (squeeze(sum(IRoct.^2 .* repmat(TsTimes,size(IRoct)))) ./ ...
        all_oct)'; % Ts

    if bpo == 3
        bandfc = [100,125,160,200,250,315,400,500,630,800,1000,1250,1600, ...
            2000,2500,3150,4000,5000];
    else
        bandfc = [125,250,500,1000,2000,4000];
    end % if bpo

    %--------------------------------------------------------------------------
    % TABLE
    %--------------------------------------------------------------------------
    
    if doplot == 1
        if chans == 1
            fig1 = figure('Name','Sound Energy');
            table1 = uitable('Data',[C50;C80;D50;D80;Ts],...
                             'ColumnName',num2cell(bandfc),...
                             'RowName',{'C50','C80','D50','D80','Ts'});
            [~,tables] = disptables(fig1,table1,{'Chan1 - Sound energy'});
            out.tables = tables;
        elseif chans == 2
            fig1 = figure('Name','Sound Energy');
            table1 = uitable('Data',[C50(:,1)';C80(:,1)';D50(:,1)';D80(:,1)';Ts(:,1)'],...
                             'ColumnName',num2cell(bandfc),...
                             'RowName',{'C50','C80','D50','D80','Ts'});
            table2 = uitable('Data',[C50(:,2)';C80(:,2)';D50(:,2)';D80(:,2)';Ts(:,2)'],...
                             'ColumnName',num2cell(bandfc),...
                             'RowName',{'C50','C80','D50','D80','Ts'});
            [~,tables] = disptables(fig1,[table1 table2],{'Chan1 - Sound energy','Chan2 - Sound energy'});
            out.tables = tables;
        end
    end
                     
    %--------------------------------------------------------------------------
    % OUTPUT STRUCTURE
    %--------------------------------------------------------------------------

    if chans == 1
        out.C50 = [bandfc;C50];
        out.C80 = [bandfc;C80];
        out.D50 = [bandfc;D50];
        out.D80 = [bandfc;D80];
        out.Ts = [bandfc;Ts];
        out.funcallback.name = 'SoundEnergy.m';
        out.funcallback.inarg = {fs,startthresh,bpo,doplot};

    elseif chans == 2
        out.C50_ch1 = [bandfc;C50(:,1)'];
        out.C50_ch2 = [bandfc;C50(:,2)'];

        out.C80_ch1 = [bandfc;C80(:,1)'];
        out.C80_ch2 = [bandfc;C80(:,2)'];

        out.D50_ch1 = [bandfc;D50(:,1)'];
        out.D50_ch2 = [bandfc;D50(:,2)'];

        out.D80_ch1 = [bandfc;D80(:,1)'];
        out.D80_ch2 = [bandfc;D80(:,2)'];

        out.Ts_ch1 = [bandfc;Ts(:,1)'];
        out.Ts_ch2 = [bandfc;Ts(:,2)'];

        out.funcallback.name = 'SoundEnergy.m';
        out.funcallback.inarg = {fs,startthresh,bpo,doplot};
    else
        out = [];
        warndlg('Function calcuates for mono or stereo audio only','AARAE info');
    end % if chans
else
    out = [];
end
% eof

