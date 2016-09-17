function out = DiffusenessHanyuISRA2013(IR, fs, startindex, endindex, doplot, threshold)
% This function calculates the degree of time series fluctuation of a room
% impulse response, or part thereof, based on the method proposed in:
% T. Hanyu, 'Analysis method for estimating diffuseness of sound fields by
% using decay-cancelled impulse responses,' International Symposium on Room
% Acoustics, Toronto, Canada, 2013.
%
% Code by Densil Cabrera
% version 1.02 (4 April 2014)
%
% INPUT ARGUMENTS
%
% IR is a room impulse response, or a set of impulse responses,
% which should (each) have a reasonably exponential decay. If IR is a
% structure, then it must include IR.audio and IR.fs (and then the other
% input arguments are ignored, with values set via dialog box).
% It is likely that an IR would be filtered (e.g. into octave bands)
% prior to calling this function.
% Truncated sections of impulse responses could be input, so long as they
% are long enough to estimate reverberation time from the section (i.e.,
% the local fluctuations in the waveform should not affect the overall
% slope, and there is adequate time to build-up the reverse-integrated IR^2).
% The function is designed to accept up to 3-dimensional matrices, where
% dimension 1 is sampled time, dimension 2 is (by convention) channel, and
% dimension 3 is band (e.g., a set of octave bands).
% The IR is initially assumed to be noise-free, and to start at sample 1.
%
% fs is the audio sampling rate in Hz.
%
% doplot: 1 for plotting, 0 for no plotting (plotting is done by default).
%
% threshold is an optional argument which can be used to change the
% z(k) threshold from the default of 0.01 (which is the threshold proposed
% in Hanyu's paper). If omitted then 0.01 is used. Threshold must be
% between 1 and 0.
%
%
% OUTPUT ARGUMENTS
%
% This function outputs a structure, which includes the following leaves:
%
% out.T: estimated reverberation time in seconds. Note that this is not
%   calculated using the standard methods, and so it may deviate from
%   reverberation time calculated in the more usual ways.
% out.result: degree of time series fluctuation.
%   A high value indicates low diffusivity, and vice versa.
% out.h2: Energy fluctuation time-series
% out.h2_sort: Energy fluctuation values sorted in descending order
%   (which is used to find the threshold point in the z(k) function).
% out.z: The z(k) function (which can be plotted against out.h2_sort).
%


%**************************************************************
% Check the input data

if isstruct(IR)
    IR = choose_from_higher_dimensions(IR,3,1); 
    % if the first input is a structure, then dialog boxes may be used for
    % user controls
    % required fields from input structure
    data = IR.audio;
    fs = IR.fs;
    starttime = 0;
    endtime = (length(data)-1)/fs;
    if isfield(IR,'bandID')
        bandID = IR.bandID;
    end
    if isfield(IR,'chanID')
        chanID = IR.chanID;
    end
else
    data = IR;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
    starttime = 0;
    endtime = (length(data)-1)/fs;
end;
if nargin < 4, threshold = 0.01; end % default z(k) threshold is 0.01
if nargin < 3
    doplot = 1; % plot by default
    settimelimits = 0;
    % dialog box to get start and finish times
    while settimelimits == 0
        prompt = {'Start time in s:', 'End time in s', 'z(k) threshold','Plot'};
        dlg_title = 'Settings';
        num_lines = 1;
        def = {num2str(starttime),num2str(endtime), num2str(0.01), num2str(1)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);

        if ~isempty(answer)
            starttime = str2num(answer{1,1});
            endtime = str2num(answer{2,1});
            threshold = str2num(answer{3,1});
            doplot = str2num(answer{4,1});
        else
            out = [];
            return
        end

        startindex = round(starttime*fs)+1;
        endindex = floor(endtime*fs);

        if startindex > 0 ...
                && startindex < endindex ...
                && endindex <= length(data)
            settimelimits = 1;
        end
    end %while
end

if ~isempty(data) && ~isempty(fs) && ~isempty(startindex) && ~isempty(endindex) && ~isempty(doplot) && ~isempty(threshold)
    data = data(startindex:endindex,:,:);
    S = size(data); % size of the IR
    ndim = length(S); % number of dimensions
    switch ndim
        case 1
            len = S(1); % number of samples in IR
            chans = 1; % number of channels
            bands = 1; % number of bands
        case 2
            len = S(1); % number of samples in IR
            chans = S(2); % number of channels
            bands = 1; % number of bands
        case 3
            len = S(1); % number of samples in IR
            chans = S(2); % number of channels
            bands = S(3); % number of bands
    end

    % We will not use the final part of the data because the initial build-up
    % of the reverse-integrated curve (in reverse time) results in a systematic
    % error, by raising the tail end of the h^2 curve. Testing with
    % exponentally decaying white noise  indicates that this error is likely to
    % be sufficiently small if the final 20% of the data is truncated after
    % reverse integration.
    truncation = 0.2; % the proportion of data to truncate
    len_trunc = round((1-truncation)*len); % truncated length

    % The following renaming is done to match Hanyu's paper
    p = data;
    clear data

    %**************************************************************
    % Apply the equations in Hanyu's paper (as identified below)

    % Reverse integration of the squared IR
    Es = flipdim(cumsum(flipdim(p.^2,1)),1); % equation 1

    % decay cancelled impulse response
    g = p ./ (Es).^0.5; % equation 2

    % squared decay-cancelled impulse response
    g2 = p.^2 ./ Es; % equation 3

    % In calculating the mean of g2, we will truncate it to avoid the high
    % values that come from the initial reverse-time growth of the reverse
    % integrated IR.
    A = mean(g2(1:len_trunc,:,:)); % equation 7

    % equation 11 (reorganised)
    h = g(1:len_trunc,:,:) ./ repmat(A,[len_trunc,1,1]).^0.5;

    h2 = h.^2;

    Rtotal = sum(h2); % equation 13

    h2_sort = sort(h2,1,'descend');

    Rk = cumsum(h2_sort);

    z = Rk ./ repmat(Rtotal, [len_trunc,1,1]); % equation 14

    % preallocate
    zval = zeros(1,chans,bands);
    result = zeros(1,chans,bands);

    % find threshold k
    for dim2=1:chans
        for dim3 = 1:bands
            % find index for threshold k
            ind = find(z(:,dim2,dim3)>= threshold, 1,'first');
            if isempty(ind)
                out.error = 'Threshold index not found';
                disp('Threshold index not found!')
                disp('Try a different start time or end time (e.g., to avoid -inf dB).')
                return
            end

            % exact z(k) value(s) that were found 
            % (it should be approximately equal to threshold)
            zval(1, dim2, dim3) = z(ind,dim2,dim3);

            % k threshold value, which is the degree of time series fluctuation
            result(1, dim2, dim3) = h2_sort(ind,dim2,dim3);
        end
    end

    %**************************************************************
    % CREATE OUTPUT STRUCTURE

    % the degree of time series fluctuation
    out.result = result;

    % z(k)
    out.z = z;

    % h^2
    out.h2 = h2;

    % h^2 sorted in descending order
    % (Although this output is redundant, it is provided to facilitate plotting)
    out.h2_sort = h2_sort;

    % Estimate of reverberation time in s
    out.T = 13.82 ./A  ./fs; % equation 9

    % Batch analysing parameters
    out.funcallback.name = 'DiffusenessHanyuISRA2013.m';
    out.funcallback.inarg = {fs,startindex,endindex,doplot,threshold};

    if ~exist('chanID','var')
        chanID = {1:chans};
    end
    if ~exist('bandID','var')
        bandID = {1:bands};
    else
        bandID = num2cell(bandID);
    end

    %**************************************************************
    % CHARTS

    if doplot

        % Generate the RGB values of colours to be used in the plots
        M = plotcolours(chans, bands);


        figure('Name','Hanyu Diffusivity Analysis 1');
        % This figure shows the h-squared values as a function of time. It has
        % two subfigures: the upper figure shows the actual values, and the
        % lower figure shows the values smoothed using a 50 ms running average.
        % Note that the average value of h^2 should be equal to 1, and so
        % examining these charts should provide an indication of
        % problems with the IR analysis (if any).

        subplot(2,1,1)
        hold on
        t = (0:(length(h2)-1))./fs;
        for dim2 = 1:chans
            for dim3 = 1:bands

                % plot h^2 time series (like Hanyu's Figure 2 right-side)
                plot(t,h2(:,dim2,dim3),'Color',squeeze(M(dim3,dim2,:))')
            end
        end    
        xlabel('Time (s)')
        ylabel('h^2')
        
        hold off
        if isstruct(IR)
            if ~ismatrix(h2)
                doresultleaf(h2,'h^2',{'time'},...
                             'time',     t',     's',           true,...
                             'channels', chanID, 'categorical', [],...
                             'bands',    bandID, 'Hz',          false,...
                             'name','h^2');
            else
                doresultleaf(h2,'h^2',{'time'},...
                             'time',     t',     's',           true,...
                             'channels', chanID, 'categorical', [],...
                             'name','h^2');
            end
        end

        %***********************
        subplot(2,1,2)
        % Running average of h^2
        windowlength = 0.05;% 50 ms running average window
        n_points = round(fs * windowlength); % number of samples in window
        b = ones(1,n_points)/n_points; % filter coefficients

        % compensation function for the fade-in of the smoothing filter
        fadeincomp = ones(len_trunc,1);
        % un-comment the following line if you wish to use this
        % fadeincomp(1:n_points-1) = 1 ./ fftfilt(b,ones(n_points-1,1));

        hold on
        for dim2 = 1:chans
            for dim3 = 1:bands
                % plot running average of h^2 (using fftfilt, with the above
                % coefficients)
                rah2(:,dim2,dim3) = fadeincomp .* fftfilt(b,h2(:,dim2,dim3));
                plot(t,rah2(:,dim2,dim3), ...
                    'Color',squeeze(M(dim3,dim2,:))')
            end
        end

        % plot a dashed grey line at h^2 = 1 (this is the overall average
        % value of h^2, by definition)
        plot([0 t(end)],[1 1], ...
            'LineStyle','--','Color',[0.5 0.5 0.5])

        % plot a vertical line to mark the point at which the running average
        % is completely filled with data (rather than zeros)
        plot([windowlength windowlength], [0 1], ...
            'LineStyle','--','Color',[0.5 0.5 0.5]);

        xlabel('Time (s)')
        ylabel('Running average of h(t)^2')

        hold off
        if isstruct(IR)
            if ~ismatrix(h2)
                doresultleaf(rah2,'Running average of h(t)^2',{'time'},...
                             'time',     t',     's',           true,...
                             'channels', chanID, 'categorical', [],...
                             'bands',    bandID, 'Hz',          false,...
                             'name','run_av_h2');
            else
                doresultleaf(rah2,'Running average of h(t)^2',{'time'},...
                             'time',     t',     's',           true,...
                             'channels', chanID, 'categorical', [],...
                             'name','run_av_h2');
            end
        end

        %***********************
        % Create a figure to follow Hanyu's Figure 4
        figure2 = figure('Name','Hanyu Diffusivity Analysis 2');

        % Create axes
        axes('Parent',figure2,'YScale','log','YMinorTick','on');
        xlabel('k threshold')
        ylabel('z(k)')

        xmax = ceil(max(max(max(h2_sort)))./20) *20;
        xlim([0 xmax]);

        ymin = 1e-3;
        %ymin = 10.^floor(log10(min(min(min(z))))); % alternative ymin
        ylim ([ymin 1]);

        hold on
        for dim2 = 1:chans
            for dim3 = 1:bands

                % plot fluctuation decay curve
                semilogy(h2_sort(:,dim2,dim3), z(:,dim2,dim3), ...
                    'Color',squeeze(M(dim3,dim2,:))')

                % plot threshold intersection point
                semilogy(result(1,dim2,dim3), zval(1,dim2,dim3), ...
                    'Marker','o','LineStyle','none',...
                    'Color',squeeze(M(dim3,dim2,:))')

                % draw line to x-axis
                semilogy([result(1,dim2,dim3) result(1,dim2,dim3)], ...
                    [ymin threshold], ...
                    'LineStyle','--','Color',squeeze(M(dim3,dim2,:))')
            end
        end

        % plot threshold line (dashed grey)
        semilogy(xlim, [threshold,threshold], ...
            'LineStyle','--','Color',[0.5 0.5 0.5])

        hold off
        if isstruct(IR)
            if ~ismatrix(h2)
                doresultleaf(z,'z(k)',{'k_threshold'},...
                             'k_threshold', h2_sort, [],            true,...
                             'channels', chanID,  'categorical', [],...
                             'bands',    bandID,  'Hz',          [],...
                             'name','z');
            else
                doresultleaf(z,'z(k)',{'k_threshold'},...
                             'k_threshold',        h2_sort, [],            true,...
                             'channels', chanID,  'categorical', [],...
                             'name','z');
            end
        end
        
        % results table
        fig3 = figure('Name','Degree of time series fluctuation (Hanyu ISRA 2013)');
        
        RowName = chanID;
        ColumnName = bandID;
        
        table1 = uitable('Data',permute(result,[2,3,1]),...
                'ColumnName',ColumnName,...
                'RowName',RowName);
        set(table1,'ColumnWidth',{52});
        disptables(fig3,table1);
    end
else
    out = [];
end % eof
end





%**************************************************************

function M = plotcolours(chans, bands)
% define plot colours in HSV colour-space
% use Value for channel
% use Hue for band

% Saturation = 1
S = ones(bands,chans);

% Value
maxV = 1;
minV = 0.4;
if chans == 1
    V = minV;
else
    V = minV:(maxV-minV)/(chans-1):maxV;
end
V = repmat(V,[bands,1]);

% Hue
if bands == 1
    H = 0;
else
    H = 0:1/(bands):(1-1/bands);
end
H = repmat(H',[1,chans]);

% convert HSV to RGB
M = hsv2rgb(cat(3,H,S,V));

end % eof