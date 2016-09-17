function [OUT, varargout] = Blind_RT_LoellmannJeub2012(IN, fs)
% This function calls Loellmann and Jeub's blind reverberation time
% estimation functions. Use this function to estimate reverberation time
% from reverberant speech (e.g. 60 seconds of speech recording).
%
% Please refer to the following directory for the license and for more
% information on the algorithm:
% Analysers\Reverberation\Release_RT_estimation_MatlabFileExchange
%
% Reference:
% Heinrich W. Löllmann, Emre Yilmaz, Marco Jeub and Peter Vary:
% "An Improved Algorithm for Blind Reverberation Time Estimation"
% International Workshop on Acoustic Echo and Noise Control (IWAENC),
% Tel Aviv, Israel, August 2010.
% (availabel at www.ind.rwth-aachen.de/~bib/loellmann10a)
%
% The algorithm allows to estimate the RT within a range of 0.2s to 1.2s
% and assumes that source and receiver are not within the critical
% distance. A denoising is not performed by this function and has to be
% done in advance.

% Calling function for integration into AARAE by Densil Cabrera
% Version 1.00 (December 2013)



if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    if isfield(IN,'bandID')
        bandID = IN.bandID;
    end
    if isfield(IN,'chanID')
        chanID = IN.chanID;
    end
%     if isfield(IN,'cal')
%         cal = IN.cal;
%     else
%         cal = zeros(1,size(audio,2));
%         disp('This audio signal has not been calibrated.')
%     end
    
    
elseif ~isempty(param) || nargin > 1
    
    audio = IN;
    
end


if ~isempty(audio) && ~isempty(fs)
    
    if size(audio,3)>1
        audio = sum(audio,3);
        disp('Multiband audio has been summed for blind RT estimation')
    end
    
    if fs<8e3 || fs>24e3
        audio = resample(audio,24000,fs);
        fs = 24000;
        disp('Audio has been resampled to fs=24000 Hz for blind RT estimation')
    end
    
    [len,chans] = size(audio);
    
    
    
%     if exist('cal','var')
%         audio = audio .* repmat(10.^(cal./20),[len,1,bands]);
%     end
    
    
    dur = len /fs;
    if dur < 15
        warndlg('This analyser is designed for reverberant speech, for example, 60 duration. (It is not designed to analyse impulse responses.) The input audio must be at least 15 s long for the analyser to run.','AARAE info')
        OUT = [];
        return
    end
    
    
    simpar.fs = fs;
    simpar.block_size = round(20e-3 * simpar.fs);  % block length
    simpar.overlap = round(simpar.block_size/2);   % overlap
    
    for ch = 1:chans
        [rt_est(ch,:),rt_est_mean(ch),rt_est_dbg] = ML_RT_estimation(audio(:,ch,1)',simpar);
        
        rt_est_median(ch) = median(rt_est(ch,:));
        
        %--------------------------------------------------------------------------
        % Plot estimated RT and 'true' RT obtained by Schroeder method
        %--------------------------------------------------------------------------
        if exist('chanID','var')
            chanstring = char(chanID(ch));
        else
            chanstring = ['chan ',num2str(ch)];
        end
        
        fr2sec_idx = linspace(1,len/simpar.fs,size(rt_est,2));
        figure('Name',['Blind Reverberation Time Estimate ',chanstring])
        clf
        hold on
        plot(fr2sec_idx,rt_est(ch,:),'-r')
        plot([0 fr2sec_idx(end)],[rt_est_mean(ch) rt_est_mean(ch)])
        plot([0 fr2sec_idx(end)],[rt_est_median(ch) rt_est_median(ch)],'Color', [0,0.5,0])
        grid on,box on
        xlabel('Time [s]'),ylabel('RT [s]');
        legend('Estimated T60',['Mean Estimate ',num2str(rt_est_mean(ch)), ' s'], ...
            ['Median Estimate ',num2str(rt_est_median(ch)), ' s'],'location','southeast');
        %--------------------------------------------------------------------------
        
    end
    
    if isstruct(IN)
        doresultleaf(rt_est','RT [s]',{'Time'},...
                     'Time',     fr2sec_idx, 's',           true,...
                     'channels', chanID,     'categorical', [],...
                     'name','Blind_reverb_time');
    end
    
    % output table
    f = figure;
    
    if exist('chanID','var')
        rownamecell = chanID;
    else
        rownamecell = num2cell(1:chans);
    end
    
    t = uitable('Data',[rt_est_mean' rt_est_median'],...
        'ColumnName',{'Mean RT estimate (s)' 'Median RT estimate (s)'},...
        'RowName',rownamecell);
    disptables(f,t);
    
    
    if isstruct(IN)
        OUT.rt_est = rt_est;
        OUT.rt_est_mean = rt_est_mean;
        OUT.rt_est_dbg = rt_est_dbg;
        OUT.funcallback.name = 'Blind_RT_LoellmannJeub2012.m';
        OUT.funcallback.inarg = {}; % nothing needed for callback
    else
        
        OUT = rt_est;
    end
    varargout{1} = rt_est_mean;
    varargout{2} = rt_est_dbg;
    
    
else
    
    OUT = [];
end