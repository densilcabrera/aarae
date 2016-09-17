function OUT = WaveComparisonPlotMatrix(in,fs,lagadjust)
% This function creates a plot matrix to compare multiple channels of
% audio. (It has no use for single channel audio.)
%
% In the case of multi-band multi-channel audio, a separate figure is
% created for each band.
%
% Optionally, the cross-correlation peak's lag can be used to re-align the
% channels prior to generating the plot matrix. You can select the
% reference channel.
%
% Code by Densil Cabrera
% Version 1.01 (5 December 2013)

if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1);
    audio = in.audio;
    fs = in.fs;
    if isfield(in,'bandID')
        bandID = in.bandID;
    end
else
    audio = in;
end
if nargin < 3
    prompt = {'Inter-channel lag reference channel number (or 0 for none, or a large number for all)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        lagadjust = str2double(answer{1,1});
    end
end


[~, chans, bands] = size(audio);

% temporary block for higher values - to be removed when the average lag
% code is working
if lagadjust > chans, lagadjust = 1; end

maxlags = ceil(fs * 0.05); % maximum xcorr lags

if ~exist('bandID','var')
    bandID = 1:bands;
end

for b = 1:bands
    figure('Name',['Wave Comparison Plot-matrix, band ', num2str(bandID(b))])
    if lagadjust > 0
        c = xcorr(audio(:,:,b), maxlags, 'coeff');
        if lagadjust <= chans
            maxc = max(abs(c(:,(lagadjust-1)*chans+1:(lagadjust-1)*chans+chans)));
        else
            maxc = max(abs(c));
        end
        chanlag = zeros(1,chans);
        for ch = 1:chans
            if lagadjust <= chans
            chanlag(ch) = find(abs(c(:,ch+(lagadjust-1)*chans)) == maxc(ch),1,'first');
            else
                % THIS PART OF THE CODE IS NOT WORKING YET %%%%%%%%%%%%%
                for ch2 = 1:chans
                    if ch2 ~= ch
                    %chanlag(ch) = chanlag(ch) + ...
                        x = find(abs(c(:,ch+(chans-1)*ch2)) == maxc(ch+(chans-1)*ch2),1,'first');
                        if ch2 < ch
                        chanlag(ch) = chanlag(ch) - x;
                        else
                            chanlag(ch) = chanlag(ch) + x;
                        end
                        disp([num2str(ch), num2str(ch2),num2str(x)])
                    end
                end
                chanlag(ch) = round(chanlag(ch) / (chans-1)); % mean interchannel lag
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            audio(:,ch,b) = circshift(audio(:,ch,b),chanlag(ch));
        end
        disp(num2str(chanlag))
    end
    
    plotmatrix(audio(:,:,b));
end
OUT.funcallback.name = 'WaveComparisonPlotMatrix.m';
OUT.funcallback.inarg = {fs,lagadjust};