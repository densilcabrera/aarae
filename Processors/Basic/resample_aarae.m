function out = resample_aarae(in,newfs,n)
% This function resamples the audio, using Matlab's resample function.
% The .fs field is changed to the new audio sampling rate.
%
% Code by Densil Cabrera & Daniel Jimenez
% version 1.01 (24 March 2014)

fs = in.fs;
if nargin < 2
    prompt = {['New sampling rate (current is ',num2str(fs), ' Hz)'];...
        'Number of samples (n) used on each side of current sample'};
    dlg_title = 'Resample';
    num_lines = 1;
    def = {num2str(fs);'50'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    newfs = str2double(answer{1,1});
    n = str2double(answer{2,1});
end
if ~isempty(newfs)
    numtracks = 0;
    tracks = strfind(fieldnames(in),'audio');
    for t = 1:length(tracks)
        if tracks{t,1} == 1
            numtracks = numtracks + 1;
        end
    end
    for i = 1:size(in.audio,2)
        for j = 1:size(in.audio,3)
            for k = 1:size(in.audio,4)
                for l = 1:size(in.audio,5)
                    for m = 1:size(in.audio,6)
            out.audio(:,i,j,k,l,m) = resample(in.audio(:,i,j,k,l,m), newfs, fs, n);
                    end
                end
            end
        end
    end
%    out.audio = resample(in.audio, newfs, fs, n); % this reshapes to 2D
    if numtracks > 1
        for tracknum = 2:numtracks
            for i = 1:size(in.(genvarname(['audio' num2str(tracknum)])),2)
                for j = 1:size(in.(genvarname(['audio' num2str(tracknum)])),3)
                    out.(genvarname(['audio' num2str(tracknum)]))(:,i,j) = resample(in.(genvarname(['audio' num2str(tracknum)]))(:,i,j), newfs, fs, n);
                end
            end
        end
    end
    out.fs = newfs;
    out.funcallback.name = 'resample_aarae.m';
    out.funcallback.inarg = {newfs,n};
else
    out = [];
end

end