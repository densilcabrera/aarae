function [OUT,varargout] = thirdoctbandfilter(IN,fs,param,method)
% Implements an 8th order 1/3-octave band filterbank (by default) using the
% Audio Toolbox's octaveFilterBank. Note that this is different to the
% original implementation of this function, which previously used
% fdesign.octave (now deprecated).
%
% This revision has removed the option of backwards-forwards double
% filtering (previously indicated by method = 0). Use AARAE's FFT-based
% filterbank instead for linear-phase filtering.
%
% The 'method' input is not available via the GUI. It has the following
% meaning:
%   A value >= 0 leads to forward time filtering
%   A value < 0 leads to reverse time filtering
%   An even value ~= 0 specifies the filter order
%   Other values lead to the default filter order of 8
% Default is 8th order forward time filtering (same as previous
% implementation).
%
% This update includes the possibility of lower frequency filters (down to
% 12.5 Hz). The highest available filter depends on the Nyquist frequency
% being above the high cutoff frequency of the filter. Hence the 20 kHz
% 1/3-octave band filter is not available for 44.1 kHz sampling rate.
%
% The filtered octave bands are always contiguous in this version of the
% function (unlike the previous version). They are from the lowest to
% highest frequency bands selected using the list dialog or the 'param'
% input of the function. Hence the 'param' input can consist of just the
% lowest and highest band centre frequencies e.g. [125 8000].
%
% Update March 2024.
ok = 0;
if isstruct(IN)
    fs = IN.fs;
end
if ~exist('fs','var'), OUT = []; return; end
if nargin < 4, method = 8; end % 8th order forward time filtering (same as method = 1)

% find maximum viable octave band for the sampling rate
maxbandindex = floor(pow2db(0.5*fs/10^0.05)-10);
nominalfreq = [12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000];
if maxbandindex<33
    nominalfreq = nominalfreq(1:maxbandindex);
end
if nargin < 3
    param = nominalfreq;
    [S,ok] = listdlg('Name','Octave band filter input parameters',...
        'PromptString','Center frequencies [Hz]',...
        'ListString',[num2str(param') repmat(' Hz',length(param),1)]);
    param = param(S);
else
    for i = 1:length(param)
        check = find(nominalfreq == param(i));
        if isempty(check), check = 0; end
        S(i) = check;
    end
    if all(S), param = sort(param,'ascend'); ok = 1; else ok = 0; end
end
if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
elseif ~isempty(param)
    audio = IN;
    % if nargin < 2
    %     fs = inputdlg({'Sampling frequency [samples/s]'},...
    %         'Fs',1,{'48000'});
    %     fs = str2num(char(fs));
    % end
end
if size(audio,3)>1
    audio = sum(audio,3);
    disp('Multiband audio has been mixed down prior to filtering')
end
if ~isempty(param) && ~isempty(fs)
    if rem(method,2) == 0 && method ~= 0
        order = abs(method);
    else
        order = 8; % default to 8th order filters
    end

    if method < 0
        audio = flip(audio,1);
    end


    bandwidth = '1/3 octave';
    range = [min(param)./10.^0.02 max(param)*10^0.02];
    ref = 1000;
    base = 10;
    octFiltBank = octaveFilterBank(bandwidth,fs,...
        'FrequencyRange',range,...
        'ReferenceFrequency',ref,...
        'FilterOrder',order,...
        'OctaveRatioBase',base);
    fc = info(octFiltBank).CenterFrequencies;
    filtered = repmat(audio,[1,1,length(fc),1,1,1]); % preallocate

    for d4 = 1:size(audio,4)
        for d5 = 1:size(audio,5)
            for d6 = 1:size(audio,6)
                filtered(:,:,:,d4,d5,d6) = permute(octFiltBank(audio(:,:,1,d4,d5,d6)),[1,3,2,4,5,6]);
            end
        end
    end


    if method < 0
        filtered = flip(filtered,1);
    end
    centerf = param;
else
    filtered = [];
    centerf = [];
end
if isstruct(IN) && ~isempty(filtered)
    OUT = IN;
    OUT.audio = filtered;
    OUT.bandID = centerf;
    OUT.funcallback.name = 'thirdoctbandfilter.m';
    OUT.funcallback.inarg = {fs,param,method};
else
    OUT = filtered;
end
varargout{1} = centerf;
end