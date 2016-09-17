function out = spectrogram_simple(in, fs, transposesubplots, winlen, NOVERLAP, dBrange, maxfreq, downsamp, wfchoice, waterfall,doresultsleaf)
% Creates a spectrogram of a waveform, or a matrix of waveforms, as a
% function of time. Matlab's spectrogram function is used.
%
% Usually a spectrogram would be used to visualise audio in 1 or more
% channels, but probably not with multiple bands. However multiband display
% is available in case it is of some use.
% By default:
%  If the waveform is multichannel, then channels are in rows.
%  If the waveform is multiband, then bands are in columns.
%  If 'Transpose subplots' is set in the dialog box, then these dimensions
%  are transposed (channels in columns, bands in rows).
%
% The window length and overlap can be specified (in samples).
%
% Values are in decibels, and the decibel range of the colormap, can be 
% specified. Note that the colormap of each subplot is individually scaled.
%
% The upper frequency to display can be chosen (note data up to the Nyquist
% frequency is calculated regardless of this). If the input upper frequency
% is higher than the Nyquist frequency (e.g. if the audio is downsampled -
% see below) then the upper frequency is changed to the Nyquist frequency.
%
% The data can be downsampled prior to analysis - for example a
% downsampling factor of 2 halves the sampling rate (and removes the top
% octave).
%
% Three window functions are available:
%   Rectangular (not usually recommended, but included to demonstrate why
%     not);
%   Hann (sometimes called Hanning);
%   Blackman-Harris.
%
% Either a conventional spectrogram display or a waterfall display can be
% generated.

if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1);
    audio = in.audio;
    fs = in.fs;
    if isfield(in, 'cal')
       cal = in.cal; % get the calibration offset from the input structure
    end
else
    audio = in;
end
Nyquist = fs/2;
if nargin < 3
    %dialog box for settings
    prompt = {'Transpose subplots (0 | 1)', ...
        'Window length (samples)', ...
        'Number of overlapping samples', ...
        'Decibel range', ...
        'Highest frequency to display (Hz)', ...
        'Downsample factor', ...
        'Rectangular, Hann, Blackman-Harris window (r|h|b)',...
        'Waterfall (0 | 1)',...
        'Output results leaf (0 | 1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0','2048','1024','90',num2str(Nyquist), '1','h','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        transposesubplots = str2num(answer{1,1});
        winlen = str2num(answer{2,1});
        NOVERLAP = str2num(answer{3,1});
        dBrange = str2num(answer{4,1});
        maxfreq = str2num(answer{5,1});
        downsamp = str2num(answer{6,1});
        wfchoice = answer{7,1};
        waterfall = str2num(answer{8,1});
        doresultsleaf = str2num(answer{9,1});
    else
        out = [];
        return
    end
end

[len,chans,bands] = size(audio); % size of the audio



% apply calibration factor if it has been provided
if exist('cal','var')
    audio = audio .* 10.^(repmat(cal,[len,1,bands])/20);
end

if wfchoice == 'h'
    wf = window(@hann,winlen);
elseif wfchoice == 'b'
    wf = window(@blackmanharris,winlen);
else
    wf = ones(winlen,1);
end

if ~(downsamp == 1)
    if downsamp < 1, downsamp = 1/downsamp; end
    audio = resample(audio, 1,downsamp);
    fs = fs / downsamp;
    Nyquist = fs/2;
end

if maxfreq > Nyquist, maxfreq = Nyquist; end

titletext = 'Spectrogram';

if isfield(in,'chanID') % Get the channel ID if it exists
    chanID = in.chanID;
else
    chanID = cellstr([repmat('Chan',size(data,2),1) num2str((1:size(data,2))')]);
end
if isfield(in,'bandID') % Get the band ID if it exists
    bandID = in.bandID;
else
    bandID = 1:size(audio,3);
end


figure('Name', titletext)

k = 1; % subplot counter

if ~transposesubplots
    for ch = 1:chans
        for b = 1:bands
            subplot(chans,bands,k)
            [~,F,T,P] = spectrogram(audio(:,ch,b),wf,NOVERLAP,[],fs);
            L = 10*log10(abs(P));
            if doresultsleaf==1
            if ~exist('L_all','var')
                L_all = zeros(size(L,1),size(L,2),chans,bands);
            end
            L_all(:,:,ch,b) = L;
            end
            m =max(max(L));
            minrange = m - abs(dBrange);
            if waterfall
                mesh(T,F,L)
                zlim([minrange m])
            else
                imagesc(T,F,L)
                set(gca,'YDir','normal');
            end
            box('on')
            
            set(gca, 'Clim', [minrange m]);
            ylim([0 maxfreq]);
            
            if ch == chans
                xlabel('Time (s)')
            end
            if b == 1
                ylabel(['Chan ',num2str(ch)])
            end
            if ch == 1
                title(['Band ',num2str(b)])
            end
            k = k+1;
        end
    end
else
    for b = 1:bands
        for ch = 1:chans
            subplot(bands,chans,k)
            [~,F,T,P] = spectrogram(audio(:,ch,b),wf,NOVERLAP,[],fs);
            L = 10*log10(abs(P));
            if doresultsleaf==1
            if ~exist('L_all','var')
                L_all = zeros(size(L,1),size(L,2),chans,bands);
            end
            L_all(:,:,ch,b) = L;
            end
            if waterfall
                mesh(T,F,L)
                zlim([minrange m])
            else
                imagesc(T,F,L)
                set(gca,'YDir','normal');
            end
            box('on')
            m =max(max(L));
            minrange = m - abs(dBrange);
            set(gca, 'Clim', [minrange m]);
            ylim([0 maxfreq]);
            
            if b == bands
                xlabel('Time (s)')
            end
            if b == 1
                title(['Chan ',num2str(ch)])
            end
            if ch == 1
                ylabel(['Band ',num2str(b)])
            end
            k = k+1;
        end
    end
end

if doresultsleaf==1
   % L_all(isinf(L_all))=min(min(min(min(L_all(~isinf(L_all))))));
doresultleaf(L_all,['level' '[dB]'],{'time'},...
                 'frequency',   F,      'Hz',          true,...
                 'time',        T,      's',           true,...
                 'channels',    chanID, 'categorical', [],...
                 'bands',       bandID, 'Hz',          false,...
                 'name',        'Spectrogram');
end
% out.transposesubplots = transposesubplots;
% out.window = winlen;
% out.noverlap = NOVERLAP;
% out.dBrange = dBrange;
% out.maxfreq = maxfreq;
% out.downsamp = downsamp;
% out.wfchoice = wfchoice;
% out.waterfall = waterfall;
out.funcallback.name = 'spectrogram_simple.m';
out.funcallback.inarg = {fs,transposesubplots,winlen,NOVERLAP,dBrange,maxfreq,downsamp,wfchoice,waterfall,doresultsleaf};

%**************************************************************************
% Copyright (c) 2014, Densil Cabrera
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