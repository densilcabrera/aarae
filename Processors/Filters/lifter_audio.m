function OUT = lifter_audio(in, fs, t1, t2, wl)
% This function lifters the cepstrum of a wave between t1 and t2 (in
% milliseconds, in the quefrency domain) for a wave with sampling frequency
% fs (in Hz). wl is the window length (in millisecods).
%
% This can be used to smooth the spectrum of the audio, or for other
% purposes. Processing is done in a single window (rather than using
% multi-window processing).
%
% Code by Densil Cabrera
% version 1.01 (30 September 2014)


if nargin < 5, wl = 1000; end
if nargin < 4, t2 = 200; end
if nargin < 3
    t1 = 0;
    if isstruct(in)
        in = choose_from_higher_dimensions(in,3,1);
        data = in.audio;
        fs = in.fs;
    else
        data = in;
        if nargin < 2
            fs = inputdlg({'Sampling frequency [samples/s]'},...
                               'Fs',1,{'48000'});
            fs = str2num(char(fs));
        end
    end
    len = size(data,1);
    %dialog box for settings
    prompt = {'Lifter start time (ms)', ...
        'Lifter end time (ms)', ...
        'Cepstrum window length (ms)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0','200',num2str(floor(1000*len/fs))};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        t1 = str2num(answer{1,1});
        t2 = str2num(answer{2,1});
        wl = str2num(answer{3,1});
    end
end
if isstruct(in)
    data = in.audio;
    fs = in.fs;
else
    data = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

[~,chans,bands] = size(data);

if ~isempty(data) && ~isempty(fs) && ~isempty(t1) && ~isempty(t2) && ~isempty(wl)
    nsamples = ceil(wl * 0.001 * fs);

    % cepstral window in samples
    s1 = round(0.001 * t1 * fs) + 1;
    s2 = round(0.001 * t2 * fs) + 1;

    % avoid errors
    if s1 < 1 || s1 > nsamples-1
        s1 = 1;
    end
    if s2 < 2 || s2 > nsamples
        s2 = nsamples;
    end
    if s1 >= s2
        s1 = 1;
        s2 = nsamples;
    end

    if nsamples < length(data)
        x= data(1:nsamples,:,:);
    else
        x = zeros(nsamples,chans,bands);
        x(1:length(data),:,:)=data;
    end

    if ~isreal(x)
        disp('Real values only have been used for liftering')
        x = real(x);
    end
    
    ND = zeros(chans,bands); % used for circular shift in cceps and icceps
    % cceps only operates on vectors (not matrices) - hence the following loop
    cepstrum = zeros(nsamples,chans,bands);
    for ch = 1:chans
        for b = 1:bands
            [cepstrum(:,ch,b),ND(ch,b)] = cceps(x(:,ch,b));
        end
    end

    lifteredcepstrum = zeros(nsamples,chans,bands);
    lifteredcepstrum(s1:s2,:,:) = cepstrum(s1:s2,:,:);
    lifteredcepstrum(nsamples-s2+1:nsamples-s1+1,:,:) = cepstrum(nsamples-s2+1:nsamples-s1+1,:,:);

    y = zeros(size(cepstrum));
    for ch = 1:chans
        for b = 1:bands
            y(:,ch,b)=icceps(lifteredcepstrum(:,ch,b),ND(ch,b));
        end
    end
    
    if isstruct(in)
        OUT.audio = y;
        OUT.funcallback.name = 'lifter_audio';
        OUT.funcallback.inarg = {fs,t1,t2,wl};
    else
        OUT = y;
    end


else
    OUT = [];
end

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