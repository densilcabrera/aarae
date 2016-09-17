function OUT = stepped_tone_log(duration, lof, hif, fstep, fs)
% Generates a sequence of pure tones (with frequencies exponentially
% spaced), separated by silence. By default each tone has a duration of 1
% s (or a little less, because the tone is present for an integer number of
% cycles to avoid clicks). The silence duration is half that of the nominal
% tone duration. This signal may be useful in the direct measurement of
% harmonic distortion, or indeed in the measurement of linear transfer
% function (although there are much more efficient ways of doing that!).
%
% The frequencies used are adjusted so that they have an integer number of
% samples in their exact period - which facilitates analsysis (e.g., using
% an fft with an analysis period equal to that of the fundamental).
% However, when frequencies are high, this may result in a substantial
% adjustment unless a high sampling rate is used (which is advisable in
% harmonic distortion analysis anyway).
%
% The highest frequency allowed is fs/4, and the lowest is 10 Hz.
%
% The tone sequence is replicated in audio2. Additional properties fields
% are generated to be used by an analyser.
%
%
% Code by Densil Cabrera
% Version 1.00 (7 August 2014)


if nargin == 0
    param = inputdlg({'Approximate duration of each tone (minimum 1 s) [s]';...
                       'Lowest tone frequency (integer) [Hz]';...
                       'Highest tone frequency (integer) [Hz]';...
                       'Fractional octave step (e.g. 1 | 3 | 12)';...
                       'Sampling rate [Hz]'},...
                       'Stepped tone sequence',1,...
                       {'1';'125';'8000';'1';'48000'});
    param = str2num(char(param));
    if length(param) < 5, param = []; end
    if ~isempty(param)
        duration = param(1);
        lof = param(2);
        hif = param(3);
        fstep = param(4);
        fs = param(5);
    else
        OUT = [];
        return
    end
end
if ~isempty(param) || nargin ~= 0
    if duration < 1, duration = 1; end
    if lof < 10, lof = 10; end
    if hif > fs/4, hif = fs/4; end
    flist = lof * 2.^(0:1/fstep:ceil(log2(hif/lof)));
    % adjust frequencies so that they have integer periods for the chosen fs
    flist = 1./(round(fs./flist)./fs);  
    t = linspace(0,duration,fs*duration);
    y = zeros(length(flist) * 1.5 * duration * fs, 1);
    startindex = 1+((1:length(flist))-1)*1.5*duration*fs;
    % make sure only an integer number of cyles is generated
    tonelen = floor(fs*duration./(fs./flist)) * fs./flist;
    endindex = startindex + tonelen-1;
    for n = 1:length(flist)
        y(startindex(n):endindex(n)) = sin(2*pi*flist(n).*t(1:tonelen(n))');
    end
    
    tag = 'SteppedToneSeqence';
    
    OUT.audio = y;
    OUT.audio2 = y; % replicate audio in 2nd output
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.properties.flist = flist;
    OUT.properties.startindex = startindex;
    OUT.properties.endindex = endindex;
    OUT.funcallback.name = 'stepped_tone_log.m';
    OUT.funcallback.inarg = {duration, lof, hif, fstep, fs};
else
    OUT = [];
end

end % End of function