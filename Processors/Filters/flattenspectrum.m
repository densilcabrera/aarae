function OUT = flattenspectrum(IN,writeaudio2)
% This function returns a flat spectrum magnitude audio signal, while
% retaining the phase of the input. DC is zeroed.
%
% Optionally, an audio2 field is written (as a kind of inverse filter, for
% educational demonstration purposes).


if isstruct(IN)
    audio = IN.audio;
    if nargin < 2
        answer = inputdlg({'Write inverse filter in audio2 [0 | 1]?'},...
            'Flatten Spectrum settings',1,{'0'});
        writeaudio2 = str2num(char(answer));
    end
else
    audio = IN;
end




if ~isempty(audio)
    audio = fft(audio);
    meanval = rms(abs(audio));
    audio = real(ifft([zeros(size(meanval));...
        repmat(meanval,[size(audio,1)-1,1])]...
        .* exp(1i .* angle(audio))));
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
        if writeaudio2 == 1
            OUT.audio2 = flipud(audio);
        end
        OUT.funcallback.name = 'flattenspectrum.m';
        OUT.funcallback.inarg = {writeaudio2};
    else
        OUT = processed;
    end
else
    OUT = [];
end