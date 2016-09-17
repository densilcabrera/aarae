function OUT = multicycleflatspectrum(cycles)
% This function loads an audio waveform, and flattens its spectrum (apart
% from DC, which is zeroed). The audio is concatenated with itself to yield
% a multicycle signal (similar to MLS or IRS). An inverse filter is created
% (which is simply the time-reversed spectrally flattened signal).
%
% Like MLS or IRS, this test signal may be used to measure impulse
% responses. Recordings should be analysed using circular (not linear)
% processing: use CircXcorrforIR (exactly the same as for IRS).


OUT = choose_audio; % call AARAE's choose_audio function
if ~isempty(OUT)
    answer = inputdlg({'Number of cycles (>=2)';...
        'Remove DC [0 | 1]'},...
            'multicycleflatspectrum settings',1,{'2','1'});
        answer = str2num(char(answer));
        cycles = round(abs(answer(1)));
        removeDC = answer(2);
        if cycles < 2, cycles = 2; end
    audio = fft(OUT.audio);
    meanval = rms(abs(audio));
    if removeDC == 1
    audio = real(ifft([zeros(size(meanval));...
        repmat(meanval,[size(audio,1)-1,1])]...
        .* exp(1i .* angle(audio))));
    else
        audio = real(ifft(...
        repmat(meanval,[size(audio,1),1])...
        .* exp(1i .* angle(audio))));
    end
    if max(max(max(max(max(abs(audio)))))) > 1
        % this normalization should be unnecessary, but just in case...
        audio = audio ./ max(max(max(max(max(abs(audio))))));
    end
    OUT.audio = [repmat(audio,[cycles,1,1,1,1,1]);zeros(size(audio))];
    OUT.audio2 = flipud(audio);
    if isfield(OUT,'name')
        OUT.tag = [OUT.name ' gen'];
    else
        OUT.tag = 'multicycleflatspectrum';
    end
    OUT.properties.cycles = cycles;
    OUT.funcallback.name = 'multicycleflatspectrum.m';
    OUT.funcallback.inarg = {cycles};
else
    OUT = [];
    return
end

