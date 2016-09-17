function OUT = loadaudiogenerator(invfiltmethod)
% The purpose of this very simple generator is to allow you to load an
% audio file while writing the properties fields that are normally
% associated with a generator. An inverse filter can be written in audio2
% if it is not present in the selected audio.





OUT = choose_audio; % call AARAE's choose_audio function
if ~isempty(OUT)
    if ~isfield(OUT,'audio2')
        answer = inputdlg({'Audio2: No inverse filter [0]?; Inverse filter excluding spectrum components <-60 dB [1]; Inverse filter excluding spectrum components <-80 dB [2]; Inverse filter excluding spectrum components <-100 dB [3]'},...
            'Audio2 inverse filter settings',1,{'1'});
        invfiltmethod = str2num(char(answer));
        switch invfiltmethod
            case 1
                % zap weak components (crude but useful) at -60 dB
                spectrum = fft(OUT.audio);
                maxval = max(max(max(max(max(abs(spectrum))))));
                badcomponents = abs(spectrum)<maxval/1e3;
                spectrum = 1./spectrum;
                spectrum(badcomponents) = 0;
                OUT.audio2 = ifft(spectrum);
            case 2
                % zap weak components (crude but useful) at -80 dB
                spectrum = fft(OUT.audio);
                maxval = max(max(max(max(max(abs(spectrum))))));
                badcomponents = abs(spectrum)<maxval/1e4;
                spectrum = 1./spectrum;
                spectrum(badcomponents) = 0;
                OUT.audio2 = ifft(spectrum);
            case 3
                % zap weak components (crude but useful) at -100 dB
                spectrum = fft(OUT.audio);
                maxval = max(max(max(max(max(abs(spectrum))))));
                badcomponents = abs(spectrum)<maxval/1e5;
                spectrum = 1./spectrum;
                spectrum(badcomponents) = 0;
                OUT.audio2 = ifft(spectrum);
                % write more methods here!
        end
    end
    if isfield(OUT,'name')
        OUT.tag = [OUT.name ' gen'];
    else
        OUT.tag = 'loadaudiogenerator';
    end
    OUT.properties.invfiltmethod = invfiltmethod;
    OUT.funcallback.name = 'loadaudiogenerator.m';
    OUT.funcallback.inarg = {invfiltmethod};
else
    OUT = [];
    return
end


