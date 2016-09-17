function out = gain_aarae(in,gain)
% This function applies gain to the audio.
% Gain is specified in dB, where:
% 0 dB means no change
% a positive value amplifies the waveform
% a negative value attenuates the waveform
%
% The default value normalizes the audio, whilst preserving the amplitude
% relationships between the channels and bands (where relevant)
%
% Inputting 'n' normalizes each data column individually
%
% If a complex valued gain is input, then it is used directly to multiply the
% waveform (instead of converting it from decibels). A purely real argument
% that includes 0i will be also be used to directly multiply the waveform.
%
% Code by Densil Cabrera
% version 1.01 (14 February 2013)

if isstruct(in)
    Lmax = 20*log10(max(max(max(abs(in.audio)))));
    prompt = {['Real gain (in dB), or ''n'' (max is ',num2str(Lmax),...
        ' dBFS), or complex gain factor (not in dB)']};
    dlg_title = 'Gain';
    num_lines = 1;
    def = {num2str(-Lmax)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
else
    in.audio = in;
    answer = gain;
    if ~ischar(answer)
        answer = num2str(answer);
    end
end


if ~isempty(answer)
    gain = char(answer{1,1});
    out.audio = in.audio;
    indices = partial_selection(in);
    if gain == 'n'
        % normalize column individually
        maxval = max(abs(in.audio(indices{:})));
        out.audio(indices{:}) = in.audio(indices{:}) ./ repmat(maxval,[length(in.audio(indices{:})),1,1,1,1,1]);
    elseif ~isempty(gain(gain == 'i'))
        gain = str2double(gain);
        out.audio(indices{:}) = in.audio(indices{:}) * gain;
    else
        gain = str2double(gain);
        out.audio(indices{:}) = in.audio(indices{:}) * 10.^(gain/20);
    end
    out.funcallback.name = 'gain_aarae.m';
    out.funcallback.inarg = {gain};
else
    out = [];
end

end