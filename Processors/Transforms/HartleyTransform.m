function out = HartleyTransform(in,complexchoice)
% This function performs the Hartley transfrom of the input audio.
% 
% The Hartley transform is closely related to the Fourier transform,
% but the result of transforming a real signal is entirely real. 
% Furthermore, the Hartley transform is its own inverse transform.

if ~isreal(in.audio) && nargin == 1
    complexchoice = questdlg('Audio data is complex - how do you wish to treat it?', ...
	'Hartley Transform', ...
	'Use Complex','Use Real','Use Real');
end
if ~isreal(in.audio)
    switch complexchoice
        case 'Use Complex'
            X = fft(in.audio);
        otherwise
            % use real
            X = fft(real(in.audio));
            
            %X = fft(abs(in.audio).*cos(angle(in.audio))); %inv-Hilbert
    end
else
X = fft(in.audio);
complexchoice = 'default';
end
%X = fft(real(in.audio)); % alternatively only transform real data
out.audio = (real(X) - imag(X)) ./ size(X,1).^0.5;
out.funcallback.name = 'HartleyTransform.m';
out.funcallback.inarg = {complexchoice};