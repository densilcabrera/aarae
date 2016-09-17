function H = minphasefreqdomain(mag,thresholddB,timeshift)
% This function caclulates the frequency domain complex coefficients
% equivalent to a minimum phase filter, using a magnitude spectrum as
% input (which is assumed to be real valued).
%
% INPUTS:
% mag - magnitude spectrum from which to derive the complex coefficients
%
% thresholddB - threshold in decibels (relative to the maximum value of the
% magnitude spectrum): values below this are assigned this value to avoid
% problems that arise with very large dynamic ranges or deep notches in the
% spectrum. The default value is -120 dB. Note that this only affects the
% phase response of the filter - the output magnitude of the filter is
% replicated from the function's input mag.
%
% timeshift (0 or 1) - apply a circular time shift so that the impulse
% response starts and ends at a zero-crossing. If this input is missing,
% then a time shift is not applied. This time shift is most useful if the
% filter is going to be used with linear convolution (rather than simply by
% spectrum multiplication equivalent to circular convolution). However, by
% shifting the impulse response, the minimum phase characteristic of the
% filter is compromised. Note that a longer filter may reduce the need for
% this operation (by reducing temporal aliasing, or filter impulse response
% wrap-around). Other approaches to fixing the ends of the impulse response
% could be to apply a small fade-in and fade-out (a window function) after
% calling this function, or to smooth the spectrum magnitude (and reduce
% its contrast) prior to calling this function.
%
%
%
% To run this filter, multiply in frequency domain, and then return to time
% domain using ifft.


% APPLY THRESHOLD FOR CEPSTRUM PROCESSING
H = mag;
%H(H==0) = 1e-200;
if ~exist('thresholddB','var'), thresholddB = 120; end
if isempty(thresholddB), thresholddB = 120; end
threshold = db2mag(-abs(thresholddB));
H(H<max(max(H))*threshold) = max(max(H))*threshold;

[n,chans] = size(mag);
% CALCULATE THE CEPSTRUM
H = real(ifft(log(H)));

% WINDOW IN CEPSTRUM DOMAIN
w = [1;2*ones(round(n/2)-1,1);ones(1-rem(n,2),1);zeros(round(n/2)-1,1)];

% GENERATE FREQUENCY DOMAIN COEFFICIENTS
H = exp(fft(repmat(w,[1,chans]).*H));

% FIX PRECISION ERRORS THAT CREATE SMALL IMAGINARY COMPONENTS AFTER IFFT
% This is done by averaging above and below the Nyquist frequency
% (This avoids the need to use real() after the ifft.)
H(2:end,:) = (H(2:end,:)+flip(conj(H(2:end,:))))./2;

% FORCE THE MAGNITUDE TO PRECISELY MATCH THE INPUT MAGNITUDE
% Note that this may contribute to the impulse response wrapping around
% (for components below thresholddB), although it is not the only cause of this.
%H = H .* mag ./ abs(H);
H = mag.*exp(1i*angle(H));

% OPTIONALLY APPLY CIRCULAR TIME-SHIFT TO START AND END AROUND A ZERO CROSSING
if exist('timeshift','var')
    if timeshift ~= 0
        H = ifft(H);
        for ch = 1:chans
            % Find the last zero-crossing
            s=find((abs(H(:,ch))<abs([H(2:end,ch);H(1,ch)]))...
                &(abs(H(:,ch))<abs([H(end,ch);H(1:end-1,ch)])),1,'last');
            H(:,ch) = circshift(H(:,ch),n-s);
        end
        H = fft(H);
    end
end
