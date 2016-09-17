function y = convolvedemo(wave1, wave2, method, fs)
% This function performs convolution using a variety of methods (for
% demonstration purposes). The various methods may have different errors
% due to precision, and some are much faster than others.
% Method 2 is used by AARAE's framework (frequency domain multiplication).

len1 = length(wave1);
len2 = length(wave2);
outputlength = len1 + len2 - 1;
%times = ((1:outputlength)-1) ./ fs;

switch method
    case 1
        y = conv(wave1, wave2);
    case 2
        y = ifft(fft(wave1, outputlength) .* fft(wave2, outputlength));
    case 3
        y = filter(wave1, 1, [wave2; zeros(len1-1,1)]);
    case 4
        if len1 < len2
            wave1 = [wave1; zeros(len2-len1,1)];
        elseif len2 < len1
            wave2 = [wave2; zeros(len1-len2,1)];
        end %if
        y = xcorr(wave1,flipud(wave2));
        y = y(1:outputlength);
    case 5
       bigmatrix = toeplitz([wave1;zeros(len2-1,1)], [wave1(1); zeros(outputlength-1,1)]);
       bigmatrix = bigmatrix .* repmat([wave2;zeros(len1-1,1)]',outputlength,1);
       y = sum(bigmatrix, 2);
end % switch