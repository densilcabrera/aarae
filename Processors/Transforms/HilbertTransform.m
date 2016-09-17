function out = HilbertTransform(in,method)
% A Hilbert transform is applied to the input audio. Note that the
% resulting analytic waveform is complex. This can be dealt with in various
% ways:
% 1. a complex wave is returned
% 2. the real wave is returned to audio, and imaginary to audio2
% 3. The absolute value (or envelope) of the complex wave is returned

% other options... probably not useful (so not implemented yet)
% 4. The instantaneous phase is returned
% 5. The instantaneous frequency is returned (not audio)



if nargin < 2
method = menu('Output Format', ...
              'Complex wave', ...
              'Real (audio) and Imag (audio2)',...
              'Absolute value');
end
% Matlab's hilbert() seems to only work on up to 2-dimensional data,
% hence the following loop to allow for 6 dimensions
[len,chans,bands,dim4,dim5,dim6] = size(in.audio);
transformed = zeros(len,chans,bands,dim4,dim5,dim6);
for b = 1:bands
    for d4 = 1:dim4
        for d5 = 1:dim5
            for d6 = 1:dim6
    transformed(:,:,b,d4,d5,d6) = hilbert(in.audio(:,:,b,d4,d5,d6));
            end
        end
    end
end

switch method
    case 1
        out.audio = transformed;
    case 2
        out.audio = real(transformed);
        out.audio2 = imag(transformed);
    otherwise
        out.audio = abs(transformed);
end
if isstruct(in)
    out.funcallback.name = 'HilbertTransform.m';
    out.funcallback.inarg = {method};
end