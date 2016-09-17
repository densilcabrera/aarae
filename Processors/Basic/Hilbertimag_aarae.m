function out = Hilbertimag_aarae(in)
% A Hilbert transform is applied to the input audio, and its imaginary value
% is returned. Hence a quadrature phase shift is applied to the input audio.
        [len,chans,bands,dim4,dim5,dim6] = size(in.audio);
    out.audio = zeros(len,chans,bands,dim4,dim5,dim6);
    for b = 1:bands
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
        out.audio(:,:,b,d4,d5,d6) = imag(hilbert(in.audio(:,:,b,d4,d5,d6)));
                end
            end
        end
    end
    out.funcallback.name = 'Hillbertmag_aarae.m';
    out.funcallback.inarg = {};
end