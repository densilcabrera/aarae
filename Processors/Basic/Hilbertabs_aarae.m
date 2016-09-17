function out = Hilbertabs_aarae(in)
% A Hilbert transform is applied to the input audio, and its absolute value
% is returned. This is a simple way to derive the envelope function
% of the audio.
    [len,chans,bands,dim4,dim5,dim6] = size(in.audio);
    out.audio = zeros(len,chans,bands,dim4,dim5,dim6);
    for b = 1:bands
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
        out.audio(:,:,b,d4,d5,d6) = abs(hilbert(in.audio(:,:,b,d4,d5,d6)));
                end
            end
        end
    end
    out.funcallback.name = 'Hilbertabs_aarae.m';
    out.funcallback.inarg = {};
end