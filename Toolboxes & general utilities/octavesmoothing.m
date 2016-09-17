function [smoothmagspectrum]= octavesmoothing(fftdbin, octsmooth,fs)

% fftdbin is the fft of the signal in dB
% this function goes a bit funny in the very low frequency range, so it's
% best to throw that data out after running this function.
% Also, it makes more sense to use this to smooth the power spectrum rather
% than level spectrum.
fftpoints = length(fftdbin);

freqeval = linspace(0,fs/2,fftpoints);

if octsmooth > 0
    octsmooth=2*octsmooth;
    % octave center freqevaluencies
    f1=1;
    i=0;
    while f1 < (fs/2)
        f1=f1*10^(3/(10*octsmooth));
        i=i+1;
        fc(i,:)=f1;
    end

    % octave edge freqevaluencies
    for i=0:length(fc)-1
        i=i+1;
        f1=10^(3/(20*octsmooth))*fc(i);
        fe(i,:)=f1;
    end

    % find nearest freqevaluency edges
    for i=1:length(fe)
        fe_p=find(freqeval>fe(i),1,'first');
        fe_m=find(freqeval<fe(i),1,'last');
        fe_0=find(freqeval==fe(i));
        if isempty(fe_0)==0
            fe(i)=fe_0;
        else
            p=fe_p-fe(i);
            m=fe(i)-fe_m;
            if p<m
                fe(i)=fe_p;
            else
               fe(i)=fe_m;
            end
        end
    end

    for i=1:length(fe)-1
        fftdbin_i=fftdbin(fe(i):fe(i+1),:);
        smoothmagspectrum(i,1:size(fftdbin,2))=mean(fftdbin_i);
    end
    fc=fc(2:end);
    if ~isreal(fftdbin)
        absfftdbin=abs(fftdbin);
        for i=1:length(fe)-1
            abssmoothspec_i=absfftdbin(fe(i):fe(i+1),:);
            abssmoothspec(i,1:size(fftdbin,2))=mean(abssmoothspec_i);
        end
        smoothmagspectrum = smoothmagspectrum.*abssmoothspec./abs(smoothmagspectrum);
    end
    smoothmagspectrum = interp1(fc,smoothmagspectrum,freqeval,'spline');
    if size(smoothmagspectrum,1) < size(smoothmagspectrum,2)
        smoothmagspectrum = smoothmagspectrum';
    end
else
    smoothmagspectrum = fftdbin;
end