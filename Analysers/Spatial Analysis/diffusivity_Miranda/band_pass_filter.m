function [hpF lpF] = band_pass_filter(low_octave_band, high_octave_band, fs)


N     = 6;      % Order
Fpasslow = high_octave_band.*sqrt(2);   % Passband Frequency
Fpasshigh = low_octave_band./sqrt(2);   % Passband Frequency
Apass = 1;      % Passband Ripple (dB)
Astop = 60;     % Stopband Attenuation (dB)

hlp = fdesign.lowpass('n,fp,ap,ast', N, Fpasslow, Apass, Astop, fs);

lpF = design(hlp, 'ellip','SOSScaleNorm', 'Linf');

hhp = fdesign.highpass('n,fp,ast,ap', N, Fpasshigh, Astop, Apass, fs);

hpF = design(hhp, 'ellip', 'SOSScaleNorm', 'Linf');
