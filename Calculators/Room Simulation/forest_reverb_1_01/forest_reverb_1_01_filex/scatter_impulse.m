function [y Y]=scatter_impulse(a,phi,F_s,N_fft,N_m,f_f, minphase,c)
%scatter_impulse(a, phi, F_s, N_fft, N_m, f_f, minphase, c)
%this calculates the impulse response and complex transfer function for scattering
%off of a hard cylinder.
% a - (m) cylinder radius
% phi - (rad) scattering angle (0 is no path change, pi is reflection
% F_s - (Hz) sample rate
% N_fft - size of dft (works best as a power of 2)
% N_m - number of terms to calculate in series. Set to 0 to automatically
%   set for a tolerance of approx 1e-6
% f_f - low pass filter cutoff frequency (applies a 4th order Butterworth
%   filter forward and backward for no phase distortion)
% minphase - this may be set to true to alter the phase response to that of
%   a minimum phase filter
% c - (m/s) speed of sound

%This is currently copyright Travis Wiens 2008 but I intend to GPL it. 
%Please contact me if you want me to hurry this process up or you want
%a commerical license. email: travis.mlfx@nutaksas.com

if nargin<1
    a=.05;%(m) tree radius
end
if nargin<2
    phi=pi*0.5;%(rad) angle (0 is no change of path, pi is full reflection)
end
if nargin<3
    F_s=44100;%(Hz) sample frequency
end
if nargin<4
    N_fft=1024;%number of points in fft
end
if nargin<5
    N_m=0;%number of terms in series
end
if nargin<6
    f_f=0;%(Hz)filter cutoff frequency
end
if nargin<7
    minphase=false;%use minimum phase response
end
if nargin<8
   c=340;%(m/s)speed of sound
end

if round(N_fft/2)==N_fft/2
    f=((1:((N_fft)/2)))*F_s/N_fft;%(Hz) frequencies in dft
else
    f=((1:((N_fft-1)/2)))*F_s/N_fft;
end


gain_P=zeros(1,numel(f));%gain at frequency f
for i=1:numel(f);
    [gain_P(i)]=scatter_cyl(a,phi,f(i),N_m,c);%calculate complex gain at each frequency
end

if round(N_fft/2)==N_fft/2
    %reconstruct dft
    Y=[0 gain_P(1:end-1) abs(gain_P(end)) conj(fliplr(gain_P(1:end-1)))];
else
    Y=[0 gain_P conj(fliplr(gain_P))];
end

y=ifft(Y);%get impulse response from dft

if f_f>0
    %lowpass filter response
    N_butter=4;%order of butterworth filter
    [b a]=butter(N_butter,f_f/(F_s/2));
    y_tmp=[y y y];%to avoid initial condition effects
    y_tmp=filtfilt(b,a,y_tmp);
    y=y_tmp(N_fft+(1:N_fft));
    Y=fft(y);%reconstruct dft
end

if minphase
    Y=mpf(abs(Y));%convert response to minimum phase (not changing magnitude)
    y=real(ifft(Y));
end