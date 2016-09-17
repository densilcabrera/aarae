function [gain_P gamma psi_m]=scatter_cyl(a, phi, f, N_m, c)
%[gain_P gamma psi_m]=scatter_cyl(a, phi, f, N_m, c)
%Calculates the complex frequency response of scattering of a wave by a hard
%cylinder at a large distance, using the theoretical method from Morse [1948]. 
% gain_P - complex gain. This is the pressure gain, not intensity gain.
% gamma - scattered wave coefficient (see Morse)
% psi_m - scattered wave coefficient (see Morse)
% a - (m) radius
% phi - (rad) scattering angle (0 is no path change, pi is reflection
% f - (Hz) frequency
% N_m - number of terms to calculate in series. Set to 0 to automatically
%   set, for a tolerance of approx 1e-6 on gamma. The automatic N_m was 
%   determined experimentally and may not be accurate for very large or
%   small f*a/c.
% c - (m/s) speed of sound
%
%Reference: Morse, Philiop M., "Vibration and Sound 2nd Ed", McGraw Hill,
%New York, pp 347, 1948.

%This is currently copyright Travis Wiens 2008 but I intend to GPL it. 
%Please contact me if you want me to hurry this process up or you want
%a commerical license. email: travis.mlfx@nutaksas.com

if nargin<1
    a=1;%(m) tree radius
end
if nargin<2
    phi=pi/3;%(rad) scattering angle
end
if nargin<3
    f=1000;%(Hz) freqency
end
if nargin<4
    N_m=0;%automatically estimate N_m
end
if nargin<5
    c=340;%(m/s) sound speed
end
r=1;%(m) radius



lambda=c/f;%(m) wavelength
k=2*pi/lambda;%(1/m)wave number 

if N_m<1;
    N_m=ceil(1.1*k*a+9.7);%use 9.7 for 10^-6, 10.7 for 10^-7, 12.50 for 10^-9
end

epsilon=2*ones(1,N_m+1);
epsilon(1)=1;

bj=besselj(0:(N_m+1),k*a);%precalculate bessel functions
by=bessely(0:(N_m+1),k*a);

gamma=zeros(1,N_m+1);

%gamma(1)=atan2(besselj(1,k*a),-bessely(1,k*a));
gamma(1)=atan2(bj(2),-by(2));

for m=1:N_m
    gamma(m+1)=atan2(-(bj(m)-bj(m+2)),-(by(m+2)-by(m)));
%    gamma(m+1)=atan2(-(besselj(m-1,k*a)-besselj(m+1,k*a)),-(bessely(m+1,k*
%    a)-bessely(m-1,k*a)));

end

psi_m=zeros(1,N_m+1);%this is psi
clear i


for m=0:N_m
    psi_m(m+1)=epsilon(m+1)*sin(gamma(m+1))*exp(-i*gamma(m+1))*cos(m*phi);
end
psi_m=psi_m/sqrt((k*a));
psi=sum(psi_m);
    
    
%gain_I=2*a/(pi*r)*abs(psi).^2;
gain_P=sqrt(2*a/(pi*r))*psi;

