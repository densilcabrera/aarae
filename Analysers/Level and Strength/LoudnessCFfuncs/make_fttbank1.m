function S = make_fttbank1(fs, cf)
% fcoefs = make_fttbank1(fs,cf);
% calculates sos-representation of a critical band filterbank with 24
% filters (bandwidth  = 1 bark) with an impulse response of 4th order.
% This implementation is based on M. Slaney "An efficient
% implementation of the Patterson-Holdsworth filter bank"
% see Apple Technical Report #35. It is an efficient way to implement an
% FTT filterbank or Gammatone filterbank.
% Note that this implementation results in filters with less stop band
% attenuation than complex modulation filterbanks.
% cf: column vector mit 24 center-frequencies (optional, default:
% Bark center-frequencies according to Zwicker)
% B: column vector mit 24 critical band widths (optional, default:
% Bark bandwidths according to Zwicker)
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Altered by flatmax is Matt Flax for the Psy-Sound project
% Jan. 2007

if nargin < 2
    cf=fliplr([50,150,250,350,450,570,700,840,1000,1170,1370,1600,1850,2150,2500,2900,3400,4000,4800,5800,7000,8500,10500,13500])';
end    
ERB=fliplr([100,100,100,100,110,120,140,150,160 ,190 ,210 ,240 ,280 ,320 ,380 ,450 ,550 ,700 ,900 ,1100,1300,1800,2500 ,3500])';
B=1.019*2*pi*ERB; %3dB-bandwidth = 0.887*ERB
T=1/fs;

A0 = T;
A2 = 0;
B0 = 1;
B1 = -2*cos(2*cf*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;
A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;
A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;
A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;

gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                         (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                          sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cf*pi*T))).* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
           (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
           2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);
	
allfilts = ones(length(cf),1);
S.fcoefs = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];

% filter save states ....
for j=1:length(gain)
  S.Zf1{j} = [];
  S.Zf2{j} = [];
  S.Zf3{j} = [];
  S.Zf4{j} = [];
end