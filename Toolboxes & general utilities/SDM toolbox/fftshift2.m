function y = fftshift2(x)
% y = fftshift2(x)
% Circularly shifts matrix x with respect to first dimension size(x,1)/2 
%
% Does the same as fftshift but faster
%

% SDM toolbox : fftshift2
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

i = [(size(x,1)/2+1):size(x,1), (1:size(x,1)/2)];

y = x(i,:);