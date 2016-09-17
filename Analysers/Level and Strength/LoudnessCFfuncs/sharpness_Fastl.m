function [sharp] = sharpness_Fastl(loudspec)

% SHARPNESS
%**************************************************************
% Method FASTL (1991)
% Expression for weighting function obtained by fitting an 
% equation to data given in 'Psychoacoustics: Facts and Models'
% using MATLAB basic fitting function
% x: time signal
% sh = sharpness [acum]
%**************************************************************
% Claire Churchill Sep 2004

n = length(loudspec);

gz(1:140) = 1;
z = 141:n;
gz(z) = 0.00012*(z/10).^4-0.0056*(z/10).^3+0.1*(z/10).^2-0.81*(z/10)+3.5;

z = 0.1:0.1:(n/10);

sharp = 0.11 * sum(loudspec.*gz.*z.*0.1) / sum(loudspec.*0.1 +eps);