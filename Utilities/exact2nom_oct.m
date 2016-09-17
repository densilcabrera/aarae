function y = exact2nom_oct(x)
% This function returns nominal 1/3-octave band (and octave band) centre 
% frequencies corresponding to the input frequencies - which are exact or
% approximate values.
%
% It assumes that the input values are within 1/6-octave of the exact
% 1/3-octave band centre frequency.
%
% Input can be up to 3 dimensional.
%
% Code by Densil Cabrera
% Version 1.00 (2 December 2013)

% list of nominal 1/3 octave frequencies (x 10^d) within any decade (d)
nomfreqtemplate = [1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10];
exactfreqtemplate = 10.^((0:10)/10);
highcutoff = exactfreqtemplate .* 10.^(0.05);

% calculate the decade of each value
d = 10.^floor(log10(x));

% allow for 3-dimensional input
[R, C, B] = size(x);
y = zeros(R,C,B);

% find nominal frequencies
for r = 1:R;
    for c = 1:C
        for b = 1:B
           ind = find((x(r,c,b)/d(r,c,b))<highcutoff,1,'first');
           y(r,c,b) = d(r,c,b) .* nomfreqtemplate(ind);
        end
    end
end






