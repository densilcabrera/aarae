function [b, a] = butter_hp(fs)
% function [b, a] = butter_hp(fs);
% Filter to approximate transmission through outer and middle ear according
% to E. Zwicker: Procedure for calculating loudness of temporally
% variable sounds. J. Acoust. Soc. Am. 62 (1977) 675-682.
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Altered by Matt Flax is flatmax for the Psy-Sound project
% Jan. 2007

if nargin < 1
  fs = 44100;
end  

fn = fs/2;

Wp = 95/fn; % lower cut-off frequency of pass band at 95 Hz
Ws = 45/fn; % upper cut-off frequency of stop band at 45 Hz
Rp = 0.5;   % Passband-Ripple
Rs = 10;    % Stopband attenuation

[n, Wn] = buttord(Wp, Ws, Rp, Rs);   
[b, a ] = butter(n, Wn, 'high');				

% end butter_hp