function [b,a]=int_tp(f_abt);
% n_t=int_tp(n,f_abt);
% models temporal integration of loudness
% with simple 1st order low pass filter [butter 8Hz]
% filter designed to match experimental data on pure tone loudness integration and
% fluctuation strength. 
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Altered by MFFM : Matt Flax is flatmax for the Psy-Sound project
% Jan. 2007

%if nargin<2
%   f_abt = 500; 
%end   

[b,a]=butter(1,8/(f_abt/2));
%n_t=filter(b,a,n);
%y=find(n_t <0);
%n_t(y) =0;
