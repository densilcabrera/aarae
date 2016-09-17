function [y, returnData] = rms_tep(fout, f_abt, fs)
% y = rms_tep(fout,f_abt,fs);
% Calculation of short term RMS levels with auditory temporal
% window (Plack&Moore)
% erd=4ms, w=-51 dB
% f_abt (envelope), fs(time signal)
% default f_abt=500ms and fs=44100
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Altered by flatmax is Matt Flax for the Psy-Sound project
% Jan. 2007

if nargin <2
   f_abt=500;
   fs=44100;
end

if size(fout,2) < size(fout,1)
	fout = fout';
end

[t_pa,w,t_sb,t_sa,t_pb] = staticParamDLM;
[h, t, erd] = tep_window(t_pb, t_pa, t_sb, t_sa, w, fs);
h = fliplr(h.^2)'; %due to convolution and intensity; 

dauer   = erd*fs;
step    = round(1/f_abt*fs);
wlen    = length(h);
n_steps = floor((length(fout)-wlen)/step)+1;

rms   = [];
power = fout'.^2;
for i=0:(n_steps-1)
  rms(i+1) = sum(power(i*step+1:i*step+wlen) .* h)/dauer;
end
returnData = fout(n_steps*step:end);

y=sqrt(rms);
j=find(y==0);  
y(j) = realmin;
