function [h,t,erd,h_orig]=tep_window(t_pb,t_pa,t_sb,t_sa,w,fs)
%[h,t,erd,h_orig]=tep_window(t_pb,t_pa,t_sb,t_sa,w);
%calculates TEP-window (h) 
%Plack & Moore: Temporal window shape as afunction of frequency and level
%t: time vector; erd: equivalent rectangular duration; h_orig: impulse response
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Altered by Matt Flax is flatmax for the Psy-Sound project
% Jan. 2007

if nargin < 5	%from Plack and Moore
	t_pb=5.4e-3;
	t_pa=2.6e-3;
	t_sb=30e-3;
	t_sa=14e-3;
	w=-39;
end

w=10^(w/10); %transform to intensity
t=[0:1/fs:2*(t_sa+t_sb)]; 
intensity_a=(1-w)*(1+(2*t)/t_pa).*exp(-(2*t)/t_pa) + w*(1+(2*t)/t_sa).*exp(-(2*t)/t_sa);
intensity_b=(1-w)*(1+(2*t)/t_pb).*exp(-(2*t)/t_pb) + w*(1+(2*t)/t_sb).*exp(-(2*t)/t_sb);
intensity=[fliplr(intensity_a) intensity_b(2:end)];
intensity=intensity(find(intensity>=1e-6)); %cut at -60dB 
h=sqrt(intensity); %amplitude window
h_orig=h/sum(h);

t=[0:1/fs:(length(h)-1)/fs];
%figure,plot(t,20*log10(h))
erd=sum(h.^2)/fs;
