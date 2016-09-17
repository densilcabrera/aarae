function kernlaut_dyn = post_maskn(kernlaut_stat,f_abt)
% kernlaut_dyn = post_maskn(kernlaut_stat,f_abt);
% kernlaut_stat is the input matrix (stationary main loudness; (N,24)-matrix),
% kernlaut_dyn is the output matrix (dynamic main loudness; (N,24)-matrix),
% f_abt (optional): sampling rate (default=500)
%   Ae	output signal of a critical band
%	Ze	state of accumulator
% Parameters of nonlinear low-pass: 
%	tau1  =  5 ms	(Zwicker:3.5 ms)
%	tau2  = 75 ms (Zwicker: 70 ms)
%	tau12 = 11..17 ms (80..40 dB)
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Initialisation
if nargin<2
   f_abt = 500; 
end   
N = length(kernlaut_stat(:,1));
kernlaut_dyn = zeros(N,size(kernlaut_stat,2));  
Ze = zeros(N,1);
if size(kernlaut_stat,1)<size(kernlaut_stat,2)
    error('not enough iterations of algorithm to calculate post_maskn');
end
for i=1:size(kernlaut_stat,2)   
   B = init_tpn(0.005,0.075,kernlaut_stat(i,1),f_abt); 
   lastcase=3; 
   Aold = 0;
   Zold = 0;
   for k=1:N 
      [kernlaut_dyn(k,i),Ze(k),B,lastcase] = nltp_v(Aold, Zold, kernlaut_stat(k,i), B,lastcase,f_abt);
      Aold = kernlaut_dyn(k,i);
      Zold = Ze(k);
   end
end   
