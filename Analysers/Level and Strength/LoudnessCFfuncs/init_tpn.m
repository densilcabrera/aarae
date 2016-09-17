function [B] = init_tpn(tau1, tau2, E, f_abt)
% [B] = init_tpn(tau1, tau2, E, f_abt)
% initialisation of nonlinear low-pass filter for modeling 
% forward masking effects on loudness
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% calculation of tau12, taking into account loudness
tau12=0.007.*exp(-E)+0.011;
% constants
wrz = sqrt(((tau12+tau2)^2-4*tau1*tau2)/(2*tau1*tau2)^2);
l1 = -(tau12+tau2)/(2*tau1*tau2)+wrz;
l2 = -(tau12+tau2)/(2*tau1*tau2)-wrz;
nenner = tau2*(l1-l2);
t2l1 = tau2*l1+1;
t2l2 = tau2*l2+1;

B=zeros(7,1);
B(1) = (exp(l1/f_abt)-exp(l2/f_abt))/nenner;			%K22ZA
B(2) = (t2l2*exp(l1/f_abt)-t2l1*exp(l2/f_abt))/nenner;		%K22ZZ
B(3) = (t2l1*exp(l1/f_abt)-t2l2*exp(l2/f_abt))/nenner; 		%K22AA
B(4) = (t2l1*t2l2*(exp(l1/f_abt)-exp(l2/f_abt)))/nenner;	%K22AZ
B(5) = exp(-1/(tau12*f_abt));					%=K21AA
B(6) = exp(-1/(tau2*f_abt));					%=K1ZZ
B(7) = 1-exp(-1/(tau2*f_abt));  				%=K1ZE

