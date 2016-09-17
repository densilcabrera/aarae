function fgp_d=damp_a0(fgrp, HL_ihc)
% fgp_d=damp_a0(fgrp, HL_ihc)
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007


a0 = [ 0 0 0 0 0 0 0 0 0 0 -.2 -.5 -1.2 -2.1 -3.2 -4.6 -5.5 -5.6 -4.3 -2.5 -0.1 2.8 6.4 20.0];
fgp_d=zeros(length(fgrp(:,1)),24);  
for i = 1:length(fgrp(:,1))
   fgp_d(i,:)=fgrp(i,:)-a0-HL_ihc;
end   