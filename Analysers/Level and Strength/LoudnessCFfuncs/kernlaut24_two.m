function out = kernlaut24_two(rms, HL_ohc)
% out=kernlaut24_two(rms, HL_ohc);
% calculates main loudness for main excitation (rms) according to DIN45631
% accounts for nonlinear component of hearing loss (HL_ohc)
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

% Data for thq and a0 from  ACUSTICA Vol.27 1972 issue 5 S.261 TableI
thq = [42 18.5 11.5 8.3 6.7 5.5 4.8 4.3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];

thq = thq+HL_ohc;

s = 10.^(0.22-0.005*[0.5:23.5])-1;
k = 0.23;
out = zeros(length(rms(:,1)),24);  

for i = 1:length(rms(:,1))
   le=rms(i,:);
   mp1= 0.04925*(1./s).^k.*10.^(0.1*k* thq);
   mp2=(1-s +s.*10.^(0.1*(le-thq))).^k -1;
   nm=mp1 .* mp2; 

   j=find(le <= thq | nm<0);   
   nm(j)=0;
   
   out(i,:)=nm;
 end
 