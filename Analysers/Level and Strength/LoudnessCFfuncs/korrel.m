function y_k=korrel(kl)
% y_k=korrel(kl);
% calculates cross-correlation coefficients of kl
% y_k: squared product of cross-correlation coefficients to adjacent channels
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.2.2007

% Divide by zero warning as par of corrcoed are unavoidable - so do
% not display them
divideByZeroWarn = warning('off', 'MATLAB:log:divideByZero');

t=size(kl);%size of for loop
kl(find(kl==0))=realmin;

kk=corrcoef(kl);%calculate cross-correlation coefficients

nd=isnan(kk);%  set undefined cross-correlation coefficients
ind=find(nd==1);%            to
kk(ind)=realmin;%        realmin

nd=isinf(kk);%  set undefined cross-correlation coefficients
ind=find(nd==1);%            to
kk(ind)=realmin;%        realmin


for i=1:(t(2)-1) %t(2): channel
   kkup(i)=kk(i,(i+1));%cross-correlation coefficient to next higher channel
   kkdown(i+1)=kkup(i);%cross-correlation coefficient to next lower channel
end

kkup(t(2))=kkup(t(2)-1); %cross-correlation coefficients of highest
                         %and lowest channel
kkdown(1)=kkdown(2);     %are set to the value of the afjacent channel

y_k=(kkup.*kkdown).^2; %according to roughness model of Daniel

% Restore warning state
warning(divideByZeroWarn);

% EOF
