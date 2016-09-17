function [ns,lautheit]=flankenlautheit24(kernlautheit)
% calculates specific loudness pattern and loudness for all time frames 
% calls flankenlautheit_t24.m to calculate specific loundess for single
% time frames
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007



ns=zeros(length(kernlautheit(:,1)),240);
for i=1:length(kernlautheit(:,1))
   [ns(i,:), lautheit(i,1)]=flankenlautheit_t24(kernlautheit(i,:));
end   