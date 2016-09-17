function lf=fluct(main_N)
% lf=fluct(main_N);
% calculates loudness fluctuation lf for normal and hearing-impaired listeners
% References:
% Chalupper, J. (2001) - in german - : Perzeptive Folgen von
% Innenohrschwerhörigkeit: 
% Modellierung, Simulation und Rehabilitation. Dissertation at the Technical
% University of Munich, Shaker Verlag.
% Chalupper, J.(2000): Modellierung der Lautstärkeschwankung für
% Normal- und Schwerhörige.  
% Tagungsband DAGA 2000, 26. Jahrestagung Akustik, Oldenburg,
% 20.-24.3.2000, S. 254-255. 
% Chalupper, J. (2007): Modeling loudness fluctuation for norm and
% hearing-impaired listeners. Proceedings of EFAS 2007,
% Heidelberg. in preparation 
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.2.2007
% required functions: korrel.m, kernl2lg.m


HL_ohc=zeros(1,24);
HL_ihc=zeros(1,24);

y_k=korrel(main_N); %cross-channel correlation
nmax=prctile(main_N,95);
nmin=prctile(main_N,5);

%maximale und minimale LE's berechnen (aus 5% und 95% Lautheits-Perzentil)
le_max=kernl2lg(nmax,HL_ohc,HL_ihc); %inverse loundess transformation for normal hearing
le_min=kernl2lg(nmin,HL_ohc,HL_ihc);
le_diff=le_max-le_min;
le_diff(find(le_diff>30))=30;  %limit to 30 dB

%account for channel correlation 
le_diff=y_k.*le_diff;

dlksum=sum(le_diff);

%transfor into categorical units
% How strong (not how fast) is loudness fluctuating?
% 0=not at all, 1=very weak, 2=weak, 3=medium, 4=strong, 5=very strong,
% 6=extremely strong
a=-0.19151712;
b=0.199654111;

lf=a+b.*(dlksum.^0.5);

lf(find(lf<0))=0;
lf(find(lf>6))=6;


