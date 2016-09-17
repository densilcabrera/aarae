function [ns_t,lautheit_t]=flankenlautheit_t24(kernlautheit_t)
% [ns_t,lautheit_t]=flankenlautheit_t24(kernlautheit_t)
% Calculates loudness based on DIN 45631 / ISO 532 B (Zwicker)
% kernlautheit_t is a row vector (1x24) representing the main loudnesses
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007


% Definition of constants and variables

% Upper limits of approximated critical bands in terms of critical 
% band rate

zup=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];

% Range of specific loudness for the determination of the steepness of 
% the upper slopes in the specific loudness critical
% band rate pattern

rns=[21.5,18.0,15.1,11.5,9.0,6.1,4.4,3.1,2.13,1.36,0.82,0.42,0.30,0.22,0.15,0.10,0.035,0.0];

% Steepness of the upper slopes in the specific loudness critical 
% band rate pattern for the ranges RNS as a function of the number of 
% the critical band

usl=[ 13.00,8.20,6.30,5.50,5.50,5.50,5.50,5.50
       9.00,7.50,6.00,5.10,4.50,4.50,4.50,4.50
       7.80,6.70,5.60,4.90,4.40,3.90,3.90,3.90
       6.20,5.40,4.60,4.00,3.50,3.20,3.20,3.20
       4.50,3.80,3.60,3.20,2.90,2.70,2.70,2.70
       3.70,3.00,2.80,2.35,2.20,2.20,2.20,2.20
       2.90,2.30,2.10,1.90,1.80,1.70,1.70,1.70
       2.40,1.70,1.50,1.35,1.30,1.30,1.30,1.30
       1.95,1.45,1.30,1.15,1.10,1.10,1.10,1.10
       1.50,1.20,0.94,0.86,0.82,0.82,0.82,0.82
       0.72,0.67,0.64,0.63,0.62,0.62,0.62,0.62
       0.59,0.53,0.51,0.50,0.42,0.42,0.42,0.42
       0.40,0.33,0.26,0.24,0.22,0.22,0.22,0.22
       0.27,0.21,0.20,0.18,0.17,0.17,0.17,0.17
       0.16,0.15,0.14,0.12,0.11,0.11,0.11,0.11
       0.12,0.11,0.10,0.08,0.08,0.08,0.08,0.08
       0.09,0.08,0.07,0.06,0.06,0.06,0.06,0.05
       0.06,0.05,0.03,0.02,0.02,0.02,0.02,0.02 ];
    
    
% Specific loudness critical band rate pattern
    
ns_t=zeros(1,240);   

% Initialisation

lautheit_t=0;
n1=0;
z=0.1;       
z1=0;        
j=0;
i_z=0;

 
for i = 0: 23 
   KernL = kernlautheit_t(i+1);
   xzup = zup(i+1);
   
   xzup = xzup + 0.0001;
   ig = i - 1;
   if ig > 7
      ig = 7;
   end   
   
   while z1 < xzup
      if n1 > KernL           
         xusl = usl(j+1,ig+1);
         
         
         n2 = rns(j+1);
         if n2 < KernL
            n2 = KernL;
         end   
         dz = (n1 - n2) / xusl;
         z2 = z1 + dz;
         
         if z2 > xzup
            z2 = xzup;
            dz = z2 - z1;
            n2 = n1 - dz * xusl;
         end
         

         
         lautheit_t = lautheit_t + (dz * (n1 + n2) / 2.);
         
         
         while k <= z2
            ns_t(i_z+1) = n1 - (k - z1) * xusl;
            i_z = i_z + 1;
            k = k + 0.1;
         end
         
         z = k;
      else 
         if n1 < KernL
           
            j = 0;
            while (j < 18) & (rns(j+1) >= KernL)
               j=j+1;
            end   
         end
         
         
         z2 = xzup;
         n2 = KernL;
         lautheit_t = lautheit_t + (n2 * (z2 - z1));
         k = z;
         while k <= z2 
            ns_t(i_z+1) = n2;
            i_z = i_z + 1;
            k = k + 0.1;
         end
         z = k;
      end
      
 
      while (n2 <= rns(j+1)) & (j < 17)
         j=j+1;
      end   
      if j > 17
         j = 17;
      end
      z1 = z2;
      n1 = n2;
      
   end % Ende while z1 < xzup
end    % Ende for i = 0: 23


if lautheit_t < 0
   lautheit_t = 0;
end