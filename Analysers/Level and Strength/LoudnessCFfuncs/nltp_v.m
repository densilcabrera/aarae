function [Ae,Ze,B,lastcase] = nltp_V(Aa,Za,E,B,lastcase,f_abt)
% applies nonlinear low pass to model forward masking
%	Aa	old output signal
%	Za	old state of accumulator
%	E	input signal
%	B	matrix with 7 constants (initialized by init_tpn.m)
%	Ae	new output signal
%	Ze	new state of accumulator
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007

if E < Aa
   % case 1: "discharge" i.e. forward masking: time constant stays
   % the same
  lastcase=1;
  if Aa > Za
     % case 1.1: select short time constant and charge accumulator
     % (tau_var)
    Ze = Aa * B(1) - Za * B(2);
    Ae = Aa * B(3) - Za * B(4);
    if E > Ae
      Ae = E;
    end
    if Ze > Ae
      Ze = Ae;
    end
  else
     % case 1.2: long time constant
    Ae = Aa * B(5);
    Ze = Ae;
    if E > Ae
      Ae = E;
      Ze = E;
    end
  end
elseif E == Aa
   
   if lastcase ==1 
      B = init_tpn(0.005, 0.075, E, f_abt);
      lastcase=2;
   end 
   
  Ae = E;
  % case 2: charge accumulator (no forward masking)
  if Aa > Za
     % case 2.1:
     Ze = Za * B(6) + E * B(7);
  else
     % case 2.2:
      Ze = E;
  end
else
   % case 3:
  lastcase=3;
  B = init_tpn(0.005, 0.075, E, f_abt);
  Ae = E;
  Ze = Za * B(6) + E * B(7);
end