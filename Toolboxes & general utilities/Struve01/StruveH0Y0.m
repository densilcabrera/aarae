function fun=StruveH0Y0(z)
%
% StruveH0Y0 calculates the function StruveH0-BesselY0 for complex argument z
%
% Author : T.P. Theodoulidis
% Date   : 11 June 2012
%
% Arguments
% z : can be scalar, vector, matrix
% 
% External routines called         : cheval, StruveH0
% Matlab intrinsic routines called : bessely, besselh
%
nom=[4,0,8816             ,0,6778206          ,...
       0,2317961250       ,0,374638650750     ,...
       0,28306147182390   ,0,937687098467700  ,...
       0,11970864124436700,0,46174193179143750,...
       0,32071055725170000,0,840808761085125];
%
den=[4,0,8820             ,0,6786990          ,...
       0,2324669760       ,0,376904178000     ,...
       0,28663562736900   ,0,963414191990250  ,...
       0,12740661151218000,0,54025303535453250,...
       0,50023429199493750,0,4092826025413125];
%
x=z(:);
%
% |x|<=16
i1=abs(x)<=16;
x1=x(i1);
if isempty(x1)==0
    fun1=StruveH0(x1)-bessely(0,x1);
else
    fun1=[];
end
%
% |x|>16 and real(x)<0 and imag(x)<0
i2=(abs(x)>16 & real(x)<0 & imag(x)<0);
x2=x(i2);
if isempty(x2)==0
    x2=-x2;
    fun2=-(2/pi./x2.*polyval(nom,x2)./polyval(den,x2))+2i*besselh(0,1,x2);
else
    fun2=[];
end
% |x|>16 and real(x)<0 and imag(x)>=0
i3=(abs(x)>16 & real(x)<0 & imag(x)>=0);
x3=x(i3);
if isempty(x3)==0
    x3=-x3;
    fun3=-(2/pi./x3.*polyval(nom,x3)./polyval(den,x3))-2i*besselh(0,2,x3);
else
    fun3=[];
end

% |x|>16 and real(x)>=0
i4=(abs(x)>16 & real(x)>=0);
x4=x(i4);
if isempty(x4)==0
    fun4=2/pi./x4.*polyval(nom,x4)./polyval(den,x4);
else
    fun4=[];
end
fun=x*0;
fun(i1)=fun1;
fun(i2)=fun2;
fun(i3)=fun3;
fun(i4)=fun4;
%
fun=reshape(fun,size(z));
%