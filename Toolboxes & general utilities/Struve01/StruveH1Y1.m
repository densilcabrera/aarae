function fun=StruveH1Y1(z)
%
% StruveH1Y1 calculates the function StruveH1-BesselY1 for complex argument z
%
% Author : T.P. Theodoulidis
% Date   : 11 June 2012
%
% Arguments
% z : can be scalar, vector, matrix
% 
% External routines called         : cheval, StruveH1
% Matlab intrinsic routines called : bessely, besselh
%
nom=[4,0,9648             ,0,8187030           ,...
       0,3120922350       ,0,568839210030      ,... 
       0,49108208584050   ,0,1884052853216100  ,...
       0,28131914180758500,0,126232526316723750,...
       0,97007862050064000,0,2246438344775625];
%
den=[4,0,9660              ,0,8215830           ,...
       0,3145141440        ,0,577919739600      ,...
       0,50712457149900    ,0,2014411492343250  ,...
       0,32559467386446000 ,0,177511711616489250,...
       0,230107774317671250,0,31378332861500625];
%
x=z(:);
%
% |x|<=16
i1=abs(x)<=16;
x1=x(i1);
if isempty(x1)==0
    fun1=StruveH1(x1)-bessely(1,x1);
else
    fun1=[];
end
%
% |x|>16 and real(x)<0 and imag(x)<0
i2=(abs(x)>16 & real(x)<0 & imag(x)<0);
x2=x(i2);
if isempty(x2)==0
    x2=-x2;
    fun2=2/pi+2/pi./x2.^2.*polyval(nom,x2)./polyval(den,x2)-2i*besselh(1,1,x2);
else
    fun2=[];
end
% |x|>16 and real(x)<0 and imag(x)>=0
i3=(abs(x)>16 & real(x)<0 & imag(x)>=0);
x3=x(i3);
if isempty(x3)==0
    x3=-x3;
    fun3=2/pi+2/pi./x3.^2.*polyval(nom,x3)./polyval(den,x3)+2i*besselh(1,2,x3);    
else
    fun3=[];
end
% |x|>16 and real(x)>=0
i4=(abs(x)>16 & real(x)>=0);
x4=x(i4);
if isempty(x4)==0
    fun4=2/pi+2/pi./x4.^2.*polyval(nom,x4)./polyval(den,x4);
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