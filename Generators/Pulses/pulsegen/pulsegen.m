function p=pulsegen(fs,T,edge,type,varargin)
%p=pulsegen(fs,T,edge,type,f,opt);
%a signal generation program
%fs is the sampling frequency
%T is the total signal length
%edge is a decay parameter for some waveforms
%   it is used in 'gaussian', 'monocycle',  'biexponential', 'mexican hat', 'sinc', 'double sinc', 'sinc squared'
%   and windowed sweep
%   it is mostly a parameter to describe how much the edge of the pulse is decayed.
%type is the type of the waveform desired
%   allowable types are 'gaussian', 'square', 'triangle', 'monocycle',
%   'biexponential', 'mexican hat', 'sinc', 'raised cosine','double sinc', 'sinc squared','sweep', and 'windowed sweep'
%f is the modulation frequency if left out it is assumed 0.
%opt is an optional argument for pulse waveforms requiring a lower and higher frequency
%    it is used in 'double sinc' ,'sweep' and 'windowed sweep' for the low and high frequency.
%  for the raised cosine function it is used as the damping factor between
%  0 and 1
% the pulses are always normalized to a peak amplitude of 1

if nargin<4 %test for optional arguments
    error('not enough input arguments');
end

p=inputParser();
p.addParamValue('start_frequency',-1,@(x) isreal(x) &&(x>=0));
p.addParamValue('stop_frequency',-1,@(x) isreal(x) &&(x>=0));
p.addParamValue('modulation',0,@isnumeric);
p.addParamValue('window',0.1,@(x) (x>=0)&&(x<=1));
p.addParamValue('dispersion',0,@isnumeric);
p.addParamValue('high_pass',0,@(x) (x>=0) && (x<=1));
p.addParamValue('low_pass',1,@(x) (x>=0) && (x<=1));
p.addParamValue('args',[]);

p.parse(varargin{:});
params=p.Results;
if (params.start_frequency==-1)
    params.start_frequency=16*edge/(5*T);
end

if (params.stop_frequency==-1)
    params.stop_frequency=64*edge/(5*T);
end


t=-T/2:1/fs:T/2;

sig=(T/8/edge)^2;

switch type
case {'gaussian','guassian'}  %generate a guassian pulse
    y=exp(-(t).^2/sig);
case {'square'}  %generate a square pulse
    y=(t>-T/edge/2)&(t<T/edge/2);
case {'triangle'} %generate a triangle pulse
    y=(t+T/edge/2).*(t<0)-(t-T/edge/2).*(t>=0);
    y(y<0)=0;
case {'monocycle'} %generate a gaussian monocycle
    if (sig==0)
        y=t;
    else
        y=2*t./sig.*exp(-(t).^2/sig); 
    end
case {'exponential'} %generate a exponential pulse
    y=exp(-t*8*edge/T);
    y(t<0)=0;
case {'biexponential'} %generate a biexponential pulse
    y=exp(-abs(t)*8*edge/T);
case {'mexican hat'}  %generate a gaussian second deriviative
    z=t./sqrt(0.75*sig);
    y=sqrt(1/2*pi).*(1-z.^2).*exp(-z.^2/2);
case {'sinc'} %generate a sinc function
    y=sinc(2*pi*edge*50.*t/(5*T));
case {'raised cosine'}
    rb=2*edge*50/(5*T);
    beta=rb/2*params.args(1);
    y=sinc(pi.*t.*rb).*(cos(2*pi*beta.*t)./(1-(4*beta.*t).^2));
case {'double sinc'} %generate a bandlimited function from two sinc functions
    f1=params.start_frequency;
    f2=params.stop_frequency;
    y=f2*sinc(2*f2*pi.*t)-f1*sinc(2*f1*pi.*t);
case {'sinc squared'} %generates sinc squared function
    y=sinc(2*pi*edge*16.*t/(5*T)).^2;
case {'sweep'} %generate frequency sweep
    f1=params.start_frequency;
    f2=params.stop_frequency;
    theta=f1.*(t+T/2)+((f2-f1)/(T)).*(t+T/2).^2;
    y=real(exp(1j*(2*pi.*theta-pi/2)));
otherwise
    error('invalid pulse type');
end


if (params.low_pass<1)||(params.high_pass>0)||(params.dispersion~=0)
    c=length(y);
s=fft(y);
sA=abs(s);
sP=angle(s);
if (params.low_pass<1)
    cP=ceil(params.low_pass*c/2);
    if (cP==0)
        sA(:)=0;
    else
    	sA(cP:end-cP+2)=0;
    end
end
if (params.high_pass>0)
    cP=floor(params.high_pass*c/2);
    if (cP~=0)
    	sA(1:cP)=0;
        sA(end-cP+2:end)=0;
    end
end
if (params.dispersion~=0)
    pp=params.dispersion.*linspace(0,2*pi,c);
    sP=sP+pp;
end
s2=sA.*cos(sP)+1j*sA.*sin(sP);
y=real(ifft(s2));

end

if (params.window>0)
    c=length(y);
    w=ones(size(y));
    p1=floor(c*params.window/2); 
    % Window is defined in three sections: taper, constant, taper
    w(1:p1+1) = (-cos((0:p1)/p1*pi)+1)/2;  
    w(end-p1:end)=(cos((0:p1)/p1*pi)+1)/2; 
    y=w.*y;
end
%apply a modulation
if params.modulation~=0
    y=y.*cos(pi*t*params.modulation);
end

%normalize the peak of the pulse to 1
p=y./max(abs(y));

end

function y=sinc(x)
y=sin(x)./x;
y(x==0)=1;
end