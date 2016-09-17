function micFreqRespHandle = MicFreqRespFunction(frq,dBs)
%
% micFreqRespHandle = MicFreqRespFunction(frq,lvlDbs)
%
% Creates a handle to a microphone frequency response (positive gain as a 
% function of the frequency). The gain is calculated by interpolating
% the gain values in dB in the frequency log. scale.
% 
% Inputs:
% - frq is a vector of frequency values
% - dBs is a vector of gain values (in dB)
%
% Remarks:
% - frq and dBs must have the same length
% - frequency values must be sorted in ascending order
% - if no value is given for frequency 0, the slope between 0 et the first
% frequency value is supposed to be the same as between the first and
% second value.
% - if no value is given for frequency Inf, the slope between the last
% frequency value and Inf is supposed to be the same as between the former
% last and the last value.
%
% N. Epain - 07/12/2009 

% A few tests
if ~isvector(frq) || ~isreal(frq) || any(frq)<0
    error('Frequencies must be passed in a vector of real positive values') ;
end
if ~isvector(dBs) || ~isreal(dBs)
    error('Gains must be passed in a vector of real values') ;
end
frq = frq(:).' ;
dBs = dBs(:).' ;
if size(frq)~=size(dBs)
    error('Frequency and gain vectors must have the same length.') ;
end    
if length(frq)<2
    error('You must pass at least two frequency values') ;
end
if any(sort(frq,1,'ascend')~=frq) || any(unique(frq)~=frq)
    error('Frequency values must be sorted in ascending order.') ;
end

% Adding 0 frequency if necessary
if frq(1)~=0
    coe = (dBs(2)-dBs(1))/(log(frq(2))-log(frq(1))) ;
    frq = [ eps , frq ] ;
    dBs = [ dBs(1)+coe*(log(frq(1))-log(frq(2))) , dBs ] ;
else
    frq(1) = eps ;
end

% Adding Inf frequency if necessary
if frq(end)~=Inf
    coe = (dBs(end)-dBs(end-1))/(log(frq(end))-log(frq(end-1))) ;
    frq = [ frq , 1/eps ] ;
    dBs = [ dBs , dBs(end)+coe*(log(frq(end))-log(frq(end-1)))] ;
else
    frq(end) = 1/eps ;
end

% Create the function handle
eval(['micFreqRespHandle = @(f) 10.^(interp1([' num2str(log(frq)) ...
    '],[' num2str(dBs) '],log(min(max(f,eps),1/eps)),''cubic'')/20) ;']) ;

