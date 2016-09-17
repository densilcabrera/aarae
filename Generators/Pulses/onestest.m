function OUT = onestest(dur,fs)
% This function generates a vector consisting of ones.
% The most likely use for this is for step function measurement (or
% simulation by convolution with an impulse response)

if nargin == 0
    
    param = inputdlg({'duration';... 
        'fs'},...
        'ones test settings',... % This is the dialog window title.
        [1 60],... 
        {'1';'48000'}); % default settings
    
    param = str2num(char(param)); 
    
    if length(param) < 2, param = []; end 
    if ~isempty(param) 
        dur = param(1);
        fs = param(2);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end

OUT.audio = ones(round(dur*fs),1);
% OUT.audio2 = OUT.audio;
OUT.fs = fs;
OUT.tag = 'ones';
OUT.funcallback.name = 'onestest.m';
OUT.funcallback.inarg = {dur,fs};