function OUT = forest_reverb_Wiens2008
if nargin == 0
    param = inputdlg({'Number of trees';...
                       'Maximum number of scatterings';...
                       'Forest width (m)';...
                       'Minimum tree radius (m)';...
                       'Maximum tree radius (m)';...
                       'Audio sampling rate (Hz)'},...
                       'Settings',1,...
                       {'10';'5';'100';'0.1';'0.4';'44100'});
    param = str2num(char(param));
    if length(param) < 6, param = []; end
    if ~isempty(param)
        N_trees = param(1);
        N_sc_max = param(2);
        forestwidth = param(3);
        minradius = param(4);
        maxradius = param(5);
        fs = param(6);
    end
end

if ~isempty(param) || nargin ~= 0
    OUT.audio = forest_reverb(N_trees,N_sc_max,forestwidth,minradius,maxradius,fs);
    OUT.fs = fs;
    % Note that the following callback is for the main function rather than
    % for this calling function. This callback is mainly done for logging.
    OUT.funcallback.name = 'forest_reverb.m';
    OUT.funcallback.inarg = {N_trees,N_sc_max,forestwidth,minradius,maxradius,fs};
else
    OUT = [];
end