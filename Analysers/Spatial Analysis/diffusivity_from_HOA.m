function OUT = diffusivity_from_HOA(IN,fs,max_order,diffusivity_setup,omni_setup,plot_setup)

diffusivity_setup = get_diffusivity_setup;
omni_setup = get_omni_setup;
plot_setup = get_plot_setup;

if nargin < 3, max_order = 4; end

if isstruct(IN)
    hoaSignals = IN.audio;
    fs = IN.fs;
else
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2double(char(fs));
    end
    hoaSignals = IN;
end

hoaFmt = GenerateHoaFmt('res2d',max_order,'res3d',max_order);

diffusivity = difussivity_calc(hoaSignals, fs, hoaFmt, diffusivity_setup);

omniparameters = omni_acoustic_parameters(hoaSignals, fs, omni_setup);

plot_diff_omni(diffusivity, omniparameters, plot_setup,[],[]);

if isstruct(IN)
    OUT.funcallback.name = 'diffusivity_from_HOA.m';
    OUT.funcallback.inarg = {fs,max_order,diffusivity_setup,omni_setup,plot_setup};
else
    OUT = diffusivity;
end