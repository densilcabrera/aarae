function OUT = units_aarae(IN,units,units_ref,units_type)
% This function edits or creates the following properties fields:
%
% units - string (e.g. 'Pa', 'm/s', 'V')
%
% units_ref - reference value, a number (e.g. 2e-5 for Pa)
%
% units_type - use 1 for amplitude units (e.g. Pa, V), use 2 for intensity
% and power units (e.g. W)


if nargin >= 3
    OUT = IN;
    OUT.properties.units = units;
    OUT.properties.units_ref = units_ref;
    if exist('unittype','var')
        OUT.properties.units_type = units_type;
    else
        if regexp(units,'W')
            OUT.properties.units_type = 2;
        else
            OUT.properties.units_type = 1;
        end
    end
    return
end

units1 = 'Pa';
units_ref1 = 2e-5; 
units_type1 = 1;
if isstruct(IN)
    if isfield(IN,'properties')
        if isfield(IN.properties,'units')
            units1 = IN.properties.units;
        end
        if isfield(IN.properties,'units_ref')
            units_ref1 = IN.properties.units_ref;
        end
        if isfield(IN.properties,'units_type')
            units_type1 = IN.properties.units_type;
        end
    end
    
    param = inputdlg({'Units (e.g. Pa, V, m/s, W/m^2, W)';... 
        'Reference value (e.g. 2e-5 for sound pressure)';...
        'Unit type: use 1 for amplitude units (Pa, V, m/s); use 2 for power units (W/m^2, W)'},...
        'Units',... % This is the dialog window title.
        [1 60],... % height & width
        {units1;num2str(units_ref1);num2str(units_type1)}); % default values
    if isempty(param)
        OUT = [];
        return
    else
        OUT = IN;
        OUT.properties.units = char(param{1});
        OUT.properties.units_ref = str2double(char(param{2}));
        OUT.properties.units_type = str2double(char(param{3}));
        OUT.funcallback.name = 'units_aarae.m';
        OUT.funcallback.inarg = {OUT.properties.units,OUT.properties.units_ref,OUT.properties.units_type};
    end
else
    OUT = [];
end