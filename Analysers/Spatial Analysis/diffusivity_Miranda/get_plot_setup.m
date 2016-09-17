function plot_setup = get_plot_setup(varargin)

% Default parameters
plot_setup.type = 'plot_setup';
plot_setup.diffusivity = 1;
plot_setup.diffusivityplot_length = 1; % in seconds
plot_setup.st = 1;
plot_setup.rt = 1;
plot_setup.cummulative = 1;
plot_setup.cumm_anLength = 10;
plot_setup.cumm_edgebands = [250 4000];
plot_setup.cumm_plotspan = [20 100];
plot_setup.pancake = 1;
plot_setup.pan_anLength = 1;
plot_setup.pan_edgebands = [250 4000];
plot_setup.pan_plotspan = [20 100];
plot_setup.range = -20;


for I = 1:2:length(varargin)-1;
    switch varargin{I}
        
        case 'diffusivity'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                plot_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end
            
        case 'st'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                plot_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end
            
        case 'rt'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                plot_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end
            
        case 'cummulative'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                plot_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end
            
        case 'cumm_anLength'
            plot_setup.cumm_anLength = varargin{I+1};
            
        case 'cumm_edgebands'
            plot_setup.cumm_edgebands = varargin{I+1};
            
        case 'pancake'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                plot_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end

            
    end
end