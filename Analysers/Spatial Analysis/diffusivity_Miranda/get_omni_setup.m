function omni_setup = get_omni_setup(varargin)

% Default parameters
omni_setup.type = 'omni_setup';
omni_setup.calculation = 1; % on or off
omni_setup.hoaOrder = 0;
omni_setup.bands_per_octave = 1; % Bands per octave
omni_setup.center_bands = [250 500 1000 2000];
omni_setup.wideband_edgebands = [250 2000];
omni_setup.analysis_length = 10; % in ms
omni_setup.overlap = 50; % in percentage
omni_setup.normalize = 1;

for I = 1:2:length(varargin)-1;
    switch varargin{I}
        
        case 'calculation'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                omni_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end
            
        case 'hoaOrder'
            omni_setup.hoaOrder = varargin{I+1};
            
        case 'bands_per_octave'
            omni_setup.bands_per_octave = varargin{I+1};
            
        case 'center_bands'
            omni_setup.center_bands = varargin{I+1};
            
        case 'wideband_edgebands'
            omni_setup.wideband_edgebands = varargin{I+1};
            
        case 'analysis_length'
            omni_setup.analysis_length = varargin{I+1};
            
        case 'overlap'
            omni_setup.overlap = varargin{I+1};
            
        case 'normalize'
            omni_setup.normalize = varargin{I+1};
            
    end
end