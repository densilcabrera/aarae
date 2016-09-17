function diffusivity_setup = get_diffusivity_setup(varargin)

% Default parameters
diffusivity_setup.type = 'diffusivity_setup';
diffusivity_setup.calculation = 1; % on or off
diffusivity_setup.hoaOrder = 5;
diffusivity_setup.bands_per_octave = 1; % Bands per octave
diffusivity_setup.center_bands = [250 500 1000 2000];
diffusivity_setup.wideband_edgebands = [250 2000];
diffusivity_setup.method = 'Gover';
diffusivity_setup.analysis_length = 10; % in ms
diffusivity_setup.overlap = 50; % in percentage

for I = 1:2:length(varargin)-1;
    switch varargin{I}
        
        case 'calculation'
            if varargin{I+1} == 1 || varargin{I+1} == 0
                diffusivity_setup.calculation = varargin{I+1};
            else
                error('Calculation on or off represented by 1 or 0')
            end
            
        case 'hoaOrder'
            diffusivity_setup.hoaOrder = varargin{I+1};
            
        case 'bands_per_octave'
            diffusivity_setup.bands_per_octave = varargin{I+1};
            
        case 'center_bands'
            diffusivity_setup.center_bands = varargin{I+1};
            
        case 'wideband_edgebands'
            diffusivity_setup.wideband_edgebands = varargin{I+1};
            
        case 'method'
            diffusivity_setup.method = varargin{I+1};
            
        case 'analysis_length'
            diffusivity_setup.analysis_length = varargin{I+1};
            
        case 'overlap'
            diffusivity_setup.overlap = varargin{I+1};
            
    end
end