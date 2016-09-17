function impedanceHandle = ImpedanceFunction(type,varargin)
%
% impedanceHandle = ImpedanceFunction(type,typeArguments)
%
% Creates a specific surface acoustic impedance function handle (SSAIFH).
%
% Three types of impedances are available:
%
%  * 'constant' -> Constant impedance function.
%    impedanceHandle = ImpedanceFunction('constant',impedanceValue)
%    Remark: impedanceValue must be scalar.
%
%  * 'delanyBazley' -> Delany-Bazley model for porous materials
%    impedanceHandle = ImpedanceFunction('delanyBazley',airFlowResistivity)
%    Remark: airFlowResistivity must be scalar.
%
%  * 'interp' -> cubic interpolation of frequency/impedance data
%    impedanceHandle = ImpedanceFunction('interp',freqValues,impeValues)
%    Remarks: - freqValues and impeValues must be vectors having the same
%               length. 
%             - Uses MATLAB interp1 function.
%
% N. Epain - 21/12/2007 


switch type
    
    case 'constant'
        if length(varargin)==1
            impe = varargin{1} ; 
            if isscalar(impe)
                eval(['impedanceHandle = @(f) ' num2str(impe) ...
                      '*ones(size(f)) ;']) ;
            else
                error('Impedance value must be scalar') ;
            end
        else
            error('Wrong number of arguments.') ;
        end
        
    case 'delanyBazley'
        if length(varargin)==1
            flowResist = varargin{1} ; 
            if isscalar(flowResist)
                eval(['impedanceHandle = @(f) ' ...
                      '1 + 0.0571*(1.2*f/' num2str(flowResist) ').^(-.754)'...
                      ' + j*0.087*(1.2*f/' num2str(flowResist) ').^(-.732)'...
                      ';']) ;
            else
                error('Static air flow resistivity value must be scalar') ;
            end
        else
            error('Wrong number of arguments.') ;
        end
        
    case 'interp'
        if length(varargin)==2
            freq = varargin{1}(:).' ;
            impe = varargin{2}(:).' ;
            if isvector(freq) & isvector(impe) & size(freq)==size(impe)
                eval(['impedanceHandle = @(f) interp1([' num2str(freq) ...
                                 '],[' num2str(impe) '],f,''cubic'') ;']) ;
            else
                error(['Frequency and impedance values must be passed ' ...
                       'in vectors having the same length.']) ;
            end
        else
            error('Wrong number of arguments.') ;
        end
        
    otherwise
        error(['Unknown argument, type "help ImpedanceFunction" ' ...
               'for options']) ;
           
end
