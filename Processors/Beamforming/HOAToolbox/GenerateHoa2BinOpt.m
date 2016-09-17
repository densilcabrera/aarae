function hoa2BinOpt = GenerateHoa2BinOpt(varargin)

% Default parameters
hoa2BinOpt.method       = 'magOptim' ;
hoa2BinOpt.filterLength = 512 ;
hoa2BinOpt.equalisation = true ;
hoa2BinOpt.magOptimFreq = 'auto' ;

% Number of arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% Parameters values are checked and assigned to the structure if possible
for I = 1 : 2 : length(varargin)-1
    switch lower(varargin{I})
        
        case 'method'
            switch lower(varargin{I+1})
                case {'magoptim','leasterr','virtualspk'}
                    hoa2BinOpt.method = varargin{I+1} ;
                otherwise
                    error('Unknown rendering method.') ;
            end            

        case 'filterlength'
            if (numel(varargin{I+1})==1) && isreal(varargin{I+1}) ...
                    && (round(varargin{I+1})==(varargin{I+1}))
                hoa2BinOpt.filterLength = varargin{I+1}(:) ;
            else
                error('Filter length must be an integer number.') ;
            end
            
        case 'equalisation'
            if (numel(varargin{I+1})==1) && islogical(varargin{I+1})
                hoa2BinOpt.equalisation = varargin{I+1}(:) ;
            else
                error('Equalisation can only be true or false.') ;
            end
            
        case 'magoptimfreq'
            if (numel(varargin{I+1})==1) && isreal(varargin{I+1}) 
                hoa2BinOpt.magOptimFreq = varargin{I+1}(:) ;
            elseif strcmpi(varargin{I+1},'auto')
                hoa2BinOpt.magOptimFreq = varargin{I+1} ;
            else
                error(['The magnitude optimisation transition frequency ' ...
                    'must be either a real number or ''auto''.']) ;
            end
            
        otherwise
            error([varargin{I} ' is not a valid parameter.']) ;
            
    end
end

% Set the magnitude optimsation transition frequency to [] if the method is
% not 'magOptim'
if ~strcmpi(hoa2BinOpt.method,'magoptim')
    hoa2BinOpt.magOptimFreq = [] ;
end
