function hoa2SpkOpt = GenerateHoa2SpkOpt(varargin)

% Note : the default transition frequency should be the kr=n one

% Default parameters
hoa2SpkOpt.decodType    = 'mixed' ;
hoa2SpkOpt.sampFreq     = 48000 ;
hoa2SpkOpt.filterLength = 256 ;
hoa2SpkOpt.transFreq    = 'kr' ;
hoa2SpkOpt.transWidth   = 200 ; 
hoa2SpkOpt.spkDistCor   = true ; 

% Number of arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% Parameters values are checked and assigned to the structure if possible
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        
        case 'decodType'
            switch varargin{I+1}
                case {'mixed','basic','maxrE','inPhase'}
                    hoa2SpkOpt.decodType = varargin{I+1} ;
                otherwise
                    error('Unknown type of decoding') ;
            end            

        case 'sampFreq'
            if isvector(varargin{I+1}) && isreal(varargin{I+1})
                hoa2SpkOpt.sampFreq = varargin{I+1}(:) ;
            else
                error(['Required sampling frequencies must be passed' ...
                       'in a vector of real values']) ;
            end
            
        case 'filterLength'
            if isvector(varargin{I+1}) && isreal(varargin{I+1})
                hoa2SpkOpt.filterLength = varargin{I+1}(:) ;
            else
                error(['Filter lengths must be passed' ...
                       'in a vector of real values']) ;
            end
            
        case 'transFreq'
            if ( isscalar(varargin{I+1}) && isreal(varargin{I+1}) ) || ...
               ( ischar(varargin{I+1}) && strcmp(varargin{I+1},'kr') ) 
                hoa2SpkOpt.transFreq = varargin{I+1} ;
            else
                error('Transition frequency must be a real value or ''kr''') ;
            end

        case 'transWidth'
            if isscalar(varargin{I+1}) && isreal(varargin{I+1})
                hoa2SpkOpt.transWidth = varargin{I+1} ;
            else
                error('Transition width must be a real value') ;
            end

            
        case 'spkDistCor'
            if islogical(varargin{I+1})
                hoa2SpkOpt.spkDistCor = varargin{I+1} ;
            else
                error(['Speaker distance correction can be ''true''' ...
                    'or ''false''']) ;
            end

            
        otherwise
            error([varargin{I} ' is not a valid parameter name']) ;
            
    end
end

% Checks if filter length and sampling frequency vectors are consistent
if ~isscalar(hoa2SpkOpt.filterLength)
    if length(hoa2SpkOpt.filterLength) ...
            ~= length(hoa2SpkOpt.sampFreq)
        error (['filterLength must be either a scalar or a vector' ...
                ' having the same size as sampFreq']) ;
    end
else
    hoa2SpkOpt.filterLength = hoa2SpkOpt.filterLength ...
                          * ones(size(hoa2SpkOpt.sampFreq)) ;
end

% If the 'decodType' is not 'mixed', some of the fields are irrelevant
% If the speaker distance correction is not used, even more are discarded
if ~strcmp(hoa2SpkOpt.decodType,'mixed')    
    hoa2SpkOpt.transFreq  = [] ;
    hoa2SpkOpt.transWidth = [] ;
    if ~hoa2SpkOpt.spkDistCor
        hoa2SpkOpt.sampFreq     = [] ;
        hoa2SpkOpt.filterLength = [] ;
    end
end
