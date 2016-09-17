function mic2HoaOpt = GenerateMic2HoaOpt(varargin)

% Default parameters
mic2HoaOpt.sampFreq      = 48000 ;
mic2HoaOpt.filterType    = 'postEqFir' ;
mic2HoaOpt.filterLength  = 256 ;
mic2HoaOpt.limiterMethod = 'tikh' ;
mic2HoaOpt.limiterLevel  = 6 ;
mic2HoaOpt.higherOrders  = false ;
mic2HoaOpt.subArrayFilt  = false ;
mic2HoaOpt.highFreqEq    = true ;
mic2HoaOpt.lowPassFreq   = 22000 ;

% Number of arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% Parameters values are checked and assigned to the structure if possible
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        
        case 'filterType'
            switch varargin{I+1}
                case {'firMatrix','postEqFir'}
                    mic2HoaOpt.filterType = varargin{I+1} ;
                otherwise
                    error('Unknown filter type') ;
            end
            
        case 'sampFreq'
            if isvector(varargin{I+1}) && isreal(varargin{I+1})
                mic2HoaOpt.sampFreq = varargin{I+1}(:) ;
            else
                error(['Required sampling frequencies must be passed' ...
                       'in a vector of real values']) ;
            end
            
        case 'filterLength'
            if isvector(varargin{I+1}) && isreal(varargin{I+1})
                mic2HoaOpt.filterLength = varargin{I+1}(:) ;
            else
                error(['Filter lengths must be passed' ...
                       'in a vector of real values']) ;
            end
            
        case 'limiterMethod'
            switch varargin{I+1}
                case {'tikh','expo','max'}
                    mic2HoaOpt.limiterMethod = varargin{I+1} ;
                otherwise
                    error('Unknown limiter method') ;
            end
            
        case 'limiterLevel'
            if isscalar(varargin{I+1}) && isreal(varargin{I+1})
                mic2HoaOpt.limiterLevel = varargin{I+1} ;
            else
                error('Limiter level must be a real scalar') ;
            end
            
        case 'higherOrders'
            if islogical(varargin{I+1})
                mic2HoaOpt.higherOrders = varargin{I+1} ;
            else
                error('Higher orders option must be a logical value') ;
            end
            
        case 'subArrayFilt'
            if islogical(varargin{I+1})
                mic2HoaOpt.subArrayFilt = varargin{I+1} ;
            else
                error('Sub-array filtering option must be a logical value') ;
            end
            
        case 'highFreqEq'
            if islogical(varargin{I+1})
                mic2HoaOpt.highFreqEq = varargin{I+1} ;
            else
                error('High frequency equalization option must be a logical value') ;
            end
            
         case 'lowPassFreq'
            if isscalar(varargin{I+1}) && isreal(varargin{I+1})
                mic2HoaOpt.lowPassFreq = varargin{I+1} ;
            else
                error('Low-pass frequency must be a real scalar') ;
            end
            
       otherwise
            error([varargin{I} ' is not a valid parameter name']) ;
            
    end
end

% Checks if filter length and sampling frequency vectors are consistent
if ~isscalar(mic2HoaOpt.filterLength)
    if length(mic2HoaOpt.filterLength) ...
            ~= length(mic2HoaOpt.sampFreq)
        error (['filterLength must be either a scalar or a vector' ...
                ' having the same size as sampFreq']) ;
    end
else
    mic2HoaOpt.filterLength = mic2HoaOpt.filterLength ...
                          * ones(size(mic2HoaOpt.sampFreq)) ;
end

% A few final tests
if mic2HoaOpt.higherOrders && strcmp(mic2HoaOpt.filterType,'postEqFir')
    fprintf(['WARNING: higher orders cannot be taken into account ' ...
        'using the ''postEqFir'' method... ' ...
        'switching to FIR matrix method \n']) ;
    mic2HoaOpt.filterType = 'firMatrix' ;
end
if mic2HoaOpt.subArrayFilt && strcmp(mic2HoaOpt.filterType,'postEqFir')
    fprintf(['WARNING: filtering cannot be done with the ''sub-array''' ...
        ' method when using ''postEqFir''... ' ...
        'switching to FIR matrix method \n']) ;
    mic2HoaOpt.filterType = 'firMatrix' ;
end
