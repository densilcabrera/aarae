function micFmt = GenerateMicFmt(arraysOptions,sphereRadius,sphereImpedance)

% General options are checked and assigned
% All are related to the possible presence of a diffracting sphere
if nargin < 3
    micFmt.sphereImpedance = ImpedanceFunction('constant',Inf) ;
else
    if strcmp(class(sphereImpedance),'function_handle')
        micFmt.sphereImpedance = sphereImpedance ;
    else
        error(['sphereImpedance must be a function handle, ' ...
            'please use ImpedanceFunction.m']) ;
    end
end
if nargin < 2
    micFmt.sphereRadius = 0 ;
else
    if isreal(sphereRadius) && isscalar(sphereRadius)
        micFmt.sphereRadius = sphereRadius ;
    else
        error('Wrong sphere radius value') ;
    end
end
micFmt.nbMic = 0 ;

% In the case where no argument is passed, the function generates
% a structure corresponding to a soundfield microphone.
if nargin < 1
    arraysOptions = { {'xyzCoord',zeros(1,3)} ; ...
        {'xyzCoord',zeros(3),'xyzDirect',eye(3),'micType','eight'} } ;
else
    if ~iscell(arraysOptions)
        error('arraysOptions must be a cell array') ;
    end
end

% Sub-arrays options are checked and assigned
micFmt.nbArrays = size(arraysOptions,1) ;
for I = 1 : micFmt.nbArrays
    if micFmt.nbArrays == 1
        opt = arraysOptions ;
    else
        opt = arraysOptions{I,:} ;
    end
    if round(size(opt,2)/2)~=size(opt,2)/2
        error(['illegal number of arguments in the options ' ...
            'of array number ' num2str(I)]) ;
    end
    micFmt.arrays(I).micType = 'omni' ;
    micFmt.arrays(I).micFreqResp = @(f) ones(size(f)) ;
    micFmt.arrays(I).measurements = [] ;
    distribType = 'RegularPolyhedron' ;
    flagCoord = [ 0 0 0 ] ;
    flagDirect = [ 0 0 ] ;
    for J = 1 : 2 : length(opt)-1
        switch opt{J}
            case 'nbMic'
                if ( (isscalar(opt{J+1})) ...
                        && (opt{J+1}==abs(round(opt{J+1}))) )
                    micFmt.arrays(I).nbMic = opt{J+1} ;
                    flagCoord(1) = 1 ;
                else
                    error(['The number of microphone ' ...
                        'must be a positive integer']) ;
                end
            case 'sphCoord'
                if size(opt{J+1},2)==3 && isreal(opt{J+1})
                    micFmt.arrays(I).sphCoord = opt{J+1} ;
                    flagCoord(2) = 1 ;
                else
                    error(['Microphone spherical coordinates must ' ...
                        'be passed in a [Nx3] real matrix']);
                end
            case 'xyzCoord'
                if size(opt{J+1},2)==3 && isreal(opt{J+1})
                    micFmt.arrays(I).xyzCoord = opt{J+1} ;
                    flagCoord(3) = 1 ;
                else
                    error(['Microphone cartesian coordinates must ' ...
                        'be passed in a [Nx3] real matrix']);
                end
            case 'micType'
                switch opt{J+1}
                    case {'omni','cardio','super', ...
                            'hyper','eight','measured'}
                        micFmt.arrays(I).micType = opt{J+1} ;
                    otherwise
                        error('Unknown microphone type') ;
                end
            case 'micFreqResp'
                if strcmp(class(opt{J+1}),'function_handle')
                    micFmt.arrays(I).micFreqResp = opt{J+1} ;
                else
                    error(['The microphone frequency response must' ...
                        'be a function handle']) ;
                end
            case 'measurements'
                if isstruct(opt{J+1})
                    micFmt.arrays(I).measurements = opt{J+1} ;
                    micFmt.arrays(I).nbMic = ...
                        size(micFmt.arrays(I).measurements.impulseResponses,2) ;
                    micFmt.arrays(I).micType = 'measured' ;
                else
                    error(['Microphone measurements must be passed' ...
                        'in a structure']) ;
                end
            case 'distribType'
                switch opt{J+1}
                    case {'RegularPolyhedron', ...
                            'SpherePacking','SphereCovering'}
                        distribType = opt{J+1} ;
                    otherwise
                        error('Unknown angular distribution type') ;
                end
            case 'sphDirect'
                if size(opt{J+1},2)==2 && isreal(opt{J+1})
                    micFmt.arrays(I).sphDirect = opt{J+1} ;
                    flagDirect(1) = 1 ;
                else
                    error(['Microphone spherical orientation must ' ...
                        'be passed in a [Nx2] real matrix']);
                end
            case 'xyzDirect'
                if size(opt{J+1},2)==3 && isreal(opt{J+1})
                    micFmt.arrays(I).xyzDirect = ...
                        opt{J+1} ./ repmat(sqrt(sum(opt{J+1}.^2,2)),1,3) ;
                    flagDirect(2) = 1 ;
                else
                    error(['Microphone cartesian orientation must ' ...
                        'be passed in a [Nx3] real matrix']);
                end
            otherwise
                error([opt{J} ' is not a valid sub-array option']) ;
        end
    end
    if strcmp(micFmt.arrays(I).micType,'measured') ...
            && isempty(micFmt.arrays(I).measurements)
        error(['A measurement structure must be passed in the \n'...
            'cas where the microphone type is ''measured''']);
    end
    if ( flagCoord(1) && flagCoord(2) ) || ( flagCoord(1) && flagCoord(3) )
        error(['Number of microphones cannot be passed ' ...
            'at the same time as spherical coordinate array ' ...
            'or cartesian coordinate array']);
    elseif flagCoord(2) && flagCoord(3)
        error(['Cartesian coordinates cannot be passed ' ...
            'at the same time as spherical coordinates']);
    elseif isempty(find(flagCoord,1)) || flagCoord(1)
        if ~any( micFmt.arrays(I).nbMic ...
                == [4 6 8 12 14 20 24 32 60 92 162 642] ) ...
                && strcmp(distribType,'RegularPolyhedron')
            distribType = 'SpherePacking' ;
            disp(['Warning: No ' num2str(micFmt.arrays(I).nbMic) '-point ' ...
                'regular polyhedron available... ' ...
                'switching to sphere packing']) ;
        end
        eval(['[micFmt.arrays(' num2str(I) ').sphCoord(:,1),' ...
            'micFmt.arrays(' num2str(I) ').sphCoord(:,2)] = ' ...
            distribType '(micFmt.arrays(' num2str(I) ').nbMic) ;']);
        micFmt.arrays(I).sphCoord(:,3) = ...
            micFmt.sphereRadius * ones(micFmt.arrays(I).nbMic,1) ;
        [micFmt.arrays(I).xyzCoord(:,1), ...
            micFmt.arrays(I).xyzCoord(:,2), ...
            micFmt.arrays(I).xyzCoord(:,3)] = ...
            sph2cart(micFmt.arrays(I).sphCoord(:,1), ...
            micFmt.arrays(I).sphCoord(:,2), ...
            micFmt.arrays(I).sphCoord(:,3)) ;
    elseif flagCoord(2)
        [micFmt.arrays(I).xyzCoord(:,1), ...
            micFmt.arrays(I).xyzCoord(:,2), ...
            micFmt.arrays(I).xyzCoord(:,3)] = ...
            sph2cart(micFmt.arrays(I).sphCoord(:,1), ...
            micFmt.arrays(I).sphCoord(:,2), ...
            micFmt.arrays(I).sphCoord(:,3)) ;
        micFmt.arrays(I).nbMic = size(micFmt.arrays(I).sphCoord,1) ;
    elseif flagCoord(3)
        [micFmt.arrays(I).sphCoord(:,1), ...
            micFmt.arrays(I).sphCoord(:,2), ...
            micFmt.arrays(I).sphCoord(:,3)] = ...
            cart2sph(micFmt.arrays(I).xyzCoord(:,1), ...
            micFmt.arrays(I).xyzCoord(:,2), ...
            micFmt.arrays(I).xyzCoord(:,3)) ;
        micFmt.arrays(I).nbMic = size(micFmt.arrays(I).xyzCoord,1) ;
    end
    if any( micFmt.arrays(I).sphCoord(:,3) < micFmt.sphereRadius )
        error('Microphones cannot be inside the sphere') ;
    end
    if ~any(flagDirect~=0)
        micFmt.arrays(I).sphDirect = micFmt.arrays(I).sphCoord(:,1:2) ;
        [micFmt.arrays(I).xyzDirect(:,1), ...
            micFmt.arrays(I).xyzDirect(:,2), ...
            micFmt.arrays(I).xyzDirect(:,3)] = ...
            sph2cart(micFmt.arrays(I).sphCoord(:,1), ...
            micFmt.arrays(I).sphCoord(:,2),1) ;
    elseif flagDirect(1)
        if size(micFmt.arrays(I).sphDirect,1) ~= micFmt.arrays(I).nbMic
            error(['Dimensions of ''sphDirect'' array are not ' ...
                'coherent with those of ''sphCoord'' or with ''nbMic''']) ;
        end
        [micFmt.arrays(I).xyzDirect(:,1), ...
            micFmt.arrays(I).xyzDirect(:,2), ...
            micFmt.arrays(I).xyzDirect(:,3)] = ...
            sph2cart(micFmt.arrays(I).sphDirect(:,1), ...
            micFmt.arrays(I).sphDirect(:,2),1) ;
    elseif flagDirect(2)
        if size(micFmt.arrays(I).xyzDirect,1) ~= micFmt.arrays(I).nbMic
            error(['Dimensions of ''xyzDirect'' array are not ' ...
                'coherent with those of ''sphCoord'' or with ''nbMic''']) ;
        end
        [micFmt.arrays(I).sphDirect(:,1), ...
            micFmt.arrays(I).sphDirect(:,2)] = ...
            cart2sph(micFmt.arrays(I).xyzDirect(:,1), ...
            micFmt.arrays(I).xyzDirect(:,2), ...
            micFmt.arrays(I).xyzDirect(:,3)) ;
    end
    micFmt.arrays(I).sphCoord(micFmt.arrays(I).sphCoord(:,3)==0,1:2) = ...
        micFmt.arrays(I).sphDirect(micFmt.arrays(I).sphCoord(:,3)==0,:) ;
    micFmt.nbMic = micFmt.nbMic + micFmt.arrays(I).nbMic ;
end

