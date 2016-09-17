function spkFmt = GenerateSpkFmt(varargin)

% Test on the number of arguments
if round(size(varargin,2)/2)~=size(varargin,2)/2
    error('illegal number of arguments') ;
end

% If no argument is passed, the function generates a structure
% corresponding to a cube of speakers. Else, the different options are
% assigned.
if nargin < 1
    setupRadius = 2 ;
    distribType = 'even3d' ;
    spkFmt.nbSpk = 8 ;
    flagCoord = [ 1 0 0 ] ;
else
    setupRadius = 2 ;
    distribType = 'even2d' ;
    flagCoord = [ 0 0 0 ] ;
    for I = 1 : 2 : length(varargin)-1
        switch varargin{I}
            case 'nbSpk'
                if ( (isscalar(varargin{I+1})) ...
                        && (varargin{I+1}==abs(round(varargin{I+1}))) )
                    spkFmt.nbSpk = varargin{I+1} ;
                    flagCoord(1) = 1 ;
                else
                    error(['The number of microphone ' ...
                        'must be a positive integer']) ;
                end
            case 'sphCoord'
                if size(varargin{I+1},2)==3 && isreal(varargin{I+1})
                    spkFmt.sphCoord = varargin{I+1} ;
                    flagCoord(2) = 1 ;
                else
                    error(['Loudspeaker spherical coordinates must ' ...
                        'be passed in a [Nx3] real matrix']);
                end
            case 'xyzCoord'
                if size(varargin{I+1},2)==3 && isreal(varargin{I+1})
                    spkFmt.xyzCoord = varargin{I+1} ;
                    flagCoord(3) = 1 ;
                else
                    error(['Loudspeaker cartesian coordinates must ' ...
                        'be passed in a [Nx3] real matrix']);
                end
            case 'distribType'
                switch varargin{I+1}
                    case {'even2d','even3d','ITU5.1'}
                        distribType = varargin{I+1} ;
                    otherwise
                        error('Unknown angular distribution type') ;
                end
            case 'setupRadius'
                if isscalar(varargin{I+1}) && isreal(varargin{I+1}) ...
                        setupRadius = abs(varargin{I+1}) ;
                else
                    error('The setup radius should be a real value') ;
                end
            otherwise
                error([varargin{I} ' is not a valid spkFmt option']) ;

        end
    end
end

% Applying the configuration option
if ( flagCoord(1) && flagCoord(2) ) || ( flagCoord(1) && flagCoord(3) )
    error(['Number of speakers cannot be passed ' ...
        'at the same time as spherical coordinate array ' ...
        'or cartesian coordinate array']);
elseif flagCoord(2) && flagCoord(3)
    error(['Cartesian coordinates cannot be passed ' ...
        'at the same time as spherical coordinates']);
elseif isempty(find(flagCoord,1)) || flagCoord(1)
    switch distribType
        case 'even3d'
            if any(spkFmt.nbSpk==[4 6 8 12 14 20 24 32 60 92 162 642])
                [spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2)] = ...
                    RegularPolyhedron(spkFmt.nbSpk) ;
            else
                [spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2)] = ...
                    SpherePacking(spkFmt.nbSpk) ;
            end
        case 'even2d'
            spkFmt.sphCoord(:,1) = ...
                (0:spkFmt.nbSpk-1)'*2*pi/spkFmt.nbSpk ;
            spkFmt.sphCoord(:,2) = zeros(spkFmt.nbSpk,1) ;
        case 'ITU5.1'
            spkFmt.nbSpk = 5 ;            
            spkFmt.sphCoord(:,1) = [30 -30 0 110 -110]'*pi/180 ;
            spkFmt.sphCoord(:,2) = zeros(5,1) ;
    end
    spkFmt.sphCoord(:,3) = setupRadius * ones(spkFmt.nbSpk,1) ;
    [spkFmt.xyzCoord(:,1),spkFmt.xyzCoord(:,2),spkFmt.xyzCoord(:,3)] = ...
        sph2cart(spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2),spkFmt.sphCoord(:,3)) ;
elseif flagCoord(2)
    [spkFmt.xyzCoord(:,1),spkFmt.xyzCoord(:,2),spkFmt.xyzCoord(:,3)] = ...
        sph2cart(spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2),spkFmt.sphCoord(:,3)) ;
    spkFmt.nbSpk = size(spkFmt.sphCoord,1) ;
elseif flagCoord(3)
    [spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2),spkFmt.sphCoord(:,3)] = ...
        cart2sph(spkFmt.xyzCoord(:,1),spkFmt.xyzCoord(:,2),spkFmt.xyzCoord(:,3)) ;
    spkFmt.nbSpk = size(spkFmt.xyzCoord,1) ;
end

end



