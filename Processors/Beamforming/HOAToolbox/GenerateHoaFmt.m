function hoaFmt = GenerateHoaFmt(varargin)
%
%
% hoaFmt = GenerateHoaFmt('fieldName',fieldValue,...)
%    
% generate a structure hoaFmt describing an Higher Order Ambisonics (HOA) 
% signal format. 
%
% DESCRIPTION:
% The format of an HOA signal is described by a set of spherical harmonic
% functions to which it refers. This set can be summarized in a structure
% with the following fields
% 
% |----------|---------------------------------------|--------------------|
% |  Field   |               Field                   |      Field         |
% |   name   |            description                |      values        |
% |          |                                       |     {defaultValue} |
% |----------|---------------------------------------|--------------------|
% |   type   | real or complex spherical harmonics ? | 'real', 'complex', |
% |          | (see explanations below)              |           {'real'} |
% |----------|---------------------------------------|--------------------|
% |  res2d   | resolution (max angular frequency) in | any real-positive  |
% |          | the horizontal plane (azimuth)        | integer            |
% |          | (see explanations below)              |                {1} |
% |----------|---------------------------------------|--------------------|
% |  res3d   | resolution (max angular frequency) in | any real-positive  |
% |          | the vertical plane (elevation)        | integer            |
% |          | (see explanations below)              |                {1} |
% |----------|---------------------------------------|--------------------|
% |   conv   | encoding convention                   | 'N2D', 'SN2D',     |
% |          |                                       | 'N3D', 'SN3D'      |
% |          |                                       |            {'N3D'} |
% |----------|---------------------------------------|--------------------|
% |  index   | Two-column array of m and n spherical | any consistant     |
% |          | harmonic indices                      | 2-column array     |
% |          | (see explanations below)              |  {depends on res.} |
% |----------|---------------------------------------|--------------------|
%
% GenerateHoaFmt generates this structure given some field values.
%
% REMARKS:
% - ('fieldName',fieldValue) pairs can be passed in any order
% - About the type of spherical harmonics:
% Real and complex spherical harmonics are not arranged in the same way, so
% that rows of the index array are not going to be in the same order for
% given res2d and res3d values (see explanations below).
% - About res2d and res3d:
% Usually, horizontal and vertical resolutions are equal. However, res2d
% and res3d can take any real positive integer value: in this case the 
% spherical harmonics chosen are those having less or as much as res3d 
% nodes in the elevation dimension and less or as much as 2*res2d nodes 
% in the azimuthal dimension.
% - About the index array:
% Index array consists in two columns describing which spherical harmonic
% functions the HOA signal refers to, including the way these functions are
% arranged in. In a "normal" case, this array is computed as a function of
% res2d, res3d, and type field values.
% Note that the order in which spherical harmonics are arranged is not the
% same in the case of real and complex harmonics. Real ones are ordered in
% the same way as in B-format. Complex ones are ordered in a way that makes
% the use of Fourier Transform easier. Please give a try to see how it
% works.
% If arbitrary values are passed simultaneously for res2d, res3d AND index
% array, then res2d and res3d values can be modified according to the index
% array. res2d/res3d then become maximum horizontal/vertical resolutions. 
% - About conventions:
% The default convention is 'N3D', which ensures that there is the same
% amount of energy in every order for a set of 3D harmonics up to order L.
% 'SN2D' is also known as the Furse-Malham convention, used in the *.amb
% file format.
% 
% EXAMPLES:
% - To describe an order (degree) 4, 3D, complex HOA signal:
%  hoaFmt = GenerateHoaFmt('res2d',4,'res3d',4,'type','complex')
% - To describe an order (degree) 3, 2D, real HOA signal:
%  hoaFmt = GenerateHoaFmt('res2d',3,'res3d',0,'conv','N2D')
% (type is 'real' by default) 
% - To describe a specific HOA signal, with only two tracks corresponding 
% to the (5,2) and (1,7) complex spherical harmonics :
%  hoaFmt = GenerateHoaFmt('type','complex','index',[5 2;1 7])
% (in this case res2d will be set to 7 and res3d will be set to 3)
%
%
% N. Epain on a J. Daniel idea - Last update: 17/12/2007

% Default values for hoaFmt structure fields
hoaFmt.type     = 'real' ;
hoaFmt.res2d    = 1 ;
hoaFmt.res3d    = 1 ;
hoaFmt.conv     = 'N3D' ;

% hoafmt.index is empty for the moment
hoaFmt.index = [] ;

% 'Plane' option is empty
planeOpt = [] ;

% Number of arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% Field values are checked and assigned to hoaFmt structure
for I = 1 : 2 : length(varargin)-1
    
    switch varargin{I}
        
        case 'type'
            
            switch varargin{I+1}
                case {'real','complex'}
                    hoaFmt.type = varargin{I+1} ;
                otherwise
                    error([varargin{I+1} ' is not a valid HOA type']) ;
            end

        case {'res2d','res3d'}
            
            if ( (isscalar(varargin{I+1})) ...
               && (varargin{I+1}==abs(round(varargin{I+1}))) )
                hoaFmt.(varargin{I}) = varargin{I+1} ;
            else
                error([varargin{I} ' cannot be equal to ' num2str(varargin{I+1})]) ;
            end

        case 'conv'
            
            switch varargin{I+1}
                case {'N2D','SN2D','N3D','SN3D'}
                    hoaFmt.conv = varargin{I+1} ;
                otherwise
                    error([varargin{I+1} ...
                           ' is not a valid HOA convention']) ;
            end

        case 'index'
            
            if varargin{I+1} == round(varargin{I+1})
                hoaFmt.index = varargin{I+1} ;
            else
                error('illegal index values') ;
            end

        case 'plane'
            
            switch varargin{I+1}
                case {'xOy','xOz','yOz'}
                    planeOpt = varargin{I+1} ;
                otherwise
                    error([varargin{I+1} ' is not a valid plane']) ;
            end

        otherwise
            
            error([varargin{I} ' is not available as an hoaFmt field']) ;
            
    end
    
end


if isempty(hoaFmt.index)
    
    % Computes the index array, based on type, res2d and res3d values
    switch hoaFmt.type
        
        case 'real'
            for I = 0 : max(hoaFmt.res2d,hoaFmt.res3d)
                J = I : -1 : 0 ;
                for J = J( (J<=hoaFmt.res2d) & (J>=(I-hoaFmt.res3d)) )             
                    if J~=0
                        hoaFmt.index = [ hoaFmt.index ; [I,J;I,-J] ] ;
                    else
                        hoaFmt.index = [ hoaFmt.index ; [I,J] ] ;
                    end
                end
            end
            
        case 'complex'
            for J = 0 : max(hoaFmt.res2d)
                I = J : max(hoaFmt.res2d,hoaFmt.res3d) ;
                for I = I(I<=(hoaFmt.res3d+J))
                    hoaFmt.index = [ hoaFmt.index ; [I,J] ] ;
                end
            end
            for J = -max(hoaFmt.res2d) : -1
                I = max(hoaFmt.res2d,hoaFmt.res3d) : -1 : -J ;
                for I = I(I<=(hoaFmt.res3d)-J)
                    hoaFmt.index = [ hoaFmt.index ; [I,J] ] ;
                end
            end
            
    end
    
end

% If asked, remove the harmonics that are null in the specified plane
if ~isempty(planeOpt)
    index = hoaFmt.index ;
    switch hoaFmt.type
        case 'real'
            switch planeOpt
                case 'xOy'
                    hoaFmt.index = ...
                        index(~mod(index(:,1)-abs(index(:,2)),2),:) ;
                case 'xOz'
                    hoaFmt.index = index(index(:,2)>=0,:) ;
                case 'yOz'
                    hoaFmt.index = index(index(:,2)<=0,:) ;
            end
        case 'complex'
            if strcmp(planeOpt,'xOy')
                hoaFmt.index = ...
                    index(~mod(index(:,1)-abs(index(:,2)),2),:) ;
            end
    end
end

% Resolution is corrected to fit with index : res2d/res3d then
% represent MAXIMUM azimuthal/elevation resolution
hoaFmt.res2d = max(abs(hoaFmt.index(:,2))) ;
hoaFmt.res3d = max( hoaFmt.index(:,1) - abs(hoaFmt.index(:,2)) ) ;

% Number of spherical harmonics selected
hoaFmt.nbComp = size(hoaFmt.index,1) ;

