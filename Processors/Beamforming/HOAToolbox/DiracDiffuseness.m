function psi = DiracDiffuseness(ambSig,cnv)
% DIRACDIFFUSENESS     Diffuseness of an order-1 HOA sound scene
%
% as used in the DIRAC - DIrectional Audio Coding - technique
%
% psi = DiracDiffuseness(ambSig,cnv)
%
% Inputs: - ambSig is a N x M matrix, where N is the numbe of time samples 
%            or frequency bins and M is the number of signals (3 or 4).
%         - cnv is the convention used for the normalisation of spherical
%            harmonics. Use 'SN2D' for B-format signals. The default value
%            is 'N3D'.
%
% Output: - psi is the diffuseness index for the input set of signals, 
%            comprised between 0 and 1. 


% Default normalisation convention
if nargin < 2
    cnv = 'N3D' ;
end

% Number of harmonics (3/4 for 2d/3d)
nmbHrm = min(size(ambSig,2),4) ;

% Compensate for the normalisation convention
switch lower(cnv)
    
    case 'n3d'
        
        ambSig(:,2:nmbHrm) = ambSig(:,2:nmbHrm) / sqrt(3) ;
        
    case 'n2d'

        ambSig(:,2:nmbHrm) = ambSig(:,2:nmbHrm) / sqrt(2) ;

    case 'sn2d'

        ambSig(:,1) = ambSig(:,1) * sqrt(2)/2 ;

end

% Average intensity vector
avgInt = mean(repmat(ambSig(:,1),1,nmbHrm-1).*ambSig(:,2:nmbHrm)) ; 

% Average energy
avgEng = mean(sum(abs(ambSig(:,1:nmbHrm)).^2,2)) ;

% Diffuseness
psi = 1 - 2 * norm(avgInt) / avgEng ;
