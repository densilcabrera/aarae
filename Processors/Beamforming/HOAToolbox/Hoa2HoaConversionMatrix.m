function Mat = Hoa2HoaConversionMatrix(hoaFmtInp,hoaFmtOut)
% 
% HOAFMT CONVERSION UTILITY
%
% If bInp is a vector of spherical harmonic components in the basis
% described by hoaFmtInp, and bOut is the vector of corresponding components
% in the basis described by hoaFmtOut, then this function calculates a
% matrix Mat that verifies
%                         bOut = Mat * bInp
% Note: in the case where some harmonics are present in the hoaFmtOut basis
% but not in the hoaFmtInp basis, the corresponding output components are 
% equal to 0.

% Harmonic index mapping from input to output
if strcmp(hoaFmtInp.type,hoaFmtOut.type)
    Mat = repmat(permute(hoaFmtOut.index,[1 3 2]),1,hoaFmtInp.nbComp) ...
        == repmat(permute(hoaFmtInp.index,[3 1 2]),hoaFmtOut.nbComp,1) ;
    Mat = Mat(:,:,1) & Mat(:,:,2) ;
else
    MatOrd = repmat(hoaFmtOut.index(:,1),1,hoaFmtInp.nbComp) ...
        == repmat(hoaFmtInp.index(:,1).',hoaFmtOut.nbComp,1) ;
    MatHrm = repmat(abs(hoaFmtOut.index(:,2)),1,hoaFmtInp.nbComp) ...
        == repmat(abs(hoaFmtInp.index(:,2)).',hoaFmtOut.nbComp,1) ;
    Mat = MatOrd & MatHrm ;
    MatInpPos = repmat(hoaFmtInp.index(:,2).',hoaFmtOut.nbComp,1) >=0 ;
    MatInpNeg = repmat(hoaFmtInp.index(:,2).',hoaFmtOut.nbComp,1) <=0 ;
    MatOutSgn = sign(repmat(hoaFmtOut.index(:,2),1,hoaFmtInp.nbComp)) ;
    if strcmp(hoaFmtInp.type,'real') && strcmp(hoaFmtOut.type,'complex')
        Mat = Mat .* ( MatInpPos + 1i*MatInpNeg.*MatOutSgn ) ;
        Mat = Mat .* ( (MatOutSgn==0) + 1/sqrt(2)*(MatOutSgn~=0) ) ;
    else
        Mat = Mat .* ( (MatOutSgn==0) ...
            + (MatOutSgn>0)/sqrt(2) + (MatOutSgn<0)/sqrt(2)/1i ) ;
        Mat(MatInpNeg&(MatOutSgn<0)) = - Mat(MatInpNeg&(MatOutSgn<0)) ;
    end
end

% Normalisation conversion
if ~strcmp(hoaFmtInp.conv,hoaFmtOut.conv)
    switch hoaFmtInp.conv
        case 'SN3D'
            convCoefInp = @(m) 1 ;
        case 'N3D'
            convCoefInp = @(m) sqrt(2*m+1) ;
        case 'SN2D'
            convCoefInp = @(m) sqrt(2*m+1) .* 2.^m .* gamma(m+1) ...
                ./ sqrt(gamma(2*m+1+1)) ./ sqrt(2) ;
        case 'N2D'
            convCoefInp = @(m) sqrt(2*m+1) .* 2.^m .* gamma(m+1) ...
                ./ sqrt(gamma(2*m+1+1)) ;
    end
    switch hoaFmtOut.conv
        case 'SN3D'
            convCoefOut = @(m) 1 ;
        case 'N3D'
            convCoefOut = @(m) sqrt(2*m+1) ;
        case 'SN2D'
            convCoefOut = @(m) sqrt(2*m+1) .* 2.^m .* gamma(m+1) ...
                ./ sqrt(gamma(2*m+1+1)) ./ sqrt(2) ;
        case 'N2D'
            convCoefOut = @(m) sqrt(2*m+1) .* 2.^m .* gamma(m+1) ...
                ./ sqrt(gamma(2*m+1+1)) ;
    end
    Mat = diag(convCoefOut(hoaFmtOut.index(:,1))) ...
        * Mat * diag(1./convCoefInp(hoaFmtInp.index(:,1))) ;
end

