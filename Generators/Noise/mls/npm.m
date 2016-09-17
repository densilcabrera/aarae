function [npm_val npm_val_dB] = npm(h, hhat)

% Normalized Projection Misalignment 
%
%   [npm_val] = npm(h, hhat)
%
%	Input Parameters [size]:
%       h       : true impulse responses [L x M]
%       hhat    : estimated impulse responses [L x M]
%
%	Output Parameter:
%       npm_val : Normalize Projection Misalignment
%
%	References:
%       [1] D. R. Morgan, J. Benesty and M. M. Sondhi, "On the evaluation of
%           estimated impulse responses," IEEE Signal Processing Lett., Vol. 5, No.
%           7, pp. 174-176 Jul 1998.
%   
%       [2] Y. Huang and J. Benesty, "Frequency-Domain adaptive approaches to
%           blind multi-channel identification," IEEE Trans. Sig. Process. Vol. 51
%           No. 1, pp/ 11-24, Jan 2003.
%
% Authors: N.D. Gaubitch 
%
% History: 2004-04-26 - Initial version
%
% Copyright (C) Imperial College London 2009
% Version: $Id$

hv = h(:);
hhatv = hhat(:);

epsilon = hv-((hv.'*hhatv)/(hhatv.'*hhatv))*hhatv;
npm_val = norm(epsilon)/norm(hv);
npm_val_dB = 20*log10(norm(hv-((hv.'*hhatv)/(hhatv.'*hhatv))*hhatv)/norm(hv));