function DOA = SDMPar(IR, p)
% LOC = SDMPar(IR, p)
% Spatial Decomposition Method [1], implemented as parallel processing
% Returns locations (or direction of arrival) at each time step of the IR. 
%
% USAGE:
% 
% IR  : A matrix of measured impulse responses from a microphone array 
%       [N numberOfMics]
% DOA : 3-D Cartesian locations of the image-sources [N 3];
% p   : struct from createSDMStruct.m or can be manually given
% 
% IR matrix includes impulse responses (pressure) from a microphone array.
% SDM assumes an omni-directional directivity pattern, however,
% the analysis will provide reasonable results with cardioid patterns
% 
% Parameters can be given as a struct:
% Radius = 0.025;      % Radius of the microphones [m]
%        p.micLocs = ... % Microphone positions in Cartesian coordinates [m]
%            [Radius 0 0;
%            -Radius 0 0;
%            0 Radius 0;
%            0 -Radius 0;
%            0 0 Radius;
%            0 0 -Radius];
% p. fs = 48000;
% p.c = 345;
% p.winLen = 0;
% p.parFrames = 8192;
% p.showArray = false;
% 
% 
% Parameters can also obtained from createSDMStruct();
% Syntax:
% p = createSDMStruct('micLocs', micLocs,'c', c, ...
% 'fs',fs,'parFrames',a_number,'winLen',b_number)
% 
% where
% 'micLocs': a matrix, Microphone Array Geometry in Cartesian coordinates, 
%          Default mic array:
% Diameter = 0.1;      % Radius of the microphones [m]
% micLocs = ...      % Microphone positions in Cartesian coordinates [m]
%     [Diameter/2 0 0;
%     -Diameter/2 0 0;
%     0 Diameter/2  0;
%     0 -Diameter/2 0;
%     0 0  Diameter/2;
%     0 0 -Diameter/2]; 
%     
% 'fs' : a single value, sampling frequency [Hz] [default
% 48e3]
% 'c' : a single value speed of sound [m/s] [default 345]
%
% 'winLen': a single value, length of the processing window, if empty, 
% minimum size is used. [default minimum size]
% minimum size is winLen = 2*ceil(max(D)/c*fs) + 8, where D are the
% distances between microphone pairs
% 
% 'parFrames' : a single value, parallel frames in the processing, 
% maximum is length(IR), and minimum is 1, usually 2^10-2^15 is a good 
% choice for modern computers. [default 2^13]
% 
% EXAMPLES
% 
% For further examples for setting the parameters, see createSDMStruct.m
% or demo*.m
% 
% References
% [1] Spatial decomposition method for room impulse responses
% S Tervo, J Pätynen, A Kuusinen, T Lokki - Journal of the Audio
% Engineering Society, vol. 61, no. 1/2, pp. 16-27, 2013

% SDM toolbox : SDMpar
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

if nargin < 1
    error(['SDM toolbox : SDMpar : Two inputs are required ' ... 
    ' SDMPar(IR, p)']);
end

disp('Started SDM processing'); tic;

% Constants for array processing
numOfMics = size(p.micLocs,1); % number of microphones
pairs = nchoosek(1:numOfMics,2); % microphone pairs
numOfPairs = size(pairs,1); % number of microphone pairs
V = p.micLocs(pairs(:,1),:) - p.micLocs(pairs(:,2),:); % microphone vector difference matrix
D = sqrt(sum(V.^2,2)); % Distances between microphone pairs

% Choose the frame size
winLenMin = 2*ceil(max(D)/p.c*p.fs) + 8;
if ~exist('p.winLen','var') || p.winLen < winLenMin
   winLen = winLenMin;
   disp(['Using frame size ' num2str(winLen)])
else
    winLen = p.winLen;
end

% Windowing
W = hanning(winLen);
maxTD = ceil(max(D)/p.c*p.fs)+1; % maximum accepted time delay difference
inv_V = pinv(V); % pseudo inverse of the vector difference matrix
offset = winLen/2 + 1; % time delay offset in the vectors
offset_v = offset+(-maxTD:maxTD+1); % Region of interest in TDOA estimation

% Variables for the frame based processing in the for-loop
startx = 1;
nforward = 1;
endx = startx + winLen - 1;
numOfFrames = floor((length(IR)-winLen)/nforward);

% This indexing is required for the parallel processing
base_inds = (startx:endx)' -1;
base_par_inds = repmat(base_inds,[1 p.parFrames]);
for nn = 1:p.parFrames
    base_par_inds(:,(nn)) = ...
        base_par_inds(:,(nn)) + nforward*(nn-1);
end

% Also this is for the parallel processing
par_pairs = zeros(numOfPairs*p.parFrames,2);
for nn = 1:p.parFrames
    par_pairs((nn-1)*numOfPairs+(1:numOfPairs),:) = ...
        pairs+(nn-1)*numOfMics;
end

% --- Frame based processing ---- 
DOA = zeros(length(IR),3);
for n = 1:p.parFrames:numOfFrames

    % Select the frames that are to be processed
    cur_frames = n:min(n+p.parFrames-1,numOfFrames);
    n_cur_frames = length(cur_frames);

    % Create overlapping indexes
    par_inds = base_par_inds(:,1:n_cur_frames) + n;

    % Index time-domain signal into parallel overlapping windows
    temp = IR(par_inds,:);
    temp = reshape(temp,[winLen n_cur_frames numOfMics]);
    temp = permute(temp,[1 3 2]);
    temp = reshape(temp,[winLen n_cur_frames*numOfMics]);
    
    % Apply the windowing function to the frames
    par_samples = bsxfun(@times, W, temp);
        
    % Run parallel fft for overlapping windows
    X = fft(par_samples,[],1);
    
    X1 = X(:,par_pairs(1:n_cur_frames*numOfPairs,1));
    X2 = X(:,par_pairs(1:n_cur_frames*numOfPairs,2));
    
    % Cross Power Spectrum between microphone pairs
    P12 = X1.*conj(X2);

    % Cross-Correlation between microphone pairs
    Rtemp = real(ifft(P12));
    R = fftshift2(Rtemp);
    
    % Get rid off too low values (only positive correlation accepted)
    R(R<(eps*2)) = eps*2;    
    
    % Remove time delays that are beyond [-t_max, t_max]
    R2 = R(offset_v,:);
    
    % Find maximum
    [~, indk] = max(R2);
    
    % Compensate the offset
    indk = indk + offset - (maxTD+1);

    % Interpolate the maximum peak by assuming exponential shape
    tau = interpolateTau(R,indk,true);
    tau_mat = reshape(tau,[numOfPairs n_cur_frames]);
    d = tau_mat/p.fs;
    
    % Solve the direction of arrival
    % This is the least squares solution assuming plane wave propagation
    % model, i.i.d. noise in microphones, and wideband reflections
    k = -inv_V*d;

    % Transfer to Cartesian coordinate system with distance
    % equal to the center of the current frame
    % [AZ,EL,~] = cart2sph(k(1,:),k(2,:),k(3,:));
    % clear x
    % [x(:,1),x(:,2),x(:,3)] = sph2cart(AZ,EL,0.5*(par_inds(1,:)+par_inds(winLen,:))/fs*c);
    
    % Normalize k to unity to obtain DOA
    k = bsxfun(@times, k, 1./sqrt(sum(k.^2)));
    
    % The distance that the sound has travelled
    d_n = (cur_frames + winLen/2)/p.fs*p.c;
    
    % Save the "image-sources", i.e, locations of the reflections
    DOA(cur_frames+winLen/2,:) = bsxfun(@times, k, d_n)';

end
% --- EOF Frame based processing ---- 
if p.showArray
    plotArray(p.micLocs);
end
disp(['Ended SDM processing in ' num2str(toc) ' seconds.' ])
end % <-- EOF SDMpar()

function plotArray(micLocs)
% plotArray(micLocs)
% Plot the geometry of microphone, defined in micLocs

% SDM toolbox : plotArray
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Copyleft

numOfMics = size(micLocs,1);
figure
hold on
for m = 1:numOfMics
    plot3([0 micLocs(m,1)],[0 micLocs(m,2)],[0 micLocs(m,3)],'color','k')
    plot3(micLocs(m,1),micLocs(m,2),micLocs(m,3),'color','k','Marker','o')
    text(micLocs(m,1),micLocs(m,2),micLocs(m,3),[ 'mic \# ' num2str(m)]);
end
hold off;grid on;
view([35 45]);
axis square;
title('Microphone Array Geometry in the Analysis');
xlabel('X-coord.');ylabel('Y-coord.');zlabel('Z-coord.');

end

function tau = interpolateTau(R,ind,offsetflag)
% tau = InterpolateTau(R,ind,offsetflag) interpolates time delay
% values for a matrix R of cross correlation vectors, ind describes the
% indices where the maxima lie in the vectors.
%
% Solve the coefficients of an exponential function 
% with linear regression.
% 
% The assumed function shape is f = a*exp(b*(t-c))
%
% R   : cross correlation matrix, it is required that R > 0, for all ind +- 1
% ind : maximum indices in the vector
% offsetflag: if true, sets the time delay to length(R)/2 + 1;

% References:
% [1] "On cross correlation based-discrete time delay estimation"
% L. Zhang, X. Wu, ICASSP 2005

% SDM toolbox : interpolateTau
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Copyleft

offset = 0;
if offsetflag
    offset = floor(size(R,1)/2)+1.0;
end
 
% index to time difference of arrival conversion
t = ind;
t2 = t - offset;

% Time indices for all the rows in the matrix
inds_nom = size(R,1)*(0:(size(R,2)-1)) + t;

% Solve the time delay from [-1 0 1] samples around the maximum

c = (log(R(inds_nom + 1)) - log(R(inds_nom - 1)))./...
    (4*log(R(inds_nom)) - 2*log(R(inds_nom - 1)) - 2*log(R(inds_nom + 1)));

% Add the initial integer delay to the interpolated delay
tau = t2 + c;

end
