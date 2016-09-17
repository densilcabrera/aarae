function OUT = FDTDdemo(roomDims,srcPos,recPos,fs,c,simLength)
% Finite difference time domain demonstrator, adapted from the following
% source:
% **************************************************
%
% Example 3D-FDTD code to simulate scattering
%
% v1.0 - July 2015
%
% Jonathan Sheaffer, Ben-Gurion University, Israel
% Trevor Cox, University of Salford, UK
% (BSD license from Matlab Central)
% ************************************************** 

% clear
% close
% clc

%% ************************************************** 
% DEFAULT SIMULATION PARAMETERS
%

if nargin == 0 % FIX
roomDims = [6,6,4];        % Room dimensions [Lx, Ly, Lz] in meters
srcPos = [3,1.5,2];        % Source position [x,y,z] in meters 

fs=10000;                  % Sample rate (Hz)
lam=sqrt(1/3);             % Courant number
c=343.5;                   % Speed of sound in air
simLength = 0.05;          % Simulation length (seconds)
ra = 1e-10;                % Admittance of rigid materials

T=1/fs;                    % (temporal) sample period
X=c*T/lam;                 % (spatial) sample period

maxc = 0.03;               % caxis limits (for visualisation only)

roomDimsD = round(roomDims/X)+2;    % Room dimensions [Nx,Ny,Nz] in nodes
srcPosD = round(srcPos/X);          % Discrete source position
maxN = round(simLength/T);

scattering = 1; % surface conditions of the room: 0 for no scattering
doplot = 1; % whether or not to provide animated slice through the room
rec


    % dialog box for settings
    prompt = {'Length, width and height of the room or box (m)', ...
        'Coordinates of the source (m)', ...
        'Coordinates of the receiver (m)', ...
        'Speed of sound (m/s)', ... 
        'Spatial jitter of reflections (m)', ...
        'Maximum order of reflections', ...
        'Ambisonics order', ...
        'Sampling rate (Hz)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {[num2str(roomDims(1)),', ',num2str(roomDims(2)),', ',num2str(roomDims(3))], ...
        [num2str(srcPos(1)),', ',num2str(srcPos(2)),', ',num2str(srcPos(3))], ...
        [num2str(xr),', ',num2str(yr),', ',num2str(zr)], ...
        num2str(c), num2str(jitter), num2str(maxorder), ...
        num2str(ambiorder), num2str(fs)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        OUT = [];
        return
    else
        roomDims = str2num(answer{1,1});
        srcPos = str2num(answer{2,1});
        r = str2num(answer{3,1});
        xr = s(1); yr = s(2); zr = s(3);
        c = str2num(answer{4,1});
        jitter = str2num(answer{5,1});
        maxorder = str2num(answer{6,1});
        ambiorder = str2num(answer{7,1});
        fs = str2num(answer{8,1});
    end
    
end


%% **************************************************
% DESIGN EXCITATION FUNCTION
%
% We design a simple Gaussian pulse, and pass it through a 
% DC blocking filter (to avoid exciting near zero frequency)
%
fc = 0.075;     % Cutoff frequency (normalised 0.5=nyquist)
n0 = 30;        % Initial delay (samples)
sigma=sqrt(2*log(2))/(2*pi*(fc/T));
n=0:maxN;
srcFn=exp(-T^2*(n-n0).^2/(2*sigma^2));
srcFn=filter([1 -1], [1 -0.995], srcFn); % DC Blocker

%% **************************************************
% ALLOCATE MATRICES
%
p = zeros([roomDimsD 3]);   % Pressure matrix (Nx x Ny x Nz x 3)
                            % Fourth dimension is discrete time
                            % (:,:,:,1) is (n+1)
                            % (:,:,:,2) is (n)
                            % (:,:,:,3) is (n-1)
                            
K = zeros(roomDimsD);       % Node classification mask
                            % 6 = air, 5 = surface, 4 = edge, 3 = corner
                            % All nodes are air by default
                            
beta = ones(roomDimsD)*ra;  % Boundary admittance mask
                            % default is rigid material ("ra")
                            
% Define shift vectors
ll=2:roomDimsD(1)-1;
mm=2:roomDimsD(2)-1;
ii=2:roomDimsD(3)-1;

%% **************************************************
% DEFINE BOUNDARY CONDITIONS
%

% Walls - absorptive
K(2,:,:)=5;         beta(2,:,:)=1;      
K(end-1,:,:)=5;     beta(end-1,:,:)=1;
K(:,2,:)=5;         beta(:,2,:)=1;
K(:,end-1,:)=5;     beta(:,end-1,:)=1; 
K(:,:,2)=5;         beta(:,:,2)=1;
K(:,:,end-1)=5;     beta(:,:,end-1)=1;

K(3:end-2,3:end-2,3:end-2)=6; 

switch scattering
    case 1
        %Scatterer - reflective
        K(30:35,2:10,2:end-1)=0;
        K(30,2:10,2:end-1)=5;
        K(35,2:10,2:end-1)=5;
        K(30:35,10,2:end-1)=5;
        K(30:35,2,2:end-1)=5;
        K(30,10,2:end-1)=6;
        K(35,10,2:end-1)=6;
        
        K(37:42,2:6,2:end-1)=0;
        K(37,2:6,2:end-1)=5;
        K(42,2:6,2:end-1)=5;
        K(37:42,6,2:end-1)=5;
        K(37:42,2,2:end-1)=5;
        K(37,6,2:end-1)=6;
        K(42,6,2:end-1)=6;
        
        K(45:50,2:8,2:end-1)=0;
        K(45,2:8,2:end-1)=5;
        K(50,2:8,2:end-1)=5;
        K(45:50,8,2:end-1)=5;
        K(45:50,8,2:end-1)=5;
        K(45,8,2:end-1)=6;
        K(50,8,2:end-1)=6;
        
        K(53:58,2:5,2:end-1)=0;
        K(53,2:5,2:end-1)=5;
        K(58,2:5,2:end-1)=5;
        K(53:58,5,2:end-1)=5;
        K(53:58,5,2:end-1)=5;
        K(53,5,2:end-1)=6;
        K(58,5,2:end-1)=6;
        
        K(61:66,2:9,2:end-1)=0;
        K(61,2:9,2:end-1)=5;
        K(66,2:9,2:end-1)=5;
        K(61:66,9,2:end-1)=5;
        K(61:66,9,2:end-1)=5;
        K(61,9,2:end-1)=6;
        K(66,9,2:end-1)=6;
        % add more cases here
        
        % otherwise no scattering
end


%% **************************************************
% PREPARE SIMULATION
%

% Intermediate boundary matrix
BK = (6-K)*lam.*beta/2;
zero_K=find(K==0);

% Visualisation infrastructure
bv=K(:,:,srcPosD(3)); bv(bv==6)=0;
pl = abs(p(:,:,srcPosD(3),1));
h1=surface(bv); hold on;
h2=pcolor(pl); colormap(flipud(bone)); hold off;
colorbar;
caxis([0 maxc]);
xlim([2 roomDimsD(1)-1]);
ylim([2 roomDimsD(2)-1]);
shading interp;

%% **************************************************
% EXECUTE SIMULATION
%
for nn=2:maxN

    % Update the grid
    p(ll,mm,ii,1) = ((2-K(ll,mm,ii)*lam^2).*p(ll,mm,ii,2) + (BK(ll,mm,ii)-1).*p(ll,mm,ii,3) ...
        + (lam^2) * (p(ll+1,mm,ii,2) + p(ll-1,mm,ii,2) + p(ll,mm+1,ii,2) ...
        + p(ll,mm-1,ii,2) + p(ll,mm,ii+1,2) + p(ll,mm,ii-1,2))) ./ (1+BK(ll,mm,ii));

    % Enforce zero pressure at 'dead spaces' (inside objects)
    p(zero_K)=0;

    % Inject the source (soft source)
    p(srcPosD(1),srcPosD(2),srcPosD(3),2) = p(srcPosD(1),srcPosD(2),srcPosD(3),2) + srcFn(nn);
    
    % Shift matrices in time
    p(:,:,:,3) = p(:,:,:,2);
    p(:,:,:,2) = p(:,:,:,1);
    
    % Visualise in 2D
    % Plane parallel to the floor at the source's height
    % Absolute value of pressure, linear scale
    pl = abs(p(:,:,srcPosD(3),1));
    set(h2,'CData',pl);
    if nn==150, maxc=maxc/2; caxis([0 maxc]); end
    title(sprintf('Sample number %d out of %d\n Note that for visual clarity, the colorscale changes at nn=150.',nn,maxN));
    drawnow;    
end