function hnd = PolarPlot(theta,Rho,scale,varargin) 
% POLARPLOT  Improved polar coordinate plot.
%
%     PolarPlot(theta,Rho,scale,linestyle)
%
%     Draw a polar plot of the data defined in theta (polar angle) and Rho
%     (corresponding radius values).
%
%     Inputs:
%     - theta is the polar angle in radiants. theta should be a vector
%     - Rho is an array containing radius values corresponding to theta.
%       Rho MUST have as many rows as the number of elements in theta.
%     - scale is the scale used in the polar plot. It can be set to 'log'
%       for a log scale with a 40dB dynamic range, or 'linear'
%       ('linear' is the default if no value is provided).
%     - linestyle is a list of linestyle options
%       (as used in the plot function)
%
%     Output:
%     - hnd = PolarPlot(theta, ...) returns a handle to the plotted object
%
%     Usage example:
%   
%       azm = (0:179)'*2*pi/180 ;
%       PolarPlot(azm,abs([cos(azm) sin(azm)]),[],':','linewidth',2) ;
%
% N.Epain 2010

fontSize = 16 ;
lineWidth = 1 ;

% Default scale: linear
if nargin<3 || isempty(scale)
    scale = 'linear' ;
end

% Check the scale
if ~strcmpi(scale,'log') && ~strcmpi(scale,'linear')
    error('Unknown scale') ;
end

% Check the size of Rho
if numel(theta) ~= size(Rho,1)
    error('Rho must have as many rows as there are elements in theta') ;
elseif length(size(Rho))>2
    error('Rho must be a vector or a matrix') ;
end

% Theta is now a column vector
theta = theta(:) ;

% Complete the data to close the circle
theta = [ theta ; theta(1) ] ;
Rho   = [ Rho ; Rho(1,:) ] ;

% A set of theta values regularly sampled over the circle
azm = (-pi:2*pi/360:pi)' ;

% Initialise the plot
% (and check if the current plot is already a "PolarPlot")

tag = get(gcf,'Tag') ;
drawAxes = ~strcmpi(tag,'PolarPlot') ;
set(gcf,'color',[1 1 1],'Tag','PolarPlot')
hold on

% Plot the data
if strcmpi(scale,'log')
    
    % Rho in dB
    Rho = 20*log10(abs(Rho)) ;
    
    % Maximum Rho
    maxRho = max(max(Rho)) ;
    
    % Maximum dB level on the plot
    maxLvl = 10*ceil(maxRho/10) ;
    
    % Apply an offset to Rho
    ofs = 40 - maxLvl ;
    Rho = Rho + ofs ;
    maxLvl = 40 ;
    
    % Crop the small values
    Rho(Rho<0) = 0 ;
    
    % Plot the circles if drawAxes is true
    if drawAxes == true
        patch(40*cos(azm),40*sin(azm),[1 1 1],'edgeColor',[0 0 0],'LineWidth',lineWidth) ;
        for I = 30:-10:10
            patch(I*cos(azm),I*sin(azm),[1 1 1], ...
                'edgeColor',[0 0 0],'linestyle',':','LineWidth',lineWidth);
        end
    end
    
    % Plot the +30deg, +60deg, +90deg etc lines and labels if drawAxes is true
    if drawAxes == true        
        for I = (0:30:150)*pi/180
            line(maxLvl*cos(I)*[1 -1],maxLvl*sin(I)*[1 -1], ...
                'color',[0 0 0],'linestyle',':','LineWidth',lineWidth) ;
        end      
        for I = 0:30:330
            text(1.1*maxLvl*cos(I*pi/180),1.1*maxLvl*sin(I*pi/180), ...
                num2str(I,'%i'),'horizontalalignment','center','FontSize',fontSize) ;
        end       
    end
    
    % Add amplitude labels
    if drawAxes == true
        text(-40*sqrt(2)/2,40*sqrt(2)/2,[num2str(40-ofs) 'dB'], ...
            'HorizontalAlignment','center', ...
            'BackgroundColor',[1 1 1],'color',[0 0 0],'FontSize',fontSize) ;
        for I = 30:-10:10
            text(-I*sqrt(2)/2,I*sqrt(2)/2,[num2str(I-ofs) 'dB'], ...
                'HorizontalAlignment','center', ...
                'BackgroundColor',[1 1 1],'color',[0 0 0],'FontSize',fontSize) ;
        end
    end

    
else
    
    % Maximum Rho
    maxRho = max(max(Rho)) ;
    
    % Number of circles (equivalent of ticks) and corresponding levels
    prvPowTen = (10^(floor(log10(maxRho)))) ;
    if maxRho / prvPowTen == 1
        tckStp = maxRho/5 ;
        maxLvl = maxRho ;
    elseif maxRho / prvPowTen <= 2.5
        tckStp = .5*prvPowTen ;
        maxLvl = tckStp*ceil(maxRho/tckStp) ;
    elseif maxRho / prvPowTen <= 5
        tckStp = 1*prvPowTen ;
        maxLvl = tckStp*ceil(maxRho/tckStp) ;
    else
        tckStp = 2*prvPowTen ;
        maxLvl = tckStp*ceil(maxRho/tckStp) ;
    end

    % Plot the circles if drawAxes is true
    if drawAxes == true
        patch(maxLvl*cos(azm),maxLvl*sin(azm),[1 1 1],'edgeColor',[0 0 0],'LineWidth',lineWidth) ;
        for I = maxLvl-tckStp:-tckStp:tckStp
            patch(I*cos(azm),I*sin(azm),[1 1 1], ...
                'edgeColor',[0 0 0],'linestyle',':','LineWidth',lineWidth);
        end
    end
    
    % Plot the +30deg, +60deg, +90deg etc lines and labels if drawAxes is true
    if drawAxes == true
        for I = (0:30:150)*pi/180
            line(maxLvl*cos(I)*[1 -1],maxLvl*sin(I)*[1 -1], ...
                'color',[0 0 0],'linestyle',':','LineWidth',lineWidth) ;
        end
        for I = 0:30:330
            text(1.1*maxLvl*cos(I*pi/180),1.1*maxLvl*sin(I*pi/180), ...
                num2str(I,'%i'),'horizontalalignment','center','FontSize',fontSize) ;
        end
    end
    
    % Add amplitude labels
    if drawAxes == true
        text(-maxLvl*sqrt(2)/2,maxLvl*sqrt(2)/2,num2str(maxLvl), ...
            'HorizontalAlignment','center', ...
            'BackgroundColor',[1 1 1],'color',[0 0 0],'FontSize',fontSize) ;
        for I = maxLvl-tckStp:-tckStp:tckStp
            text(-I*sqrt(2)/2,I*sqrt(2)/2,num2str(I), ...
                'HorizontalAlignment','center', ...
                'BackgroundColor',[1 1 1],'color',[0 0 0],'FontSize',fontSize) ;
        end
    end
    
end

% Plot Rho
if nargout < 1
    plot(Rho.*repmat(cos(theta),1,size(Rho,2)), ...
         Rho.*repmat(sin(theta),1,size(Rho,2)),varargin{:}) ;
else
    hnd = plot(Rho.*repmat(cos(theta),1,size(Rho,2)), ...
               Rho.*repmat(sin(theta),1,size(Rho,2)),varargin{:}) ;
end

% Last commands
axis equal off ;
hold off ;
