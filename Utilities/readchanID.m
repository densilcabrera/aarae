function [out,format] = readchanID(chanID)
% This aarae utility function extracts the numerical data from chanIDs,
% assuming that they are in one of the formats supported by aarae's 
% makechanID. 
% 
% The following formats are supported:
% 0: numbered channels in the form Chan1;Chan2;...
% 1: spherical harmonic order and degree, Y 0 0; Y 1 1; Y 1 -1; Y 1 0;...
% 2: polar coordinates using degrees, e.g., 90 deg,  45 deg, 1.4142 m
% 3: polar coordinates using radians, e.g., 0 rad,  1.3 rad, 1.4142 m
% 4: Cartesian coordinates using metres, e.g., 1 m, 3 m, -4 m
%
% Note that the chanID formats are cell arrays of strings, in a column.
%
% See the code for more information on the imputs and outputs.


if ~isempty(regexp(char(chanID{1}),'InChan','once'))
    % generic chanID in the form of 'InChan1', 'InChan2', etc
    out = zeros(length(chanID),1);
    for n = 1:length(chanID)
         cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
         if size(cellval,2)~=1
            out = [];
            return
        end
         out(n) = str2double(cellval{1});
    end
    format = 11;
    
elseif ~isempty(regexp(char(chanID{1}),'Chan','once'))
    % generic chanID in the form of 'Chan1', 'Chan2', etc
    out = zeros(length(chanID),1);
    for n = 1:length(chanID)
         cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
         if size(cellval,2)~=1
            out = [];
            return
        end
         out(n) = str2double(cellval{1});
    end
    format = 0;
    
elseif ~isempty(regexp(char(chanID{1}),'OutCh','once'))
    % generic dimension 5 (output channel) ID in the form of 'OutCh1', 'OutCh2', etc
    out = zeros(length(chanID),1);
    for n = 1:length(chanID)
         cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
         if size(cellval,2)~=1
            out = [];
            return
        end
         out(n) = str2double(cellval{1});
    end
    format = 10;
    
elseif ~isempty(regexp(char(chanID{1}),'Dim6_','once'))
    % generic dimension 6 ID in the form of 'Dim6_1', 'Dim6_2', etc
    out = zeros(length(chanID),1);
    for n = 1:length(chanID)
         cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
         if size(cellval,2)~=1
            out = [];
            return
        end
         out(n) = str2double(cellval{1});
    end
    format = 20;
    
elseif ~isempty(regexp(char(chanID{1}),'Y','once'))
    % Spherical harmonic chanIDs, showing order and degree (e.g., for 
        % HOA signals)
    out = zeros(length(chanID),2);
    for n = 1:length(chanID)
        cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
        if size(cellval,2)~=2
            out = [];
            return
        end
        out(n,:) = [str2double(cellval{1,1}),str2double(cellval{1,2})];
    end
    format = 1;
    
elseif ~isempty(regexp(char(chanID{1}),'deg','once'))
    % chanID using spherical coordinates in degrees (az,elev,radius)
    % but output of this function is always Cartesian
    out = zeros(length(chanID),3);
    for n = 1:length(chanID)
        cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
        if size(cellval,2)~=3
            out = [];
            return
        end
        out(n,:) = [str2double(cellval{1,1}),str2double(cellval{1,2}) ...
            ,str2double(cellval{1,3})];
    end
    out(:,1:2) = out(:,1:2) * pi/180;
    [out(:,1),out(:,2),out(:,3)] = sph2cart(out(:,1),out(:,2),out(:,3));
    format = 2;
    
elseif ~isempty(regexp(char(chanID{1}),'rad','once'))
    % chanID using spherical coordinates in radians (az,elev,radius)
    % but output of this function is always Cartesian
    out = zeros(length(chanID),3);
    for n = 1:length(chanID)
        cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
        if size(cellval,2)~=3
            out = [];
            return
        end
        out(n,:) = [str2double(cellval{1,1}),str2double(cellval{1,2}) ...
            ,str2double(cellval{1,3})];
    end    
    [out(:,1),out(:,2),out(:,3)] = sph2cart(out(:,1),out(:,2),out(:,3));
    format = 3;
    
elseif ~isempty(regexp(char(chanID{1}),'m','once'))
    % chanID using Cartesian coordinates (x,y,z)
    out = zeros(length(chanID),3);
    for n = 1:length(chanID)
        cellval = regexp(chanID{n},'-?\d+\.?\d*|-?\d*\.?\d+','match');
        if size(cellval,2)~=3
            out = [];
            return
        end
        out(n,:) = [str2double(cellval{1,1}),str2double(cellval{1,2}) ...
            ,str2double(cellval{1,3})];
    end
    format = 4;
    
else
    out = [];
    format = [];  
end