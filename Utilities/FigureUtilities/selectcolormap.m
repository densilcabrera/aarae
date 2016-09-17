function c = selectcolormap(choice, parameter)
% A library of custom colormaps available for use within the AARAE project.
%
% The second input argument (parameter) is defined in relation to
% particular colormap choices (and is not necessarily used).
%
% You can use set(fig,'Colormap', c) to apply a colormap to the current
% figure.

switch choice
    case 'bipolar1'
        % use this colormap to plot bipolar data (where the two extremes
        % are the most interesting values)
        grey = 0.2;
        c = zeros(256,3);
        c(130:256,1) = (1:127)./127;
        c(:,2) = (-1:1/127.5:1).^2;
        c(1:127,3) = flipud((1:127)'./127);
        c = c * (1-grey) + grey;
    case 'red'
        % a monochrome red colormap
        ncolors = 256;
        c = zeros(ncolors,3);
        c(:,1) = 1-(((1:ncolors)-1)/(ncolors*8));
        c(:,2) = 1-(((1:ncolors)-1)/(ncolors*1.1));
        c(:,3) = 1-(((1:ncolors)-1)/(ncolors*1.1));

    case 'blue'
        % a monochrome blue colormap
        ncolors = 256;
        c = zeros(ncolors,3);
        c(:,1) = 1-(((1:ncolors)-1)/(ncolors*1.1));
        c(:,2) = 1-(((1:ncolors)-1)/(ncolors*1.1));
        c(:,3) = 1-(((1:ncolors)-1)/(ncolors*8));
    case 'rainbow7'        
        % a 7 discrete level rainbow colormap (might be useful for 7 octave
        % bands, for example). Yellow, green and cyan have been darkened
        % relative to the other colors (otherwise they look too bright).
        c = [255, 0, 0; ... % red
        255, 128, 0; ... % orange
        204, 204, 0; ... % dark yellow
        0, 204, 0; ... % mid green
        0, 204, 204; ... % dark cyan
        0, 0, 255; ... % blue
        127, 0, 255]; % violet
        c = c / 255; % rescale to 0-1 range
    case 'rainbow256'
        % 8-bit pseudo-continuous RGB rainbow using the same colors as
        % rainbow7. It could be useful for signifying the frequency scale.
        r = [255*ones(43,1); ...
            linspace(255,204,43)'; ...
            linspace(199,0,43)'; ...
            zeros(84,1); ...
            linspace(0,127,43)'];
        g = [linspace(0,125,43)'; ...
            linspace(128,204,43)'; ...
            204*ones(84,1); ...
            linspace(204,0,44)'; ...
            zeros(42,1)];
        b = [zeros(128,1); ...
            linspace(0,199,42)'; ...
            linspace(204,255,44)'; ...
            255*ones(42,1)];
        c = [r, g, b] / 255;
    otherwise
        % default Matlab colormap
        c = 'jet';
          
end
        