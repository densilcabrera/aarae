function pixels = get_axes_width(h)

% pixels = get_axes_width(h)
% 
% Returns the width of the axes object, h, in pixels.
%
% Tucker McClure
% Copyright 2013, The MathWorks, Inc.


    % Record the current axes units setting.
    axes_units = get(h, 'Units');

    % Change axes units to pixels.
    set(h, 'Units', 'pixels');

    % Get axes width in pixels.
    axes_position = get(h, 'Position');
    pixels = round(axes_position(3));

    % Return the units.
    set(h, 'Units', axes_units);
    
end
