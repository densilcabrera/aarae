function closefigures(category)
% Closes figures except for the AARAE GUI
% Use no input argument to close all figures.
% Use an input argument of 'T' to only close tables.
% Use an input argiment of 'C' to only close non-tables (charts)
%
% This function may be useful in writing workflows in cases where you do
% not want any previous figures to be preserved.
h = findobj('type','figure','-not','tag','aarae');
if nargin == 0 % delete all figures
    delete(h)
    return
end

dodelete = false(length(h),1);

% delete tables only
if strcmp(category,'T') || strcmp(category,'t')
    for i = 1:length(h)
        Figtag = get(h(i),'Tag');
        if ~isempty(Figtag)
            if strcmp(Figtag,'AARAE table')
                 dodelete(i) = true;
            end
        end
    end
% delete non-tables (charts) only
elseif strcmp(category,'C') || strcmp(category,'c')
    for i = 1:length(h)
        Figtag = get(h(i),'Tag');
        if ~isempty(Figtag)
            if ~strcmp(Figtag,'AARAE table')
                 dodelete(i) = true;
            end
        else
            dodelete(i) = true;
        end
    end
end

delete(h(dodelete));
    
   