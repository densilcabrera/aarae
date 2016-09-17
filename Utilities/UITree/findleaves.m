function leaves = findleaves
% This function mirrors part of AARAE's choose_audio, so that selection
% of audio from AARAE's workspace can be done without a dialog box.
handles = guidata(findobj('Tag','aarae'));
root = handles.root; % Get selected leaf
root = root(1);
first = root.getFirstChild;
nbranches = root.getChildCount;
branches = cell(nbranches,1);
branches{1,1} = char(first.getValue);
nleaves = 0;
nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{1,1}))(1).getChildCount;
next = first.getNextSibling;
for n = 2:nbranches
    branches{n,1} = char(next.getValue);
    nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{n,1}))(1).getChildCount;
    next = next.getNextSibling;
end
if nleaves == 0, return; end
leaves = cell(nleaves,1);
i = 0;
for n = 1:size(branches,1)
    currentbranch = handles.(matlab.lang.makeValidName(branches{n,1}));
    if currentbranch.getChildCount ~= 0
        first = currentbranch.getFirstChild;
        data = first(1).handle.UserData;
        if isfield(data,'audio')
            i = i + 1;
            %leafnames(i,:) = first.getName;
            leaves{i,:} = char(first.getValue);
        end
        next = first.getNextSibling;
        if ~isempty(next)
            data = next(1).handle.UserData;
            for m = 1:currentbranch.getChildCount-1
                if isfield(data,'audio')
                    i = i + 1;
                    %leafnames(i,:) = next.getName;
                    leaves{i,:} = char(next.getValue);
                    next = next.getNextSibling;
                    if ~isempty(next)
                        data = next(1).handle.UserData;
                    end
                end
            end
        end
    end
end
if nleaves ~=0
    leaves = leaves(~cellfun(@isempty,leaves));
else
    leaves = [];
end