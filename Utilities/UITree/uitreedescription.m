function uitreedescription
% this function identifies and displays (in the command window) the parts
% of the uitree, for debug

st = dbstack;
[~,filename]=st.file;
[~,funname] = st.name;
[~,line] = st.line;
disp(['uitreedescription from ' filename ', running ' funname ', line ' num2str(line)])

handles = guidata(findobj('Tag','aarae'));

root = handles.root; % Get selected leaf
root = root(1);
disp(['root: ' char(root.getName)])

first = root.getFirstChild;
% disp(['first child: ' char(root.getFirstChild.getName)])

nbranches = root.getChildCount;
branches = cell(nbranches,1);
branches{1,1} = char(first.getValue);
nleaves = 0;
nbranchleaves = zeros(nbranches,1);
nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{1,1}))(1).getChildCount;
nbranchleaves(1) = nleaves;
next = first.getNextSibling;
for n = 2:nbranches
    branches{n,1} = char(next.getValue);
    nleaves = nleaves + handles.(matlab.lang.makeValidName(branches{n,1}))(1).getChildCount;
    nbranchleaves(n) = handles.(matlab.lang.makeValidName(branches{n,1}))(1).getChildCount;
    next = next.getNextSibling;
end
disp(['Number of branches: ' num2str(nbranches)])
disp(['Number of leaves: ' num2str(nleaves)])
leaves = cell(nleaves,1);
i = 0;
for n = 1:size(branches,1)
    currentbranch = handles.(matlab.lang.makeValidName(branches{n,1}));
    if currentbranch.getChildCount ~= 0
        i = i + 1;
        first = currentbranch.getFirstChild;
        disp(['Leaves in branch ' branches{n,1} ':'])
        
        %leafnames(i,:) = first.getName;
        leaves{i,:} = char(first.getValue);
        disp(['   Leaf ' num2str(i) ': ' leaves{i,:}])
        next = first.getNextSibling;
        if ~isempty(next)
            for m = 1:currentbranch.getChildCount-1
                i = i + 1;
                %leafnames(i,:) = next.getName;
                leaves{i,:} = char(next.getValue);
                disp(['   Leaf ' num2str(i) ': ' leaves{i,:}])
                next = next.getNextSibling;
            end
        end
    end
end


selectedNodes = handles.mytree.getSelectedNodes;
for i = 1:length(selectedNodes);
    disp(['selected node ' num2str(i) ': ' matlab.lang.makeValidName(char(selectedNodes(i).getName))])
end