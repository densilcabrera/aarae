function identifyselectednodes
% function for debugging in relation to uitree selection

handles = guidata(findobj('Tag','aarae'));
selectedNodes = handles.mytree.getSelectedNodes;
st = dbstack;
[~,filename]=st.file;
[~,funname] = st.name;
[~,line] = st.line;
disp(['identifyselectednodes from ' filename ', running ' funname ', line ' num2str(line)])
for i = 1:length(selectedNodes);
    disp(['selected node ' num2str(i) ': ' matlab.lang.makeValidName(char(selectedNodes(i).getName))])
end