function logtext(inputstring)
% This function writes text to AARAE's logfile
handles = guidata(findobj('Tag','aarae')); 
fprintf(handles.fid,inputstring);