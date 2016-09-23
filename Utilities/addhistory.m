function signaldata = addhistory(signaldata,str,varargin)
% adds a time stamp, string and function callback to the history field of
% the input structure
options = varargin;
if nargin == 1, str = ''; end
   row = cell(1,4);
   row{1,1} = datestr(now);
    row{1,2} = str;
    if isfield(signaldata,'funcallback') && (checkinput(options,'callback') || checkinput(options,'all'))
        if isfield(signaldata.funcallback,'name') && isfield(signaldata.funcallback,'inarg')
            row{1,3} = signaldata.funcallback.name;
            row{1,4} = signaldata.funcallback.inarg;
        end
    end

    if isfield(signaldata, 'audio') && (checkinput(options,'audio') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'audio size';
        newrow{1,4} = [num2str(size(signaldata.audio,1)) 'x' num2str(size(signaldata.audio,2)) 'x' ...
            num2str(size(signaldata.audio,3)) 'x' num2str(size(signaldata.audio,4)) 'x' ...
            num2str(size(signaldata.audio,5)) 'x' num2str(size(signaldata.audio,6))];
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'fs') && (checkinput(options,'fs') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'fs';
        newrow{1,4} = signaldata.fs;
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'audio2') && (checkinput(options,'audio2') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'audio2 size';
        newrow{1,4} = [num2str(size(signaldata.audio,1)) 'x' num2str(size(signaldata.audio,2)) 'x' ...
            num2str(size(signaldata.audio,3)) 'x' num2str(size(signaldata.audio,4)) 'x' ...
            num2str(size(signaldata.audio,5)) 'x' num2str(size(signaldata.audio,6))];
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'chanID') && (checkinput(options,'chanID') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'chanID';
        newrow{1,4} = signaldata.chanID;
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'bandID') && (checkinput(options,'bandID') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'bandID';
        newrow{1,4} = signaldata.bandID;
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'cal') && (checkinput(options,'cal') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'cal';
        newrow{1,4} = signaldata.cal;
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'datatype') && (checkinput(options,'datatype') || checkinput(options,'all'))
        newrow = cell(1,4);
        newrow{1,3} = 'datatype';
        newrow{1,4} = signaldata.datatype;
        row = [row;newrow];
    end
    
    if isfield(signaldata, 'properties') && (checkinput(options,'properties') || checkinput(options,'all'))
        names = fieldnames(signaldata.properties);
        nfields = size(names,1);
        newrow = cell(nfields,4);
        for n = 1:nfields
            newrow{n,3} = ['properties.' names{n}];
            [~, newrow{n,4}] = evalc(['signaldata.properties.' names{n}]);
        end
        row = [row;newrow];
    end
    
    
    if ~isfield(signaldata,'history')
        signaldata.history = row;
    else
        signaldata.history = [signaldata.history;row];
    end
    
end

function out = checkinput(varcell,str)
% maybe this is more straightforward that inputParser
out = false;
    for n = 1:length(varcell)
        if strcmp(varcell{n},str)
            out = true;
            return
        end
    end
end