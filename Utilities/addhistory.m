function signaldata = addhistory(signaldata,str)
% adds a time stamp, string and function callback to the history field of
% the input structure
if nargin == 1, str = ''; end
   row = cell(1,4);
   row{1,1} = datestr(now);
    row{1,2} = str;
    if isfield(signaldata,'funcallback')
        if isfield(signaldata.funcallback,'name') && isfield(signaldata.funcallback,'inarg')
            row{1,3} = signaldata.funcallback.name;
            row{1,4} = signaldata.funcallback.inarg;
        end
    end
if ~isfield(signaldata,'history')
    signaldata.history = row;
else
    signaldata.history = [signaldata.history;row];
end
        
