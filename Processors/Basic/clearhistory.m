function OUT = clearhistory(IN)
% This function replaces the entire history field with a row indicating
% 'Cleared history'
OUT = IN;
OUT.history = cell(1,4);
OUT.history{1,1} = datestr(now);
OUT.history{1,2} = 'Cleared history';
OUT.funcallback.name = 'clearhistory.m';
OUT.funcallback.inarg = {};
end
    