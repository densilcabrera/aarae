function data = historytable(audiodata,nofigure)
% This function parses the history field of AARAE structures.
% The history field is a cell array consisting of 4 columns and any number
% of rows.
% It is a record of events that occured such as generation, loading
% processing, recording, etc. The generic row written at events such as
% these has the datetime string in column 1, the general description of the
% action in column 2 (string), the name of the function (or the name of the
% function in the current callback) in column 3 (string), and specific data
% in column 4. In many cases, the fourth column is a cell array, which may
% or may not follow this generic format.
% Some detail may not be converted in the historytable function.
%
% The output is done via a figure (formatted using AARAE's disptable
% function) which has 'copy-to-clipboard' functionality. The nofigure
% argument can be included (with any content whatsoever) to prevent the
% figure from being generated, and the output argument (data) used instead.
% Data is a cell array of characters.
%
% code by Densil Cabrera (Sept 2016)

if ~isfield(audiodata,'history')
    return
end
mainrows = size(audiodata.history,1);

data = {};
for r = 1:mainrows
    data1 = audiodata.history(r,:);
    data2 = {};
    data3 = {};
    data4 = {};
    if ~isempty(data1)
        if size(data1,1)==1 && ischar(data1{1,1}) && ischar(data1{1,2}) ...
                && (iscell(data1{1,4}) || ischar((data1{1,4})))
            [data1, data2] = expandcell(data1);
            if ~isempty(data2)
                if size(data2,1)==1 && ischar(data2{1,1}) && ischar(data2{1,2}) ...
                        && (iscell(data2{1,4}) || ischar((data2{1,4})))
                    [data2, data3] = expandcell(data2);
                    if ~isempty(data3)
                        if size(data3,1)==1 && ischar(data3{1,1}) && ischar(data3{1,2}) ...
                                && (iscell(data3{1,4}) || ischar((data3{1,4})))
                            [data3, data4] = expandcell(data3);
                        end
                    end
                end
            end
        end
    end
data = [data;data1;data2;data3;data4];
end
tablesize = size(data);
% only strings allowed in the table
for r = 1:tablesize(1)
    for c = 1:tablesize(2)
        if isempty(data{r,c})
            data{r,c} = '';
        elseif iscell(data{r,c})
            data{r,c} = ['{', num2str(size(data{r,c},1)), 'x',...
            num2str(size(data{r,c},2)), ' cell}'];
        elseif isstruct(data{r,c})
            data{r,c} = 'struct';
        elseif isnumeric(data{r,c}) && numel(data{r,c})>1
            if numel(data{r,c})<101
                data{r,c} = num2str(data{r,c}(:)');
            else
                data{r,c} = ['[', num2str(size(data{r,c},1)), 'x',...
                    num2str(size(data{r,c},2)), ' num]'];
            end
        elseif isnumeric(data{r,c})
            data{r,c} = num2str(data{r,c});
        elseif ~ischar(data{r,c})
            data{r,c} = '?';
        elseif size(data{r,c},1) > 1
            data{r,c} = ['[', num2str(size(data{r,c},1)), 'x',...
            num2str(size(data{r,c},2)), ' char]'];
        end
    end
end
if ~exist('nofigure','var')
f = figure('Name','AARAE history table');
t = uitable(f,'Data',data,'ColumnWidth',{200},...
    'ColumnName',{'A','B','C','D'},...
    'RowName',cellstr(num2str((1:size(data,1))')));
disptables(f,t);
end
end %eof



function [outdata, outcell] = expandcell(in)
outdata = in;
outcell = {};
if ischar(in{1,4}) && strcmp(in{1,2},'COMMENT') % special case: comment
    cellsize = size(in{1,4},1);
    outdata = [outdata;cell(cellsize-1,4)];
    data = outdata{1,4};
    for n = 1:cellsize
        outdata{n,4} = char(data(n,:));
    end
elseif iscell(in{1,4})
    cellsize = size(in{1,4});
    outdata{1,4} = ['{', num2str(cellsize(1)), 'x',...
        num2str(cellsize(2)), ' cell}'];
    if cellsize(2) == 4
        outcell = in{1,4};
    elseif cellsize(1) == 1 % convert row of arbitrary length to a column
        data = cell(cellsize(2),4);
        data(:,4) = in{1,4}';
        outdata = [outdata;data];
    elseif cellsize(2)==2
        data = cell(cellsize(1),4);
        data(:,3:4) = in{1,4};
        outdata = [outdata;data];
    elseif cellsize(2)==3
        data = cell(cellsize(1),4);
        data(:,2:4) = in{1,4};
        outdata = [outdata;data];
    elseif cellsize(2)==4
        outdata = [outdata;in{1,4}];
    end
end
end

