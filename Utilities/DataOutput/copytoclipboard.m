function copytoclipboard
% copytoclipboard
%   copytoclipboard is a function called as an extension of disptables.m
%   that enables data export of the data displayed in the uitables
%   contained in the figures adjusted with the aforementioned function.
%
%   See also disptables
%
%   Code by Daniel R. Jimenez for the AARAE project. (28 March 2014)

fighandle = gcf;
children = findobj(fighandle,'Type','uitable');

bigone = [];
for i = 1:length(children)
    clear header
    clear table
    child = get(children(i));
    RowName = child.RowName;
    ColName = child.ColumnName;
    Data = child.Data;

    ColName = reshape(ColName',1,[]);
    [m n] = size(ColName);
    header(1:m,1:2:2*n) = {sprintf('\t')};
    header(1:m,2:2:2*n) = ColName;
    Data = num2cell(Data);
    [m n]  = size(Data);
    for k = 1:m
        for j = 1:n
            Data{k,j} = num2str(Data{k,j});
        end
    end
    RowandData = [RowName Data];
    [m n] = size(RowandData);
    table(1:m,1:2:2*n) = RowandData;
    table(1:m,2:2:2*n) = {sprintf('\t')};
    table(1:m,end) = {sprintf('\n')};
    table = reshape(table',1,[]);
    table = [header {sprintf('\n')} table {sprintf('\n')}];
    bigone = [bigone table];
end
t = sprintf('%s', bigone{:});
clipboard('Copy',t)
msgbox('Figure data copied to clipboard!','AARAE info');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2013,2014, Daniel Jimenez
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the University of Sydney nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%