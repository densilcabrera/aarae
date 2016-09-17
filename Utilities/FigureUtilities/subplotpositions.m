function [r, c] = subplotpositions(number_of_subplots, aspect)
% returns the number of rows and colums of subplots to lay out the
% required number of subplots
%
% Code by Densil Cabrera 2013
%
% number_of_subplots is the number of subplots to be generated
%
% aspect is a number that controls the aspect ratio:
%  0.5 yields an approximately square layout (more columns than rows when a
%    square is not possible).
% A smaller value (e.g. 0.3) yields more columns than rows.
% A value of 0 always yields 1 row.
% A value larger than 0.5 (e.g. 0.8) yields more rows than columns.
% A value of 1 always yields 1 column.
%

r = floor(number_of_subplots^aspect); % number of subplot rows
c = ceil(number_of_subplots/r); % number of subplot columns