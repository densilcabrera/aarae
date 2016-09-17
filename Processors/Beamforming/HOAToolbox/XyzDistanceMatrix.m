function D = XyzDistanceMatrix(X,Y)
%XYZDISTANCEMATRIX Matrix of the distances between two sets of points
%
%  D = XyzDistanceMatrix(X,Y) returns the matrix of the distances between
%  every point defined by matrix X and every point defined by matrix Y.
%  
%  The sets of points must be defined in the same space, thus X and Y must
%  have the same number of columns, for instance,
%  X = [ 0 0 0 ] and Y = [ 1 0 0 ; 0 1 0 ; 0 0 1 ]
%
%  D is a M x N matrix, where M is the number of points in X and N is the 
%  number of points in Y.

%  N.Epain, 2010

% Number of columns in X and Y
P = size(X,2) ;
Q = size(Y,2) ;

% X and Y MUST have the same number of columns
if P ~= Q
    error('X and Y must have the same number of columns') ;
end

% Number of points in X and Y
M = size(X,1) ;
N = size(Y,1) ;

% Compute the distance matrix D
D = sqrt(sum(bsxfun(@minus,reshape(X,M,1,P),reshape(Y,1,N,Q)).^2,3)) ;

end