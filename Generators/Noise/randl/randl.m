function y = randl(varargin)
%RANDL Laplacian distributed pseudorandom numbers.
%   R = RANDL(N) returns an N-by-N matrix containing pseudorandom values drawn
%   from the laplacian distribution.  RANDL(M,N) or RANDL([M,N]) returns
%   an M-by-N matrix. RANDL(M,N,P,...) or RANDN([M,N,P,...]) returns an
%   M-by-N-by-P-by-... array. RANDL returns a scalar.  RANDL(SIZE(A)) returns
%   an array the same size as A.
%
%   Note: The size inputs M, N, P, ... should be nonnegative integers.
%   Negative integers are treated as 0.
%
%
%   Examples:
%
%      Example 1: Generate values from a laplacian distribution with mean 1
%       and standard deviation 2.
%         r = 1 + 2.*randl(100,1);
%
%      Example 2: Generate values from a bivariate laplacian distribution with
%      specified mean vector and covariance matrix.
%         mu = [1 2];
%         Sigma = [1 .5; .5 2]; R = chol(Sigma);
%         z = repmat(mu,100,1) + randl(100,2)*R;
%
x = rand(varargin{:});
y = sign(0.5-x).*(1/sqrt(2)).*log(2*min(x,1-x));


