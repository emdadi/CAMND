function [K, B] = findCircuitBase( S )
%FINDCIRCUITBASE Given a matrix S, computes a base of null(S) of circuits.
%
% Input:
%		S				a real matrix
%
% Output:
%		K				a null space matrix of S where each column is a fundamental circuit of the basis B
%		B				a basis of S (maximal set of linearly independent columns of S). The columns are specified by a logical index array.


% Copyright (C) 2013  Arne Müller <arne.mueller@fu-berlin.de>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.




[m,n] = size(S); % m = # metabolites, n = # reactions

% initialize empty nullspace matrix (we add columns when we find them)
K = sparse(n,0);

% if the matrix is too big, don't check if the rank computed by matlab coincides with our result
if size(S,1) > 3000
    nulldim = -1;
else
    nulldim = size(S,2) - rank(full(S));
    fprintf('Computing %d short linearly independent flux vectors\n     ',nulldim);
end

% initially, there is no element in the basis
B = false(n,1);
% i counts the number of columns added to K
i = 1;

for r=1:n
    fprintf('\b\b\b\b\b%5d',r);
    C = S(:,B) \ S(:,r);
    qual = norm(S(:,B)*C - S(:,r),inf);
    if qual < 1e-8
        % consistent
        K(B,i) = -C;
        K(r,i) = 1;
        i = i+1;
    else
        % infeasible
        B(r) = true;
    end
end

% each element either entered the basis or induced a column in K
assert(n-sum(B) == size(K,2));

% warn if the results don't coincide with the expected result (in terms of matrix rank).
if nulldim == -1
		% For big matricies, we don't want to rely on the matlab rank function. It may produce false negative results.
    warning('input matrix was too big, rank of null space is not checked');
elseif size(K,2) > nulldim
		% we found more columns. This is not fatal. It probably just means that some columns of K are nearly linearly dependent.
    warning('computed solution has not same dimension as null space: %d expected: %d!',size(K,2), nulldim);
elseif size(K,2) < nulldim
		% we found less columns. This can be fatal, because the null-space is smaller. We may have missed an important flux mode.
    error('computed solution has not same dimension as null space: %d expected: %d!',size(K,2), nulldim);
end

end

