function C = findConnectedComponents(K)
%FINDCONNECTEDCOMPONENTS interprets K as the adjacency matrix of a
%bipartite graph and finds the connected components in that one
%   returns only the indices of the rows
%
% The bipartite graph contains a node for each columns and row. There is an edge between column i and row j if K(j,i) != 0. 
%
% Input:
%		K		a matrix (values != 0 are considered true, others are considered false). Zero-Test is done using an epsilon of 1e-8
%
% Output:
%		C		a logical matrix encoding the connected components. Each column encodes a connected component. 
%				C(i,j) = true if connected component i contains row j (of K) and false otherwise.


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




% epsilon to test if entries of K are equal to zero
eps = 1e-8;

[m,n] = size(K);

% build logical matrix from K
K2 = K < -eps | K > eps;

% build ordinary adjacency matrix from K2 (remember that K2 encodes a bipartite graph)
A = [zeros(m), K2; K2', zeros(n)];
% observe that A is symmetric.

% f stores the rows that we have already found
f = zeros(m,1);

% C stores the components that we found so far
C = false(m,0);

% j counts the connected components that we found
j=0;

for i=1:m
    if f(i) == 0
    		% we did not find row i yet
        if any(A(i,:)) 
        		% we actually have incidencies
            j = j+1;
            found = findConnectedComponent(A,i); % start search from node i
            C(:,j) = found(1:m) ~= 0; 
            f = f + found(1:m); % update the list of nodes that we have already seen
        else
        		% trivial connected component
            f(i) = 1;
            j = j+1;
            C(:,j) = false(m,1);
            C(i,j) = true;
        end
    end
end

end

function found = findConnectedComponent(A,i)
%FINDCONNECTEDCOMPONENT interprets K as the adjacency matrix of a
%bipartite graph and finds the connected component containing i in that one


m = size(A,1);

% we run a DFS using a stack (recursion is a bad idea in matlab)
% we assume that A is symmetric.

found = zeros(m, 1);

% todo stores a stack, todos is the current size
todo = [i];
todos = 1;
found(i) = 1;

count = 1;

while todos > 0
  k = todo(todos);
  todos = todos-1;
  for j = find(A(k,:) ~= 0)
    if found(j) == 0
      todos = todos+1;
      todo(todos) = j;
      count = count+1;
      found(j) = count;
    end
  end
end

end
