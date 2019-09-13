function hasCycle  = util_hasCycle (matrix)

%
%This function takes a matrix, and returns true if the graph representing
%the matrix has a cycle.

% The a simple way to test whether a directed graph is cyclic is to attempt a topological traversal of its vertices.
% If all the vertices are not visited, the graph must be cyclic.

% At each step of the traversal, a vertex with in-degree of zero is visited.
% After a vertex is visited, the vertex and all the edges emanating from that vertex are removed from the graph.


%n = size(matrix,2);

hasCycle = 0;

while (size(matrix,2) ~= 1 && ~hasCycle)
    
    % find the incoming degrees of all vertices:
    
    noIncomingEdges = sum(matrix);
    [minValue,minValueIndex] = min(noIncomingEdges); % the index has the first of such minValues
    
    if (minValue >0)
        hasCycle = 1;
    else
        % remove row/column minValueIndex:
        matrix(minValueIndex,:) = [];  % remove row
        matrix(:,minValueIndex) = [];  % remove column
    end
    
end
end
