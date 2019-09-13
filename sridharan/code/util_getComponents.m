function[ componentNum, componentLengths, nodeComponents] = util_getComponents(adjMatrix)

% Takes in an adjacency matrix and returns the number of connected components in the graph, a
% vector of the lengths of all components , and a vector identifying the component
% of each node (i.e. if nodeComponents(i) == 5, node i is in component 5)

% this function does not find strongly connected compnents -- simply
% connected ones!



    
    numNodes = length(adjMatrix);

    nodeComponents(1:numNodes) = 0;
    visited(1:numNodes) = 0;

    componentNum = 0;
    componentLengths = [];
    
    while min(visited) == 0
        thisComponentLength = 0;
        componentNum = componentNum + 1;
        [value, node] = min(visited);
        to_visit = node;
        
        while numel(to_visit) ~= 0
            currentNode = to_visit(1);
            visited(currentNode) = 1;
            to_visit(1) = [];
            for i = 1:numNodes
                if (adjMatrix(currentNode, i) >= 1 || adjMatrix(i, currentNode) >= 1) && visited(i) == 0
                    to_visit = [to_visit i];
                    visited(i) = 1;
                end
            end
        end

        for i = 1:numNodes
            if visited(i) == 1
                nodeComponents(i) = componentNum;
                thisComponentLength = thisComponentLength + 1;
                visited(i) = 2;
            end
        end
        
        %if thisComponentLength > 1
            componentLengths = [componentLengths thisComponentLength];
        %end
    end
%     disp componentNum
%     disp (componentNum)
%     disp currentNode
%     disp (currentNode)
%     disp componentLengths
%     disp (componentLengths)
%     disp nodeComponents
%     disp (nodeComponents)
%     disp ('finished GetComponents');
end 