function [path,root,d] = find_path_dg(DG, node)
%
% INPUTS:
%   DG    - Directed graph object representing the rooted tree
%   node - Starting node of the path
%
% OUTPUTS:
%   path - Vector containing the nodes in the path from the starting node to the root
%   root - Final node in the path (the root of the tree)
%   d    - Weighted distance of the path
%
% DESCRIPTION:
%   This function follows a path in a rooted tree (directed graph) until the end is reached.
%   It records the path, the final node (root), and the weighted distance (d).
%
% Author: Zachary T. Berndsen

    % Initialize the path and weighted distance
    path = [];
    d = 0;

    % Get the edges and their weights
    edges = DG.Edges;
    edgeNodes = edges{:, 1};
    edgeWeights = edges{:, 2};
    
    % Find the successors of the starting node
    nextNode = successors(DG, node);

    % Follow the path through the graph
    while ~isempty(nextNode)
        % Add the current node to the path
        path = [path, node];

        % Update the weighted distance
        d = d + edgeWeights(edgeNodes(:, 1) == node);

        % Move to the next node
        node = nextNode;
        nextNode = successors(DG, node);
    end

    % Add the root node to the path
    path = [path, node];
    root = node;

end