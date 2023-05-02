function [branches, BD] = branch_decomposition(DG, THRESH, toplot)
% BRANCH_DECOMPOSITION Decomposes a merge tree into branches using elder rule
%
% SYNTAX:
%   [branches, BD] = branchDecomposition(DG, THRESH, toplot)
%
% INPUTS:
%   DG      - Directed graph object representing the merge tree
%   THRESH  - Vector of threshold values for every node in the merge tree
%   toplot  - Boolean flag to create a plot of the branch decomposition (optional)
%
% OUTPUTS:
%   branches - Cell array containing the branches of the branch decomposition
%   BD       - Directed graph object representing the merge tree after branch decomposition
%
% DESCRIPTION:
%   This function assigns pairs of nodes as persistence pairs using the elder rule.
%   It starts by sorting all leaf nodes in descending order based on their threshold
%   values, then loops through every leaf node and finds the path from that node to the
%   current root of the tree and assigns those two nodes as persistence pairs. In the
%   case of degenerate merge nodes (merge nodes with ID > 2), the merge node is replicated
%   so that every input branch gets paired.
%
% Author: Zachary T. Berndsen

    % Check if toplot is provided, otherwise set it to 0 (do not plot)
    if nargin < 3
        toplot = 0;
    end
    
    BD = DG;
    ID = indegree(DG);
    leafNodes = find(ID == 0);
    
    % Sort the leaf nodes by the function value in descending order
    [~, sortedIndex] = sort(THRESH(ID == 0), 'descend');
    sortedLeafNodes = leafNodes(sortedIndex);
    
    % Initialize the branches cell array
    branches = cell(length(leafNodes), 2);
    
    if toplot == 1
        % Initiate plot
        figure();
        LWidths = 5 * BD.Edges.Weight / max(BD.Edges.Weight);
        p = plot(BD, 'Layout', 'layered', 'Direction', 'up', 'LineWidth', LWidths, 'MarkerSize', 5);
    end
    
    % Iterate over the sorted leaf nodes
    for i = 1:length(sortedLeafNodes)
        
        % Current birth node
        node = sortedLeafNodes(i);
        
        % Find the path from the root to the current leaf node
        [path, ~, pl] = find_path_dg(BD, node);
    
        if toplot == 1
            highlight(p, path, 'EdgeColor', 'g');
            pause
            highlight(p, path, 'EdgeColor', 'r');
        end
        
        % Add the path to the branches cell array
        branches{i, 1} = path;
        branches{i, 2} = pl;
    
        % Remove all edges in the path
        for j = 1:length(path) - 1
            BD = rmedge(BD, path(j), path(j + 1));
        end
    end
end
