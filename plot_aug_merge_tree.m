function plot_aug_merge_tree(DG)
% PLOT_AUG_MERGE_TREE - Plot an augmented merge tree as a directed graph
%
% INPUTS:
%   - DG: Directed graph object representing the merge tree
%
% OUTPUTS:
%   None
%
% This function plots the merge tree as a directed graph, highlighting
% birth nodes, mergers, and path nodes with different colors. Nodes with ID
% > 3 are made larger. The highest intensity nodes are on the left and the
% lowest on the right.
%
% Author: Zachary T. Berndsen


    % Compute in-degree and out-degree of nodes
    ID = indegree(DG);
    OD = outdegree(DG);
    
    % Identify birth, merger, and path nodes
    Births = find(ID == 0);
    Mergers = find(ID > 1);
    Paths = find(OD == 1 & ID == 1);

    % Create plot and highlight nodes
    figure()
    p = plot(DG, 'layout', 'layered', 'Direction', 'right');
    highlight(p, Births, 'Marker', 'o', 'MarkerSize', 4, 'NodeColor', [0.9 0 0])
    highlight(p, Mergers, 'Marker', 'o', 'MarkerSize', 4, 'NodeColor', [0.9 0 0.9])
    highlight(p, Paths, 'Marker', 'o', 'MarkerSize', 4, 'NodeColor', [0 0.9 0])
    highlight(p, find(ID > 2), 'Marker', 'o', 'MarkerSize', 6)

end
