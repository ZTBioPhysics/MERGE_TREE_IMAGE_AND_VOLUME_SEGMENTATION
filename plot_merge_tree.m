function plot_merge_tree(DG, THRESH2)
% PLOT_MERGE_TREE Plots a merge tree as a directed graph in two differnt 
% layouts.
%
% INPUTS:
%   - DG: Directed graph object representing the merge tree
%   - THRESH2: Vector of thresholds for the merge tree
%
% OUTPUTS:
%   None
%
% This function plots a merge tree as a directed graph. It colors birth
% nodes, mergers, and path nodes with different colors. Nodes with ID > 3
% are made larger. The highest intensity nodes are on the left, and the
% lowest intensity nodes are on the right.
%
% Author: Zachary T. Berndsen

    % Remove root node for force layout
    DG2 = rmnode(DG, numnodes(DG));
    
    % Create force layout plot
    figure();
    LWidths = 5 * DG2.Edges.Weight / max(DG2.Edges.Weight);
    p = plot(DG2, 'Layout', 'force', 'UseGravity', false, 'LineWidth', LWidths, 'MarkerSize', 5);
    p.NodeCData = THRESH2(1:end-1);
    cbh = colorbar;
    cbh.Label.String = 'Threshold';
    path1 = shortestpath(DG2, 1, numnodes(DG2));
    highlight(p, path1, 'EdgeColor', 'g');
    
    % Create layered layout plot
    figure();
    LWidths = 5 * DG.Edges.Weight / max(DG.Edges.Weight);
    p = plot(DG, 'Layout', 'layered', 'Direction', 'up', 'LineWidth', LWidths, 'MarkerSize', 5);
    p.NodeCData = THRESH2;
    cbh = colorbar;
    cbh.Label.String = 'Threshold';
    path1 = shortestpath(DG, 1, numnodes(DG));
    highlight(p, path1, 'EdgeColor', 'g');
end
