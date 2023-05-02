function [DG3, THRESH3, PIXID3] = simplify_merge_tree(DG,branches,PLs,THRESH,PIXID,plco,toplot)
% SIMPLIFY_MERGE_TREE Simplify the merge tree by removing short branches
%
% SYNTAX:
%   [DG3, THRESH3, PIXID3] = simplify_merge_tree(DG, branches, PLs, THRESH, PIXID, plco, toplot)
%
% INPUTS:
%   DG      - Directed graph object representing the merge tree
%   branches - Cell array containing the branches of the merge tree
%   PLs      - Vector of persistence lengths for every branch
%   THRESH   - Vector of threshold values for every node in the merge tree
%   PIXID    - Vector of pixel indices for every node in the merge tree
%   plco     - Persistence length cutoff value; branches with persistence length
%              below this value will be removed
%   toplot   - (Optional) Boolean flag to create a plot of the simplified merge tree;
%              if not provided, defaults to 0 (do not plot)
%
% OUTPUTS:
%   DG3     - Directed graph object representing the simplified merge tree
%   THRESH3 - Vector of threshold values for every node in the simplified merge tree
%   PIXID3  - Vector of pixel indices for every node in the simplified merge tree
%
% DESCRIPTION:
%   This function simplifies the input merge tree by removing branches with persistence
%   length below the specified cutoff value. The simplified merge tree is returned
%   along with the updated threshold values and pixel indices.
%
% Author: Zachary T. Berndsen

    % Set default value for toplot if not provided
    if nargin < 7
        toplot = 0;
    end

    if toplot == 1
        figure();
        LWidths = 5 * DG.Edges.Weight / max(DG.Edges.Weight);
        plot(DG, 'Layout', 'layered', 'Direction', 'up', 'LineWidth', LWidths, 'MarkerSize', 5);
    end

    % initialize arrays
    ShortBranches = find(PLs <= plco);
    N2R = [];

    for n = 1:length(ShortBranches)
        path = branches{ShortBranches(n), 1};
        N2R = [N2R, path(1:end-1)];
    end

    N2R = sort(N2R, 'descend');
    temp = 1:numnodes(DG);
    N2K = setdiff(temp, N2R);
    THRESH2 = THRESH(N2K);
    PIXID2 = PIXID(N2K);

    % Remove all branches below the plco
    DG2 = rmnode(DG, N2R);

    if toplot == 1
        figure(100);
        LWidths = 5 * DG2.Edges.Weight / max(DG2.Edges.Weight);
        plot(DG2, 'Layout', 'layered', 'Direction', 'up', 'LineWidth', LWidths, 'MarkerSize', 5);
    end

    [DG3, THRESH3, PIXID3] = compress_tree(DG2, THRESH2, PIXID2);
    plot_merge_tree(DG3, THRESH3);
    
end