function [DG2, THRESH2, PIXID] = compress_aug_merge_tree(DG, THRESH, NUMCC, CC, thresh)
% COMPRESS_AUG_MERGE_TREE - Compress an augmented merge tree by removing 
% nodes with ID=1 and OD=1. This will create the true merge tree
%
% INPUTS:
%   - DG: Directed graph object representing the merge tree
%   - THRESH: Vector of thresholds used to construct the merge tree
%   - NUMCC: Matrix of connected components for each threshold level
%   - CC: Cell array of connected components for each threshold level
%   - thresh: Vector of threshold values
%
% OUTPUTS:
%   - DG2: Directed graph object representing the compressed merge tree
%   - THRESH2: Vector of thresholds for the compressed merge tree
%   - PIXID: Cell array of pixel indices for each unique connected component
%
% This function compresses an augmented merge tree by removing nodes with
% in-degree and out-degree equal to 1. The resulting merge tree is 
% represented as a directed graph object. The function also calculates
% the pixel indices for each unique connected component and the threshold
% values for the merge tree.
%
% Author: Zachary T. Berndsen
    
    % calculate the size of the discretization step
    dt = thresh(1)-thresh(2);

    % Find all nodes of each type.
    ID = indegree(DG);
    OD = outdegree(DG);
    Births = find(ID == 0);
    Mergers = find(ID > 1);
    
    % Calculate the number of unique connected components and their IDs.
    NumUCCs = length(Births) + length(Mergers);
    UCCIDs = sort([Births; Mergers]);
    
    % Initialize THRESH2 and PIXID.
    THRESH2 = zeros(NumUCCs + 1, 1);
    PIXID = cell(NumUCCs, 3);
    
    % Preallocate reduced adjacency matrix.
    reduced_tree = zeros(length(UCCIDs) + 1);
    
    % Iterate through each unique connected component.
    for n = 1:length(UCCIDs)
        nid = UCCIDs(n);
        [path, d] = find_branch_dg(DG, nid);

        % Assign a representative connected component to each branch.
        largest_contour = path(end - 1);
        CCs = CC{thresh == THRESH(largest_contour)};
        numcc = CCs.NumObjects;
        numcc2 = NUMCC(find(NUMCC(:, 3) >= largest_contour, 1, 'first'), 3);
        ccid = numcc - (numcc2 - largest_contour);
        PIXID{n} = CCs.PixelIdxList{ccid};
        
        % Update adjacency matrix with edge weight as the function span.
        PL = d * dt;
        if OD(path(end)) == 0
            reduced_tree(n, end) = PL;
            THRESH2(n) = THRESH(nid);
        else
            reduced_tree(n, (UCCIDs == path(end))') = PL;
            THRESH2(n) = THRESH(nid);
        end
    end
    
    % Create the compressed directed graph.
    DG2 = digraph(reduced_tree);
end
