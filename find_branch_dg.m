function [path,d] = find_branch_dg(DG,node)
% FIND_BRANCH_DG - Find the path and length of a branch in a directed graph
%
% INPUTS:
%   - DG: Directed graph object
%   - node: Starting node of the branch
%
% OUTPUTS:
%   - path: Vector of node IDs representing the branch
%   - d: Length of the branch
%
% DESCRIPTION: This function finds the path and length of a branch in a directed graph
% starting from a given node. The branch is defined as a sequence of nodes
% that have in-degree 1 and out-degree 1, and the path includes all nodes
% in the branch. The function returns the path and length of the branch as
% outputs.
%
% Author: Zachary T. Berndsen
    
    ID = indegree(DG);
    parent = successors(DG,node);
    PD = ID(parent);
    edgs = DG.Edges;
    edgn = edgs{:,1};
    edgw = edgs{:,2};
    d = 0;
    d = d + edgw(edgn(:,1)==node);

    % while loop runs until it reaches a node with input degree >= 1
    % record path of every branch (all nodes in the branch)
    path = [node parent];
    while PD == 1
        newid = parent;
        parent = successors(DG,newid);
        path = [path parent];
        if ~isempty(parent)
            d = d + edgw(edgn(:,1)==newid);
        end
        try
            PD = ID(parent);
        catch
            PD = 0;
        end
    end
end
