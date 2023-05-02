function [CLUSTERS, MIN_CLUSTERS] = expand_merge_tree_clusters(DG, leaf_clusters, toplot)
% EXPAND_MERGE_TREE_CLUSTERS attempts to assigne cluster ids to internal
% nodes of a merge tree based on the cluster ids of the leaf nodes
%
% Syntax:
%   CLUSTERS = expand_merge_tree_clusters(DG, leaf_clusters)
%
% Inputs:
%   DG: A directed graph
%   leaf_clusters: A vector of cluster IDs for each leaf node in the directed graph
%   toplot: boolean flag for ploting intermediate results
%
% Outputs:
%   CLUSTERS: A vector of cluster IDs for each node in the directed graph
%
% Description:
%   This function assigns a cluster ID to each node in the directed graph
%   based on finding subtrees that contain leaves that all have the same 
%   cluster id. Starting at each leaf node, the algorithm searches upward
%   through the tree until it reaches a merge node whos subtree contains
%   leaf nodes with more than one cluster id. If a subtree is found that
%   contains leaves that all have the same cluster id, those leaves are
%   removed from the list of nodes that need to be processed and the merge
%   node is assigned as the representative of that subtree cluster
%
% Author: Zachary T. Berndsen

    % Flip the direction of the digraph
    DGf = flipedge(DG);
    
    % Find all the leaf nodes
    leaves = find(outdegree(DGf) == 0)';
    
    % Initialize cluster IDs to 0 for all nodes
    CLUSTERS = zeros(numnodes(DG),1);
    CLUSTERS(leaves) = leaf_clusters;
    
    % initialize MIN_CLUSTERS cell array
    unique_clusters = unique(leaf_clusters);
    MIN_CLUSTERS = cell(length(unique_clusters),2);
    for i = 1:length(unique_clusters)
        MIN_CLUSTERS{i,1} = unique_clusters(i);
        MIN_CLUSTERS{i,2} = leaves(leaf_clusters==unique_clusters(i));
    end

    all_vectors = [];
    for i = 1:size(MIN_CLUSTERS, 1)
        all_vectors = [all_vectors MIN_CLUSTERS{i, 2}];
    end
    unseen_leaf_nodes = sort(intersect(all_vectors,leaves));

    while ~isempty(unseen_leaf_nodes)

        % Find the predecessor of the current leaf node
        node = unseen_leaf_nodes(1);
        pred = predecessors(DGf,node);
        
        % Find all the other leaf nodes in the subtree with the predecessor as root
        subtree = dfsearch(DGf,pred);  % find all nodes in the subtree
        subtree = subtree(2:end); % remove parent
        sub_leaves = intersect(leaves,subtree); % remove merge nodes
        unique_clusters = unique(CLUSTERS(sub_leaves));
        
        % If all the leaf nodes have the same cluster ID, assign that cluster ID to the predecessor node
        if length(unique_clusters)==1
            % assign the pred node cluster id its leaf nodes
            CLUSTERS(pred) = unique_clusters;

            % update the MIN_CLUSTERS array by removing subtree nodes and
            % adding the predecessor node
            current_nodes = MIN_CLUSTERS{unique_clusters,2};
            keep_nodes = setxor(sub_leaves',current_nodes);
            new_set = sort([pred keep_nodes]);
            MIN_CLUSTERS{unique_clusters,2} = new_set;

            % remove leaf nodes in subtree from list of leaf nodes still
            % needing to be processed
            unseen_leaf_nodes = setdiff(unseen_leaf_nodes,sub_leaves);

            if toplot==1
                figure(90)
                LWidths = 5*DG.Edges.Weight/max(DG.Edges.Weight);
                p=plot(DG,'Layout','layered','Direction','up','LineWidth',LWidths,'MarkerSize',5);
                highlight(p,subtree,'NodeColor',[0 0.9 0])
                pause
                highlight(p,subtree,'NodeColor',[0 0.0 0])
            end 

            % Continue moving up the tree to the next predecessor until we encounter a node whose subtree does not contain leaf nodes with the same cluster ID
            while length(unique_clusters)==1
                pred = predecessors(DGf,pred);
                subtree = dfsearch(DGf,pred);  % find all nodes in the subtree
                subtree = subtree(2:end); % remove parent
                sub_leaves = intersect(leaves,subtree); % remove merge nodes
                unique_clusters = unique(CLUSTERS(sub_leaves));
                if length(unique_clusters)==1
                    CLUSTERS(pred) = unique_clusters;
                    current_nodes = MIN_CLUSTERS{unique_clusters,2};
                    keep_nodes = setdiff(current_nodes,subtree);
                    new_set = sort([pred keep_nodes]);
                    MIN_CLUSTERS{unique_clusters,2} = new_set;
                    unseen_leaf_nodes = setdiff(unseen_leaf_nodes,sub_leaves);
                    
                    if toplot==1
                        figure(90)
                        LWidths = 5*DG.Edges.Weight/max(DG.Edges.Weight);
                        p=plot(DG,'Layout','layered','Direction','up','LineWidth',LWidths,'MarkerSize',5);
                        highlight(p,subtree,'NodeColor',[0 0.9 0])
                        pause
                        highlight(p,subtree,'NodeColor',[0 0.0 0])
                    end
                end
                
            end
            
        else
            unseen_leaf_nodes = setdiff(unseen_leaf_nodes,node);
        end
    end
end
