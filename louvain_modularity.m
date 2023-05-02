function [clusters] = louvain_modularity(adj_matrix)
% LOUVAIN_MODULARITY Perform modularity maximization using the Louvain method
%
% SYNTAX:
%   clusters = louvain_modularity(adj_matrix)
%
% INPUTS:
%   adj_matrix - Square adjacency matrix representing the graph
%
% OUTPUTS:
%   clusters - Vector with cluster IDs for each node
%
% DESCRIPTION:
%   This function performs modularity maximization using the Louvain
%   method, taking an adjacency matrix as input and returning a vector
%   with cluster IDs for each node.
%
% Authors: Zachary T. Berndsen, ChatGPT

    % Initialize clusters
    num_nodes = size(adj_matrix, 1);
    clusters = 1:num_nodes;

    % Calculate degrees and total edge weight
    degrees = sum(adj_matrix, 2);
    total_edge_weight = sum(degrees) / 2;

    % Initialize modularity and improvement flag
    Q = -inf;
    improved = true;

    while improved
        % Perform local optimization
        [new_clusters, new_Q] = local_optimization(adj_matrix, clusters, degrees, total_edge_weight);

        % Check for improvement
        if new_Q - Q <= 1e-10
            improved = false;
        else
            % Update clusters and modularity
            clusters = new_clusters;
            Q = new_Q;

            % Perform community aggregation
            adj_matrix = aggregate_communities(adj_matrix, clusters);
            num_nodes = size(adj_matrix, 1);
            degrees = sum(adj_matrix, 2);
            clusters = 1:num_nodes;
        end
    end
end

function [new_clusters, Q] = local_optimization(adj_matrix, clusters, degrees, total_edge_weight)
    num_nodes = size(adj_matrix, 1);
    new_clusters = clusters;
    improved = true;

    while improved
        improved = false;

        for i = 1:num_nodes
            % Remove node i from its current community
            curr_community = new_clusters(i);
            dQ_remove = delta_Q(adj_matrix, degrees, total_edge_weight, i, curr_community, -1);
            new_clusters(i) = -1;

            % Calculate modularity change for moving node i to each neighboring community
            neighbors = find(adj_matrix(i, :) > 0);
            neighbor_communities = unique(new_clusters(neighbors));
            neighbor_communities(neighbor_communities == -1) = [];

            dQ = arrayfun(@(c) delta_Q(adj_matrix, degrees, total_edge_weight, i, c, 1), neighbor_communities);
            [max_dQ, max_dQ_idx] = max(dQ);

            % Move node i to the community with the largest modularity gain
            if max_dQ > dQ_remove
                new_clusters(i) = neighbor_communities(max_dQ_idx);
                improved = true;
            else
                new_clusters(i) = curr_community;
            end
        end
    end

    Q = modularity(adj_matrix, new_clusters, total_edge_weight);
end

function new_adj_matrix = aggregate_communities(adj_matrix, clusters)
    unique_clusters = unique(clusters);
    num_communities = length(unique_clusters);
    new_adj_matrix = zeros(num_communities);

    for i = 1:num_communities
        for j = i:num_communities
            community_i_nodes = (clusters == unique_clusters(i));
            community_j_nodes = (clusters == unique_clusters(j));
            total_weight = sum(sum(adj_matrix(community_i_nodes, community_j_nodes)));
            new_adj_matrix(i, j) = total_weight;
            if i ~= j
                new_adj_matrix(j, i) = total_weight;
            end
        end
    end
end

function dQ = delta_Q(adj_matrix, degrees, total_edge_weight, node, community, add_or_remove)
    community_nodes = find(community == community);
    A_iC = sum(adj_matrix(node, community_nodes));
    k_i = degrees(node);
    K_C = sum(degrees(community_nodes));
    dQ = (A_iC - k_i * K_C / (2 * total_edge_weight)) * add_or_remove / total_edge_weight;
end

function Q = modularity(adj_matrix, clusters, total_edge_weight)
    num_nodes = size(adj_matrix, 1);
    Q = 0;

    for i = 1:num_nodes
        for j = 1:num_nodes
            if clusters(i) == clusters(j)
                Q = Q + adj_matrix(i, j) - degrees(i) * degrees(j) / (2 * total_edge_weight);
            end
        end
    end

    Q = Q / (2 * total_edge_weight);
end
