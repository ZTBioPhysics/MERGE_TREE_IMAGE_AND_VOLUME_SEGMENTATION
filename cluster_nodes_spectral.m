function [cluster_labels] = cluster_nodes_spectral(affinity_matrix)
% CLUSTER_NODES_SPECTRAL  Spectral clustering of graph nodes using an affinity matrix
% 
% Syntax:
%   [cluster_labels] = CLUSTER_NODES_SPECTRAL(affinity_matrix)
%
% Inputs:
%   affinity_matrix - An NxN affinity matrix representing the similarity between nodes, where N is the number of nodes
% 
% Outputs:
%   cluster_labels - A vector of length N containing the cluster labels of the nodes
%
% Description:
%   This function takes an NxN affinity_matrix as input, where N is the number of nodes, and returns a vector of length N containing the cluster labels of the nodes. The function uses spectral clustering to partition the nodes into clusters based on their similarity in the affinity matrix. The number of clusters is determined automatically using the eigengap heuristic.
%
% Example:
%   % Generate a synthetic dataset with 3 clusters
%   data = [randn(50,2); randn(50,2)+4; randn(50,2)+8];
%   % Compute the affinity matrix using Gaussian kernel
%   affinity_matrix = exp(-squareform(pdist(data)).^2 / (2 * 1^2));
%   % Perform spectral clustering
%   cluster_labels = cluster_nodes_spectral(affinity_matrix);
%   % Visualize the results
%   scatter(data(:, 1), data(:, 2), 30, cluster_labels, 'filled');
%   title('Spectral Clustering Results');
%   xlabel('Feature 1');
%   ylabel('Feature 2');
%
% Author: Zachary T. Berndsen, ChatGPT


    % Input validation
    if nargin < 1
        error('Not enough input arguments. Please provide an affinity matrix.');
    end
    if size(affinity_matrix, 1) ~= size(affinity_matrix, 2)
        error('Affinity matrix must be a square matrix.');
    end

    % Compute the Laplacian matrix
    D = diag(sum(affinity_matrix, 2));
    L = D - affinity_matrix;
    
    % Compute the normalized Laplacian
    D_inv_sqrt = diag(1 ./ sqrt(diag(D)));
    L_norm = D_inv_sqrt * L * D_inv_sqrt;
    
    % Compute eigenvectors and eigenvalues
    num_eig = min(10, size(affinity_matrix, 1)); % Number of eigenvectors and eigenvalues to compute
    [eig_vec, eig_val] = eigs(L_norm, num_eig, 'smallestabs');
    eig_val = diag(eig_val); % Extract eigenvalues
    
    % Find the eigengap
    eig_gaps = diff(eig_val);
    [~, num_clusters] = max(eig_gaps);
    
    % Normalize the rows of the eigenvectors
    eig_vec_norm = eig_vec(:, 1:num_clusters) ./ vecnorm(eig_vec(:, 1:num_clusters), 2, 2);
    
    % Apply k-means clustering
    cluster_labels = kmeans(eig_vec_norm, num_clusters);
    
    % Visualize the eigenvalues
    figure;
    plot(1:num_eig, eig_val, 'o-');
    hold on;
    plot(num_clusters, eig_val(num_clusters), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    title('Eigenvalues and Eigengap');
    xlabel('Index');
    ylabel('Eigenvalue');
    legend('Eigenvalues', 'Eigengap', 'Location', 'best');

end
