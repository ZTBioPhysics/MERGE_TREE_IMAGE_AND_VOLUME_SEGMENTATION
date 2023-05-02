function [permuted_matrix] = rcm_permuted_matrix(matrix, toplot)
% RCM_PERMUTED_MATRIX  Computes the Reverse Cuthill-McKee permuted matrix and visualizes the original and permuted matrices using heatmap plots
%
% Syntax:
%   [permuted_matrix] = RCM_PERMUTED_MATRIX(matrix)
%
% Inputs:
%   matrix - A square matrix
%   toplot - boolean flag for plotting results
%
% Outputs:
%   permuted_matrix - The Reverse Cuthill-McKee permuted matrix
%
% Description:
%   This function takes a square matrix as input, computes the Reverse Cuthill-McKee (RCM) permutation, applies the permutation to the matrix, and returns the permuted matrix. It also creates heatmap plots of both the original matrix and the RCM permuted matrix.
%
% Example:
%   A = gallery('tridiag', 100); % Create a tridiagonal matrix
%   A_perm = rcm_permuted_matrix(A);

% Author: ChatGPT

    % Input validation
    if nargin < 1
        error('Not enough input arguments. Please provide a square matrix.');
    end
    if size(matrix, 1) ~= size(matrix, 2)
        error('The input matrix must be a square matrix.');
    end

    % Compute the Reverse Cuthill-McKee (RCM) permutation
    p = symrcm(matrix);

    % Apply the permutation to the matrix
    permuted_matrix = matrix(p, p);

    if toplot == 1
        % Create heatmap plots
        figure;
        subplot(1, 2, 1);
        heatmap(matrix);
        title('Original Matrix');
    
        subplot(1, 2, 2);
        heatmap(permuted_matrix);
        title('RCM Permuted Matrix');
    end
end
