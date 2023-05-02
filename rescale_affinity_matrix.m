function Rescaled_AffMat = rescale_affinity_matrix(AffMat, sigma, kernel_type)
% RESCALE_AFFINITY_MATRIX Rescale the affinity matrix with a specified kernel.
%
% Syntax:
%   RescaledAffMat = rescale_affinity_matrix(AffMat, sigma, kernel_type)
%
% Inputs:
%   AffMat: Affinity matrix
%   sigma: Kernel width parameter
%   kernel_type: Kernel type ('gaussian' or 'alternative')
%
% Outputs:
%   RescaledAffMat: Rescaled affinity matrix
%
% Description:
%   This function rescales the affinity matrix using either the Gaussian kernel or an alternative kernel. The choice of kernel affects the rescaled values in the affinity matrix.
%
% Author: Zachary T. Berndsen

    % Get the size of the affinity matrix
    [nRows, nCols] = size(AffMat);

    % Initialize the rescaled affinity matrix
    Rescaled_AffMat = zeros(nRows, nCols);

    % Set default kernel type to 'gaussian' if not provided
    if nargin < 3
        kernel_type = 'gaussian';
    end


    % Rescale the affinity matrix using the chosen kernel
    for i = 1:nRows
        for j = 1:nCols
            if i ~= j
                if strcmp(kernel_type, 'gaussian')
                    Rescaled_AffMat(i, j) = exp(-(AffMat(i, j)^2) / (2 * sigma^2));
                elseif strcmp(kernel_type, 'alternative')
                    Rescaled_AffMat(i, j) = exp((-(AffMat(i, j)^2)) / (sigma^2));
                else
                    error('Invalid kernel type. Choose "gaussian" or "alternative".');
                end
            end
        end
    end

end
