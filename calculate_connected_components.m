function [num_cc, connected_components] = calculate_connected_components(input_image, thresholds)
% CALCULATE_CONNECTED_COMPONENTS  Calculates the number of connected components as a function of threshold for a binary or grayscale image
%
% SYNTAX:
%   [num_cc, connected_components] = calculate_connected_components(input_image, thresholds)
%
% INPUTS:
%   input_image - A grayscale image or volume
%   thresholds - A vector of intensity or threshold values
%
% OUTPUTS:
%   num_cc - A matrix that stores the number of connected components and the threshold value for each threshold
%   connected_components - A cell array of structure arrays that store the connected components at each threshold level
%
% DESCRIPTION:
%   This function takes a binary or grayscale image or volume and a vector of intensity or threshold values and calculates the number of connected components as a function of threshold. The output of the function is a matrix called "num_cc" that stores the number of connected components and the threshold value for each threshold, and a cell array of structure arrays called "connected_components" that stores the connected components at each threshold level. The function can be used for various applications such as image segmentation, feature extraction, and object tracking.
%
% EXAMPLE:
%   img = imread('cameraman.tif');
%   thresholds = [50, 100, 150, 200];
%   [num_cc, connected_components] = calculate_connected_components(img, thresholds);
%   figure;
%   plot(num_cc(:, 2), num_cc(:, 1));
%   xlabel('Threshold');
%   ylabel('Number of connected components');
%   title('Connected components as a function of threshold');
%
% Author: Zachary T. Berndsen

    num_levels = length(thresholds);
    connected_components = cell(num_levels, 1);
    num_cc = zeros(num_levels, 3);

    for i = 1:num_levels
        level = thresholds(i);
        bw_img = input_image > level;
        connected_components{i} = bwconncomp(bw_img);
        
        num_cc(i, 1) = connected_components{i}.NumObjects;
        num_cc(i, 2) = level;
        
        if i == 1
            num_cc(i, 3) = num_cc(i, 1);
        else
            num_cc(i, 3) = num_cc(i - 1, 3) + num_cc(i, 1);
        end
    end
end
