function [num_cc, connected_components, THRESH, aug_merge_tree, directed_graph] = create_aug_merge_tree(input_image, thresholds)
% CREATE_AUG_MERGE_TREE  Generates a merge tree of a grayscale image
% 
% Syntax:
%   [num_cc, thresholds, merge_tree, directed_graph] = create_merge_tree(input_image, thresholds)
%
% Inputs:
%   input_image - A grayscale image to generate the merge tree for
%   thresholds - A vector of intensity thresholds used to segment the image
% 
% Outputs:
%   num_cc - a matrix containing the total number of connected components found across all threshold levels
%   connected_components - A cell array of structure arrays that store the connected components at each threshold level
%   THRESH - a vector of threshold values associated with every connected component
%   aug_merge_tree - A matrix representing the augmented merge tree of the image, where rows and columns correspond to indexed connected components and entries are 1 if the component is contained within the other and 0 otherwise
%   directed_graph - A directed graph structure representing the augmented merge tree of the image
%
% Description:
%   This function takes a grayscale image and a vector of intensity thresholds and generates an augmented merge tree that represents the containment relationships between connected components at different threshold levels. The output of the function is a matrix called "merge_tree" and a directed graph structure called "directed_graph" that can be used for various applications such as image segmentation, feature extraction, and object tracking.
%
% Example:
%   img = imread('cameraman.tif');
%   input_image = mat2gray(img);
%   thresholds = fliplr(linspace(min(input_image(:)),max(input_image(:)),100));
%   [num_cc, connected_components, THRESH, aug_merge_tree, directed_graph] = create_aug_merge_tree(input_image, thresholds);
%   figure;
%   plot_aug_merge_tree(directed_graph);
%   title('Augmented Merge Tree of Image')
%
% Author: Zachary T. Berndsen

    % Get connected components at each intensity level and store them in an
    % indexed structure array

    num_levels = length(thresholds);
    [num_cc,connected_components] = calculate_connected_components(input_image,thresholds);
    
    % Initialize merge_tree matrix as sparse
    numcc = num_cc(:,3);
    aug_merge_tree = zeros(numcc(end), numcc(end));
    THRESH = zeros(numcc(end),1);
    
    % For each pair of adjacent intensity levels loop through each connected component and
    % find which connected component at the next level it is contained within. Record a 1 in
    % the merge_tree matrix at (i,j)
    count = 0;
    for i = 1:(num_levels-1)
        % Get connected components at this level and the next
        connected_components1 = connected_components{i};
        connected_components2 = connected_components{i+1};
        num_objects1 = connected_components1.NumObjects;
        num_objects2 = connected_components2.NumObjects;
    
        % For each connected component at this level
        for j = 1:num_objects1
            count = count+1;
            THRESH(count) = num_cc(i,2);
            % Get the pixels for this connected component
            pixels1 = connected_components1.PixelIdxList{j};
    
            % For each connected component at the next level
            for k = 1:num_objects2
                % Get the pixels for this connected component
                pixels2 = connected_components2.PixelIdxList{k};
                
                % Check if the connected component at the this level is contained
                % within the connected component at the next level
                if all(ismember(pixels1, pixels2))
                    % add at 1 to the matrix if yes
                    aug_merge_tree(count,count+(num_objects1-j)+k) = 1;
                end
            end
        end
    end

    % Convert to a directed graph
    directed_graph = digraph(aug_merge_tree);

end
