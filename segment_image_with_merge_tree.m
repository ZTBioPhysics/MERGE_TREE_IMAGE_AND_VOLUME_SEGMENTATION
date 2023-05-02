function [Lab] = segment_image_with_merge_tree(I, DG3, PIXID2, toplot)
% SEGMENT_IMAGE_WITH_MERGE_TREE Performs 2D image segmentation based on a merge tree
%
% INPUTS:
%   - I: 2D input image
%   - DG3: Directed graph object representing the simplified merge tree
%   - PIXID2: Cell array containing the indices of pixels in each connected component
%   - toplot: Boolean flag to create a plot of the segmented image (optional)
%
% OUTPUTS:
%   - Lab: Label matrix of the segmented image
%
% DESCRIPTION: This function assigns labels to the pixels in the input image based
% on the simplified merge tree (DG3) and the pixel indices (PIXID2). If toplot
% is set to 1, the segmented image will be plotted during the process.
%
% Author: Zachary T. Berndsen

    % Set default value for toplot if not provided
    if nargin < 4
        toplot = 0;
    end

    % Remove the global minimum from the merge tree
    DGf = rmnode(DG3, numnodes(DG3));
    numNodes = numnodes(DGf);

    % Initialize label matrix
    Lab = zeros(size(I));

    % Assign labels to pixels based on merge tree
    for nodeIdx = fliplr(1:numNodes)
        pix = PIXID2{nodeIdx};
        Lab(pix) = nodeIdx;

        % Plot segmented image if toplot is set to 1
        if toplot == 1
            figure(100)
            imshow(imresize(label2rgb(Lab), 10));
            pause
        end
    end
end


