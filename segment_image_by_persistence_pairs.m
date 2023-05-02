function [Lab] = segment_image_by_persistence_pairs(I, DG3, DG3_branches, PIXID2, toplot)
% PERSISTENCEPAIR_SEGMENTATION_2D Performs 2D image segmentation based on persistence pairs
%
% INPUTS:
%   - I: 2D input image
%   - DG3_branches: Cell array containing the branches of the merge tree
%   - PIXID2: Cell array containing the indices of pixels in each connected component
%   - toplot: Boolean flag to create a plot of the segmented image (optional)
%
% OUTPUTS:
%   - Lab: Label matrix of the segmented image
%
% DESCRIPTION: This function assigns labels to the pixels in the input image based
% on the persistence pairs derived from the merge tree. If toplot
% is set to 1, the segmented image will be plotted during the process.
%
% Author: Zachary T. Berndsen

    % Set default value for toplot if not provided
    if nargin < 5
        toplot = 0;
    end

    % Get the number of rows
    num_rows = size(DG3_branches, 1);
    
    % Iterate through the rows and print the last value of each vector
    last_value = zeros(num_rows,1);
    for i = 1:num_rows
        last_value(i) = DG3_branches{i, 1}(end);
    end

    % Sort branches by their terminal node
    [~, idx] = sort(last_value, 'descend');

    % Initialize label matrix
    mask = zeros(size(I));
    Lab = double(mask);
    clrmp = fliplr(1:length(idx));

    % Assign labels to pixels based on merge tree
    for n = 1:length(idx)
        path = [DG3_branches{idx(n),1}];
        if n==1
            path=path(1:end-1);
        end
        node = path(end);
        pix = PIXID2{node};

        % Update segmented image
        Lab(pix) = clrmp(n);

        % Plot segmented image if toplot is set to 1
        if toplot == 1

            figure(99)
            LWidths = 5*DG3.Edges.Weight/max(DG3.Edges.Weight);
            p=plot(DG3,'Layout','layered','Direction','up','LineWidth',LWidths,'MarkerSize',5);
            highlight(p,path,'NodeColor',[0 0.9 0])

            figure(100)
            imshow(label2rgb(Lab));
            pause
        end
    end
end

