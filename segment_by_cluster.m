function [Lab] = segment_by_cluster(I, PIXID, CLUSTERS, DG, toplot)
% SEGMENT_BY_CLUSTER Segment and image based on the clustered merge tree.
%
% Inputs:
%   I            -  Input image
%   PIXID        -  Cell array with pixel indices of connected components associated with each node in the merge tree
%   CLUSTERS     -  Cell array with the minimum number of nodes necessary to represent each cluster of leaf nodes in the merge tree or simply a vector with cluster IDs 
%   DG           -  Directed graph representing the merge tree
%   toplot       -  Boolean flag designating if you want to plot intermediate results
%
% Output:
%   Lab          - Segmented image
%
% Author: Zachary T. Berndsen

    % Initialize variables
    mask = zeros(size(I));
    Lab = double(mask);

    if iscell(CLUSTERS)
        clrmp = fliplr(1:size(CLUSTERS, 1));
        num_clusters = size(CLUSTERS, 1);
    else
        clrmp = fliplr(1:max(CLUSTERS));
        num_clusters = max(CLUSTERS);
    end

    % Set default kernel type to 'gaussian' if not provided
    if nargin < 5
        toplot = 0;
    end

    for i = 1:num_clusters
        
        if iscell(CLUSTERS)
            % Extract nodes for each cluster
            nodes = CLUSTERS{i, 2};
        else
            nodes = find(CLUSTERS==i);
        end

        % Combine all connected components
        pix = unique(vertcat(PIXID{nodes}));

        % Update segmented image
        Lab_pix = find(Lab);
        pix2 = setdiff(pix, Lab_pix);
        Lab(pix) = clrmp(i);

        if toplot == 1
            figure(99)
            LWidths = 5*DG.Edges.Weight/max(DG.Edges.Weight);
            p=plot(DG,'Layout','layered','Direction','up','LineWidth',LWidths,'MarkerSize',5);
            highlight(p,nodes,'NodeColor',[0 0.9 0])
            
            figure(100)
            imshow(label2rgb(Lab));
            pause
        end

    end
end

    