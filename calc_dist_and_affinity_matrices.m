function [AffMat, DistMat] = calc_dist_and_affinity_matrices(DG3, THRESH, PIXID, I, AP, type, pixelsize)
% CALC_DIST_AND_AFFINITY_MATRICES Calculate affinity and distance matrices for a merge tree.
%
% Syntax:
%   [AffMat, DistMat] = calc_dist_and_affinity_matrices(DG3, THRESH, PIXID, I, AP_idx, type, pixelsize)
%
% Inputs:
%   DG3: Directed graph representing the merge tree
%   THRESH: Vector of threshold values for each node
%   PIXID: Cell array containing pixel indices for each node
%   I: Image matrix
%   AP: Affinity parameter (1 for inverse isovalue, 2 for volume/area, 3 for hypervolume, 4 for intrinsic tree distance)
%   type: Node type (1 for birth nodes, 2 for birth and merge nodes)
%   pixelsize (optional): The size of a pixel, defaults to 1
%
% Outputs:
%   AffMat: Affinity matrix
%   DistMat: Distance matrix
%
% See also: REGIONPROPS, SHORTESTPATH
%
% Authors: Zachary T. Berndsen, ChatGPT


    % Set default value for toplot if not provided
    if nargin < 7
        pixelsize = 1;
    end
    
    % calculate the adjacency matrix of the merge tree and convert to an
    % undirected graph structure
    A = adjacency(DG3,'weighted');
    UDG = graph(A,'upper');

    % Remove the global minimum from the merge tree
    UDG = rmnode(UDG, numnodes(UDG));

    if type == 1
        nodes = find(degree(UDG) == 1)';
    elseif type == 2
        nodes = 1:numnodes(UDG);
    end

    % Initialize matrices
    AffMat = zeros(length(nodes));
    DistMat = zeros(length(nodes));

    % Compute affinity and distance matrices
    for i = 1:length(nodes)
        n1 = nodes(i);
        cent1 = get_centroid(PIXID, n1, I);

        for j = 1:length(nodes)
            if i == j
                AffMat(i, j) = 0;
                DistMat(i, j) = 0;
            else
                n2 = nodes(j);
                cent2 = get_centroid(PIXID, n2, I);
                DistMat(i, j) = pdist([cent1; cent2], 'euclidean')*pixelsize;

                [path, d] = shortestpath(UDG, n1, n2);
                lca = max(path); % least common ancestor node

                switch AP
                    case 1
                        AffMat(i, j) = 1 / THRESH(lca);
                    case 2
                        ConMat(i, j) = 1 / length(PIXID{lca});
                    case 3
                        AffMat(i, j) = 1 / sum(I(PIXID{lca}));
                    case 4
                        AffMat(i, j) = 1 / d;
                end
            end
        end
    end
end

function centroid = get_centroid(PIXID, node, img)
    img(:) = 0;
    pix = PIXID{node};
    img(pix) = 1;
    rp = regionprops(img);
    centroid = rp.Centroid;
end
