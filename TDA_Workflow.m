%% TDA Workflow %%
close all
clearvars

%% Select image or volume
[file,path] = uigetfile({'*.mrc;*.tif;*.tiff;*.png;*.jpg;*.jpeg;*.gif', 'Image Files';...
                          '*.*', 'All Files'}, 'Select Image or Volume');
if isequal(file,0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(path,file)]);
    [I, dim2] = load_image_or_volume(fullfile(path,file));
end

%% or load test image

dim2 = 1;
I = imread('pears.png');
I = rgb2gray(I);
I = mat2gray(I);

% gaussian filter
sigma = 2;
I = imgaussfilt(I,2);

%% Display image or volume slices

% Display image or volume
figure();
if dim2 == 1
    imshow(I);
else
    montage(I, 'DisplayRange', [min(I(:)) max(I(:))]);
end

%% Generate and plot the augmented merge tree

% Set starting and ending threshold of the filtration as the max and min of image
NT = min(min(min(I)));
ST = max(max(max(I)));

% determine the number of steps in the filtration
steps = 100;

% % Automatically determine the final threshold based on connected components (0th Betti Number) vs Threshold plot
% NT = find_noise_threshold(I,steps); 

% Create threshold vector
thresh = fliplr(linspace(NT, ST, steps));

if dim2 == 1
    plot_contours_and_extrema2(I,NT)
end

% Generate and plot the augmented merge tree
[NUMCC, CC, THRESH, aug_merge_tree, DG] = create_aug_merge_tree(I, thresh);
plot_aug_merge_tree(DG)

% Plot # of connected components vs threshold
numCC = NUMCC(:,1);
figure();
plot(thresh, numCC);
set(gca, 'xdir', 'reverse')
xlabel('Threshold')
ylabel('0th Betti Number')
set(gca, 'YScale', 'log')

%% Generate and plot the merge tree
% Compresses branches of the augmented merge tree by removing all nodes 
% with ID and OD = 1. These nodes represent topological features that
% are not changing. Each branch is replaced by a weighted edge equal to the
% function span between the birth and merge node.

[DG2, THRESH2, PIXID] = compress_aug_merge_tree(DG, THRESH, NUMCC, CC, thresh);
plot_merge_tree(DG2, THRESH2)

%% Perform branch decomposition then calculate and plot persistence diagram
% Using the elder rule, assign critical point pairs of birth nodes with
% merge nodes and record their persistence length

% assign persistence pairs by performing branch decomposition
toplot = 0; % boolean flag indicating if you want to interatively plot
[branches, ~] = branch_decomposition(DG2, THRESH2, toplot);

% Calculate persistence lengths and plot the distribution
toplot = 1; % boolean flag if you want to plot the distribution
[PLs] = calculate_persistence_length_distribution(branches, toplot);

% Calculate and plot the persistence curve
toplot = 1; % boolean flag if you want to plot the persistence curve
[~] = calculate_and_plot_persistence_curve(PLs, toplot);

% Perform "stability analysis"
[~, ~] = stability_analysis(branches, THRESH2, PLs);

% Calculate and plot the persistence diagram
toplot = 1; % boolean flag if you want to plot the persistence diagram
[~] = calculate_and_plot_persistence_diagram(branches, THRESH2, toplot);

%% Persistence-based simplification

% Ask user if they want to perform persistence-based simplification
simplify_button = questdlg('Do you want to perform persistence-based simplification?', ...
    'Persistence-based Simplification', 'Yes', 'No', 'No');

if strcmpi(simplify_button, 'yes')
    % Prompt user for persistence length cutoff
    prompt = {'Enter persistence length cutoff:'};
    dlgtitle = 'Persistence Length Cutoff';
    dims = [1 35];
    definput = {'0.1'}; 
    persistence_cutoff = str2double(inputdlg(prompt, dlgtitle, dims, definput));

    % Perform persistence-based simplification using the specified cutoff
    cutoffline = linspace(1, 0, 100);
    plot(cutoffline, cutoffline-persistence_cutoff, 'k--')
    
    % Simplify merge tree by pruning all edges < plco
    toplot = 0; % boolean flag if you want to plot
    [DG3, THRESH3, PIXID2] = simplify_merge_tree(DG2, branches, PLs, THRESH2, PIXID, persistence_cutoff, toplot);
    
    % recalculate the persistence pairs
    toplot = 0; % boolean flag if you want to plot
    [DG3_branches, ~] = branch_decomposition(DG3, THRESH3, toplot);

    % Calculate and plot the persistence diagram
    toplot = 1; % boolean flag if you want to plot
    [~] = calculate_and_plot_persistence_diagram(DG3_branches, THRESH3, toplot);
    plot(cutoffline, cutoffline-persistence_cutoff, 'k--')
else
    DG3 = DG2;
    THRESH3 = THRESH2;
    PIXID2 = PIXID;
end

%% Hierarchical segmentation of 2D image based on simplified merge tree

% segment the image based on the nodes of the simplified merge tree
toplot = 1; % boolean flag if you want to iteratively plot
[Lab] = segment_image_with_merge_tree(I, DG3, PIXID2, toplot);

if toplot == 0

    figure();
    if dim2 == 1
        imshow(label2rgb(Lab));
    else
        montage(label2rgb(Lab), 'DisplayRange', [0 max(Lab(:))]);
    end
end

%% segment 2D image by persistence pairs

% segment the image by persistence pairs
toplot = 1; % boolean flag if you want to iteratively plot
[Lab] = segment_image_by_persistence_pairs(I, DG3, DG3_branches, PIXID2, toplot);
    
if toplot == 0
    if dim2 == 1
        imshow(label2rgb(Lab));
    else
        montage(label2rgb(Lab), 'DisplayRange', [0 max(Lab(:))]);
    end
end

%% calculate the merge tree affinity and distance matrices
% The affinity matrix is a square matrix with rows and columns corresponding
% to nodes in the merge tree, and entries representing the affinity parameter.
% The affinity parameter is designed to capture a measure of similarity or
% connectivity between two nodes in the merge tree, and it can be one of four
% different measures. The type of matrix specifies whether to calculate the
% affinity between just leaf nodes, or both leaf and merge nodes. The distance
% matrix, on the other hand, contains the Euclidean distances between the centers
% of mass of the connected components corresponding to each node in the merge
% tree. This provides an additional metric for understanding the spatial
% relationships between the nodes in the merge tree.

% Prompt user to select the type of affinity parameter
AP_options = {'Isovalue', 'Volume/Area', 'Total Intensity', 'Intrinsic Tree Distance'};
[selected_AP, AP_idx] = listdlg('PromptString', 'Select the type of affinity parameter:', ...
                                'SelectionMode', 'single', ...
                                'ListString', AP_options, ...
                                'ListSize', [250, 100]);
                            
if isempty(AP_idx)
    error('No affinity parameter selected. Aborting.');
end

% Prompt user to select node type
type_options = {'Birth nodes only', 'All nodes'};
[selected_type, type_idx] = listdlg('PromptString', 'Select the node type:', ...
                                    'SelectionMode', 'single', ...
                                    'ListString', type_options, ...
                                    'ListSize', [200, 100]);
                                
if isempty(type_idx)
    error('No node type selected. Aborting.');
end

% Prompt user for pixel size
prompt = {'Enter pixel size:'};
dlgtitle = 'Pixel Size';
dims = [1 35];
definput = {'1'};
pixel_size_input = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(pixel_size_input)
    error('No pixel size value provided. Aborting.');
else
    pixelsize = str2double(pixel_size_input{1});
end

% calculate the affinity and distance matricies
[AffMat, DistMat] = calc_dist_and_affinity_matrices(DG3, THRESH3, PIXID2, I, AP_idx, type_idx, pixelsize);

figure()
heatmap(AffMat)
xlabel('Node Index');
ylabel('Node Index');
title('Affinity Matrix');
figure()
histogram(AffMat,20)
xlabel('Affinity');
ylabel('Count');
title('Affinity Histogram');

figure()
heatmap(DistMat)
xlabel('Node Index');
ylabel('Node Index');
title('Euclidean Distance Matrix');
figure()
histogram(DistMat,20)
xlabel('Euclidean Distance');
ylabel('Count');
title('Distance Histogram');

% % repermute matrix to block diagonal
% toplot =1;
% [DistMat] = rcm_permuted_matrix(DistMat, toplot);

% % calculate average and total connectivity of each column
% calculate_and_plot_avg_and_total_connectivity(AffMat)

% % hard cutoff on Distmat?
% DistMat(DistMat>200)=0;

%% rescale the affinity and distance matrices

% Prompt the user to select if they want to rescale the matrices
rescale = questdlg('Do you want to rescale the affinity matrix:', ... 
    'Rescale?', 'no', 'yes', 'no');

if strcmpi(rescale, 'yes')

    % Prompt the user to select a kernel type
    kernel_type = questdlg('Select the kernel type:', ...
        'Kernel Type', 'gaussian', 'alternative', 'gaussian');
    
    % Provide the kernel width parameter sigma
    prompt = {'Enter the kernel width parameter sigma:'};
    dlgtitle = 'Kernel Width Parameter';
    dims = [1 35];
    definput = {'1'};
    sigma = str2double(inputdlg(prompt, dlgtitle, dims, definput));
    
    % Calculate the rescaled affinity matrix using the chosen kernel
    Rescaled_AffMat = rescale_affinity_matrix(AffMat, sigma, kernel_type);

else
    Rescaled_AffMat = AffMat;
end

% Prompt the user to select if they want to rescale 
rescale = questdlg('Do you want to rescale the distance matrix:', ... 
    'Rescale?', 'no', 'yes', 'no');

if strcmpi(rescale, 'yes')
    % Prompt the user to select a kernel type
    kernel_type = questdlg('Select the kernel type:', ...
        'Kernel Type', 'gaussian', 'alternative', 'gaussian');
    
    % Provide the kernel width parameter sigma
    prompt = {'Enter the kernel width parameter sigma:'};
    dlgtitle = 'Kernel Width Parameter';
    dims = [1 35];
    definput = {'1'};
    sigma = str2double(inputdlg(prompt, dlgtitle, dims, definput));
    
    % Calculate the rescaled distance matrix using the chosen kernel
    Rescaled_DistMat = rescale_affinity_matrix(DistMat, sigma, kernel_type);

else
    Rescaled_DistMat = DistMat;
end

%% Create the aggregate matrix
% multiply the DistMat by the AffMat

AggMat = Rescaled_DistMat.*Rescaled_AffMat;

% % use just the AffMat or DistMat
% AggMat = Rescaled_AffMat;
% AggMat = Rescaled_DistMat;

%% Cluster nodes in the merge tree using the aggregate matrix
% Cluster nodes of the merge tree using graph community detection algorithms
% such as modularity maximization or spectral clustering.

% Prompt the user to select a clustering algorithm
    algo_type = questdlg('Select the clustering algorithm:', ...
        'Algorithm', 'ModMax', 'Spectral', 'Spectral');
if strcmp(algo_type, 'Spectral')
    % spectral clustering.
    clusters = cluster_nodes_spectral(AggMat);
    max(clusters) % number of clusters
else
    clusters = GCModulMax1(AggMat);
    max(clusters) % number of clusters
end

% Plot merge tree with nodes colored by cluster id.
f=figure();
LWidths = 5*DG3.Edges.Weight/max(DG3.Edges.Weight);
p=plot(DG3,'Layout','layered','Direction','up','LineWidth',LWidths,'MarkerSize',5);

% Color nodes of merge tree by cluster ID.
ID = indegree(DG3);
if type_idx == 1
    Births = find(ID==0); 
    CID = zeros(1,numnodes(DG3)); % cluster ID for each leaf node
    for n = 1:length(Births)
        CID(Births(n)) = clusters(n);
    end
elseif type_idx == 2
    UCCS = find(or(ID==0,ID>1));    
    CID = zeros(1,numnodes(DG3)); % cluster ID for each leaf node
    for n = 1:length(UCCS)
        CID(UCCS(n)) = clusters(n);
    end
end

CLUSTERS1 = CID;
p.NodeCData = CLUSTERS1;
cbh = colorbar;

% Plot the AggMat as a graph/network and color nodes by cluster_id
network = graph(AggMat);
figure();
LWidths = 5*network.Edges.Weight/max(network.Edges.Weight);
p=plot(network,'Layout','force','LineWidth',LWidths,'MarkerSize',5);
p.NodeCData = clusters;
colorbar;

%% segment 2D image by cluster

% expand each cluster upwards through the tree such that merge nodes whos
% subtree has all the same cluster_id are assinged to the same cluster

if type_idx == 1
    toplot = 0;
    [CLUSTERS, MIN_CLUSTERS] = expand_merge_tree_clusters(DG3, clusters, toplot);
end

% input can take MIN_CLUSTERS if expand_merge_tree was used, or CLUSTERS1
% from previous section
toplot = 1; % 1 if you want to plot results iteratively, 0 if not
[Lab] = segment_by_cluster(I, PIXID2, MIN_CLUSTERS, DG3, toplot);

if toplot == 0
    figure();
    imshow(imresize(label2rgb(Lab),1));

    figure();
    LWidths = 5*DG3.Edges.Weight/max(DG3.Edges.Weight);
    p=plot(DG3,'Layout','layered','Direction','up','LineWidth',LWidths,'MarkerSize',5);
    p.NodeCData = CLUSTERS;
    cbh = colorbar;
end




