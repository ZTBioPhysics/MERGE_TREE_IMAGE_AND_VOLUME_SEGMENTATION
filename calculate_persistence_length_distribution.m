function [PLs] = calculate_persistence_length_distribution(branches,toplot)
% CALCULATE_PERSISTENCE_LENGTH_DISTRIBUTION Computes and plots the distribution of persistence pair lengths for a given set of branches
%
% Syntax:
%   PLs = PersistenceDist(branches, plot)
%
% Inputs:
%   branches - A cell array containing the persistence pairs, with each row representing a branch and the second column containing the persistence length
%   toplot - A boolean flag indicating whether to plot the distribution of persistence lengths (optional; default is false)
%
% Outputs:
%   PLs - A vector containing the persistence lengths for each branch
%
% Description:
%   This function computes the persistence lengths for a given set of 
%   branches and optionally plots the distribution of persistence lengths 
%   as a histogram with a logarithmic scale on the y-axis.
%
% Author: Zachary T. Berndsen

    % Check if toplot is provided, otherwise set it to 0 (do not plot)
    if nargin < 2
        toplot = 0;
    end

    % plot the distribution of persistence pair persistence lengths
    PLs = zeros(length(branches),1);
    for i = 1:length(PLs)
        PLs(i) = branches{i,2};
    end
    
    if toplot == 1
        figure();
        [n, xout] = hist(PLs,20);
        bar(xout, n, 'barwidth', 1, 'basevalue', 0.1);
        set(gca,'YScale','log')
        ylabel('counts')
        xlabel('Persistence')
    end

end