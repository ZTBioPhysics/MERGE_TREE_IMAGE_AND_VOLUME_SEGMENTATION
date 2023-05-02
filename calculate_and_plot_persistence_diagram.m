function [PD] = calculate_and_plot_persistence_diagram(branches,THRESH2, toplot)
% PERSISTENCEDIAGRAM Calculate and optionally plot the persistence diagram
%
% SYNTAX:
%   PD = PersistenceDiagram(branches, THRESH2, toplot)
%
% INPUTS:
%   branches - Cell array containing the branches of the merge tree
%   THRESH2  - Vector of threshold values for every node in the merge tree
%   toplot   - (Optional) Boolean flag to create a plot of the persistence diagram; 
%              if not provided, defaults to 0 (do not plot)
%
% OUTPUTS:
%   PD - Matrix containing the birth and death thresholds for each
%   persistence pair
%
% DESCRIPTION:
%   This function calculates the persistence diagram for the input branches of the merge tree, 
%   and, if requested, plots the persistence diagram.
%
% Author: Zachary T. Berndsen

    % Set default value for toplot if not provided
    if nargin < 3
        toplot = 0;
    end

    numBranches = length(branches);
    B = zeros(1, numBranches);
    D = zeros(1, numBranches);

    for n = 1:numBranches
        p = branches{n, 1};
        B(n) = THRESH2(p(1));
        D(n) = THRESH2(p(end));
    end

    if toplot == 1
        figure();
        set(gca, 'XDir', 'reverse');
        set(gca, 'YDir', 'reverse');
        ylabel('Death');
        xlabel('Birth');
        title('Persistence Diagram');
        hold on;
        scatter(B, D, 'bo');
        thresh = linspace(1, 0, 100);
        plot(thresh, thresh, '--');
        xlim([thresh(end), thresh(1)]);
    end

    PD = [B', D'];
end