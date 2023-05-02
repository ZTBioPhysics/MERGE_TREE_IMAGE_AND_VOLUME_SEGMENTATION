function [PC] = calculate_and_plot_persistence_curve(PLs,toplot)
% CALCULATE_AND_PLOT_PERSISTENCE_CURVE Calculate and optionally plot the persistence curve
%
% SYNTAX:
%   PC = calculate_and_plot_persistence_curve(PLs, toplot)
%
% INPUTS:
%   PLs    - Vector containing the persistence pair lengths for each branch
%   toplot - (Optional) Boolean flag to create a plot of the persistence curve; 
%            if not provided, defaults to 0 (do not plot)
%
% OUTPUTS:
%   PC - Vector containing the calculated persistence curve
%
% DESCRIPTION:
%   This function calculates the persistence curve for the input persistence pair lengths
%   and, if requested, plots the persistence curve.
%
% Author: Zachary T. Berndsen

    % Set default value for toplot if not provided
    if nargin < 2
        toplot = 0;
    end

    numIndices = 1000;
    ind = linspace(min(PLs), max(PLs), numIndices);
    PC = zeros(1, numIndices);

    for n = 1:numIndices
        PC(n) = sum(PLs >= ind(n));
    end
    
    if toplot == 1
        figure();
        plot(ind, PC);
        set(gca, 'YScale', 'log');
        set(gca, 'XScale', 'log');
        ylabel('counts');
        xlabel('Persistence');
    end

    find_plateaus(ind, PC, 3)

end