function [NumPairs, TotalPL] = stability_analysis(branches, THRESH2, PLs)
%% STABILITY_ANALYSIS Perform stability analysis and plot the number of pairs and total persistence against threshold levels
%
% SYNTAX:
%   [NumPairs, TotalPL] = stability_analysis(branches, THRESH2, PLs)
%
% INPUTS:
%   branches - Cell array containing branch information
%   THRESH2  - Vector with threshold values corresponding to each node in the merge tree
%   PLs      - Vector containing the persistence pair lengths for each branch
%
% OUTPUTS:
%   NumPairs - Vector containing the number of pairs at each threshold level
%   TotalPL  - Vector containing the total persistence at each threshold level
%
% DESCRIPTION:
%   This function performs stability analysis on the input branches data and
%   calculates the number of pairs and total persistence at different threshold
%   levels. It also plots the number of pairs and total persistence against the
%   threshold levels.
%
% Authors: Zachary T. Berndsen, ChatGPT

    % Calculate B and D values for each branch
    B = arrayfun(@(n) THRESH2(branches{n, 1}(1)), 1:length(branches));
    D = arrayfun(@(n) THRESH2(branches{n, 1}(end)), 1:length(branches));
    
    % Initialize output vectors
    NumPairs = zeros(1, 1000);
    TotalPL = zeros(1, 1000);
    
    % Calculate number of pairs and total persistence for each threshold level
    thresh = linspace(1, 0, 1000);
    for n = 1:length(thresh)
        id3 = find(B >= thresh(n) & D < thresh(n));
        NumPairs(n) = length(id3);
        TotalPL(n) = sum(PLs(id3));
    end
    
    % Plot the number of pairs and total persistence against threshold levels
    figure();
    yyaxis right
    plot(thresh, NumPairs)
    set(gca, 'XDir', 'reverse')
    set(gca, 'YScale', 'log')
    ylabel('number of pairs');
    ylim([0, max(NumPairs) + 2])
    xlabel('threshold');
    yyaxis left
    plot(thresh, TotalPL)
    ylabel('total persistence');
    ylim([0, max(TotalPL) + 2])
    
    % Plateau information
    find_plateaus(thresh, TotalPL, 3)

end
