function calculate_and_plot_avg_and_total_connectivity(AffMat)
% CALCULATE_AND_PLOT_AVG_AND_TOTAL_CONNECTIVITY Calculate and plot the average and total connectivity of each column in the affinity matrix.
%
% Syntax:
%   calculate_and_plot_avg_and_total_connectivity(AffMat)
%
% Inputs:
%   AffMat: Affinity matrix
%
% Outputs:
%   None
%
% Description:
%   This function calculates the average and total connectivity of each column
%   in the input affinity matrix (excluding the diagonal) and plots the histograms
%   for both measures on the same plot with two y-axes on the right and left side.
%
% Author: Zachary T. Berndsen

    lc = length(AffMat);
    AVG = zeros(lc, 1);
    TCON = zeros(lc, 1);
    
    for n = 1:lc
        AVG(n) = (sum(AffMat(:, n)) - AffMat(n, n)) / (lc - 1);
        TCON(n) = sum(AffMat(:, n)) - AffMat(n, n);
    end
    
%     % Normalize the values
%     AVG_normalized = AVG / max(AVG);
%     TCON_normalized = TCON / max(TCON);

    % Plot average connectivity
    figure();
    bar(AVG);
    xlabel('Node Index');
    ylabel('Average Connectivity');
    title('Average Connectivity');

    % Plot total connectivity
    figure();
    bar(TCON);
    xlabel('Node Index');
    ylabel('Total Connectivity');
    title('Total Connectivity');

end

