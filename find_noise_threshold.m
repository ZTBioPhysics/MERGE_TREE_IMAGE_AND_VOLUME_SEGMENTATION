function noiseThresh = find_noise_threshold(I, numSteps)
% FIND_NOISE_THRESHOLD finds the optimal threshold for a given image or
% volume that best separates the signal from the noise.
% The function computes the 0th Betti number (number of connected components) for
% a range of threshold values and then identifies two key peaks in the
% threshold vs. connected components curve. The first peak corresponds to
% the high-intensity signal, while the second peak corresponds to the
% low-intensity noise. The function then finds the global minimum between these
% two peaks and returns the threshold value corresponding to this minimum.
% In cases where the minimum is degenerate, the function returns the
% threshold value corresponding to the last index in the degenerate minimum
% valley.
%
% INPUTS:
%   I        - input image or volume (2D or 3D matrix)
%   numSteps - (optional) number of threshold steps to be used in the analysis
%              Default value is 100.
%
% OUTPUT:
%   noiseThresh - the threshold value that separates signal from noise
%
% USAGE:
%   noiseThresh = find_noise_threshold(I)
%   noiseThresh = find_noise_threshold(I, numSteps)
% 
% Author: Zachary T. Berndsen


    % Set the default number of threshold steps
    default_value = 100;

    % Check if the user has provided a value for the numSteps parameter
    if nargin < 2
        numSteps = default_value;
    end

    % Create a vector of threshold values from high to low
    thresh = fliplr(linspace(min(I(:)), max(I(:)), numSteps));

    % Compute connected components for each threshold
    [NUMCC, ~] = CC_Thresh(I, thresh);
    numCC = NUMCC(:, 1);

    % Find the local maxima (peaks) in the numCC data
    [~, locs] = findpeaks(numCC);

    % Return early if there are less than 2 peaks
    if length(locs) < 2
        minIndex = [];
        return;
    end

    % Find the first peak and the last peak
    firstPeakIdx = locs(1);
    lastPeakIdx = locs(end);

    % Find the global minimum between the first peak and the last peak
    [~, minIndex] = min(numCC(firstPeakIdx:lastPeakIdx));
    minIndex = minIndex + firstPeakIdx - 1;

    % Split the numCC data into two parts
    numCC_part1 = numCC(1:minIndex);
    numCC_part2 = numCC(minIndex+1:end);

    % Find the highest peak in each part
    [~, peakIdx1] = max(numCC_part1);
    [~, peakIdx2] = max(numCC_part2);

    % Add the offset for the second part
    peakIdx2 = peakIdx2 + minIndex;

    % Find the global minimum between the two highest peaks
    [~, minIndex] = min(numCC(peakIdx1:peakIdx2));
    minIndex = minIndex + peakIdx1 - 1;

    % Find the last index in the degenerate minimum valley
    minValue = numCC(minIndex);
    for i = minIndex+1:peakIdx2
        if numCC(i) == minValue
            minIndex = i;
        else
            break;
        end
    end

    % Plot the threshold vs connected components with the minIndex highlighted
    figure();
    plot(thresh, numCC);
    hold on;
    plot(thresh(minIndex), numCC(minIndex), 'ro', 'MarkerSize', 10);
    set(gca, 'xdir', 'reverse');
    xlabel('Threshold');
    ylabel('0th Betti number');
    set(gca, 'YScale', 'log');

    % Return the threshold at minIndex as the noise threshold
    noiseThresh = thresh(minIndex);
end
