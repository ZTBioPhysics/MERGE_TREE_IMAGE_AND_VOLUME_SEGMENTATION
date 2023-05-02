function plot_contours_and_extrema(I, NT)
% PLOT_CONTOURS_AND_EXTREMA Plots an input image as a surface, displays a 
% contour plot with regional extrema, and plots histograms and line plots 
% of the intensity values at these extrema
%
% Syntax:
%   plot_contours_and_extrema(I, NT)
%
% Inputs:
%   I - A grayscale image to plot and analyze for regional extrema
%   NT - A threshold value to plot as a horizontal dashed line in the line plot of regional extrema intensities
%
% Outputs:
%   None. This function generates a series of plots.
%
% Description:
%   This function takes a grayscale image and plots it as a surface, 
%   displays a contour plot with regional maximum and minimum intensities 
%   indicated by red and blue stars, respectively, and plots histograms and 
%   line plots of the intensity values at these extrema.
%
% Authors: Zachary T. Berndsen, ChatGPT

    % Plot the image as a surface and show the contour plot with regional
    % maximum and minimum indicated by red and blue stars, respectively

    % Plot the surface
    plot_surface(I);

    % Plot the contour with regional extrema
    plot_contour_with_extrema(I);

    % Plot histograms of the regional extrema intensities
    plot_histograms(I);

    % Plot line plots of the regional extrema intensities
    plot_line_plots(I, NT);
end

function plot_surface(I)
    figure();
    surf(I);
end

function plot_contour_with_extrema(I)
    figure();
    [~, h] = imcontour(I, 20);
    h.LineWidth = 1;
    hold on;

    plot_extrema(I, @imregionalmin, 'b*');
    plot_extrema(I, @imregionalmax, 'r*');
end

function plot_extrema(I, extrema_fun, plot_style)
    extrema_mask = extrema_fun(I);
    [x, y] = find(extrema_mask);
    plot(y, x, plot_style);
end

function plot_histograms(I)
    [regminInts, regmaxInts] = get_regional_intensities(I);

    figure();
    hold on;
    histogram(regmaxInts, 70, 'FaceColor', 'blue');
    histogram(regminInts, 30, 'FaceColor', 'red');
    xlabel('Intensity at regional extrema');
    ylim([-inf 130]);
end

function plot_line_plots(I, NT)
    [regminInts, regmaxInts] = get_regional_intensities(I);

    figure();
    hold on;
    plot(regminInts);
    plot(regmaxInts);
    xlabel('regional extrema');
    ylabel('intensity');
    xlim([0 length(regmaxInts)]);
    plot(ones(size(regmaxInts)) * NT, 'k--');
end

function [regminInts, regmaxInts] = get_regional_intensities(I)
    regminInts = I(imregionalmin(I));
    regmaxInts = I(imregionalmax(I));
end
