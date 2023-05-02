function Plateaus = find_plateaus(x_data, y_data, sort_col)
% FIND_PLATEAUS Find plateaus in the given data and return a table with XValue, YValue, and Length of all plateaus
%
% SYNTAX:
%   Plateaus = find_plateaus(x_data, y_data, sort_col)
%
% INPUTS:
%   x_data   - Vector containing the x values of the data
%   y_data   - Vector containing the y values of the data
%   sort_col - Column to sort the data on (3 corresponds to plateau length)
%
% OUTPUTS:
%   Plateaus - Table with columns XValue, YValue, and PlateauLength of all
%              plateaus in the input y_data, sorted based on the specified column
%
% DESCRIPTION:
%   This function identifies plateaus in the input y_data and returns a table with
%   columns XValue, YValue, and PlateauLength for all plateaus found. The data is
%   sorted based on the column specified by sort_col. The function uses vectorized
%   operations and built-in MATLAB functions to simplify the code and improve efficiency.
%
% Author: Zachary T. Berndsen, ChatGPT

    % Calculate differences between adjacent y_data values
    y_diff = diff(y_data);
    
    % Find start and end indices of plateaus
    start_indices = find([1; y_diff(:) ~= 0]);
    end_indices = find([y_diff(:) ~= 0; 1]) - 1;
    
    % Calculate lengths, y_values, and x_values for plateaus
    plateau_lengths = end_indices - start_indices + 1;
    plateau_y_values = y_data(start_indices);
    plateau_x_values = x_data(start_indices);
    
    % Combine the data into a single matrix
    plateau_data = [plateau_x_values(:), plateau_y_values(:), plateau_lengths(:)];
    
    % Sort the data based on the specified column
    sorted_data = sortrows(plateau_data, sort_col, 'descend');
    
    % Create the output table
    Plateaus = table(sorted_data(:, 1), sorted_data(:, 2), sorted_data(:, 3), 'VariableNames', {'XValue', 'YValue', 'PlateauLength'});
end
