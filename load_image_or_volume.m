function [I, dim2] = load_image_or_volume(filepath)
% LOAD_IMAGE_OR_VOLUME Loads an image or volume from a file, converts it to 
% grayscale, and optionally applies a Gaussian filter
%
% Syntax:
%   [I, dim2] = load_image_or_volume(filepath)
%
% Inputs:
%   filepath - A string specifying the path to the image or volume file to load
%
% Outputs:
%   I - The loaded image or volume converted to grayscale
%   dim2 - A flag indicating if the loaded data is a 2D image (1) or a 3D volume (0)
%
% Description:
%   This function loads an image or volume from the specified file path, 
%   converts it to grayscale, and optionally applies a Gaussian filter to 
%   the data. The function supports loading MRC files as well as common 
%   image formats such as JPEG, PNG, and TIFF. If an MRC file is loaded, 
%   the user is prompted to select the type of electron microscopy volume 
%   (Cryo-EM or NS-EM). After loading the data, the user is asked if they 
%   want to apply a Gaussian filter to the image or volume.
%
% Authors: Zachary T. Berndsen, ChatGPT

    % Load image or volume and convert to grayscale
    [~,~,ext] = fileparts(filepath);
    
    if strcmpi(ext, '.mrc')
        % Ask user to select the type of volume being loaded
        button = questdlg('Select the type of electron microscopy volume:', ...
            'Volume Type', 'NS', 'Cryo', 'Cryo');
        
        % Load volume based on the selected type
        [~, I] = read_mrc(filepath);
        dim2 = 0;
        if ndims(I) == 3
            if strcmpi(button, 'cryo')
                % Cryo-EM volume
                I = I - min(I(:));
                I = I / max(I(:));
            elseif strcmpi(button, 'ns')
                % NS-EM volume
                I = mat2gray(I);
            else
                error('Invalid volume type specified. Must be "Cryo" or "NS".');
            end
        else
            % 2D image
            I = imadjust(I);
            if ndims(I) == 3
                % RGB image
                I = rgb2gray(I);
                dim2 = 1;
            else
                dim2 = 0;
            end
        end
    else
        % Load non-MRC image
        I = imread(filepath);
        dim2 = 1;
        if ndims(I) == 3
            % RGB image
            I = rgb2gray(I);          
        end
        I = mat2gray(I);
    end

    % Ask user if they want to perform a Gaussian filter
    filter_button = questdlg('Do you want to perform a Gaussian filter?', ...
        'Gaussian Filter', 'Yes', 'No', 'No');
    if strcmpi(filter_button, 'yes')
        % Prompt user for sigma value
        prompt = {'Enter sigma value:'};
        dlgtitle = 'Gaussian Filter Sigma';
        dims = [1 35];
        definput = {'2'};
        sigma = str2double(inputdlg(prompt, dlgtitle, dims, definput));

        % Apply Gaussian filter
        I = imgaussfilt(I, sigma);
    end

end