function [header, data] = read_mrc(filename)
% READ_MRC Read an MRC file and output the header information and the image/volume.
%
% Inputs:
%   filename - Name of the MRC file
%
% Outputs:
%   header - Header information of the MRC file
%   data - Image/volume data
%
% Authors: Zachary T. Berndsen, ChatGPT4

    % Open the file in 'little-endian' mode for reading
    fid = fopen(filename, 'r', 'ieee-le');

    % Read the header information
    header.NX = fread(fid, 1, 'int32');
    header.NY = fread(fid, 1, 'int32');
    header.NZ = fread(fid, 1, 'int32');
    header.MODE = fread(fid, 1, 'int32');
    header.NXSTART = fread(fid, 1, 'int32');
    header.NYSTART = fread(fid, 1, 'int32');
    header.NZSTART = fread(fid, 1, 'int32');
    header.MX = fread(fid, 1, 'int32');
    header.MY = fread(fid, 1, 'int32');
    header.MZ = fread(fid, 1, 'int32');
    header.CELLA = fread(fid, 3, 'float32');
    header.CELLB = fread(fid, 3, 'float32');
    header.MAPC = fread(fid, 1, 'int32');
    header.MAPR = fread(fid, 1, 'int32');
    header.MAPS = fread(fid, 1, 'int32');
    header.DMIN = fread(fid, 1, 'float32');
    header.DMAX = fread(fid, 1, 'float32');
    header.DMEAN = fread(fid, 1, 'float32');
    header.ISPG = fread(fid, 1, 'int32');
    header.NSYMBT = fread(fid, 1, 'int32');
    header.EXTRA = fread(fid, 25, 'int32');
    header.ORIGIN = fread(fid, 3, 'float32');
    header.MAP = fread(fid, 4, 'char=>char')';
    header.MACHST = fread(fid, 1, 'int32');
    header.RMS = fread(fid, 1, 'float32');
    header.NLABL = fread(fid, 1, 'int32');
    header.LABELS = fread(fid, [80, 10], 'char=>char')';

    % Determine the data type based on the MODE
    switch header.MODE
        case 0
            data_type = 'int8';
        case 1
            data_type = 'int16';
        case 2
            data_type = 'float32';
        otherwise
            error('Unsupported data mode');
    end

    % Read the image/volume data
    data = fread(fid, [header.NX, header.NY * header.NZ], data_type);
    data = reshape(data, [header.NX, header.NY, header.NZ]);

    % Close the file
    fclose(fid);
end
