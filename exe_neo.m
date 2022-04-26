% This function is used to calculate single-channel neo value and write out result values.
% Input:        result_dir  - output directory
%               step        - including raw, median, filter, re-filter, normal
%               idx         - channel index
%               sigma       - neo parameter (neo_value = x[n]^2 - x[n-sigma]*x[n+sigma])
%               data        - single-channel extracellular data
%               thS_factor  - threshold parameter (neo_thS = neo_mean * thS_factor)
%               thH_factor  - threshold parameter (neo_thH = neo_mean + thH_factor*neo_std)
% Output:       NULL

function exe_neo(result_dir, step, idx, sigma, data, thS_factor, thH_factor)

    % initialization
    SNEO_thH = zeros(1, size(data, 2));
    SNEO_thS = zeros(1, size(data, 2));
    dataSNEO = zeros(1, size(data, 2));
    NEO_thH = zeros(1, size(data, 2));
    NEO_thS = zeros(1, size(data, 2));
    dataNEO = zeros(1, size(data, 2));
    
    % performa neo and sneo calculation
    disp('    -NEO start');
    [dataNEO, NEO_thS, NEO_thH] = neo(512, sigma, data, thS_factor, thH_factor); disp('    -NEO done');
    [dataSNEO, SNEO_thS, SNEO_thH] = sneo(512, dataNEO, thS_factor, thH_factor); disp('    -SNEO done');
    
    % output single-channel data including raw/pre-processed extracellular data, 
    % neo_value, neo typeS threshold, neo typeH threshold, sneo_value, sneo typeS threshold and sneo typeH threshold
    disp('    -write out data');
    fileID = fopen(strcat(result_dir, step, '_', num2str(idx), '.txt'), 'w');
    fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'RAW', 'NEO', 'THS', 'THH', 'SNEO', 'STHS', 'STHH');
    fprintf(fileID, '%6.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n', [data; dataNEO; NEO_thS; NEO_thH; dataSNEO; SNEO_thS; SNEO_thH]);
    fclose(fileID);
end    