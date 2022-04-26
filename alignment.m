% This function is used to detect where are the spike based on single-channel raw data, neo_value and neo_threshold, modelled only in sw.
% Input:        data(ndataLen, 1)   - raw data
%               neo(1, dataLen)     - neo_value
%               th(1, dataLen)      - neo_threshold
% Output:       train(1, ndataLen)  - time-domain labels which are the end points of each spike
%               loc(1, dataLen)     - time-domain labels which are the most negative pont of each spike
% Reference:    https://ieeexplore.ieee.org/document/1608524

function [train, loc] = alignment(data, neo, th)

    i = 1;                                      % time stamp index
    max = 0;                                    % neo max value
    max_idx = 0;                                % neo max value index
    count = 0;                                  % 0 means no spike happens, !0 means spike happens
    L = size(data, 1);                          % ndataLen
    train = zeros(1, L);                        % labeling where a spike ends
    loc = zeros(1, L);                          % labeling where is the spike's the most negative point
    while i <= L                                % traverse all time stamp (1 ~ ndataLen)
        if count ~= 0                               % after detect neo > th
            if neo(i) > th(i)                           % if still neo > th
                if neo(i) > max                         % record max neo value nad max time index
                    max = neo(i);
                    max_idx = i;
                end
                count = count + 1;
            else                                        % if neo < th, which means spike ends 
                train(max_idx) = 1;                     % label the spike ends time stamp
                count = 0;                              % reset count
                max = 0;                                % reset max value
                [temp_pks, temp_loc] = findpeaks(-1*data(max_idx:max_idx+23));  
                % find the most negative peak from here to 23 time stamp after
                loc(max_idx + min(temp_loc) - 1) = 1;   % label the most negative peak
            end
        else
            if neo(i) > th(i)                       % if encounter a new spike, new neo > th
                if neo(i) > max                     % record max neo value nad max time index
                    max = neo(i);
                    max_idx = i;
                end
                count = count + 1;
            end
        end
        i = i + 1;                              % move on to the next time stamp
    end
end
        
    