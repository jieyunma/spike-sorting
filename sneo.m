% This function is used to transform neo value to smoothed neo value and calculate its adaptive threshold, modelled in sw and hw.
% Input:        win         - window length to calculate neo_mean and neo_std
%               neo         - singal-channel neo value
%               thS_factor  - threshold parameter (sneo_thS = sneo_mean * thS_factor)
%               thH_factor  - threshold parameter (sneo_thH = sneo_mean + thH_factor*neo_std)
% Output:       sneo_value  - smoothed nonlinear energy operator 
%               sneo_thS    - sneo threshold typeS
%               sneo_thH    - sneo threshold typeH
% Reference:    sneo_value: https://ieeexplore.ieee.org/document/1608524
%               sneo_thS:   https://ieeexplore.ieee.org/document/5740375
%               sneo_thH:   https://ieeexplore.ieee.org/document/6070974

function [sneo_value, sneo_thS, sneo_thH] = sneo(win, neo, thS_factor, thH_factor)

    % sneo_value
    N = size(neo, 2);
    w = [0.08, 0.54, 1, 0.54, 0.08];
    sneo_value = filter(w, 1, neo);
    
    acc = 0;                    % moving accumulation for neo_value
    acc2 = 0;                   % moving accumulation for neo_value^2
    
    sneo_thS = zeros(1, N);
    sneo_thH = zeros(1, N);
    
    for i = 1:1:N
        
        if i <= win                                             % initialization
            acc2 = acc2 + sneo_value(1, i)*sneo_value(1, i);    % acc2(n) = sum(sneo(1:n)*sneo(1:n))
            acc = acc + sneo_value(1, i);                       % acc(n)  = sum(sneo(1:n))
        else
            acc2 = acc2 + sneo_value(1, i)*sneo_value(1, i) - sneo_value(1, i-win)*sneo_value(1, i-win);    % acc2(n) = acc2(n-1) - sneo(n-win)^2 + sneo(n)^2
            acc = acc + sneo_value(1, i) - sneo_value(1, i-win);                                            % acc(n)  = acc(n-1) - sneo(n-win) + sneo(n)
        end
        
        % hardware
        mov_mean = acc/win;                                                             % moving sneo mean
        sneo_thS(1, i) = mov_mean * thS_factor;                                         % adaptive sneo threshold typeS
        sneo_thH(1, i) = mov_mean + sqrt(acc2/win - mov_mean*mov_mean) * thH_factor;    % adaptive sneo threshold typeH
        
        % software
%         sneo_thH(1, i) = mean(neo_value(1, 1:i)) + thH_factor*std(neo_value(1, 1:i));
    end    
end

