% This function is used to calculate neo value and its adaptive threshold, modelled in sw and hw.
% Input:        win         - window length to calculate neo_mean and neo_std
%               sigma       - neo_value = x[n]^2 - x[n-sigma]*x[n+sigma]
%               x           - single-channel extracellular data
%               thS_factor  - threshold parameter (neo_thS = neo_mean * thS_factor)
%               thH_factor  - threshold parameter (neo_thH = neo_mean + thH_factor*neo_std)
% Output:       neo_value   - nonlinear energy operator 
%               neo_thS     - neo threshold typeS
%               neo_thH     - neo threshold typeH
% Reference:    neo_value:  https://ieeexplore.ieee.org/document/1608524
%               neo_thS:    https://ieeexplore.ieee.org/document/5740375
%               neo_thH:    https://ieeexplore.ieee.org/document/6070974

function [neo_value, neo_thS, neo_thH] = neo(win, sigma, x, thS_factor, thH_factor)

    % neo_value
    N = size(x, 2);
    valid_value = x(1, sigma+1:N-sigma) .* x(1, sigma+1:N-sigma) - x(1, 1:N-2*sigma) .* x(1, 2*sigma+1:N);
    
    neo_value = zeros(1, N);                        % exclude the last few point, and assign to zeros
    neo_value(1, sigma+1:N-sigma) = valid_value;
    
    acc = 0;                    % moving accumulation for neo_value
    acc2 = 0;                   % moving accumulation for neo_value^2
        
    neo_thS = zeros(1, N);
    neo_thH = zeros(1, N);
    
    for i = 1:1:N
        
        if i <= win                                         % initialization
            acc2 = acc2 + neo_value(1, i)*neo_value(1, i);  % acc2(n) = sum(neo(1:n)*neo(1:n))
            acc = acc + neo_value(1, i);                    % acc(n)  = sum(neo(1:n))
        else
            acc2 = acc2 + neo_value(1, i)*neo_value(1, i) - neo_value(1, i-win)*neo_value(1, i-win);    % acc2(n) = acc2(n-1) - neo(n-win)^2 + neo(n)^2
            acc = acc + neo_value(1, i) - neo_value(1, i-win);                                          % acc(n)  = acc(n-1) - neo(n-win) + neo(n)
        end
        
        % hardware
        mov_mean = acc/win;                                                         % moving neo mean
        neo_thS(1, i) = mov_mean * thS_factor;                                      % adaptive neo threshold typeS
        neo_thH(1, i) = mov_mean + sqrt(acc2/win - mov_mean*mov_mean) * thH_factor; % adaptive neo threshold typeH
        
        % software
%         neo_thH(1, i) = mean(neo_value(1, 1:i)) + thH_factor*std(neo_value(1, 1:i));
    end
end