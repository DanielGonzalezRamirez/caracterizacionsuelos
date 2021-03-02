function [S_loc] = S_point_location(Vref_prime_smooth)
%S_POINT_LOCATION Finds the location of the S point for the TDR signal
%   Using the first derivative of the reflected voltage this function
%   locates the inflection point (S) for the first reflection of the TDR
%   signal.
[~,max_loc] = max(Vref_prime_smooth);
zeros_length = ceil(max_loc + 0.4*max_loc); % Magic number = 0.4
aux_Vref_prime = [zeros(1,zeros_length) Vref_prime_smooth(zeros_length + 1:end)];
aux_Vref_prime = smooth(aux_Vref_prime);
[~,S_loc] = max(aux_Vref_prime);
end

