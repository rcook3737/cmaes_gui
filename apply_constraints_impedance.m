    %% Apply contraints to step frequency during optimization

% This function is used during human-in-the-loop optimization of step
% length. The high and low bounds are measured before optimization and used
% to limit the requested step frequencies. This function thresholds the
% CMA-ES generated frequencies and outputs new step frequencies to be
% tested and relayed back into the algorithm.

% INPUTS
% original parameters: set of step freqs from CMA-ES
% Freq. bounds: the highest and lowest step freqs for that subject

% OUTPUTS
% new parameters: thresholded parameters based on bounding

function params = apply_constraints_impedance(original_params, param_bounds ,N)
params = original_params;
% param_bounds = n-parameters x 2 matrix (low and high limits for each param)
if N ~= size(original_params,1)
    N = input('N does not match parameter bounds matrix. Enter appropriate N');
end

for b = 1:size(original_params,1)
% Stiffness
    % All parameters are asserted to be bigger or equal to low bound.
    params(b,params(b,:)<param_bounds(b,1))=param_bounds(b,1);

    % All parameters are asserted to be lesser or equal to high bound.
    params(b,params(b,:)>param_bounds(b,2))=param_bounds(b,2);
    
% % Theta0
%     % All parameters are asserted to be bigger or equal to low bound.
%     params(2,params(2,:)<param_bounds(2,1))=param_bounds(2,1);
% 
%     % All parameters are asserted to be lesser or equal to high bound.
%     params(2,params(2,:)>param_bounds(2,2))=param_bounds(2,2);
end    
return