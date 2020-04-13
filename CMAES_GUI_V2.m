%% Sample code for the CMA-ES process used in the main study.

% This function emulates the CMA-ES optimization used in the system. In
% actual experiments, the following process was implemented using simulink
% state machine.

% USER DEFINED PARAMETERS
% N: Number of objective variables
% xmean: Initial guess of objective variable distribution mean
%        xmean(1): peak torque
%        xmean(2): 2 * percent stride at which peak torque occurs
%        xmean(3): rise time, as a percentage of stride period
%        xmean(4): 2 * fall time, as a percentage of stride period
% sigma: Standard deviation of variable distribution. Initial value needs
% selecting
% stopeval: Number of evaluations before stopping the optimization
% lambda: Population size of each generation
% pc, ps   % Evolution paths for C and sigma. Initial values needs
% selecting
% B: Matrix that defines the coordinate system. Initial values needs
% selecting
% D: Vector that defines the scaling. Initial values needs selecting
% C: Covariance matrix C of variable distribution. Initial values needs
% selecting


% OUTPUTS
% xmin: The estimated optimal objective variable value that minimizes the
% objective function.

%%% Adapated from Kirby Witte's code adapted from ...
% This function was adapted from  Nikolaus Hansen's code located at
% URL: http://www.lri.fr/~hansen/purecmaes.m and references at the end. It
% demonstrates a basic CMA-ES process as used in our main study.



function [xmin,CMAES_struct]= CMAES_GUI_V2(filename,sigma,generations,N,paramNames,paramMeans,paramMins,paramMaxes)   % (mu/mu_w, lambda)-CMA-ES
% --------------------  Initialization --------------------------------
load(filename)

%% Initialize Variables
if counteval == 0
    paramRange_size = paramMaxes-paramMins;
    paramMean_norm = (paramMeans)./paramRange_size; % normalized means (ensures equal weighting of parameters) 
    paramRange = [paramMins, paramMaxes];
    paramBounds = (paramRange)./paramRange_size;
    %paramBounds = (paramRange)./paramRange_size; % aka the normalized boundaries

    % Strategy parameter setting: Selection
    lambda = 4+floor(3*log(N));  % Population size, offspring number
    mu = lambda/2;               % Number of parents/points for recombination
    weights = log(mu+1/2)-log(1:mu)'; % mu*1 array for weighted recombination
    mu = floor(mu);
    weights = weights/sum(weights);     % Normalize recombination weights array
    mueff = sum(weights)^2/sum(weights.^2); % Variance-effectiveness of sum w_i x_i
    stopeval = lambda*generations;   % Stop after stopeval number of function evaluations
    x = zeros(N,lambda); % cols are parameter sets (child), matrix for whole gen

    % Data Structure for data allocation
    CMAES_struct = struct('gen_idx',[],'x',zeros(lambda,N),'metrates',zeros(lambda,1),...
        'min_metrate',[],'x_mean',[],'x_std',[],'New_x_mean',[],'Bounds',paramRange);

    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(N,1); ps = zeros(N,1);   % Evolution paths for C and sigma
    B = eye(N);                         % B defines the coordinate system
    D = ones(N,1);                      % Diagonal D defines the scaling
    C = eye(N);                         % Covariance matrix C
    invsqrtC = eye(N);                  % C^-1/2
    eigeneval = 0;                      % Track update of B and D
    chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % Expectation of
                                        %   ||N(0,I)|| == norm(randn(N,1))

    % Strategy parameter setting: Adaptation
    cc = (4+mueff/N) / (N+4 + 2*mueff/N);  % Time constant for cumulation for C
    cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
    c1 = 2 / ((N+1.3)^2+mueff);    % Learning rate for rank-one update of C
    cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
    damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % Damping for sigma
    % usually close to 1
end
    
% -------------------- Generation Loop --------------------------------
% %%%% Plotting for validation 
% k_mesh = linspace(paramMins(1),paramMaxes(1),500);
% theta_mesh = linspace(paramMins(2),paramMaxes(2),500);
% [K,Theta] = meshgrid(k_mesh,theta_mesh);
% x1 = K;
% x2 = Theta;
% 
% term1 = x1.^2;
% term2 = 2*x2.^2;
% term3 = -0.3 .* cos(3*pi.*x1);
% term4 = -0.4 .* cos(4*pi*x2);
% 
% y = term1 + term2 + term3 + term4 + 0.7;
% 
% figure(2)
% surf(k_mesh,theta_mesh,y)
% brighten(.5)
% colormap('parula')
% hold on
% contour3(k_mesh,theta_mesh,y,10)
% initial_y = boha1(paramMeans);
% colors = cool(generations+1);
% scatter3(paramMeans(1),paramMeans(2),initial_y,40,colors(1,:),'filled')
% view(0,90)
% set(gca,'clim',[0,.5*max(max(y))])
% xlabel 'k'
% ylabel '\theta'
% Parameter progression


    CMA_params = figure;
    % met_log = [];
    colormap('hsv')
    colors = hsv(generations+1);
    i=1; %initialize each time called
for i = 1:N
    subplot(N,1,i)
    scatter(0,paramMeans(i),40,colors(1,:))
    hold on
    ylabel(paramNames{i})
%     ylim(variable1_range_norm)
%     ylim(variable2_range_norm)
    xlabel 'Iteration'
    i=i+1;
end
% counteval = 0;  % Number of evalution finished
% x_log = zeros(N,lambda*generations);
while (counteval <= stopeval-1)
    
    
        %     metrates = zeros(lambda,1);
        % display to command
        sprintf('Sigma: %.3f',sigma)
        % Make sure paramMean_norm is a column vector
        if isrow(paramMean_norm)
            paramMean_norm = paramMean_norm';
        end
        % Generate and evaluate lambda offspring one-by-one
        for k = 1:lambda
            counteval = counteval + 1;

    %%%%%%%%%%%% Need to use each variable's sigma to generate means
            % x holds all normalized candidate parameter sets as columns with lambda
            % number of columns
            x(:,k) = paramMean_norm + (sigma * B * (D .* randn(N,1))); % m + sig * Normal(0,C)
            % Apply hard constraints to the random generated parameters.
            x(:,k) = apply_constraints_impedance(x(:,k),paramBounds,N);      
    %         % Re-roll parameter sets of there are repeats
    %             test = 1; rep_count = 0;
    %             while rep_count > 0 || test == 1
    %                 test = 0;
    %                 rep_count = 0;
    %                 for l = 1:length(x_log)
    %                     if x(:,k) == x_log(:,l)
    %                         rep_count = rep_count + 1;
    %                     end
    %                 end
    %                 if rep_count > 0
    %                     sprintf('Has repeat')
    %                     x(:,k) = xmean_norm + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
    %                     x(:,k) = apply_constraints_V2(x(:,k),param_bounds,N,scale_norm);
    %                 end
    %             end
    %         x_log(:,counteval) = x(:,k);
        end   


        % Manually evaluate lambda offspring



        % Collect all metabolic data in one run
        gen_num = ceil(counteval/lambda); 
        [metrates] = metabolic_rate_estimation_GUI_v1(x,gen_num,paramNames,fileInfo,fileInfo.speed,N);
    %     metrates= [10,20,30,40,50,60];
    %     metrates = zeros(lambda,1);
    %     for l = 1:lambda
    %         metrates(l) = boha1(x(:,l));
    %     end
        % Sort by fitness and compute weighted mean into xmean
        [sorted_metrates, idx] = sort(metrates); % Minimization

        % Storing variables
        CMAES_struct(gen_num).gen_idx = gen_num;
        CMAES_struct(gen_num).x_std = sigma;
        CMAES_struct(gen_num).min_metrate = sorted_metrates(1);
        CMAES_struct(gen_num).x_mean = paramMean_norm.*paramRange_size;
        CMAES_struct(gen_num).xmean_norm = paramMean_norm;
        CMAES_struct(gen_num).sigma = sigma;
        CMAES_struct(gen_num).x_std = std(paramRange_size.*x+paramRange(:,1),0,2);
        CMAES_struct(gen_num).x = paramRange_size.*x+paramRange(:,1);       
        CMAES_struct(gen_num).metrates = metrates;
        CMAES_struct(gen_num).pc = pc;
        CMAES_struct(gen_num).ps = ps;
        % Candidate Plotting
            i=1;
            bounded = 'Bounded Parameters:\n';
            figure(CMA_params.Number)
        for i = 1:N
            subplot(N,1,i)
            hold on
            scatter(repmat(gen_num,1,lambda),x(i,:),40,colors(gen_num,:))
            scatter(gen_num,CMAES_struct(gen_num).x_mean(i),50,'k')        
            ylim(paramRange(i,:))
    %     sprintf('Bounded Parameters:\n %.3f,%.3f\n',x)
            bounded = [bounded '%.3f,'];  %iteratively create an appropriately long string
            i=i+1;
        end

        sprintf([bounded '\n'],paramRange_size.*x)   %this was the statement works for all param #'s 

        xold = paramMean_norm;
        paramMean_norm = x(:,idx(1:mu))*weights;   % Recombination, new mean value
    %%%%%%% Evolution paths need to update each variable correctly
        % Cumulation: Update evolution paths
        ps = (1-cs)*ps ...
            + sqrt(cs*(2-cs)*mueff) * invsqrtC * (paramMean_norm-xold) / sigma;
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
        pc = (1-cc)*pc ...
            + hsig * sqrt(cc*(2-cc)*mueff) * (paramMean_norm-xold) / sigma;

        % Adapt covariance matrix C
        artmp = (1/sigma) * (x(:,idx(1:mu))-repmat(xold,1,mu));
        C = (1-c1-cmu) * C ...                  % Regard old matrix
            + c1 * (pc*pc' ...                 % Plus rank one update
            + (1-hsig) * cc*(2-cc) * C) ... % Minor correction if hsig==0
            + cmu * artmp * diag(weights) * artmp'; % Plus rank mu update

        % Adapt step size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));

        % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
        if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
            eigeneval = counteval;
            C = triu(C) + triu(C,1)'; % Enforce symmetry
            [B,D] = eig(C);           % Eigen decomposition, B==normalized eigenvectors
            D = sqrt(diag(D));        % D is a vector of standard deviations now
            invsqrtC = B * diag(D.^-1) * B';
        end
        CMAES_struct(gen_num).New_x_mean = paramMean_norm.*paramRange_size;
        CMAES_struct(gen_num).New_x_mean_norm = paramMean_norm;
        CMAES_struct(gen_num).New_sigma = sigma;
        CMAES_struct(gen_num).New_pc = pc;
        CMAES_struct(gen_num).New_ps = ps;
        save(filename) % Save mat file for each trial to review data
        % Stop CMA-ES optimization if new mean is within 5% of the old mean 
    %     if mean(abs(xmean-xold)./xold) <= .05
    %         break
    %     end

        %pause
%         sprintf('You have completed %d generations.\nStop and rest for 1-minute\nHit ENTER when ready.',gen_num);
         pause(3);
% 
%         Intext = input('\nType stop to stop experiment:\n','s');
%         if strcmpi(strtrim(Intext),'stop')
%             error('Trial Stopped'); 
%         end



    % stats = table(gen_idx,gen_xmean,min_metrates,gen_stds,...
    %     lambda1,lambda2,lambda3,lambda4,metrate1,metrate2,metrate3,metrate4);
   
    
end
xmin = paramMean_norm; % Return the last variable distribution as an estimate of the optimal value
sprintf('You have completed the trial after %d generations!',gen_num)
end
% catch 
%     sprintf('Error occurred, saved data')
%     save(filename,'stats_struct')
% end
% % Uniqueness test for repeated values
% A = x_log(1,:);
% B = x_log(2,:);
% 
% [vA, wA] = unique( A, 'stable' );
% duplicate_indices_A = setdiff( 1:numel(A), wA )
% [vB, wB] = unique( B, 'stable' );
% duplicate_indices_B = setdiff( 1:numel(B), wB )

% ---------------------------------------------------------------
% function f=getmetabolicrates(x)
% % This is a pseudo function used to reprent the process of metabolic rate
% % evalution in actual human-in-the-loop optimization. In actual operation,
% % this process was accomplished by measuring the human subject's resparitory
% % flow rates of oxygen and carbon dioxide, converting them to raw metabolic
% % rate measurements and getting a single estimate of the actual metabolic
% % rate while the subject is walking with exoskeleton under one fixed
% % assistance condition. The assistance condition was defined by parameter
% % x, with a corresponding torque curve (see torquecurve_example.m) defined
% % by HLCParams = [x(1) 0.5*x(2) x(3) 0.5*x(4)];
% 
% f = 250 + 100*randn(1); % just random noise, for illustration purposes only



%% REFERENCES
%
% Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% Strategy on Multimodal Test Functions.  Eighth International
% Conference on Parallel Problem Solving from Nature PPSN VIII,
% Proceedings, pp. 282-291, Berlin: Springer.
% (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
%
% Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% 9(2), pp. 159-195.
% (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
%
% Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% Time Complexity of the Derandomized Evolution Strategy with
% Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
%