%%  Constant Metabolic Cost Estimation
% This function is a least-squares system identification of a constant
% input to a first-order dynamic system .
% This script finds a constant actual metabolic rate 'estimated_dot_E' that result in
% the least squared error between the predicted system response y_bar and
% the series of measurements y_meas.

%INPUTS
%x: n-parameters by lambda matrix (or vector if lambda = 1) of parameter inputs, each input will be used for a
%2-min data collection set to estimate steady state metabolic cost

%OUTPUTS
% estimated_dot_E (length(x) element vector): the estimated metabolic rate of the period.
% y_bar (length(x) by n-measurements matrix): the best fit system outputs predicted by the polynomial
% relationship

% The system is modeled as a first order discrete dynamical system such that
% y(i) = (1-dt/tau)*y(i-1) + (dt/tau)*u(c,p)
% dt = t(i) - t(i-1)

%The system is equivalent to the forward Euler integration of a first-order
%continuos system of the form
% y_dot = (1/tau)*(u-y)

%The function in this script identifies the vector x_star that results in
%the least squared-error between y_bar and y_meas.

%x_star is the optimal solution to a specially formulated matrix equation
% A*x = y_meas, where x = [y_bar(1) u]'

% x_star can be found using the pseudo-inverse of A
% x_star = pinv(A)*y_meas .

% y_bar can be found by using x_star
% A*x_star = y_bar

% Adapted from the function of Wyatt Felt and C. David Remy at
% https://cn.mathworks.com/matlabcentral/fileexchange/51328-instantaneuos-cost-mapping
% Copyright@Juanjuan Zhang, Steven H Collins 11/30/2016

function [estimated_dot_E] = metabolic_rate_estimation_GUI_v1(x,gen_num,variables,fileInfo,speed,N)
    %     x = 1;
    % Initialize parameters
    %meta_est_plot = figure
    global continueFlag pauseFlag stopFlag;
    filename = sprintf('metabolic_rate_est_%s_%s_%s_gen%d_speed%s.mat',fileInfo.controller, fileInfo.subjectID,fileInfo.date,gen_num,fileInfo.speed);
    if isfile(filename)
        load(filename);
    end
    lambda = length(x);
    tau = 42; % time constant for fit in seconds (From Zhang, JJ... Collins, S)
    estimated_dot_E = zeros(lambda,1);
    mean_squared_error_int = estimated_dot_E;
    y_bar_int = cell(1,lambda);
    y_meas_int = y_bar_int;
    time = y_bar_int;
    %webcam text recognition params
    test_mode = 2; mins = 2; test_method = 0; widthBOX = 300; heightBOX = 115;

    if ~exist('trialsDone') %Creates Trial counter if first time through
       trialsDone = 0; 
    end

    for c = 1:lambda % testing for each candidate set
        if ~pauseFlag
            if c > trialsDone %Skips existing trials
                % Display testing parameters
        %         if N == 4
        %             sprintf('Parameter Set number: %d\nSet %s to %.3f\nSet %s to %.3f.\nSet %s to %.3f\nSet %s to %.3f.\nPress any key when ready'...
        %                 ,c,variables{1},x(1,c),variables{2},x(2,c),variables{3},x(3,c),variables{4},x(4,c))
        %             pause
        %         else 
        %             sprintf('Parameter Set number: %d\nSet %s to %.3f\nSet %s to %.3f.\nPress any key when ready'...
        %                 ,c,variables{1},x(1,c),variables{2},x(2,c))
        %             pause
        %         end

        %         % Send New Control Laws
        %         if N == 4
        %             setparam(slrt,'Subcontroller/Desired Torques/Ext_FFFB_Ratio','Value',x(1))
        %             setparam(slrt,'Subcontroller/Desired Torques/Ext_Gain','Value',x(2))
        %             setparam(slrt,'Subcontroller/Desired Torques/Flex_FFFB_Ratio','Value',x(3))
        %             setparam(slrt,'Subcontroller/Desired Torques/Flex_Gain','Value',x(4))
        %         else
        %             setparam(slrt,'Subcontroller/Desired Torques/Stiffness','Value',x(1))
        %             setparam(slrt,'Subcontroller/Desired Torques/RefAng','Value',x(2))
        %         end
                % Test
                        tests = 1:lambda;
                        estimated_dot_E(c) = griewank(x(:,c));
                        y_bar(c) = tests(c);
                        mean_squared_error(c) = tests(c);
                        t(c)= tests(c);
                        vo2_rate(c)= tests(c);
                        vco2_rate(c)= tests(c);
                        x_star(c)= tests(c);
                        mean_squared_error(c)= tests(c);
                        y_bar(c)= tests(c);
                        y_meas(c)= tests(c);

                %----------------------------------------------------------
        %         meta_test = 1;
        %         while meta_test == 1
        %             % Set metronome speed and collect metabolic data
        %             [t,vo2_rate,vco2_rate] = getmetabolics_V7(test_mode,test_method,mins,widthBOX,heightBOX,gen_num*lambda+c);
        %             time{1,c} = t;
        %             % Calculate respiratory output
        %             y_meas_int{1,c} = .278*vo2_rate + .075*vco2_rate;
        %             y_meas = y_meas_int{1,c};
        %             
        %             %Reshape the measurement vector if needed
        %             if isrow(y_meas)
        %                 y_meas = y_meas';
        %             elseif ~isrow(y_meas) & ~iscolumn(y_meas)
        %                 error('Measurements are not in a single column vector')
        %             end
        %             
        %             % Generate the matrix A
        %             n_samp = length(t);
        %             A = zeros(n_samp,2);
        %             A(1,:) = [1,0];
        %             for i = 2:length(t)
        %                 for j = 1:2
        %                     dt = t(i)-t(i-1);
        %                     if j == 1
        %                         A(i ,j) = A(i-1,j)*(1-dt/tau);
        %                     else
        %                         A(i ,j) = A(i-1,j)*(1-dt/tau) + (dt/tau);
        %                     end
        %                 end
        %             end
        %             
        %             %solve for the optimal parameters
        %             x_star = pinv(A)*y_meas; % should be a 2x1 column vector
        %             %solve for the best-fit predicted response (respiratory output)
        %             y_bar_int{1,c} = A*x_star;
        %             y_bar = y_bar_int{1,c};
        %             %find the error between the best-fit predicted response and the
        %             %measurement vector
        %             mean_squared_error = ((y_bar-y_meas)'*(y_bar-y_meas))/n_samp;
        %             mean_squared_error_int(c) = mean_squared_error;
        %             %solve for the optimal parameters
        %             %measurement vector
        %             estimated_dot_E(c) = x_star(2);
        %             sprintf('Metabolic Rate for NMM %.3f: %.3f',x(:,c),estimated_dot_E(c))
        %             figure(meta_est_plot)
        %             plot(t,y_meas,'ko')
        %             hold on
        %             plot(t,y_bar,'r-')
        %             hold off
        %             reply = input('Good Estimation?');
        %             if isempty(reply)
        %                 meta_test = 0;
        %             end
        %         end
                %----------------------------------------------------------
            % Candidate metabolic storage
                trialsDone = trialsDone + 1;
                sprintf('Metabolic Cost: %.2f',estimated_dot_E(c))
                save(filename,'t','vco2_rate','vo2_rate','x_star','mean_squared_error','y_bar','y_meas','trialsDone');
            %   save(sprintf('metabolic_rate_est_NMM_4D_Habit_Emily_11222019_gen1_1_%s_ZI',speed),'t','vco2_rate','vo2_rate','x_star','mean_squared_error','y_bar','y_meas')
            %   save(sprintf('metabolic_rate_est_NMM_4D_Habit_Emily_11222019_gen1_1_%s_Opt',speed),'t','vco2_rate','vo2_rate','x_star','mean_squared_error','y_bar','y_meas')
            end
            
        else
            pauseFlag = false;
            paused = true;
            while paused
               pause(0.1);
               if continueFlag
                   paused = false;
                   continueFlag = false;
                   stopFlag = false;
               elseif stopFlag
                   paused = false;
                   continueFlag = false;
                   stopFlag = false;
                   error('Experiment Stopped');
               end
            end
%             sprintf('Experiment Paused.\nYou have completed %d generations.',gen_num)
%             Intext = input('\nType STOP to stop experiment\nOtherwise, hit ENTER to continue\n','s');
%             if strcmpi(strtrim(Intext),'stop')
%                 pauseFlag = false;
%                 break;
%                 %error('Experiment Stopped'); 
%             else
%                 pauseFlag = false;
%             end
        end

    end
    save(filename,'t','vco2_rate','vo2_rate','x_star','mean_squared_error','y_bar','y_meas','trialsDone')
    %     % Optimal metabolics storage
    %     save (['test2_meta_est_Imp_RJ03072019_' num2str(gen_num) '.mat'], 'y_meas_int' ,'y_bar_int' ,'mean_squared_error_int' ,'time' ,'estimated_dot_E')
end

   
