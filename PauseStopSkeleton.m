 %While loop starts
     if ~pauseFlag
        %Code that does what you want in while loop here

     else
                pauseFlag = false;
                paused = true;
                while paused
                   pause(0.1); %small time delay to not crash matlab
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
     end
 %While loop ends