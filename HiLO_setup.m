
%% Initialize Parameters and Update File Information
clearvars fileInfo
fileInfo.controller = 'control';
fileInfo.subjectID = 'subject';
fileInfo.speed = 'speed';
fileInfo.grade = '0deg';

formatOut = 'mmddyyyy';
fileInfo.date = datestr(now,formatOut);

N=2; % number of parameters
sigma = 0.4;
generations = 10;
trial = '1';
%% Run CMAES
dateStart = input('What date did you start testing? Please format as ''mmddyyy''.]\n');
if ~ischar(dateStart)
    dateStart = num2str(dateStart);
end
if strcmp(fileInfo.date,dateStart)
    counteval = 0;  % Number of evalution finished
    filename = sprintf('%s_%s_%s_%s_%s_test%s.mat',fileInfo.controller,...
        fileInfo.subjectID,fileInfo.speed,fileInfo.grade,fileInfo.date,trial);
    save(filename);
else
    fileInfo.date = dateStart;
    filename = sprintf('%s_%s_%s_%s_%s_test%s.mat',fileInfo.controller,...
        fileInfo.subjectID,fileInfo.speed,fileInfo.grade,fileInfo.date,trial);
    sigma = load(filename,'sigma');

end

[xmin,stats_struct]=CMAES_norm_v3(filename,sigma,generations,N);
title(sprintf('Sigma: %0.1f',sigma))