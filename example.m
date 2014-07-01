clear;
clc;

% Preparation work
% add the class in default path
STARTUP % remember to change the path name in the .m file to where you put the S_Race class
Delta = 0.9;
load('results.mat')
[max_step,M,no_obj] = size(results);

% Initial a S_Race object and begin racing
obj = S_Race_v2(M,max_step - 9,Delta,1);

% we won't start racing until we have collected 10 samples
retained = obj.Racing(results(1:10,:,:));

% then each step we input one test instance
for i = 11:max_step
    retained = obj.Racing(results(i,retained,:));
end
