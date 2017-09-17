%% Convert data to adapt to dpc code

clear all
close all
clc

load('Data_PreProccessing/cetemps_interp_10min.mat')
load('Data_PreProccessing/et_interp_10min.mat')
load('Data_PreProccessing/ti_interp_10min.mat')

sTrain = 67; % Training start index
eTrain = 6734; % Training end index
sTest = 7346; % Testing start index
eTest = 9505; % Testing end index

%% Creating "tod - time of day" and "dow - day of week" for the training and testing processes
tTrain = datenum(ti_interp_10min(sTrain:eTrain,1)); % Time vector for the train dataset
temp = datevec(tTrain);
todTrain = temp(:,4)';
[dowTrain ~] = weekday(tTrain);
dowTrain = dowTrain';

tTest = datenum(ti_interp_10min(sTest:eTest,1));
temp = datevec(tTest);
todTest = temp(:,4)';
[dowTest ~] = weekday(tTest);
dowTest = dowTest';

clear temp

%% Creates data sets

% Creates training data

x = cell2mat(ti_interp_10min(sTrain:eTrain,2:end))';
d = cell2mat(cetemps_interp_10min(sTrain:eTrain,2:end))';
d = [d;todTrain;dowTrain];
u = cell2mat(et_interp_10min(sTrain:eTrain,[5]))';
y = cell2mat(et_interp_10min(sTrain:eTrain,6))';

save('train.mat','d','u','x','y');
save('tTrain','tTrain')
% Creates testing data

x = cell2mat(ti_interp_10min(sTest:eTest,2:end))';
d = cell2mat(cetemps_interp_10min(sTest:eTest,2:end))';
d = [d;todTest;dowTest];
u = cell2mat(et_interp_10min(sTest:eTest,[5]))';
y = cell2mat(et_interp_10min(sTest:eTest,6))';

save('test.mat','d','u','x','y');
save('tTest','tTest')