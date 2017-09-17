%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to do before compiling: change path in the load_data_new_states and power - Choose between 2 inputs folder and 3 inputs folder
%                                check orderAR and ctrlHzn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orderAR = 2;
ctrlHzn = 1;
nc = 1;
modelIndex = 3; % 1 for roomT1,...,4 for roomT4, 5 for power
if ((modelIndex > 0) && (modelIndex < 5))
    [xTest_d, xTest_c, yTest] = load_data_new_states('test', orderAR, ctrlHzn,1,modelIndex,nc);
elseif modelIndex == 5
    [xTest_d, xTest_c, yTest] = load_data_new_power('test', orderAR, ctrlHzn,1,nc);
end
nc = size(xTest_c,2)/ctrlHzn;

n = 1000;
xTest_d = xTest_d(1:n,:);
xTest_c = xTest_c(1:n,:);
yTest = yTest(1:n,:);

% calculate root mean square error with our definition
yTrue = yTest;
yPred = zeros(size(xTest_d,1),ctrlHzn);

disp('Loading ramdom forest models')
load(['data/',num2str(nc),' inputs/models_random_forest_N=1_AR=2.mat']);
disp('Loading linear models in the leaves')
load(['data/',num2str(nc),' inputs/lin_random_forest_N=1_AR=2.mat']);
disp('Loading done!')

% find linear models
for idn = 1:size(xTest_d,1)
    for idm = 1
        coeff = find_linearmodel_in_leaves(model{modelIndex,idm}, leafmodels{modelIndex,idm}, xTest_d(idn,:));
        yPred(idn,idm) = [1, xTest_c(idn,1:nc*idm)]*coeff;
    end
    disp(idn);
end

yPredNaN = isnan(yPred);
yPred(yPredNaN) = [];
yTrue(yPredNaN) = [];
yTrueNaN = isnan(yTrue);
yPred(yTrueNaN) = [];
yTrue(yTrueNaN) = [];
% nrmseTest = sqrt(mean((yPred-yTrue).^2,1))./mean(yTrue,1);
% fprintf('nrmse on test data using our definition = %f\n', min(nrmseTest));

[a,b]=rsquare(yTrue,yPred); % Error computation
fprintf('Random Forests 2013(Testing) RMSE(W): %.2f, R2: %.3f, RMSE/peak %0.4f, NRMSD: %0.2f \n\n'...
    ,b,a,(b/max(yTrue)),(100*b/(max(yTrue)-min(yTrue))));

plot_ml_results(yTrue, yPred', 'true-predicted', ctrlHzn);
plot_ml_results(yTrue, yPred', 'predictive', ctrlHzn);

%% Validation on Energy consumption
Energy = zeros(size(yPred));
Energy(1) = yPred(1)/6;
for ii = 2:length(yPred)
    Energy(ii) = yPred(ii)/6 + Energy(ii-1);
end
sTest = 7346; % Testing start index
eTest = 9505; % Testing end index
load('et_interp_10min')
YtestEnergy = cell2mat(et_interp_10min(sTest+1:eTest,2));
tTest = datenum(et_interp_10min(sTest:eTest-1,1));
figure(4);
YtestEnergy = YtestEnergy(1:n);
YtestEnergy(yPredNaN) = [];
YtestEnergy(yTrueNaN) = [];
tTest = tTest(1:n);
tTest(yPredNaN) = [];
tTest(yTrueNaN) = [];
YtestEnergy = YtestEnergy-YtestEnergy(1);
plot(tTest,YtestEnergy,'k',tTest,Energy,'r','LineWidth',1.5);
datetick('x','dd/mm');
xlabel 'Time' ;
ylabel 'Thermal Energy [kWh]';
legend('true','predicted')

[a,b]=rsquare(YtestEnergy,Energy); % Error computation
fprintf('Random Forests 2013(Testing) RMSE(W): %.2f, R2: %.3f, RMSE/peak %0.4f, NRMSD: %0.2f \n\n'...
    ,b,a,(b/max(YtestEnergy)),(100*b/(max(YtestEnergy)-min(YtestEnergy))));