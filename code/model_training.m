%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to do before compiling: change nc - Choose between 2 inputs and 3 inputs
%                                check orderAR and ctrlHzn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

% set params
orderAR = 2; % Autoregression order
ctrlHzn = 1; % Control horizon
nc = 1;

%% Creating random forest models
numRooms = 4;
model = cell(numRooms+1,ctrlHzn); % The +1 is for the power component

warning('off','all')
leafmodels = cell(numRooms+1,ctrlHzn);

% Creating models for the rooms
for ii = 1:numRooms %% 5 is the number of rooms in ti_intermp_10min
   
    % Create random forests
    catcol = sort([7*(1:ctrlHzn)-1 7*(1:ctrlHzn)]);
    leafsize = 20;
    disp('Room models computation');
    for idm = 1:ctrlHzn % separate regression tree for each output
        disp(['Loading data for room: ',num2str(ii)]);
        [xTrain_d, xTrain_c, yTrain] = load_data_new_states('train', orderAR, ctrlHzn,idm,ii,nc);
        nc = size(xTrain_c,2)/ctrlHzn; % Number of control inputs
    
            % Remove NaN values
            total = [yTrain,xTrain_d,xTrain_c];
            [m,n] = size(xTrain_d);
            [m,p] = size(total);
            total(any(isnan(total),2),:) = [];
            yTrain = total(:,1:ctrlHzn);
            xTrain_d = total(:,ctrlHzn + 1:ctrlHzn + n);
            xTrain_c = total(:,ctrlHzn + n + 1:p);
        
        % Create random forest
        disp(['Random forests creation for room: ',num2str(ii),' - Horizon: ',num2str(idm)]);
        model{ii,idm} = TreeBagger(50, xTrain_d, yTrain(:,idm),...
        'Method', 'regression', 'OOBPred', 'On', 'OOBVarImp', 'on',...
        'CategoricalPredictors', catcol, 'MinLeaf', leafsize+idm);
        
        % Create linear models in leaves
        disp(['Training linear models in leaves for room ',num2str(ii),' - Horizon ',num2str(idm)]);
        leafmodels{ii,idm} = train_linearmodel_in_leaves(model{ii,idm}, xTrain_c, yTrain(:,idm), idm, nc);
    end
end

% creating model for the power

% Create random forest
catcol = sort([7*(1:ctrlHzn)-1 7*(1:ctrlHzn)]);
leafsize = 20;

% separate regression tree for each output
disp('Power model computation');
for idm = 1:ctrlHzn
    disp('Loading data for power');
    [xTrain_d, xTrain_c, yTrain] = load_data_new_power('train', orderAR, ctrlHzn,idm,nc);
    nc = size(xTrain_c,2)/ctrlHzn; % Number of control inputs
    
            % Remove NaN values
            total = [yTrain,xTrain_d,xTrain_c];
            [m,n] = size(xTrain_d);
            [m,p] = size(total);
            total(any(isnan(total),2),:) = [];
            yTrain = total(:,1:ctrlHzn);
            xTrain_d = total(:,ctrlHzn + 1:ctrlHzn + n);
            xTrain_c = total(:,ctrlHzn + n + 1:p);
    
    % Create random forest for power
    disp(['Random forest creation for power',' - Horizon: ',num2str(idm)]);
    model{numRooms+1,idm} = TreeBagger(50, xTrain_d, yTrain(:,idm),...
    'Method', 'regression', 'OOBPred', 'On', 'OOBVarImp', 'on',...
    'CategoricalPredictors', catcol, 'MinLeaf', leafsize+idm);
    
    % Create linear models in leaves for power
    disp(['Training linear models in leaves for power - Horizon ',num2str(idm)]);
    leafmodels{numRooms+1,idm} = train_linearmodel_in_leaves(model{numRooms+1,idm}, xTrain_c, yTrain(:,idm), idm, nc);
end

clearvars -except model leafmodels ctrlHzn orderAR

disp('Saving models');
save(['data/models_random_forest_N=' num2str(ctrlHzn) '_AR=' num2str(orderAR) '.mat'],'model')

disp('Saving leaves models');
save(['data/lin_random_forest_N=' num2str(ctrlHzn) '_AR=' num2str(orderAR) '.mat'],'leafmodels')
