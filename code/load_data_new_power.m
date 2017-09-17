% idm is the tree number on the horizon we want to load disturbance for

function [features_d, features_c, outputs] = load_data_new_power(argStr, orderAR, ctrlHzn, idm, nc)

if strcmp(argStr, 'train')
    load(['Progetto Tullio/DPC/data/',num2str(nc),' inputs/train.mat']);
elseif strcmp(argStr, 'test')
    load(['Progetto Tullio/DPC/data/',num2str(nc),' inputs/test.mat']);
else
    error('input string should be either ''train'' or ''test.');
end

features_d = [lagmatrix(d', 0-idm+1:1:ctrlHzn-1-idm+1),...
              lagmatrix(x(1,:), 0:1:orderAR-1),...
              lagmatrix(x(2,:), 0:1:orderAR-1),...
              lagmatrix(x(3,:), 0:1:orderAR-1),...
              lagmatrix(x(4,:), 0:1:orderAR-1)]; %#ok<*NODEF>

features_c = lagmatrix(u', -(0:1:ctrlHzn-1));

outputs = lagmatrix(y', -(1:1:ctrlHzn));

features_d(1:orderAR-1,:) = [];
features_c(1:orderAR-1,:) = [];
outputs(1:orderAR-1,:) = [];
features_d(end-ctrlHzn:end,:) = [];
features_c(end-ctrlHzn:end,:) = [];
outputs(end-ctrlHzn:end,:) = [];