function x = updateModel(power,roomT1,roomT2,roomT3,roomT4,d,u,xtot,model,leafmodels)

Xtest = [d' u' xtot'];
powerUpdate  = predict(power,Xtest);
roomT1Update = predict(roomT1,Xtest);
roomT2Update = predict(roomT2,Xtest);
roomT3Update = predict(roomT3,Xtest);
roomT4Update = predict(roomT4,Xtest);

x = [roomT1Update;roomT2Update;roomT3Update;roomT4Update;powerUpdate];
% idm = 1;
% Xtest1 = [d' xtot(1) xtot(1)];
% coeff = find_linearmodel_in_leaves(model{1,idm}, leafmodels{1,idm}, Xtest1);
% roomT1Update = [1, u]*coeff;
% 
% Xtest2 = [d' xtot(2) xtot(2)];
% coeff = find_linearmodel_in_leaves(model{2,idm}, leafmodels{2,idm}, Xtest2);
% roomT2Update = [1, u]*coeff;
% 
% Xtest3 = [d' xtot(3) xtot(3)];
% coeff = find_linearmodel_in_leaves(model{3,idm}, leafmodels{3,idm}, Xtest3);
% roomT3Update = [1, u]*coeff;
% 
% Xtest4 = [d' xtot(4) xtot(4)];
% coeff = find_linearmodel_in_leaves(model{4,idm}, leafmodels{4,idm}, Xtest4);
% roomT4Update = [1, u]*coeff;
% 
% XtestP = [d' xtot(1) xtot(1) xtot(2) xtot(2) xtot(3) xtot(3) xtot(4) xtot(4)];
% coeff = find_linearmodel_in_leaves(model{5,idm}, leafmodels{5,idm}, XtestP);
% powerUpdate = [1, u]*coeff;
% 
% x = [roomT1Update;roomT2Update;roomT3Update;roomT4Update;powerUpdate];