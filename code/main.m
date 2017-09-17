%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to do before compiling: comment or add input constraints for u3 in dpcCVX
%                                change nc = 1 or nc = 2 or nc = 3
%                                check orderAR and ctrlHzn
%                                reload models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Before running this function, make sure you run model_training
clearvars -except model leafmodels power roomT1 roomT2 roomT3 roomT4
close all
clc
addpath(genpath('D:\PostDoc\Matlab\cvx'))
cvx_solver gurobi

% addpath(genpath('C:\Program Files\IBM'))

orderAR = 2; % This has to be the same used in "model_training.m"
ctrlHzn = 4; % This has to be the same used in "model_training.m"
index = max(orderAR,ctrlHzn); % Used to create the initial state and input vectors of appropriate dimensions
start = 1;
nc = 1;
numRooms = 4;
simTime = 2150;
xtot = zeros(numRooms+1,simTime);
utot = zeros(nc,simTime);
epsilontotm = zeros(4,simTime);
epsilontotM = zeros(4,simTime);
temptot = zeros(4,simTime);
%% Load models
% disp('Loading ramdom forest models')
% load(['data/',num2str(nc),' inputs/models_random_forest_N=4_AR=2.mat']);
% disp('Loading linear models in the leaves')
% load(['data/',num2str(nc),' inputs/lin_random_forest_N=4_AR=2.mat']);
% disp('Loading testing data')
load(['data/',num2str(nc),' inputs/test.mat']);
load(['data/',num2str(nc),' inputs/tTest.mat']);

% disp('Loading real plant data')
% load(['data/',num2str(nc),' inputs/real plant/power.mat']);
% load(['data/',num2str(nc),' inputs/real plant/roomT1.mat']);
% load(['data/',num2str(nc),' inputs/real plant/roomT2.mat']);
% load(['data/',num2str(nc),' inputs/real plant/roomT3.mat']);
% load(['data/',num2str(nc),' inputs/real plant/roomT4.mat']);
% disp('Loading done!')
x = x(:,start:end);
d = d(:,start:end);
u = u(:,start:end);
y = y(:,start:end);
tTest = tTest(start:end);
tempTot = [];

%% Run closed loop simulation

% Define optimization parameters

% tref = [22;22;22;22];
% tref = x;
tref = zeros(2,simTime + index);
% for k = 1:simTime + index
%     if ( (d(6,k) >= 8) && (d(6,k) <= 23) )
%         tref(1,k) = 21;
%         tref(2,k) = 22.5;
%     else
%         tref(1,k) = 20;
%         tref(2,k) = 22.5;
%     end
% end
for k = 1:simTime + index
    if ( (d(6,k) >= 7) && (d(6,k) <= 9) || (d(6,k) >= 18) && (d(6,k) <= 23) )
        tref(1,k) = 21;
        tref(2,k) = 22.5;
    else
        tref(1,k) = 20;
        tref(2,k) = 22.5;
    end
end

% tref = 21.5*ones(1,simTime + index);

if nc == 1
    umin = 0;
    umax = 0.35;
elseif nc == 2
    umin = [24;0];
    umax = [65;0.35];
elseif nc == 3
    umin = [24;24;0];
    umax = [65;65;0.35];
else
    disp('Error with inputs bounds')
    return
end

P = 1e2; % Power weight
Q = [3e3 1e2 1e5 1e5]; % Temperatures weigths for T1,T2,T3,T4

% Define initial state
xtot(:,1:index) = [x(1:4,1:index)+1; y(1:index)]; % Power consumption and temperature of the 4 rooms
                                                % Needed to start the dpc with the autoregressive part
                                                % 1:4 is because we don't want the last room   
% xtot(:,1:index) = [tref(1,1:index)+0.5;tref(1,1:index)+0.5;tref(1,1:index)+0.5;tref(1,1:index)+0.5; zeros(1,index)]; % Power consumption and temperature of the 4 rooms
%                                                 % Needed to start the dpc with the autoregressive part
%                                                 % 1:4 is because we don't want the last room    
utot(:,1:index-1) = u(:,1:index-1);

uDPCtot(:,1:index-1) = u(:,1:index-1);
% u(268:284) = 0;
% Closed loop simulation
for k = index:simTime
    k
    if ( (utot(1,k-2) == 0) && (utot(1,k-1) == 0.35) )
        uprev = 1;
    elseif ( (utot(1,k-3) == 0) && (utot(1,k-2) == 0.35) )
        uprev = 2;
    else
        uprev = 0;
    end
    % Compute control
    [uDPC,epsilonm,epsilonM] = dpcCVX(model, leafmodels, k, xtot, uprev, d, orderAR, ctrlHzn, umin, umax, nc, numRooms, tref, P, Q);
% uDPC = u(:,k); % with this we generate xsim, i.e. the simulation made using real inputs
epsilontotm(:,k) = epsilonm;
epsilontotM(:,k) = epsilonM;
%     uDPCtot(:,k) = uDPC;
%     if uDPC > 0.05
%         utot(:,k) = umax(1);
%     end
    utot(:,k) = uDPC;
%     utot(:,k) = u(:,k);
    % Update the model
    xtot(:,k+1) = updateModel(power,roomT1,roomT2,roomT3,roomT4,d(:,k),utot(:,k),xtot(1:4,k),model,leafmodels);
%     xtot(:,k+1) = updateModel(power,roomT1,roomT2,roomT3,roomT4,d(:,k),u(:,k),x(1:4,k),model,leafmodels);
    
%     tempTot = [tempTot temp];
end
xOF(:,1:index) = xtot(:,1:index);
uOF(:,1:index-1) = u(:,1:index-1);
% Closed loop simulation
ON = 0;
for k = index:simTime
    k
    % Compute ON/OFF control
    if ON
        if xOF(3,k) > tref(2,k)
            uOFk = umin(1);
            ON = 0;
        else
            uOFk = umax(1);
        end
    else
        if xOF(3,k) < tref(1,k)
            uOFk = umax(1);
            ON = 1;
        else
            uOFk = umin(1);
        end
    end

    uOF(:,k) = uOFk;
    
    % Update the model
    xOF(:,k+1) = updateModel(power,roomT1,roomT2,roomT3,roomT4,d(:,k),uOF(:,k),xOF(1:4,k),model,leafmodels);
end
%% Plots


lt = size(utot,2);
% lt = 2160;
t = tTest(1:lt);

figure(1)
    plot(t,d(1,1:lt),'m',t,xOF(3,1:lt),'r',t,xtot(3,1:lt),'k',t,50*utot(1,1:lt),'k--',t,50*uOF(1,1:lt),'r*',t,tref(1,1:lt),'g',t,tref(2,1:lt),'g');
    legend('T1','T2','T3','T4','Tf','T1sim','T2sim','T3sim','T4sim')
    datetick('x','HH')
%     plot(t,xtot(5,1:lt),'k',t,xOF(5,1:lt),'k--',t,50*utot(1,1:lt),'r',t,50*uOF(1,1:lt),'r--');
%     legend('PoptInput','PrealInput')
%     datetick('x','HH')

%% Validation on Energy consumption
EnergyDPC = zeros(size(xtot(5,1:lt)));
EnergyDPC(1) = xtot(5,1)/6;
for ii = 2:lt
    EnergyDPC(ii) = xtot(5,ii)/6 + EnergyDPC(ii-1);
end
EnergyOF = zeros(size(xOF(5,1:lt)));
EnergyOF(1) = xOF(5,1)/6;
for ii = 2:lt
    EnergyOF(ii) = xOF(5,ii)/6 + EnergyOF(ii-1);
end
figure(2)
plot(t,EnergyDPC,'k',t,EnergyOF,'b')
% 
% datetick('x','HH')

% figure(1)
% plot(t,xtot(1,1:lt),'r',t,x(1,1:lt),'k');
% figure(2)
% plot(t,xtot(2,1:lt),'r',t,x(2,1:lt),'k');
% figure(3)
% plot(t,xtot(3,1:lt),'r',t,x(3,1:lt),'k');
% figure(4)
% plot(t,xtot(4,1:lt),'r',t,x(4,1:lt),'k');
% figure(5)
% plot(t,xtot(5,1:lt),'r',t,y(1:lt),'k');
