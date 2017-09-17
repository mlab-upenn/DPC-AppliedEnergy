function [input] = dpcYalmip(model, leafmodels, k, xtot, uprev, d, orderAR, N, umin, umax, xmin, xmax, nc, numRooms, tref,P,Q)

%% calculate leaf coeffs
nx = numRooms+1;
d_cur = d(:,k-N+1:k); % features corr to disturbances

% Creating empty cell for coefficients: 1 for the power and one for each room
coeff = cell(numRooms+1,1);
for ii = 1:numRooms+1
    coeff{ii} = zeros(nc*N+1,N);
end

% Computing coefficients for room models
for ii = 1:numRooms
    % prepare features_d to find an appropriate leaf for rooms
    x_cur = fliplr(xtot(ii,k-orderAR+1:k)); % features corr to indoor temp.
    features_d = [d_cur(:)', x_cur];
    for idm = 1:N % Because the first column contains a0+a1*u(0), the second b0+b1*u(0)+b2*u(1) and so on
        coeff{ii}(1:nc*idm+1,idm) = find_linearmodel_in_leaves(model{ii,idm}, leafmodels{ii,idm}, features_d);
    end
end
% prepare features_d to find an appropriate leaf for power
d_cur = d(:,k:k+N-1); % features corr to disturbances
x_cur = [fliplr(xtot(1,k-orderAR+1:k)),...
         fliplr(xtot(2,k-orderAR+1:k)),...
         fliplr(xtot(3,k-orderAR+1:k)),...
         fliplr(xtot(4,k-orderAR+1:k))]; % features corr to indoor temp.
features_d = [d_cur(:)', x_cur];

% Computing coefficients for power
for idm = 1:N % Because the first column contains a0+a1*u(0), the second b0+b1*u(0)+b2*u(1) and so on
   coeff{numRooms+1}(1:nc*idm+1,idm) = find_linearmodel_in_leaves(model{numRooms+1,idm}, leafmodels{numRooms+1,idm}, features_d);
end

tref1 = tref(1,k+1:k+N);
tref2 = tref(2,k+1:k+N);
tref3 = tref(3,k+1:k+N);
tref4 = tref(4,k+1:k+N);
% tref1 = tref(1);
% tref2 = tref(2);
% tref3 = tref(3);
% tref4 = tref(4);
%% solve optimization problem
addpath(genpath('D:\PostDoc\Matlab\cvx'))
nu = nc*N;
U     = sdpvar(nu,1);
beta0 = sdpvar(1,1);
constraints = [];
objective = 0;

M = blkdiag(P*eye(N),Q(1)*eye(N),Q(2)*eye(N),Q(3)*eye(N),Q(4)*eye(N));
objective = [([beta0, U']*coeff{5}),...     
            ([beta0, U']*coeff{1}-tref1),...
            ([beta0, U']*coeff{2}-tref2),... 
            ([beta0, U']*coeff{3}-tref3),... 
            ([beta0, U']*coeff{4}-tref4)]*M*[([beta0, U']*coeff{5})';...
                                            ([beta0, U']*coeff{1}-tref1)';...
                                            ([beta0, U']*coeff{2}-tref2)';...
                                            ([beta0, U']*coeff{3}-tref3)';...
                                            ([beta0, U']*coeff{4}-tref4)'] ;
% objective = ([beta0, U']*coeff{5})      *P*eye(N)   *([beta0, U']*coeff{5})' +...
%             ([beta0, U']*coeff{1}-tref1)*Q(1)*eye(N)*([beta0, U']*coeff{1}-tref1)' +...
%             ([beta0, U']*coeff{2}-tref2)*Q(2)*eye(N)*([beta0, U']*coeff{2}-tref2)' +...
%             ([beta0, U']*coeff{3}-tref3)*Q(3)*eye(N)*([beta0, U']*coeff{3}-tref3)' +...
%             ([beta0, U']*coeff{4}-tref4)*Q(4)*eye(N)*([beta0, U']*coeff{4}-tref4)' ;
constraints = [constraints, umin(1) <= U(1:3:nu) <= umax(1),...
                            umin(2) <= U(2:3:nu) <= umax(2),...
                            umin(3) <= U(3:3:nu) <= umax(3),...
                            beta0 == 1];



ops = sdpsettings('solver','cplex','verbose',1);  
% ops = sdpsettings('verbose',1); 
solution=solvesdp([],objective,ops);
if ~strcmp(solution.info,'Successfully solved (CPLEX)')
    display('Problem is: ')
    solution.info
end

% Convert input variables

input = zeros(nc,1);
for ii = 1:nc
    input(ii) = double(u{ii});
end

% Eliminate NaN variable

% for ii=1:N
%     nan_inputs=isnan(input{ii});
%     elements=find(nan_inputs);
%     for kk=1:length(elements)
%         input{ii}(elements(kk))=0;
%     end
% end


% % Convert state variables
% 
% state1=cell(1,room);
% 
% for jj=1:room
%     for ii=1:N+1
%             state1{jj}{ii}=double(x{jj}{ii});
%     end
% end
% 
% % Take the 1st state
% 
% real_state=cell(1,room);
% for jj=1:room
%     real_state{jj}=state1{jj}{1};
% end