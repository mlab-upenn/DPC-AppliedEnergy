function [input] = dpc(model, leafmodels, k, xtot, uprev, d, orderAR, N, umin, umax, xmin, xmax, nc, numRooms, tref,P,Q)

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

nu = nc*N;
H = 2*(coeff{5}*P*eye(N)*coeff{5}' + coeff{1}*Q(1)*eye(N)*coeff{1}' + coeff{2}*Q(2)*eye(N)*coeff{2}' + coeff{3}*Q(3)*eye(N)*coeff{3}' + coeff{4}*Q(4)*eye(N)*coeff{4}');
H = (H + H')/2;
f = -2*(tref1*Q(1)*eye(N)*coeff{1}' + tref2*Q(2)*eye(N)*coeff{2}' + tref3*Q(3)*eye(N)*coeff{3}' + tref4*Q(4)*eye(N)*coeff{4}');
       
lb = [1 repmat(umin,N,1)']';
ub = [1 repmat(umax,N,1)']';
A = [0 -1 80 0 0 0 0 0 0;
     0 0 0  0 0 0 0 0 0;
     0 0 0  0 0 0 0 0 0;
     0 0 0  0 0 0 0 0 0;
     0 -1 0  0 0 0 0 0 0 0 0 0 0;
     0 -1 0  0 -1 0 0 -1 0 0 -1 0 0;
     0 -1 0  0 -1 0 0 -1 0 0 -1 0 0;
     0 -1 0  0 -1 0 0 -1 0 0 -1 0 0;];
b = -umin(1);
options=optimset;
% options.Algorithm='interior-point-convex'; % to get rid of some warnings.
options.Display='off';
[U, fval, exitflag] = quadprog(H,f,A,b,[],[],lb,ub,[],options);

if ne(exitflag,1)
    disp('Problem is:');
    exitflag
    return
end

input = U(2:nc+1);