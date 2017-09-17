function [U,epsilonm,epsilonM] = dpcCVX(model, leafmodels, k, xtot, uprev, d, orderAR, N, umin, umax, nc, numRooms, tref,P,Q)

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

% tref1 = tref(1,k+1:k+N);
% tref2 = tref(2,k+1:k+N);
% tref3 = tref(3,k+1:k+N);
% tref4 = tref(4,k+1:k+N);
% tref1 = tref(1);
% tref2 = tref(2);
% tref3 = tref(3);
% tref4 = tref(4);
tref3 = tref(:,k+1:k+N);
%% solve optimization problem
        
nu = nc*N;
cvx_begin quiet
    variables epsilonm(nu) epsilonM(nu)
    variable u1 binary
    variable u2 binary
    variable u3 binary
    variable u4 binary


%                 minimize( ([1, U']*coeff{5})*diag(P*ones(N,1))*([1, U']*coeff{5})' +...
%                           ([1, U']*coeff{1}-tref1)*diag(Q(1)*ones(N,1))*([1, U']*coeff{1}-tref1)' +...
%                           ([1, U']*coeff{2}-tref2)*diag(Q(2)*ones(N,1))*([1, U']*coeff{2}-tref2)' +...
%                           ([1, U']*coeff{3}-tref3)*diag(Q(3)*ones(N,1))*([1, U']*coeff{3}-tref3)' +...
%                           ([1, U']*coeff{4}-tref4)*diag(Q(4)*ones(N,1))*([1, U']*coeff{4}-tref4)' )
% 
%                 minimize( ([1, U']*coeff{5})*diag(P*ones(N,1))*([1, U']*coeff{5})' + ([1, U']*coeff{3}-tref3)*diag(Q(3)*ones(N,1))*([1, U']*coeff{3}-tref3)' )
                U = [u1;u2;u3;u4];
                minimize( ([1, 0.35*U']*coeff{5})*diag(P*ones(N,1))*([1, 0.35*U']*coeff{5})'  + Q(1)*norm(epsilonm,2) + Q(2)*norm(epsilonM,2) )
%                 minimize( P*(u1^2 + u2^2 + u3^2 + u4^2) + 1e10*epsilonm + 1e10*epsilonM)
                                                   
    subject to
                if uprev == 1
                    u1 == 1;
%                     u2 == u1;
%                 elseif uprev == 2
%                     u1 == 1;
                end

%                 [1, 0.35*U']*coeff{3} <= tref3 + 1*ones(1,N) + epsilon;
%                 [1, 0.35*U']*coeff{3} >= tref3 - 0.5*ones(1,N) - epsilon;
                [1, 0.35*U']*coeff{3} <= tref3(2,:) + epsilonM(1:nu)';
                [1, 0.35*U']*coeff{3} >= tref3(1,:) - epsilonm(1:nu)';

                epsilonm(1:nu) >= zeros(4,1);
                epsilonM(1:nu) >= zeros(4,1);
%                 % 1 input constraints
%                 U(1:nc:nu) >= ones(N,1)*umin(1);
%                 U(1:nc:nu) <= ones(N,1)*umax(1);
                
    
%                 % 2 inputs constraints
%                 U(1:nc:nu) >= ones(N,1)*umin(1) + U(2:nc:nu)*80;
%                 U(1:nc:nu) <= ones(N,1)*umax(1);
%                 
%                 U(2:nc:nu) >= ones(N,1)*umin(2) + (U(1:nc:nu) - ones(N,1)*umin(1))/100;
%                 U(2:nc:nu) <= ones(N,1)*umax(2);
                
                % 3 inputs constraints TBD if needed
%                 U(1:nc:nu) - U(2:nc:nu) <= ones(N,1)*25;
%                 U(3:nc:nu) >= ones(N,1)*umin(3) + U(1:nc:nu)./300;
%                 U(3:nc:nu) >= (U(1:nc:nu) - U(2:nc:nu))./100 - ones(N,1)*0.05;
%                 U(3:nc:nu) <= ones(N,1)*umax(3);
        
cvx_end
if ( (strcmp(cvx_status,'Solved') == 0) && (strcmp(cvx_status,'Inaccurate/Solved') == 0) )
    disp('Problem is:');
    cvx_status
    return
end
temp = [ [1, U']*coeff{1}(:,1);
         [1, U']*coeff{2}(:,1);
         [1, U']*coeff{2}(:,1);
         [1, U']*coeff{2}(:,1)];
U = U(1:nc)*0.35;
negIndex = find(U < 0);
if sum(U(negIndex)) > -1e-05
    U(negIndex) = 0;
elseif sum(U(negIndex)) < -1e-05
    disp('Negative u')
    return;
end