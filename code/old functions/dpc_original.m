function [u,fval,iter,time] = mydpc(model, leafmodels, k, x, d, Q, R, orderAR, N, umin, umax, uprev, udotmax, solver, xref)

%% calculate leaf coeffs

% prepare features_d to find an appropriate leaf
d_cur = d(:,k:k+N-1); % features corr to disturbances
x_cur = fliplr(x(4,k-orderAR+1:k)); % features corr to indoor temp.
features_d = [d_cur(:)', x_cur];

coeff = zeros(N+1,N);
for idm = 1:N
    coeff(1:idm+1,idm) = find_linearmodel_in_leaves(model{idm}, leafmodels{idm}, features_d);
end

% objective function
QQ = diag(Q(4,4)*ones(N,1));
RR=kron(eye(N),R);
xref = xref(4)*ones(N,1);

% input rate change constraint
UB = uprev+udotmax;
LB = uprev-udotmax;
umin = max(umin,LB);
umax = min(UB,umax);

%% solve optimization problem

switch solver
    
    case 'cvx' % too slow
        
        tic;
        cvx_begin quiet
            variables U(N)
            minimize( U'*RR*U + ([1, U']*coeff-xref')*QQ*([1, U']*coeff-xref')' )
            subject to
                U>=umin; %#ok<*VUNUS>
                U<=umax;
        cvx_end
        time = toc;
        fval = 0;
        iter = 0;
        
    case 'quadprog'
        
        C1 = coeff(1,:);
        C2 = coeff(2:end,:);
        H = 2*( RR +  C2*QQ*C2' );
        f = 2*C2*QQ*C1'-2*C2*QQ*xref;
        c = C1*QQ*C1'-C1*QQ*xref-xref'*QQ*C1'+xref'*QQ*xref;
        lb = repmat(umin,N,1);
        ub = repmat(umax,N,1);
        options=optimset;
        options.Algorithm='interior-point-convex'; % to get rid of some warnings.
        options.Display='off';
        tic;
        [U, fval] = quadprog(H,f,[],[],[],[],lb,ub,[],options);
        time=toc;
        iter=nan;
        
end
    
u = U(1);