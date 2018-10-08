m = 100;
n = 10;

A = randn(m,n);
xtrue = ones(n,1);
b = A*xtrue;

S = diag([1:m]);
lambda = 1e1;

xp = (A'*S*A + lambda*eye(n))\(A'*S*b);

x0 = xt + randn(n,1);
niter = 500;
rho = 1e1;


%% steepest descent
% f(x) = 0.5*|A*x - b|_S^2 + 0.5*lambda*|x|^2
M = 1/norm(A'*S*A);
xk = x0;
e1 = zeros(niter,1);
fprintf(1,'k, |g|\n')
for k = 1:niter
    % gradient step
    gk = A'*S*(A*xk - b) + lambda*xk;
    xt = xk - M*gk;
    
    % print some info
    fprintf(1,'%d,%f\n',k,norm(gk));
    
    % update
    xk = xt;
    
    % compute error
    e1(k) = norm(xk - xp);
end
x1 = xk;

%% ADMM
% L(x,z,u) = 0.5*|z - b|^2_S + 0.5*rho*|A*x - z|^2 + u'*(A*x - z) +
% 0.5*lambda*|x|

xk = zeros(n,1);
zk = zeros(m,1);
uk = zeros(m,1);
e2 = zeros(niter,1);

fprintf(1,'k, |r|,|s|\n');
for k = 1:niter
    % min_z L(xk,z,uk): S*(z - b) - rho*(A*x - z) - u = 0
    zt = (eye(m) + (1/rho)*S)\((1/rho)*S*b + A*xk + uk/rho); 
    
    % min_x L(x,zk,uk): rho*A'*(A*x - z) + A'*u + lambda*x = 0
    xt = (A'*A + (lambda/rho)*eye(n))\(A'*(zk - uk/rho));
    
    % uk := uk + rho*L_u(xk,xk,uk)
    ut = uk + rho*(A*xk - zk); 
    
    % update rho
    rho = 1*rho;
    
    % print some optimality info
    fprintf(1,'%d,%f,%f\n',k,norm(A*xk - zk),rho*norm(A'*(zt - zk)));
    
    % update
    zk = zt; xk = xt; uk = ut;
    
    % compute error
    e2(k) = norm(xk - xp);
end
x2 = xk;

%% plot
semilogy(1:niter,e1,1:niter,e2);
