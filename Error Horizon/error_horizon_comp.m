% Add functions path
addpath('Functions\')

% Algorithm parameters
MAX_ITERATIONS  = 20000;
q1              = 4/5;
beta            = 0.05;
q0              = 3/5;

[m,n,A,x,b] = SystemSetup(5000,500,"gaussian",beta,"Matrix\ash958.mat",0,100,1,true,true);

%------------------------------------------------------
% Data collection values
err     = zeros(MAX_ITERATIONS+1, 1);
qerr    = zeros(MAX_ITERATIONS+1, 1);
derr    = zeros(MAX_ITERATIONS+1, 1);

graphstep       = 10;
xpoints         = 1:graphstep:MAX_ITERATIONS+1;

% Initialization
xk  = x+100*randn(n,1);

qxk = xk;
dxk = xk;

err(1)      = norm(x - xk)^2;
qerr(1)     = norm(x - qxk)^2;
derr(1)     = norm(x - dxk)^2;

% Run RK, qRK, dqRK for MAX_ITERATIONS
for iter = 1:MAX_ITERATIONS
    xk     = RKalgstep(A,xk,b);
    qxk    = qRKalgstep(A,qxk,b,q1);
    dxk    = dqRKalgstep(A,dxk,b,q0,q1);
       
    %   Data collection
    err(iter + 1)       = norm(x - xk)^2;
    qerr(iter + 1)      = norm(x - qxk)^2;
    derr(iter + 1)      = norm(x - dxk)^2;
end

% Error graphs
figure
errgraph = semilogy(xpoints, err(xpoints),'DisplayName','RK','LineWidth', 2);
hold on
semilogy(xpoints, qerr(xpoints),'DisplayName','qRK', 'LineWidth', 2);
semilogy(xpoints, derr(xpoints),'DisplayName','dqRK', 'LineWidth', 2);

xlabel("Iteration")
ylabel("Squared Approximation Error")
legend
hold off

% Save tikz plot
filename = "Figs\error_horizon_comparisons";
saveas(errgraph, filename, 'fig');
filename = convertStringsToChars(filename+".tex");
matlab2tikz(filename)

