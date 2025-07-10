% Add functions path
addpath('Functions\')

% Algorithm parameters
MAX_ITERATIONS = 15000;
q1 = 4/5;
beta = 0.05;
q0 = 3/5;

% Matrix parameters
mattype = "gaussian";

% Corruption scaling parameters (and count)
start_corruption_scale = 0;
increment_size = 0.5;
end_corruption_scale = 15;
corruption_scales = start_corruption_scale:increment_size:end_corruption_scale;
num_corruption_scales = length(corruption_scales);

% Number of points at the end of a sample run to use to compute error horizon
num_points = 100;

% Number of samples per corruption scale
sample_size = 10;

% Points to plot ratio of largest and (1-q)m+1 largest entries (in
% magnitude) v.s. error horizon
RK_points = zeros(2,num_corruption_scales*sample_size);
dqRK_points = zeros(2,num_corruption_scales*sample_size);

corruption_count = 1;
bar = waitbar(1/(num_corruption_scales*sample_size),"Loading sample 1 of corruption scale " + string(start_corruption_scale));
for corruption_scale = corruption_scales
    for run = 1:sample_size
        waitbar(((corruption_count-1)*sample_size+run)/(num_corruption_scales*sample_size),bar,"Loading sample " + string(run) + " of corruption scale " + string(corruption_scale))
        
        % Initialize the System
        [m,n,A,x,b] = SystemSetup(2500,500,mattype,beta,"Matrix\ash958.mat",0,corruption_scale,1,true,true);
        
        % Data collection values and parameters
        possible_error_horizon = zeros(num_points,1);
        possible_d_error_horizon = zeros(num_points,1);
        
        % Initial point
        xk = ones(n, 1)*100;
        dxk = xk;
        
        % Run RK and dqRK for MAX_ITERATIONS
        for iter = 1:MAX_ITERATIONS
            xk = RKalgstep(A,xk,b);
            dxk = dqRKalgstep(A,dxk,b,q0,q1);
            % Error Horizon Data collection
            if iter >= MAX_ITERATIONS - num_points
                possible_error_horizon(iter - MAX_ITERATIONS + num_points + 1) = norm(x - xk)^2;
                possible_d_error_horizon(iter - MAX_ITERATIONS + num_points + 1) = norm(x - dxk)^2;
            end
        end
        % Record entries for RK_points and dqRK_points
        largest_corruptions = maxk(abs(b-A*x),(1-q1)*m+1);
        ratio = largest_corruptions(1)/largest_corruptions(end);
        RK_points(1,(corruption_count-1)*sample_size+run) = ratio;
        RK_points(2,(corruption_count-1)*sample_size+run) = max(possible_error_horizon);
        dqRK_points(1,(corruption_count-1)*sample_size+run) = ratio;
        dqRK_points(2,(corruption_count-1)*sample_size+run) = max(possible_d_error_horizon);
    end
    corruption_count = corruption_count + 1;
end

% Plot RK_points and dqRK_points
f = figure;
set(gca,'Yscale','log')
hold on
error_horizon_graph = scatter(RK_points(1,:),RK_points(2,:),'filled','DisplayName','RK');
scatter(dqRK_points(1,:),dqRK_points(2,:),'filled','DisplayName','dqRK')
legend
hold off

% Save tikz plot
filename = "Figs\error_horizon_vs_sparse_scale";
saveas(error_horizon_graph, filename, 'fig');
filename = convertStringsToChars(filename+".tex");
matlab2tikz(filename)