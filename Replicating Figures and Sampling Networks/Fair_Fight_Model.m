clear
clf

%% set scale parameter
lambda = 1;

%% loop over dimension
Ts = (1:1:25);
for j = 1:length(Ts)
    T = Ts(j);
    
    %% sample for variance and correlation
    n_real = 10^6;
    f_sample.exp = zeros([1,2*n_real]);
    f_prod_sample.exp = zeros([1,n_real]);
    f_sample.g = zeros([1,2*n_real]);
    f_prod_sample.g = zeros([1,n_real]);
    f_sample.u = zeros([1,2*n_real]);
    f_prod_sample.u = zeros([1,n_real]);
    for k = 1:n_real
        %% sample three competitors
        x = -log(rand([T,1]))*lambda;
        y = -log(rand([T,1]))*lambda;
        z = -log(rand([T,1]))*lambda;
        
        %% compute f
        xy = x - y;
        xz = x - z;
        
        [~,J_xy] = min(abs(xy));
        [~,J_xz] = min(abs(xz));
        
        f_xy = x(J_xy) - y(J_xy);
        f_xz = x(J_xz) - z(J_xz);
        
        f_sample.exp(2*(k-1)+1) = f_xy;
        f_sample.exp(2*k) = f_xz;
        f_prod_sample.exp(k) = f_xy*f_xz;
        
        %% sample three competitors Gaussian
        x = randn([T,1])*lambda;
        y = randn([T,1])*lambda;
        z = randn([T,1])*lambda;
        
        %% compute f
        xy = x - y;
        xz = x - z;
        
        [~,J_xy] = min(abs(xy));
        [~,J_xz] = min(abs(xz));
        
        f_xy = x(J_xy) - y(J_xy);
        f_xz = x(J_xz) - z(J_xz);
        
        f_sample.g(2*(k-1)+1) = f_xy;
        f_sample.g(2*k) = f_xz;
        f_prod_sample.g(k) = f_xy*f_xz;
        
        
        %% sample three competitors uniform
        x = rand([T,1])*lambda;
        y = rand([T,1])*lambda;
        z = rand([T,1])*lambda;
        
        %% compute f
        xy = x - y;
        xz = x - z;
        
        [~,J_xy] = min(abs(xy));
        [~,J_xz] = min(abs(xz));
        
        f_xy = x(J_xy) - y(J_xy);
        f_xz = x(J_xz) - z(J_xz);
        
        f_sample.u(2*(k-1)+1) = f_xy;
        f_sample.u(2*k) = f_xz;
        f_prod_sample.u(k) = f_xy*f_xz;
        
    end
    
    %% get variance and correlation
    sigma_sqr.exp(j) = mean(f_sample.exp.^2);
    sigma_sqr.exp_std(j) = std(f_sample.exp.^2)/sqrt(n_real);
    rho.exp(j) = mean(f_prod_sample.exp)/sigma_sqr.exp(j);
    rho.exp_std(j) = (std(f_prod_sample.exp)/sigma_sqr.exp(j))/sqrt(n_real);
    
    sigma_sqr.g(j) = mean(f_sample.g.^2);
    sigma_sqr.g_std(j) = std(f_sample.g.^2)/sqrt(n_real);
    rho.g(j) = mean(f_prod_sample.g)/sigma_sqr.g(j);
    rho.g_std(j) = (std(f_prod_sample.g)/sigma_sqr.g(j))/sqrt(n_real);
    
    sigma_sqr.u(j) = mean(f_sample.u.^2);
    sigma_sqr.u_std(j) = std(f_sample.u.^2)/sqrt(n_real);
    rho.u(j) = mean(f_prod_sample.u)/sigma_sqr.u(j);
    rho.u_std(j) = (std(f_prod_sample.u)/sigma_sqr.u(j))/sqrt(n_real);
    
    
    %% Display
    figure(1)
    clf
    hold on
    errorbar(Ts(1:j),rho.exp(1:j),3*rho.exp_std,'Color',[0.25,0,0.75],'Linewidth',2)
    errorbar(Ts(1:j),rho.g(1:j),3*rho.g_std,'Color',[0.5,0,0.5],'Linewidth',2)
    errorbar(Ts(1:j),rho.u(1:j),3*rho.u_std,'Color',[0.75,0,0.25],'Linewidth',2)
    grid on
    set(gca,'FontSize',16)
    xlabel('$T$ ','FontSize',20,'interpreter','latex')
    ylabel('$\rho$ ','FontSize',20,'interpreter','latex')
    l = legend('Exponential ','Gaussian ','Uniform ');
    set(l,'FontSize',18,'Location','northeast','interpreter','latex');
    title('Fair Fight Performance Model','FontSize',20,'interpreter','latex')
    axis([0,Ts(j),0,0.5]);
    drawnow
end