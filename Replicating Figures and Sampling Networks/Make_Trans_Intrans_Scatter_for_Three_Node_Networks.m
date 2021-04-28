clear

% ps = (0.01:0.01:0.99);
% qs = ps;

n_real = 10^5;
% ps = rand([3,n_real]);
ps = betarnd(0.2,0.2,[3,n_real]);

% [p_grid,q_grid] = meshgrid(ps,qs);
% p_list = reshape(p_grid,[1,length(ps)^2]);
% q_list = reshape(q_grid,[1,length(ps)^2]);

logit = @(x) log(x./(1-x));
% f = @(p,q) logit([p;p;q]);
f = @(p) logit(p);


G = [1,-1,0;0,1,-1;-1,0,1];
[U,~,~] = svd(G);
Q = U(:,(1:2));
P_trans = Q*Q';
P_cyc = eye([3,3]) - P_trans;

% transitivities = sqrt(sum((P_trans*f(p_list,q_list)).^2));
% intransitivities = sqrt(sum((P_cyc*f(p_list,q_list)).^2));
transitivities = sqrt(sum((P_trans*f(ps)).^2));
intransitivities = sqrt(sum((P_cyc*f(ps)).^2));

transitive_set =[];
intransitive_set = [];
for j = 1:n_real
    %     p = p_list(j);
    %     q = q_list(j);
    
    %     if sign(p-0.5) ~= sign(q-0.5)
    %         transitive_set = [transitive_set,j];
    %     else
    %         intransitive_set = [intransitive_set,j];
    %     end
    
    %% check transitivity
    p = ps(:,j);
    if sign(p(1)-0.5) == sign(p(2)-0.5) && sign(p(2)-0.5) == sign(p(3)-0.5)
        intransitive_set = [intransitive_set,j];
    else
        transitive_set = [transitive_set,j];
    end
    
    
end

%% check bounds
n_intransitive_failures = length(find(intransitivities(intransitive_set) < transitivities(intransitive_set)/sqrt(2)));
n_transitive_failures = length(find(intransitivities(transitive_set) > sqrt(2)*transitivities(transitive_set)));


%% estimate bounds numerically
intransitive_bound_numerical = min(intransitivities(intransitive_set)./transitivities(intransitive_set));
transitive_bound_numerical = max(intransitivities(transitive_set)./transitivities(transitive_set));

%% display

figure(1)
clf
hold on
scatter(transitivities(transitive_set),intransitivities(transitive_set),1,'b','fill')
scatter(transitivities(intransitive_set),intransitivities(intransitive_set),1,'r')
set(gca,'FontSize',16)
shift = 0.15;
plot((0:12),(1/sqrt(2))*(0:12),'m','Linewidth',1.5)
plot((0:12),sqrt(2)*(0:12),'m','Linewidth',1.5)
scatter(0,7.959,40,'k','o','Linewidth',2)
text(0+shift,7.959+shift,'0.99,0.99,0.99','FontSize',14,'BackgroundColor','white')
scatter(3.7192,5.3291,40,'k','o','Linewidth',2)
text(3.7192+shift,5.3291+2*shift,'0.99,0.99,0.51','FontSize',14,'BackgroundColor','white')
scatter(3.7846,5.2829,40,'k','o','Linewidth',2)
text(3.7846+shift,5.2829-2*shift,'0.99,0.99,0.49','FontSize',14,'BackgroundColor','white')
scatter(11.2557,0,40,'k','o','Linewidth',2)
text(11.2557+shift,0+shift,'0.99,0.99,0.0001','FontSize',14,'BackgroundColor','white')
scatter(0,0.0693,40,'k','o','Linewidth',2)
text(0+shift,0.0693+4*shift,'0.51,0.51,0.51','FontSize',14,'BackgroundColor','white')
scatter(0.098,0,40,'k','o','Linewidth',2)
text(0.098+shift,0+shift,'0.51,0.51,0.48','FontSize',14,'BackgroundColor','white')
grid on
xlabel('Size of Transitive Component: $||f_t||_2$','FontSize',20,'interpreter','latex')
ylabel('Size of Cyclic Component: $||f_c||_2$','FontSize',20,'interpreter','latex')
axis square
xlim([0,12])
ylim([0,12])
drawnow