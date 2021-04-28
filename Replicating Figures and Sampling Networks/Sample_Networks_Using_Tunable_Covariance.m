%% Make random networks figure for SIREV paper
clear
figure(1)
clf
set(gcf,'color','w')

%% pick network sizes
Vs = 20*[1,5,10,15];

%% pick (V-1)/E
densities = [0.2,0.4,0.6]; 

%% pick variance in performance
sigma = 1;

%% pick correlations
rhos = 0.5*linspace(0,1,5);

%% pick number of sampled networks
n_networks = 20;

%% pick number of samples per network
n_flows = 50;

%% preallocate
transitivity = nan([length(Vs),length(densities),length(rhos),n_networks,n_flows]);
intransitivity = nan([length(Vs),length(densities),length(rhos),n_networks,n_flows]);
expected_trans = nan([length(densities),length(rhos)]);
expected_intrans = nan([length(densities),length(rhos)]);

%% set up axes for subplots
[handles, ~] = tight_subplot(length(densities),length(rhos), [.03 0],[.1 .05],[.1 .1]);

%% loop over network sizes
for V_index = 1:length(Vs)
    V = Vs(V_index);
    E_complete = V*(V-1)/2;
    
    %% generate complete endpoints list
    endpoints_complete = nan([E_complete,2]);
    f = nan([E_complete,1]);
    k = 0;
    for i = 1:V
        for j = i+1:V
            % update edge count
            k = k+1;
            
            % record endpoints
            endpoints_complete(k,:) = [i,j]; % records endpoints of each edge
        end
    end
    
    %% loop over sparsity
    for p_index = 1:length(densities)
        % E = p*V*(V-1)/2, so (V - 1)/E = 2/(p*V), if want to equal density
        % set p*V = 2/density, so p = 2/(V*density);
        density = densities(p_index);
        p = 2/(V*density);
        
        %% loop over correlations
        for rho_index = 1:length(rhos)
            rho = rhos(rho_index);
            
            %% compute predicted sizes of components
            E_expected = p*E_complete;
            L_expected = E_expected - (V-1);
            expected_trans(p_index,rho_index) = sigma^2*((V - 1)/(E_expected) + 2*rho*L_expected/E_expected);
            expected_intrans(p_index,rho_index) = sigma^2*(1 - 2*rho)*L_expected/E_expected;
            
            if E_expected >= (1.25)*(V-1)
                
                %% loop over networks
                for realization = 1:n_networks
                    
                    %% loop over realizations
                    stop = 0;
                    max_it = 100;
                    it = 0;
                    while stop == 0
                        
                        %% sample which endpoints to keep
                        z = rand([E_complete,1]);
                        endpoints = endpoints_complete(z < p,:);
                        
                        %% build adjacency
                        A = sparse(endpoints(:,1),endpoints(:,2),1,V,V) +...
                            sparse(endpoints(:,2),endpoints(:,1),1,V,V);
                        
                        
                        %% check connected
                        [subgraphs] = conncomp(graph(A));
                        n_components = max(subgraphs);
                        
                        if n_components == 1
                            stop = 1;
                        end
                        
                        %%don't repeat too many times
                        it = it + 1;
                        if it >= max_it
                            stop = 1;
                        end
                        
                    end
                    
                    %% get dimensions
                    [E,~] = size(endpoints);
                    
                    %% Build gradient
                    G = sparse((1:E),endpoints(:,2),1,E,V) - sparse((1:E),endpoints(:,1),1,E,V);
                    
                    %% Build Laplacian
                    degrees = sum(A);
                    L = sparse((1:V),(1:V),degrees,V,V) - A;
                    
                    %% build covariance
                    I = sparse((1:E),(1:E),1,E,E);
                    C = sigma^2*(I + rho*(G*G' - 2*I));
                    
                    %% sample f
                    F = mvnrnd(zeros([E,1]),C,n_flows)';
                    
                    %% apply HHD
                    D = G'*F; % computes all the divergences at once
                    %             R0 = sparse((1:V),(1:V),1/degrees,V,V)*D; % guess at ratings
                    L_small = L((1:V-1),(1:V-1));
                    D_small = D((1:V-1),:);
                    R = L_small\D_small;
                    R = [R;zeros(1,n_flows)];
                    
                    F_trans = G*R;
                    F_cyclic = F - F_trans;
                    
                    %% compute sizes of components
                    for j = 1:n_flows
                        transitivity(V_index,p_index,rho_index,realization,j) = norm(F_trans(:,j))/sqrt(E);
                        intransitivity(V_index,p_index,rho_index,realization,j) = norm(F_cyclic(:,j))/sqrt(E);
                    end
                    
                    %% display
                    figure(1)
                    example_index = length(rhos)*(p_index - 1) + rho_index;
%                     subplot(length(densities),length(rhos),example_index)
                    axes(handles(example_index))
                    hold on
                    if V_index == 1
                        if realization == 1
                            plot([1,0],[0,1],'k--','Linewidth',0.5)
%                             plot([0,expected_trans(p_index,rho_index)],[0,expected_intrans(p_index,rho_index)],'k-','Linewidth',2)
                            plot([0,10*expected_trans(p_index,rho_index)],[0,10*expected_intrans(p_index,rho_index)],'k-','Linewidth',0.5)
                        end
                    end
                    c = (V  - min(Vs))/(max(Vs) - min(Vs));
                    color = [c,0,1-c];
                    if rho == 0.5
                        scatter(transitivity(V_index,p_index,rho_index,realization,:).^2,intransitivity(V_index,p_index,rho_index,realization,:).^2,...
                            1,color,'fill')
                    end
                    if sum(isnan(squeeze(transitivity(V_index,p_index,rho_index,realization,:).^2))) == 0
                        boundary_samples = convhull(squeeze(transitivity(V_index,p_index,rho_index,realization,:).^2),squeeze(intransitivity(V_index,p_index,rho_index,realization,:).^2));
                        fill(squeeze(transitivity(V_index,p_index,rho_index,realization,boundary_samples).^2),squeeze(intransitivity(V_index,p_index,rho_index,realization,boundary_samples).^2),color,...
                            'facealpha',0.025,'LineStyle','none')
                    end
                    grid on
                    titlestring = strcat('$\frac{V-1}{E} = $',{' '},num2str(density),', $\rho = $',{' '},num2str(rho));
                    %                 title(titlestring,'FontSize',14,'interpreter','latex')
                    xlim([0,1.5])
                    ylim([0,1.5])
                    axis square
                    if rho_index ~=1
                        set(gca,'YTickLabel',{''})
                    else
                        set(gca,'YTickLabel',{'0','0.5','1.0','1.5'})
                        ylabel('$||f_c||^2$','FontSize',12,'interpreter','latex')
                    end
                    if p_index ~= 3
                        set(gca,'XTickLabel',{''})
                    else
                        set(gca,'XTickLabel',{'0','0.5','1.0','1.5'})
                        xlabel('$||f_t||^2$','FontSize',12,'interpreter','latex')
                    end
                    set(gca,'FontSize',12)
                    drawnow
                    
                    
                end
            end
            
   
        end
        
    end
end
