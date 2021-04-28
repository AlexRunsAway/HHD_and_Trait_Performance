function [flow,ratings,vorticities,measures] = Perform_HHD(f,G,complete_flag,C)

%% Performs the HHD given an edge flow and the operators
% Inputs: 
% 1. f, the edge flow, should be an E by 1 array (E = number of edges)
% 2. G, the gradient, should be E by V (V = number of vertices) and is more
% efficient if sparse
% 3. complete_flag = 0 or 1, if = 0 then the graph is not complete, if = 1
% then it is complete
% 4. C, the curl operator, if no vorticities are needed can be entered as nan,
% for the vorticities on a cycle basis C should be L by E, where L is the number of basis
% loops L = E - (V - 1), if C includes more loops than a cycle basis than
% the code solves a linear programming problem to estimate the most
% parsimonious representation. WARNING: If C does not have rank L then the
% solver will fail.

% Outputs:
% 1. flow, a struct containing the conservative (transitive) and rotational
% (cyclic) components
% 2. ratings, a V by 1 array containing the competitor ratings
% 3. vorticities, a L by 1 array containing the vorticities
% 4. measures, a struct containing the absolute and relative measures of
% transitivity and intransitivity, as well as the empirical estimate for
% the correlation coefficient rho


%% extract dimensions
[E,V] = size(G);

%% compute divergence of f
div = G'*f;

%% check if complete
if complete_flag == 0
    %% form Laplacian
    L = G'*G;

    %% remove a row and column (set rating of competitor 1 to zero)
    L_trunc = L((2:end),(2:end));
    
    %% solve linear system for ratings
    u = L_trunc\div(2:end); % NOTE: if this is too slow, replace with an iterative solver
    ratings = [0;u];
    ratings = ratings - sum(ratings)/V;

    
    %% compute components of flow
    flow.con = G*ratings;
    flow.rot = f - flow.con;
    
    %% set vorticities to nan
    vorticities = nan;

else
    
    %% compute ratings
    ratings = (1/V)*div;

    %% compute components of flow
    flow.con = G*ratings;
    flow.rot = f - flow.con;
    
end

%% check if C has been supplied
if isnan(C) == 0
        
        %% check if cycle basis
        [num_loops,~] = size(C);
        if num_loops == E - (V - 1)
            vorticities = C\flow.rot; % NOTE: if too slow, replace with an iterative solver
            
        elseif num_loops >= E - (V - 1)
            %% compute (most parsimonius) vorticities
            options = optimoptions('linprog','Display','none');
            theta_absolute = linprog(ones([2*num_loops,1]),[],[],[C',-C'],flow.rot,...
                zeros([2*num_loops,1]),[],[],options); % uses linear programming to solve for solution to C' theta = f.rot that minimizes ||theta||_1
            
            vorticities = theta_absolute(1:num_loops) - theta_absolute(num_loops + 1:end);
        else
            vorticities = nan;
        end
    else
        vorticities = nan;
end

%% compute measures
measures.total = norm(f);

measures.trans.abs = norm(flow.con);
measures.intrans.abs = norm(flow.rot);

measures.trans.rel = norm(flow.con)/norm(f);
measures.intrans.rel = norm(flow.rot)/norm(f);

%% compute rho
measures.rho = (E/sum(sum(abs(G*G' - 2*eye(E)))))*((norm(div)/norm(f))^2 - 2);

end