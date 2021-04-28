function [dimensions,A,adjacency_list,complete_flag,edge_to_endpoints,...
    edge_indices] = Get_Topology(Events)
%% Builds a set of data structures storing the graph topology based on observed events
% Inputs:
% 1. Events, a V by V matrix storing the number of events observed between
% each pair of competitors

% Outputs:
% 1. dimensions, a struct containing the number of competitors, edges, and
% basis loops
% 2. A: the adjacency matrix
% 3. adjacency_list: a cell array with V cells, each containing a list of
% the neighbors of the associated competitor (lists the neighborhoods)
% 4. edge_to_endpoints, a E by 2 array containing the endpoints of each
% edge
% 5. edge_indices: a V by V matrix whose i,j entry is the index of the edge
% k with endpoints i,j


%% get number of competitors
[V,~] = size(Events);
dimensions.V = V;

%% build adjacency matrix
A = Events;
A(A > 0) = 1;
for i = 1:V
    A(i,i) = 0;
end

%% build adjacency structure, edge to endpoints and endpoints to edge mappings
edge_count = 0;
edge_indices = nan(V,V);
for j = 1:V
    adjacency_list{j} = find(A(j,:) == 1);
    for k = 1:length(adjacency_list{j})
        i = adjacency_list{j}(k);
        if i > j
            edge_count = edge_count + 1;
            edge_to_endpoints(edge_count,:) = [j,i];
            
            edge_indices(i,j) = edge_count;
            edge_indices(j,i) = edge_count;
        end
    end
end

%% get number of edges
E = edge_count;
dimensions.E = E;

%% get number of loops
L = E - (V - 1);
dimensions.L = L;


%% test if complete
if E == V*(V-1)/2
    complete_flag = 1;
else
    complete_flag = 0;
end




end