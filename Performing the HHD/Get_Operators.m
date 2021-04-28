function [G,C,cycles] = Get_Operators(E,V,edge_to_endpoints,edge_indices,complete_flag)

%% generate gradient (by convention edges point from lower index to higher indexed competitor)
G = sparse((1:E),edge_to_endpoints(:,1),1,E,V) - ...
    sparse((1:E),edge_to_endpoints(:,2),1,E,V);

%% if complete generate the overdetermined curl (all triangles)
if complete_flag == 1
    k = 0;
    n_triangles = V*(V-1)*(V-2)/6;
    C = sparse([],[],[],n_triangles,E);
    for node_1 = 1:V
        for node_2 = node_1+1:V
            for node_3 = node_2+1:V
                k = k+1;
                
                edge_1 = edge_indices(node_1,node_2);
                edge_2 = edge_indices(node_2,node_3);
                edge_3 = edge_indices(node_3,node_1);
                
                cycles(k,:) = [edge_1,edge_2,edge_3];
                
                C(k,cycles(k,:)) = [1,1,-1];
                
            end
        end
    end
else
    C = nan;
    cycles = nan;
end

end