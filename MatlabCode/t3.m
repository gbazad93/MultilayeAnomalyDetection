clc
clear all
close all

%% Reading the Input File
[Graph,txt,~] = xlsread('NoordinSubset.xlsx');

%% Calculate the Feature vector for each node

NOF = size(Graph, 1); % keep the number of nodes
Matrix_Features = zeros(NOF,2); % first cloumn for number of node and second cloumn for number of edges

Matrix_Features_secondorder = zeros(NOF,2);

for idx=1:NOF
    
    Matrix_Features(idx, 1) = sum(Graph(idx,:)); % Number of nodes
    candidates = [];
    for idy=1:NOF
        
        if Graph(idx,idy)==1
            
           candidates = [candidates, idy]; 
            
        end
        
    end
    Matrix_Features(idx, 2) = Get_edges_count(Graph, candidates); % Number of edges
    
   
end

% Calculate for second order
for idx=1:NOF
    
    neighbour_nodes = [];
    for idy=1:NOF
        if Graph(idx,idy)==1
            neighbour_nodes = [neighbour_nodes, idy];
        end
    end
    temp = Matrix_Features(neighbour_nodes,:);
    Matrix_Features_secondorder(idx,1)= Matrix_Features(idx,1)/(sum(temp(:, 1)));
    Matrix_Features_secondorder(idx,2)= Matrix_Features(idx,2)/(sum(temp(:, 2)));
    
    
end
       
Matrix_Features_secondorder(isnan(Matrix_Features_secondorder))=0;
    
Matrix_bothorder = zeros(size(Matrix_Features,1),4);
Matrix_bothorder(:, 1:2) = Matrix_Features;
Matrix_bothorder(:, 3:4) = Matrix_Features_secondorder;
%% Calculate the Anomoaly score for each node

% For third layer
[suspicious_index1, lof1] = LOF(Matrix_Features, 1);
[suspicious_index2, lof2] = LOF(Matrix_Features_secondorder, 1);
[suspicious_index3, lof3] = LOF(Matrix_bothorder, 1);

a = 1;
%%%%%%%%%%%%########## Functions ##########%%%%%%%%%%%%%
function [Edge_count] = Get_edges_count(Graph, candidates)

Edge_count = 0;
nodes  = size(candidates,2);

if nodes == 0
    Edge_count = 0; 
elseif nodes ==1
    Edge_count = 1;
else    
    for idx=1:nodes
        for idy=idx+1:nodes
            if Graph(candidates(idx), candidates(idy))==1
               Edge_count = Edge_count + 1;
            end
        end
    end
end

end









