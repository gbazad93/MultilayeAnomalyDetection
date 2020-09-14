clc
clear all
close all

%% Reading the Input File
[Graph,txt,~] = xlsread('NoordinSubset.xlsx');

%% Calculate the Feature vector for each node

NOF = size(Graph, 1); % keep the number of nodes
Matrix_Features = zeros(NOF,2); % first cloumn for number of node and second cloumn for number of edges



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

%% Calculate the Anomoaly score for each node

% For third layer
[suspicious_index, lof] = LOF(Matrix_Features, 1);







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









