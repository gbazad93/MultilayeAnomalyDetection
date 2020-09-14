clc
clear all
close all

%% Reading the Input File
[Graph1,txt,~] = xlsread('01Trustlayer.xlsx','a1');
[Graph2,txt,~] = xlsread('01Trustlayer.xlsx','a2');
[Graph3,txt,~] = xlsread('01Trustlayer.xlsx','a3');
[Graph4,txt,~] = xlsread('01Trustlayer.xlsx','a4');

Graph = Graph1+Graph2+Graph3+Graph4;
Graph(Graph>1)=1;




NOF = size(Graph, 1); % keep the number of nodes

connections = zeros(size(Graph,2),3); % shows total connections, common, unique
for idx = 1: size(Graph, 2)
    
   connections(idx,1) = sum(Graph(:,idx));
    

    
end
    
    



