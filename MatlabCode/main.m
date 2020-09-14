%# Authors: Babak Azad
clc
clear all
close all

%% First Layer

[Grapht1,~,~] = xlsread('01Trustlayer.xlsx','a1');
[Grapht2,~,~] = xlsread('01Trustlayer.xlsx','a2');
[Grapht3,~,~] = xlsread('01Trustlayer.xlsx','a3');
[Grapht4,~,~] = xlsread('01Trustlayer.xlsx','a4');

Graph = Grapht1+Grapht2+Grapht3+Grapht4;
Graph(Graph>1)=1;

[lof1_1, lof2_1, lof3_1, edge1_1, edge2_1, edge3_1] = Layer_Anomaly(Graph);

%% Second layer

%% Third layer
[Graph3,~,~] = xlsread('03ComunicationLayer.xlsx');

[lof1_3, lof2_3, lof3_3, edge1_3, edge2_3, edge3_3] = Layer_Anomaly(Graph3);



%% Fourth layer
[temp1,~,~] = xlsread('04BusinessLayer.xlsx');
Graph4 = zeros(79,79);
for idx=1:size(temp1, 2)
    list = [];
    for idy=1:size(Graph4,2)
     if temp1(idy, idx)==1
        list = [list idy];
     end
    end
    if size(list,2)>1
    for idz=1:size(list,2)
        for idy=1:size(list,2)
            Graph4(list(idz),list(idy))=1;
        end
    end
    end

end
for idx=1:size(Graph4,1)
    Graph4(idx,idx)=0;  
end
Graph4(Graph4>1)= 1;
[lof1_4, lof2_4, lof3_4, edge1_4, edge2_4, edge3_4] = Layer_Anomaly(Graph4);

%% Overall anomaly detection (two layer) with first order information

Total_edg1 = edge1_1+edge1_3;
Multiscore_firstorder = zeros(79,1);

for idx=1:79
    
    Multiscore_firstorder(idx) = ((edge1_1(idx)/Total_edg1(idx))*lof1_1(idx))+((edge1_3(idx)/Total_edg1(idx))*lof1_3(idx));
    
end

Multiscore_firstorder(isnan(Multiscore_firstorder))=0;
[~,suspicious_index]=sort(Multiscore_firstorder,'descend');
suspicious_index = uint8(suspicious_index);
Final_anomaly = zeros(79,2);
for idx = 1:79
    
    Final_anomaly(idx,1) =  suspicious_index(idx);
    Final_anomaly(idx,2) =  Multiscore_firstorder(suspicious_index(idx));
    
end


%% Overall anomaly detection (two layer) with Second order information

Total_edg2 = edge2_1+edge2_3;
Multiscore_secondorder = zeros(79,1);

for idx=1:79
    
    Multiscore_secondorder(idx) = ((edge2_1(idx)/Total_edg2(idx))*lof2_1(idx))+((edge2_3(idx)/Total_edg2(idx))*lof2_3(idx));
    
end

Multiscore_secondorder(isnan(Multiscore_secondorder))=0;
[~,suspicious_index2]=sort(Multiscore_secondorder,'descend');
suspicious_index2 = uint8(suspicious_index2);
Final_anomaly2 = zeros(79,2);
for idx = 1:79
    
    Final_anomaly2(idx,1) =  suspicious_index2(idx);
    Final_anomaly2(idx,2) =  Multiscore_secondorder(suspicious_index2(idx));
    
end


%% Overall anomaly detection (two layer) with both order information

CF = [ 1 1 1]; % Coefficients
Multiscore_bothorder = zeros(79,1);

for idx=1:79
    
    LR1 = ((edge1_1(idx)/Total_edg1(idx))*lof1_1(idx))+((edge1_3(idx)/Total_edg1(idx))*lof1_3(idx));
    LR2 = ((edge2_1(idx)/Total_edg2(idx))*lof2_1(idx))+((edge2_3(idx)/Total_edg2(idx))*lof2_3(idx));
    LR3 = (((edge2_1(idx)+edge1_1(idx))/(Total_edg1(idx)+Total_edg2(idx)))*lof3_1(idx))+ (((edge1_3(idx)+edge2_3(idx)) / (Total_edg1(idx)+Total_edg2(idx))*lof3_3(idx)));
    
    Multiscore_bothorder(idx) = LR1*CF(1)+LR2*CF(2)+LR3*CF(3);
    
end

Multiscore_bothorder(isnan(Multiscore_bothorder))=0;
[~,suspicious_index3]=sort(Multiscore_bothorder,'descend');
suspicious_index3 = uint8(suspicious_index3);
Final_anomaly3 = zeros(79,2);
for idx = 1:79
    
    Final_anomaly3(idx,1) =  suspicious_index3(idx);
    Final_anomaly3(idx,2) =  Multiscore_bothorder(suspicious_index3(idx));
    
end



%% Representation Stage
figure
subplot(1,3,1)
bar(1:79,(Multiscore_firstorder/max(Multiscore_firstorder)),0.6)
xlabel('Person index');
ylabel('Anomaly score');
ylim([0,1.1])
title('Anomaly score First-Order')

subplot(1,3,2)
bar(1:79, (Multiscore_secondorder/max(Multiscore_secondorder)),0.6)
xlabel('Person index');
ylabel('Anomaly score');
ylim([0,1.1])
title('Anomaly score Second-Order')

subplot(1,3,3)
bar(1:79,(Multiscore_bothorder/max(Multiscore_bothorder)),0.6)
xlabel('Person index');
ylabel('Anomaly score');
ylim([0,1.1])
title('Anomaly score Both-Order')

%% Effect of coefficient

CF = [ 3 2 1]; % Coefficients
Multiscore_bothorder1 = zeros(79,1);

for idx=1:79
    
    LR1 = ((edge1_1(idx)/Total_edg1(idx))*lof1_1(idx))+((edge1_3(idx)/Total_edg1(idx))*lof1_3(idx));
    LR2 = ((edge2_1(idx)/Total_edg2(idx))*lof2_1(idx))+((edge2_3(idx)/Total_edg2(idx))*lof2_3(idx));
    LR3 = (((edge2_1(idx)+edge1_1(idx))/(Total_edg1(idx)+Total_edg2(idx)))*lof3_1(idx))+ (((edge1_3(idx)+edge2_3(idx)) / (Total_edg1(idx)+Total_edg2(idx))*lof3_3(idx)));
    
    Multiscore_bothorder1(idx) = LR1*CF(1)+LR2*CF(2)+LR3*CF(3);
    
end

Multiscore_bothorder1(isnan(Multiscore_bothorder1))=0;
%
CF = [ 3 1 2]; % Coefficients
Multiscore_bothorder2 = zeros(79,1);

for idx=1:79
    
    LR1 = ((edge1_1(idx)/Total_edg1(idx))*lof1_1(idx))+((edge1_3(idx)/Total_edg1(idx))*lof1_3(idx));
    LR2 = ((edge2_1(idx)/Total_edg2(idx))*lof2_1(idx))+((edge2_3(idx)/Total_edg2(idx))*lof2_3(idx));
    LR3 = (((edge2_1(idx)+edge1_1(idx))/(Total_edg1(idx)+Total_edg2(idx)))*lof3_1(idx))+ (((edge1_3(idx)+edge2_3(idx)) / (Total_edg1(idx)+Total_edg2(idx))*lof3_3(idx)));
    
    Multiscore_bothorder2(idx) = LR1*CF(1)+LR2*CF(2)+LR3*CF(3);
    
end

Multiscore_bothorder2(isnan(Multiscore_bothorder2))=0;
%

figure
subplot(1,3,1)
bar(1:79,(Multiscore_bothorder/max(Multiscore_bothorder)),0.6)
xlabel('Person index');
ylabel('Anomaly score');
ylim([0,1.1])
title('Anomaly score Both-Order Coef[1 1 1]')

subplot(1,3,2)
bar(1:79,(Multiscore_bothorder1/max(Multiscore_bothorder1)),0.6)
xlabel('Person index');
ylabel('Anomaly score');
ylim([0,1.1])
title('Anomaly score Both-Order Coef[3 2 1]')

subplot(1,3,3)
bar(1:79,(Multiscore_bothorder2/max(Multiscore_bothorder2)),0.6)
xlabel('Person index');
ylabel('Anomaly score');
ylim([0,1.1])
title('Anomaly score Both-Order Coef[3 1 2]')





