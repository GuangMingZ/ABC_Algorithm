%%%%%%
% 本函数用于在同一坐标系中显示不同NP数目下函数求解的收敛性
% 通过一系列NP值的变化，得到相应的收敛过程，得出在一定范围内增加NP的值
% 确实可以提高函数收敛速度和最优解的性能，而当达到一定数目时，继续增大NP
% 反而会出现收敛性速度下降、最优解性能降低的现象，由此可以得出NP存在一个
% 最优值，盲目地增加NP不会持续提高算法收敛速度和最优解性能！
clear;
clc;
maxCycle = 1000;
S1 = load('F:\ABC算法参数研究\NP组\Sphere\04result.txt');
S2 = load('F:\ABC算法参数研究\NP组\Sphere\30result.txt');
%S3 = load('F:\ABC算法参数研究\NP组\Sphere\30result.txt');
%S4 = load('F:\ABC算法参数研究\NP组\Sphere\04result.txt');
Y = 1:maxCycle;
[a,c] = size(S1);
b = a./1000-1;
average1 = zeros(maxCycle,1);
average2 = zeros(maxCycle,1);
average3 = zeros(maxCycle,1);
%%%%%%%%%%%% NP = 4 %%%%%%%%%%%%
for i=0:b
    for j=1:maxCycle
        temp = j+b*1000;
        average1(j) =average1(j) + S1(temp);
    end;
end;
for a=1:maxCycle
 average1(j) = average1(j)./(b+1);
end;
plot(Y,average1,'m-.');
hold on;
%%%%%%%%%%%%%%% NP = 6 %%%%%%%%%%%%%%%%%%%
% S4(1:end,1:end) = 0;
% S4 = load('E:\06result.txt');
% average1(1:maxCycle) = 0;
% for i=0:b
%     for j=1:maxCycle
%         temp = j+b*1000;
%         average1(j) =average1(j) + S4(temp);
%     end;
% end;
% for a=1:maxCycle
%  average1(j) = average1(j)./(b+1);
% end;
% plot(Y,average1,'k-.');
% hold on;
%%%%%%%%%%%%%% NP = 10 %%%%%%%%%%%%%%%%%%%

average1(1:maxCycle) = 0;
for i=0:b
    for j=1:maxCycle
        temp = j+b*1000;
        average1(j) =average1(j) + S2(temp);
    end;
end;
for a=1:maxCycle
 average1(j) = average1(j)./(b+1);
end;
plot(Y,average1,'-b');
hold on;
% clear;
% clc;
%%%%%%%%%%%%%% NP = 20 %%%%%%%%%%%%%%%%%%%
for i=0:b
    for j=1:maxCycle
        temp = j+b*1000;
        average2(j) =average2(j) + S3(temp);
    end;
end;
for a=1:maxCycle
 average2(j) = average2(j)./(b+1);
end;
plot(Y,average2,':r');
hold on;
%%%%%%%%%%%% NP = 30 %%%%%%%%%%%%%%
for i=0:b
    for j=1:maxCycle
        temp = j+b*1000;
        average3(j) =average3(j) + S4(temp);
    end;
end;
for a=1:maxCycle
 average3(j) = average3(j)./(b+1);
end;
plot(Y,average3,'--k');
hold on;
legend('NP=10','NP=20','NP=30','NP=40');