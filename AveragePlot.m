%%%%%%
% ������������ͬһ����ϵ����ʾ��ͬNP��Ŀ�º�������������
% ͨ��һϵ��NPֵ�ı仯���õ���Ӧ���������̣��ó���һ����Χ������NP��ֵ
% ȷʵ������ߺ��������ٶȺ����Ž�����ܣ������ﵽһ����Ŀʱ����������NP
% ����������������ٶ��½������Ž����ܽ��͵������ɴ˿��Եó�NP����һ��
% ����ֵ��äĿ������NP�����������㷨�����ٶȺ����Ž����ܣ�
clear;
clc;
maxCycle = 1000;
S1 = load('F:\ABC�㷨�����о�\NP��\Sphere\04result.txt');
S2 = load('F:\ABC�㷨�����о�\NP��\Sphere\30result.txt');
%S3 = load('F:\ABC�㷨�����о�\NP��\Sphere\30result.txt');
%S4 = load('F:\ABC�㷨�����о�\NP��\Sphere\04result.txt');
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