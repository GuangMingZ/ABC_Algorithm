clear;
clc;
maxCycle = 10;
Y = 1:maxCycle;
DataPath='F:\ABC����ʵ����\NP��maxCycle��\Sphere\';
%FileName = {'04','06','08','10','12','14','16','18','20','22'};%Rastrigin
%FileName ={'10','12','14','16','18','20','22','26','28','30'};%Rosenbrock
%FileName = {'10','20','30','40','50','60','70','80','90','100'};%Schaffer
FileName = {'04','06','08','10','16','20','24','30','34','40'};%Sphere
average = zeros(10,1);

for i=1:10
    path = strcat(FileName(i),'result.txt');
    url = strcat(DataPath,path);
    url=url{1};%��cell����ת��Ϊstring����
    S1 = load(url);
    a = size(S1,1);
    average(i) = sum(S1)/a;
end;

plot(Y,average,'r-','LineWidth',2);
set(gca,'xticklabel','4|6|8|10|12|14|16|18|20|22');%����X���ע
xlabel('��Ⱥ��ģNP');
ylabel('����ֵ');
%text(2.5,6,'Griewank����');%Ϊ���߼��ϱ�ע
%set(gca,'xtick',[0:1:10],'ytick',[0:0.1:1]);

for n=1:maxCycle
% str=['(' num2str(Y(n)) ',' num2str(average(n)) ')'];
% text(Y(n)+0.2,average(n)+1,str);
str=[ num2str(average(n)) ];
text(n,average(n),str);
end;
