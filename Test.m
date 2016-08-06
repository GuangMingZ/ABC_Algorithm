clear;
clc;
maxCycle = 1000;
avecycle = 1000;

average1 = zeros(avecycle,1);
average2 = zeros(avecycle,1);
average3 = zeros(avecycle,1);
average4 = zeros(avecycle,1);
average5 = zeros(avecycle,1);
Y = 1:avecycle;

S1 = load('F:\ABC����ʵ����\NP��\Sphere\04result.txt');
[a,c] = size(S1);

for i=1:a
    for j=1:avecycle
        average1(j) =average1(j) + S1(i,j);
    end;
end;
for i=1:avecycle
 average1(i) = average1(i)./a;
end;
plot(Y,average1,'m-');
hold on;

S1 = load('F:\ABC����ʵ����\NP��\Sphere\06result.txt');
for i=1:a
    for j=1:avecycle
        average2(j) =average2(j) + S1(i,j);
    end;
end;
for i=1:avecycle
 average2(i) = average2(i)./a;
end;
plot(Y,average2,'k--');
hold on;

S1 = load('F:\ABC����ʵ����\NP��\Sphere\10result.txt');
for i=1:a
    for j=1:avecycle
        average3(j) =average3(j) + S1(i,j);
    end;
end;
for i=1:avecycle
 average3(i) = average3(i)./a;
end;
plot(Y,average3,'b-');
hold on;

S1 = load('F:\ABC����ʵ����\NP��\Sphere\30result.txt');
for i=1:a
    for j=1:avecycle
        average4(j) =average4(j) + S1(i,j);
    end;
end;
for i=1:avecycle
 average4(i) = average4(i)./a;
end;
plot(Y,average4,'r--');
hold on;

S1 = load('F:\ABC����ʵ����\NP��\Sphere\50result.txt');
for i=1:a
    for j=1:avecycle
        average5(j) =average5(j) + S1(i,j);
    end;
end;
for i=1:avecycle
 average5(i) = average5(i)./a;
end;
plot(Y,average5,'g-');
hold on;

legend('NP=10','NP=20','NP=50','NP=80','NP=100');