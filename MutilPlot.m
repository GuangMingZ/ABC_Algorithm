clear;
clc;
maxCycle = 10000;
avecycle = 5000;

average1 = zeros(avecycle,1);
average2 = zeros(avecycle,1);
average3 = zeros(avecycle,1);
average4 = zeros(avecycle,1);
average5 = zeros(avecycle,1);

S1 = load('E:\Schaffer\10result.txt');
[a,c] = size(S1);
Y = 1:avecycle;

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

S2 = load('E:\Schaffer\20result.txt');
[a,c] = size(S2);
%Y = 1:5000;
for i=1:a
    for j=1:avecycle
        %temp = j+b*1000;
        if(j<=5000)
        average2(j) =average2(j) + S2(i,j);
        else
            average2(j) = average2(5000);
        end;
    end;
end;
for i=1:avecycle
 average2(i) = average2(i)./a;
end;
plot(Y,average2,'k--');
hold on;

S2 = load('E:\Schaffer\50result.txt');
%Y = 1:1000;
for i=1:a
    for j=1:2000
        average3(j) =average3(j) + S2(i,j);
    end;
end;
for i=1:2000
 average3(i) = average3(i)./a;
end;
for i=2001:avecycle
    average3(i) = average3(2000);
end;
plot(Y,average3,'b-');
hold on;

S2 = load('E:\Schaffer\80result.txt');
%Y = 1:1250;
for i=1:a
    for j=1:1250
        %temp = j+b*1000;
        average4(j) =average4(j) + S2(i,j);
    end;
end;
for i=1:1250
 average4(i) = average4(i)./a;
end;
for i=1251:avecycle
    average4(i) = average4(1250);
end;
plot(Y,average4,'r--');
hold on;

S2 = load('E:\Schaffer\100result.txt');
%Y = 1:1000;

for i=1:a
    for j=1:1000
        %temp = j+b*1000;
        average5(j) =average5(j) + S2(i,j);
    end;
end;
for i=1:1000
 average5(i) = average5(i)./2;
end;
for i=1001:avecycle
    average5(i) = average5(1000);
end;
plot(Y,average5,'g-');
hold on;
legend('NP=10','NP=20','NP=50','NP=80','NP=100');
%     for j=1:avecycle
%         average1(j) = S(j);
%     end;
%      for j=1:avecycle
%         average2(j) = S(j+1*1000);
%     end;
%      for j=1:avecycle
%         average3(j) = S(j+2*1000);
%     end;
%      for j=1:avecycle
%         average4(j) = S(j+3*1000);
%     end;
%         plot(Y,average1,'b',Y,average2,'g',Y,average3,'r',Y,average4,'c');
%end;
