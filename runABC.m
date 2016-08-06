%/* ABC algorithm coded using MATLAB language */

%/* Artificial Bee Colony (ABC) is one of the most recently defined algorithms by Dervis Karaboga in 2005, motivated by the intelligent behavior of honey bees. */

%/* Referance Papers*/

%/*D. Karaboga, AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL OPTIMIZATION,TECHNICAL REPORT-TR06, Erciyes University, Engineering Faculty, Computer Engineering Department 2005.*/

%/*D. Karaboga, B. Basturk, A powerful and Efficient Algorithm for Numerical Function Optimization: Artificial Bee Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39, Issue:3,pp:459-171, November 2007,ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x */

%/*D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony (ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, Pages 687-697. */

%/*D. Karaboga, B. Akay, A Comparative Study of Artificial Bee Colony Algorithm,  Applied Mathematics and Computation, 214, 108-132, 2009. */

%/*Copyright ?2009 Erciyes University, Intelligent Systems Research Group, The Dept. of Computer Engineering*/

%/*Contact:
%Dervis Karaboga (karaboga@erciyes.edu.tr )
%Bahriye Basturk Akay (bahriye@erciyes.edu.tr)
%*/


clear all
close all
clc

%/* Control Parameters of ABC algorithm*/
global NP; %/* The number of colony size (employed bees+onlooker bees)*/
global FoodNumber; %/*The number of food sources equals the half of the colony size*/
global limit; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
global maxCycle;
global D;
global ub;
global lb;

NP = 100;
FoodNumber = NP/2;
limit = 200;
maxCycle=1000; %/*The number of cycles for foraging {a stopping criteria}*/


%/* Problem specific variables*/
objfun='Schaffer'; %cost function to be optimized
D=15; %/*The number o f parameters of the problem to be optimized*//*当前问题的参数数目*/
ub=ones(1,D)*100; %ones产生一个1×D元素全是1的矩阵/*lower bounds of the parameters. */参数的上限
lb=ones(1,D)*(-100);%/*upper bound of the parameters.*/参数的下限

runtime=1;%/*Algorithm can be run many times in order to see its robustness*/



%Foods [FoodNumber][D]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
%ObjVal[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
%Fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
%trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
%prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
%solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
%ObjValSol; /*Objective function value of new solution*/
%FitnessSol; /*Fitness value of new solution*/
%neighbour, param2change; /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
%GlobalMin; /*Optimum solution obtained by ABC algorithm*/
%GlobalParams[D]; /*Parameters of the optimum solution*/
%GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/

GlobalMins=zeros(1,runtime);%一个 1×runtime 全零矩阵.存储当前全局函数最小值
values=cell(1,maxCycle); %生成一个1×maxCycle空矩阵

for r=1:runtime
  
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */

Range = repmat((ub-lb),[FoodNumber 1]);%FoodNumber是食物源的数量，等于蜂群种群数目NP的一半
Lower = repmat(lb, [FoodNumber 1]);

Foods = rand(FoodNumber,D) .* Range + Lower;%rand（M,N）返回一个（0,1）之间数值的M×N维矩阵
%返回一个FoodNumber×D 矩阵，生成初始解

%[y1,..,yn] = Feval_r(F,x1,...,xn) F是需要使用函数的函数名，或者句柄;xi是函数的参数，yi是函数的返回值，feavl可以指向任何自定义的函数 
%[a1,b1] = size(Foods)
ObjVal=feval('Sphere',Foods);%feval函数是matlab自定义的函数，用于函数间接调用，Foods作为objfun函数的参数
Fitness=calculateFitness(ObjVal);

%reset trial counters
trial=zeros(1,FoodNumber);%一个1×FoodNumber全零向量

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));%找到ObjVal向量中的最小值的序列号，最小值有个能有多个，可能会返回一个序列。
%min(ObjVal)是找到ObjVal向量中的最小值。
BestInd=BestInd(end);%BestInd(end)获取BestInd序列中的最后一个元素值
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);%获取矩阵的第BestInd行，也就是当前最小值
iter=1;

while ((iter <= maxCycle)),

%%%%%%%%% 雇佣蜂阶段 %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;%fix（）截尾取整，向零靠拢取整.
        % rand函数产生在(0, 1)之间均匀分布的随机数组成的数组。
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour == i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);%获取第i个蜜源解
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       % sol(Param2Change)，Param2Change是要改变的参数的位置，下面的语句是蜜蜂的邻域搜索过程，及对某一个参数在一个可控的邻域内对其进行改变。
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       % 当邻域搜索生成的新解出现越界时，对其进行规范化，直接取边界值作为变异后的新值。
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValEva=feval(objfun,sol);%计算新解的收益度值
        FitnessEva=calculateFitness(ObjValEva);
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessEva>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessEva;
            ObjVal(i)=ObjValEva;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
         
         
    end;

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities 计算轮盘赌概率 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/* A food source is chosen with the probability which is proportioal to its quality*/
%/*Different schemes can be used to calculate the probability values*/
%/*For example prob(i)=fitness(i)/sum(fitness)*/
%/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
%/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/

prob=(0.9.*Fitness./max(Fitness))+0.1;%Fitness是一个向量
  
%%%%%%%%%%%%%%%%%%%%%%%% 观察蜂阶段 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%观察蜂阶段
i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))%产生一个0-1的随机数，当此随机数小于prob（i）时，认为当前蜜源被选中，需要进行邻域搜索
        t=t+1;
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValEva=feval(objfun,sol);%这是函数值，越小越好
        FitnessEva=calculateFitness(ObjValEva);%这是收益度值，越大越好
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessEva>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessEva;
            ObjVal(i)=ObjValEva;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
    end;
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end;   
end; 

%/*The best food source is memorized*/
         ind=find(ObjVal==min(ObjVal));
         ind=ind(end);
       if (ObjVal(ind)<GlobalMin)%更新全局函数最小值
         GlobalMin=ObjVal(ind);
         GlobalParams=Foods(ind,:);%更新参数解
       end;
         
%%%%%%%%%%%% 侦查蜂阶段 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/*determine the food sources whose trial counter exceeds the "limit" value. 
%In Basic ABC, only one scout is allowed to occur in each cycle*/
%在基本蜂群算法中，只有每次循环中只有一个侦查蜂
ind=find(trial==max(trial));%找到目前最大的搜索次数
ind=ind(end);
if (trial(ind)>limit)%如果邻域搜索次数大于limit值，则放弃该解，重新生成一个新解
    %Bas(ind)=0;
    trial(ind) = 0;
    sol=(ub-lb).*rand(1,D)+lb;%产生新解
    ObjValSol=feval(objfun,sol);
    FitnessSol=calculateFitness(ObjValSol);%计算产生的新解的收益度值
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end;

fprintf('当er=%d ObjVal=%g\n',iter,GlobalMin);
values{iter} = GlobalMin;
iter=iter+1;
end % End of ABC 
%ABC算法结束，MatLab代码是典型的过程型代码，该代码中没有将过程放在一个个单独的函数中，而是整个一个流程下来

GlobalMins(r) = GlobalMin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 基本蜂群算法中的寻找过程  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial();//初始化解
% MemorizeBestSource();//记住最佳解 
% for (iter=1;iter<=maxCycle;iter++)//maxCycle为觅食的周期数
%     {
% 		SendEmployedBees();//雇佣蜂找寻解,采蜜蜂阶段
% 		CalculateProbabilities();//计算新解适应度
% 		SendOnlookerBees();//观察蜂阶段
% 		MemorizeBestSource();//记忆此时的最佳解
% 		SendScoutBees();//侦查蜂阶段
% 		/**读取结果到文本文件中**/
% 		outfile.open("result.txt",ios_base::app);//打开文本文件
% 		outfile<<"GlobalMax="<<GlobalMax<<' ' << run*10+iter<<"  ";
% 		for(int a=0;a<D;a++)
% 		{
% 			outfile<<'('<<LocParam[a].x<<','<<LocParam[a].y<<')'<<"  ";
% 		}
% 		if(iter%2==0) outfile<<endl;
% 		outfile.close( );
%     }
    
end; %算法运行结束
%a = (1:2500);
% A=cell(1,2500);
% for a=1:2500
%     A{a}=a;
% end;
%scatter(a,values,'k');
B = transpose((cell2mat(values')));
dlmwrite('result.txt',B,'-append','delimiter', ' ');%通过此方法写入的数据将会自动分隔开
%下次读取时就会是4×1000的矩阵，而不再是4000×1的列向量！

% fp = fopen('result.txt','a');
% %fprintf(fp,'程序运行的结果如下：\r\n');
% for i =1 : maxCycle
%     fprintf(fp, '%d\r ', B(i));
% end
% fprintf(fp,'\n');
% fclose(fp);
plot(B);
save all;

