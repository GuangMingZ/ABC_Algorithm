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
D=15; %/*The number o f parameters of the problem to be optimized*//*��ǰ����Ĳ�����Ŀ*/
ub=ones(1,D)*100; %ones����һ��1��DԪ��ȫ��1�ľ���/*lower bounds of the parameters. */����������
lb=ones(1,D)*(-100);%/*upper bound of the parameters.*/����������

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

GlobalMins=zeros(1,runtime);%һ�� 1��runtime ȫ�����.�洢��ǰȫ�ֺ�����Сֵ
values=cell(1,maxCycle); %����һ��1��maxCycle�վ���

for r=1:runtime
  
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */

Range = repmat((ub-lb),[FoodNumber 1]);%FoodNumber��ʳ��Դ�����������ڷ�Ⱥ��Ⱥ��ĿNP��һ��
Lower = repmat(lb, [FoodNumber 1]);

Foods = rand(FoodNumber,D) .* Range + Lower;%rand��M,N������һ����0,1��֮����ֵ��M��Nά����
%����һ��FoodNumber��D �������ɳ�ʼ��

%[y1,..,yn] = Feval_r(F,x1,...,xn) F����Ҫʹ�ú����ĺ����������߾��;xi�Ǻ����Ĳ�����yi�Ǻ����ķ���ֵ��feavl����ָ���κ��Զ���ĺ��� 
%[a1,b1] = size(Foods)
ObjVal=feval('Sphere',Foods);%feval������matlab�Զ���ĺ��������ں�����ӵ��ã�Foods��Ϊobjfun�����Ĳ���
Fitness=calculateFitness(ObjVal);

%reset trial counters
trial=zeros(1,FoodNumber);%һ��1��FoodNumberȫ������

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));%�ҵ�ObjVal�����е���Сֵ�����кţ���Сֵ�и����ж�������ܻ᷵��һ�����С�
%min(ObjVal)���ҵ�ObjVal�����е���Сֵ��
BestInd=BestInd(end);%BestInd(end)��ȡBestInd�����е����һ��Ԫ��ֵ
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);%��ȡ����ĵ�BestInd�У�Ҳ���ǵ�ǰ��Сֵ
iter=1;

while ((iter <= maxCycle)),

%%%%%%%%% ��Ӷ��׶� %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;%fix������βȡ�������㿿£ȡ��.
        % rand����������(0, 1)֮����ȷֲ����������ɵ����顣
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour == i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);%��ȡ��i����Դ��
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       % sol(Param2Change)��Param2Change��Ҫ�ı�Ĳ�����λ�ã������������۷�������������̣�����ĳһ��������һ���ɿص������ڶ�����иı䡣
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       % �������������ɵ��½����Խ��ʱ��������й淶����ֱ��ȡ�߽�ֵ��Ϊ��������ֵ��
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValEva=feval(objfun,sol);%�����½�������ֵ
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

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities �������̶ĸ��� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/* A food source is chosen with the probability which is proportioal to its quality*/
%/*Different schemes can be used to calculate the probability values*/
%/*For example prob(i)=fitness(i)/sum(fitness)*/
%/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
%/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/

prob=(0.9.*Fitness./max(Fitness))+0.1;%Fitness��һ������
  
%%%%%%%%%%%%%%%%%%%%%%%% �۲��׶� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�۲��׶�
i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))%����һ��0-1������������������С��prob��i��ʱ����Ϊ��ǰ��Դ��ѡ�У���Ҫ������������
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
        ObjValEva=feval(objfun,sol);%���Ǻ���ֵ��ԽСԽ��
        FitnessEva=calculateFitness(ObjValEva);%���������ֵ��Խ��Խ��
        
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
       if (ObjVal(ind)<GlobalMin)%����ȫ�ֺ�����Сֵ
         GlobalMin=ObjVal(ind);
         GlobalParams=Foods(ind,:);%���²�����
       end;
         
%%%%%%%%%%%% ����׶� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/*determine the food sources whose trial counter exceeds the "limit" value. 
%In Basic ABC, only one scout is allowed to occur in each cycle*/
%�ڻ�����Ⱥ�㷨�У�ֻ��ÿ��ѭ����ֻ��һ������
ind=find(trial==max(trial));%�ҵ�Ŀǰ������������
ind=ind(end);
if (trial(ind)>limit)%�������������������limitֵ��������ý⣬��������һ���½�
    %Bas(ind)=0;
    trial(ind) = 0;
    sol=(ub-lb).*rand(1,D)+lb;%�����½�
    ObjValSol=feval(objfun,sol);
    FitnessSol=calculateFitness(ObjValSol);%����������½�������ֵ
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end;

fprintf('��er=%d ObjVal=%g\n',iter,GlobalMin);
values{iter} = GlobalMin;
iter=iter+1;
end % End of ABC 
%ABC�㷨������MatLab�����ǵ��͵Ĺ����ʹ��룬�ô�����û�н����̷���һ���������ĺ����У���������һ����������

GlobalMins(r) = GlobalMin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ������Ⱥ�㷨�е�Ѱ�ҹ���  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial();//��ʼ����
% MemorizeBestSource();//��ס��ѽ� 
% for (iter=1;iter<=maxCycle;iter++)//maxCycleΪ��ʳ��������
%     {
% 		SendEmployedBees();//��Ӷ����Ѱ��,���۷�׶�
% 		CalculateProbabilities();//�����½���Ӧ��
% 		SendOnlookerBees();//�۲��׶�
% 		MemorizeBestSource();//�����ʱ����ѽ�
% 		SendScoutBees();//����׶�
% 		/**��ȡ������ı��ļ���**/
% 		outfile.open("result.txt",ios_base::app);//���ı��ļ�
% 		outfile<<"GlobalMax="<<GlobalMax<<' ' << run*10+iter<<"  ";
% 		for(int a=0;a<D;a++)
% 		{
% 			outfile<<'('<<LocParam[a].x<<','<<LocParam[a].y<<')'<<"  ";
% 		}
% 		if(iter%2==0) outfile<<endl;
% 		outfile.close( );
%     }
    
end; %�㷨���н���
%a = (1:2500);
% A=cell(1,2500);
% for a=1:2500
%     A{a}=a;
% end;
%scatter(a,values,'k');
B = transpose((cell2mat(values')));
dlmwrite('result.txt',B,'-append','delimiter', ' ');%ͨ���˷���д������ݽ����Զ��ָ���
%�´ζ�ȡʱ�ͻ���4��1000�ľ��󣬶�������4000��1����������

% fp = fopen('result.txt','a');
% %fprintf(fp,'�������еĽ�����£�\r\n');
% for i =1 : maxCycle
%     fprintf(fp, '%d\r ', B(i));
% end
% fprintf(fp,'\n');
% fclose(fp);
plot(B);
save all;

