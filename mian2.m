%% 标准粒子群算法求解手套排产
%% I. 清空环境
clc
clear
%% II. 导入数据
load F.mat
load A.mat
load V.mat
load dk.mat
load D.mat
load H.mat
%% III. 参数初始化
c1 = 2.0;
c2 = 2.0;
M = 4;  % 装配线条数
N = 6;    % 生产产品种数
T = 30;     % 计划周期
L = 8;       %订单数
b= 100;     %最小生产批量

maxgen = 500;   % 进化次数  
sizepop = 30;   %种群规模

a = zeros(L,T) ;%订单Ok在第t天是否交货矩阵
for k=1:L
    for t=1:T
        if(t==dk(k,:))
            a(k,t)=1;
        end
    end
end

R = zeros(N,T) ;%到第t天时需要交货的j种产品的数量 
for j=1:N
    s=0;
 for q=1:t
    for k=1:L
        s =s+D(k,j)*a(k,q);
    end
        R(j,q)=s;
 end
end

%% IV. 产生初始粒子和速度
pop=zeros(sizepop,M*N*T);
v=zeros(sizepop,M*N*T);
fitness=zeros(sizepop,1);
ISIF=zeros(sizepop,1);
for i = 1:sizepop
    %初始种群为最小批量b
    for j=1:M*N*T
        pop(i,j)=b;
        v(i,j)=10*rand;%初始化速度
    end
    % 计算适应度
    fitness(i,1) = funcm(pop(i,:),F,T,M,N)+func1(pop(i,:),R,V,T,N,M);   %染色体的适应度
    ISIF(i)=1; %记录迭代过程中是否满足约束
end


%% V. 个体极值和群体极值
[bestfitness, bestindex] = max(fitness);
zbest = pop(bestindex,:);   %全局最佳
gbest = pop;    %个体最佳
fitnessgbest = fitness;   %个体最佳适应度值
fitnesszbest = bestfitness;   %全局最佳适应度值

%% VI. 迭代寻优
Fg=zeros(1,maxgen);
Cm=zeros(sizepop,1);
C1=zeros(sizepop,1);
for i = 1:maxgen
     Fg(i) = fitnesszbest; 

    for j = 1:sizepop

        % 速度更新
        v(j,:) = round(v(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        
        % 种群更新
        pop(j,:) = pop(j,:) + v(j,:); 

        G=funG(pop(j,:),T,N,M,1./A);
        
        % 修复种群，使满足约束
        if(max(G-H)>0)
            pop(j,:)= pop(j,:) - v(j,:);
        end
        pop(j,pop(j,:)<b) = 0;
        
        % 适应度值更新
        Cm(j,1)=funcm(pop(j,:),F,T,M,N);
        C1(j,1)=func1(pop(j,:),R,V,T,N,M);
        fitness(j,1) = Cm(j,1)+ C1(j,1);
    end

    
    for j = 1:sizepop  
        % 个体最优更新
        if fitness(j,1) > fitnessgbest(j,1)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j,1) = fitness(j,1);
        end
        
        % 群体最优更新
        if fitness(j,1) > fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j,1);
            Cmbest=Cm(j,1);
            C1best=C1(j,1);
        end
    end 
           
end

%% VII.输出结果

figure
plot(Fg)
grid on;
title('最优个体适应度','fontsize',12);
xlabel('进化代数','fontsize',12);ylabel('适应度','fontsize',12);

H1=[8 8 8;8 8 8;8 8 8]    