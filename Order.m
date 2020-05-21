%% 改进粒子群算法求解手套模型
%% I. 清空环境
clc
clear
%% II. 导入数据
load F1.mat
load A1.mat
load V1.mat
load dk1.mat
load D1.mat
load H1.mat
%% III. 参数初始化
% D1(1,1)=1500;
c1 = 2.0;
c2 = 2.0;
winit =1.2;
wend = 0.4;
M = 3;  % 装配线条数
N = 3;    % 生产产品种数
T = 3;     % 计划周期
L = 3;       %订单数
b= 10;     %最小生产批量

maxgen = 500;   % 进化次数  
sizepop = 100;   %种群规模
popmax = 480;
popmin = 0;

a = zeros(L,T) ;%订单Ok在第t天是否交货矩阵
for k=1:L
    for t=1:T
        if(t==dk1(k,:))
            a(k,t)=1;
        end
    end
end

R = zeros(N,T) ;%到第t天时需要交货的j种产品的数量 
for j=1:N
    s=0;
 for q=1:t
    for k=1:L
        s =s+D1(k,j)*a(k,q);
    end
        R(j,q)=s;
 end
end

%% IV. 产生初始粒子和速度
for r=1:10
pop=zeros(sizepop,M*N*T);
v=zeros(sizepop,M*N*T);
fitness=zeros(sizepop,1);
for i = 1:sizepop
    %初始种群为最小批量b
    for j=1:M*N*T
        pop(i,j)=10;
        v(i,j)=10*rand(1);%初始化速度
    end
    % 计算适应度
    fitness(i,1) = funcm(pop(i,:),F1,T,M,N)+func1(pop(i,:),R,V1,T,N,M);   %染色体的适应度
   
end


%% V. 个体极值和群体极值
[bestfitness, bestindex] = min(fitness);
zbest = pop(bestindex,:);   %全局最佳
gbest = pop;    %个体最佳
fitnessgbest = fitness;   %个体最佳适应度值
fitnesszbest = bestfitness;   %全局最佳适应度值

%% VI. 迭代寻优
w=winit;
Et=0;
Fg=zeros(1,maxgen);
Cm=zeros(sizepop,1);
C1=zeros(sizepop,1);
for i = 1:maxgen
     Fg(i) = fitnesszbest; 
         % 动态惯性权重
          w=funcw(gbest,zbest,M*N*T,sizepop,Fg,i,winit);
%         w=winit-(winit-wend)*(i-1)/(maxgen-1);
    for j = 1:sizepop
        % 速度更新
        v(j,:) = round(1*v(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        % 种群更新
        pop(j,:) = pop(j,:) + v(j,:);        
         
        G=funG(pop(j,:),T,N,M,1./A1);
        uu=G-H1;

        
        % 修复种群，使满足约束1
%         if(max(max(G-H))>0)
            for t=1:T
             for j2=1:N
               for i2=1:M
                   if(uu(i2,t)>0)
                  pop(j,(t-1)*M*N+(j2-1)*M+i2)=min(pop(j,(t-1)*M*N+(j2-1)*M+i2)-v(j,(t-1)*M*N+(j2-1)*M+i2),round(H1(i2,t)*A1(i2,j2)/N));   
                   else
                        if(pop(j,(t-1)*M*N+(j2-1)*M+i2)<b)
                        pop(j,(t-1)*M*N+(j2-1)*M+i2)=0;
                        end 
                   end
               end
            end
            end
         pop(j,find(pop(j,:)>popmax)) = popmax;
        pop(j,find(pop(j,:)<popmin)) = popmin;

        % 适应度值更新
        fitness(j,1) = funcm(pop(j,:),F1,T,M,N)+func1(pop(j,:),R,V1,T,N,M);
        % 个体最优更新
        if fitness(j,1) < fitnessgbest(j,1)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j,1) = fitness(j,1);
        end
             % 群体最优更新
        if fitness(j,1) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j,1);            
        end
    end
   
           
end

%% VII.输出结果
display(fitnesszbest);
zbest
z(r)=fitnesszbest;
plot(1,z(r),'o');
hold on
% o=[77 78 78 77 78 78 133 133 134 77 78 78 77 78 78 133 133 134 78 78 78 78 78 78 133 133 134];
%         G1=funG(o,T,N,M,1./A1);
%         uu1=G1-H1;
cmm=funcm(zbest,F1,T,M,N);
c11=func1(zbest,R,V1,T,N,M);
% figure
% plot(Fg)
% grid on;
% title('最优个体适应度','fontsize',12);
% xlabel('进化代数','fontsize',12);ylabel('适应度','fontsize',12);
end
%         G1=funG(zbest,T,N,M,1./A);
%         uu1=G1-H;


% B=zeros(N,T);
% 
%     for t=1:T
%         for j=1:N
%            B(j,t)=zbest(1,(t-1)*M*N+(j-1)*M+1); 
%         end
%         
%     end

