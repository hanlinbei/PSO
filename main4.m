%% 改进粒子群求解shuffet模型
%% I. 清空环境
tic
clc
clear
[x1,x2]=meshgrid(-10:0.1:10);

figure(1);mesh(x1,x2,fitness1(x1,x2));%画出Shubert函数图象
hold on

%% III. 参数初始化
c1 = 2;
c2 = 2;
winit = 1.2;
wend = 0.4;
maxgen = 100;   % 进化次数  
sizepop = 30;   %种群规模

Vmax = 1;
Vmin = -1;
popmax = 10;
popmin = -10;

for r=1:20
    
%% IV. 产生初始粒子和速度
for i = 1:sizepop
    % 随机产生一个种群
    pop(i,:) = 10*rands(1,2);    %初始种群
    V(i,:) = rands(1,2);  %初始化速度
    % 计算适应度
    fitness(i) = fitness1(pop(i,1),pop(i,2));   %染色体的适应度

end

%% V. 个体极值和群体极值
[bestfitness bestindex] = min(fitness);
zbest = pop(bestindex,:);   %全局最佳
gbest = pop;    %个体最佳
fitnessgbest = fitness;   %个体最佳适应度值
fitnesszbest = bestfitness;   %全局最佳适应度值

%% VI. 迭代寻优


    w=winit;
for i = 1:maxgen
        yy(i) = fitnesszbest; 
         w=winit-(winit-wend)*i/maxgen;
%          w=funcw(pop,zbest,2,sizepop,yy,i,winit);
    for j = 1:sizepop
        % 速度更新
%          plot(i, w,'*');
%          hold on;
        V(j,:) = V(j,:)*w + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
%         V(j,find(V(j,:)>Vmax)) = Vmax;
%         V(j,find(V(j,:)<Vmin)) = Vmin;
        
        % 种群更新
        pop(j,:) = pop(j,:) + V(j,:);
         pop(j,find(pop(j,:)>popmax)) = popmax;
         pop(j,find(pop(j,:)<popmin)) = popmin;
        
        % 适应度值更新
        fitness(j) = fitness1(pop(j,1),pop(j,2)); 




        % 个体最优更新
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        % 群体最优更新
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end

end
% [fitnesszbest, zbest]
% z(r)=fitnesszbest;
% plot(1,z(r),'o');
% hold on
end
%% VII.输出结果

 plot3(zbest(1), zbest(2), fitnesszbest,'bo','linewidth',1.5)

% figure
% plot(yy)
% title('最优个体适应度','fontsize',12);
% xlabel('进化代数','fontsize',12);ylabel('适应度','fontsize',12);
toc



