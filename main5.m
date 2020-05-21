%% 改进粒子群求解shuffet模型
%% I. 清空环境
clc
clear
ai=[20;18;22;16];
tij=[2 1 3;2.5 2 1.5;1 3 2.5;1.5 2 1];
bj=[20;16;18];
cij=[3 6 2;2 4 4;5 2 2;3 3 5];

%% III. 参数初始化
c1 = 2;
c2 = 2;
winit = 1.2;
wend = 0.4;
maxgen = 100;   % 进化次数  
sizepop = 30;   %种群规模

Vmax = 1;
Vmin = -1;
popmax = 30;
popmin = 0;


    
%% IV. 产生初始粒子和速度
for i = 1:1
    % 随机产生一个种群
    flag=1;
    while(flag==1)
    pop(i,:) = round(20*rand(1,12));    %初始种群\
      for j3=1:3
            s3=0;
            for i3=1:4
                s3=s3+pop(i,(i3-1)*3+j3);
                flag=0;
            end
            if(s3~=bj(j3))
                flag=1;
                break;
            end
      end
      
      for i4=1:4
    ss=0;
        for j4=1:3
            ss=ss+tij(i4,j4)*pop(i,(i4-1)*3+j4);
        end
        if(ss>ai(i4))
            flag=1;
            break;
        end
     end
      
    end
    
    V(i,:) = 1*rands(1,12);  %初始化速度
    % 计算适应度
%      fitness(i) = fitness2(cij,pop(i,:));   %染色体的适应度

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
%         yy(i) = fitnesszbest; 
%          w=winit-(winit-wend)*i/maxgen;
%          w=funcw(pop,zbest,2,sizepop,yy,i,winit);
    for j = 1:sizepop
        % 速度更新
        V(j,:) = round(V(j,:)*1 + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
%         V(j,find(V(j,:)>Vmax)) = Vmax;
%         V(j,find(V(j,:)<Vmin)) = Vmin;
        
        % 种群更新
        pop(j,:) = pop(j,:) + V(j,:);
      
          
for i1=1:4
    ss=0;
    for j1=1:3
        ss=ss+tij(i1,j1)*pop(j,(i1-1)*3+j1);
    end
    if(ss>ai(i1))
       for j1=1:3
        pop(j,(i1-1)*3+j1)=pop(j,(i1-1)*3+j1)- V(j,(i1-1)*3+j1);
       end 
    end
end
       pop(j,find(pop(j,:)>popmax)) = popmax;
          pop(j,find(pop(j,:)<popmin)) = popmin;

     
        fitness(j) = fitness2(cij,pop(j,:));  
        % 适应度值更新
        for j2=1:3
            s2=0;
            for i2=1:4
                s2=s2+pop(j,(i2-1)*3+j2);
            end
            if(s2<bj(j2))
                fitness(j)=10000;
            end
        end


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
[fitnesszbest, zbest]
% z(r)=fitnesszbest;
% plot(1,z(r),'o');
% hold on

%% VII.输出结果

%  plot3(zbest(1), zbest(2), fitnesszbest,'bo','linewidth',1.5)

% figure
% plot(yy)
% title('最优个体适应度','fontsize',12);
% xlabel('进化代数','fontsize',12);ylabel('适应度','fontsize',12);
A=[1,1,1,1;1,1,1,1;1,1,1,1];
B=[20;16;18];
X=A\B;
syms x1 x2 x3 x4 x5 x6
[x,params,conds]=solve(x1+x2+x3+x4+x5+x6==6,'ReturnConditions', true); 
