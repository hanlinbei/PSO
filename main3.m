%% �Ľ�����Ⱥ�����ͨģ��
%% I. ��ջ���
clc
clear

%% II. ����Ŀ�꺯������
% figure
% [x,y] = meshgrid(-5:0.1:5,-5:0.1:5);
% z = x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y) + 20;
% mesh(x,y,z)
% hold on

%% III. ������ʼ��
c1 = 2;
c2 = 2;
winit = 1.2;
wend = 0.4;
maxgen = 100;   % ��������  
sizepop = 100;   %��Ⱥ��ģ

Vmax = 1;
Vmin = -1;
popmax = 5;
popmin = -5;

for r=1:20
    
%% IV. ������ʼ���Ӻ��ٶ�
for i = 1:sizepop
    % �������һ����Ⱥ
    pop(i,:) = 5*rands(1,2);    %��ʼ��Ⱥ
    V(i,:) = rands(1,2);  %��ʼ���ٶ�
    % ������Ӧ��
    fitness(i) = -fun(pop(i,:));   %Ⱦɫ�����Ӧ��

end

%% V. ���弫ֵ��Ⱥ�弫ֵ
[bestfitness bestindex] = min(fitness);
zbest = pop(bestindex,:);   %ȫ�����
gbest = pop;    %�������
fitnessgbest = fitness;   %���������Ӧ��ֵ
fitnesszbest = bestfitness;   %ȫ�������Ӧ��ֵ

%% VI. ����Ѱ��


    w=winit;
for i = 1:maxgen
        yy(i) = fitnesszbest; 
%         w=winit-(winit-wend)*i/maxgen;
        w=funcw(pop,zbest,2,sizepop,yy,i,winit);
    for j = 1:sizepop
        % �ٶȸ���
        V(j,:) = V(j,:)*w + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
%         V(j,find(V(j,:)>Vmax)) = Vmax;
%         V(j,find(V(j,:)<Vmin)) = Vmin;
        
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + V(j,:);
        pop(j,find(pop(j,:)>popmax)) = popmax;
        pop(j,find(pop(j,:)<popmin)) = popmin;
        
        % ��Ӧ��ֵ����
        fitness(j) = -fun(pop(j,:)); 




        % �������Ÿ���
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        % Ⱥ�����Ÿ���
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end

end
[fitnesszbest, zbest]
z(r)=fitnesszbest;
plot(1,z(r),'o');
hold on
end
%% VII.������

% plot3(zbest(1), zbest(2), fitnesszbest,'bo','linewidth',1.5)
% 
% figure
% plot(yy)
% title('���Ÿ�����Ӧ��','fontsize',12);
% xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);



