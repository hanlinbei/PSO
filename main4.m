%% �Ľ�����Ⱥ���shuffetģ��
%% I. ��ջ���
tic
clc
clear
[x1,x2]=meshgrid(-10:0.1:10);

figure(1);mesh(x1,x2,fitness1(x1,x2));%����Shubert����ͼ��
hold on

%% III. ������ʼ��
c1 = 2;
c2 = 2;
winit = 1.2;
wend = 0.4;
maxgen = 100;   % ��������  
sizepop = 30;   %��Ⱥ��ģ

Vmax = 1;
Vmin = -1;
popmax = 10;
popmin = -10;

for r=1:20
    
%% IV. ������ʼ���Ӻ��ٶ�
for i = 1:sizepop
    % �������һ����Ⱥ
    pop(i,:) = 10*rands(1,2);    %��ʼ��Ⱥ
    V(i,:) = rands(1,2);  %��ʼ���ٶ�
    % ������Ӧ��
    fitness(i) = fitness1(pop(i,1),pop(i,2));   %Ⱦɫ�����Ӧ��

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
         w=winit-(winit-wend)*i/maxgen;
%          w=funcw(pop,zbest,2,sizepop,yy,i,winit);
    for j = 1:sizepop
        % �ٶȸ���
%          plot(i, w,'*');
%          hold on;
        V(j,:) = V(j,:)*w + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
%         V(j,find(V(j,:)>Vmax)) = Vmax;
%         V(j,find(V(j,:)<Vmin)) = Vmin;
        
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + V(j,:);
         pop(j,find(pop(j,:)>popmax)) = popmax;
         pop(j,find(pop(j,:)<popmin)) = popmin;
        
        % ��Ӧ��ֵ����
        fitness(j) = fitness1(pop(j,1),pop(j,2)); 




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
% [fitnesszbest, zbest]
% z(r)=fitnesszbest;
% plot(1,z(r),'o');
% hold on
end
%% VII.������

 plot3(zbest(1), zbest(2), fitnesszbest,'bo','linewidth',1.5)

% figure
% plot(yy)
% title('���Ÿ�����Ӧ��','fontsize',12);
% xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);
toc



