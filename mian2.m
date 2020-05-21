%% ��׼����Ⱥ�㷨��������Ų�
%% I. ��ջ���
clc
clear
%% II. ��������
load F.mat
load A.mat
load V.mat
load dk.mat
load D.mat
load H.mat
%% III. ������ʼ��
c1 = 2.0;
c2 = 2.0;
M = 4;  % װ��������
N = 6;    % ������Ʒ����
T = 30;     % �ƻ�����
L = 8;       %������
b= 100;     %��С��������

maxgen = 500;   % ��������  
sizepop = 30;   %��Ⱥ��ģ

a = zeros(L,T) ;%����Ok�ڵ�t���Ƿ񽻻�����
for k=1:L
    for t=1:T
        if(t==dk(k,:))
            a(k,t)=1;
        end
    end
end

R = zeros(N,T) ;%����t��ʱ��Ҫ������j�ֲ�Ʒ������ 
for j=1:N
    s=0;
 for q=1:t
    for k=1:L
        s =s+D(k,j)*a(k,q);
    end
        R(j,q)=s;
 end
end

%% IV. ������ʼ���Ӻ��ٶ�
pop=zeros(sizepop,M*N*T);
v=zeros(sizepop,M*N*T);
fitness=zeros(sizepop,1);
ISIF=zeros(sizepop,1);
for i = 1:sizepop
    %��ʼ��ȺΪ��С����b
    for j=1:M*N*T
        pop(i,j)=b;
        v(i,j)=10*rand;%��ʼ���ٶ�
    end
    % ������Ӧ��
    fitness(i,1) = funcm(pop(i,:),F,T,M,N)+func1(pop(i,:),R,V,T,N,M);   %Ⱦɫ�����Ӧ��
    ISIF(i)=1; %��¼�����������Ƿ�����Լ��
end


%% V. ���弫ֵ��Ⱥ�弫ֵ
[bestfitness, bestindex] = max(fitness);
zbest = pop(bestindex,:);   %ȫ�����
gbest = pop;    %�������
fitnessgbest = fitness;   %���������Ӧ��ֵ
fitnesszbest = bestfitness;   %ȫ�������Ӧ��ֵ

%% VI. ����Ѱ��
Fg=zeros(1,maxgen);
Cm=zeros(sizepop,1);
C1=zeros(sizepop,1);
for i = 1:maxgen
     Fg(i) = fitnesszbest; 

    for j = 1:sizepop

        % �ٶȸ���
        v(j,:) = round(v(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + v(j,:); 

        G=funG(pop(j,:),T,N,M,1./A);
        
        % �޸���Ⱥ��ʹ����Լ��
        if(max(G-H)>0)
            pop(j,:)= pop(j,:) - v(j,:);
        end
        pop(j,pop(j,:)<b) = 0;
        
        % ��Ӧ��ֵ����
        Cm(j,1)=funcm(pop(j,:),F,T,M,N);
        C1(j,1)=func1(pop(j,:),R,V,T,N,M);
        fitness(j,1) = Cm(j,1)+ C1(j,1);
    end

    
    for j = 1:sizepop  
        % �������Ÿ���
        if fitness(j,1) > fitnessgbest(j,1)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j,1) = fitness(j,1);
        end
        
        % Ⱥ�����Ÿ���
        if fitness(j,1) > fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j,1);
            Cmbest=Cm(j,1);
            C1best=C1(j,1);
        end
    end 
           
end

%% VII.������

figure
plot(Fg)
grid on;
title('���Ÿ�����Ӧ��','fontsize',12);
xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);

H1=[8 8 8;8 8 8;8 8 8]    