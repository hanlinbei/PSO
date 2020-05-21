%% �Ľ�����Ⱥ�㷨�������ģ��
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
winit =1.2;
wend = 0.4;
M = 4;  % װ��������
N = 6;    % ������Ʒ����
T = 30;     % �ƻ�����
L = 8;       %������
b= 100;     %��С��������

maxgen = 500;   % ��������  
sizepop = 30;   %��Ⱥ��ģ
popmax = 200;
popmin = 100;

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
for r=1:1
pop=zeros(sizepop,M*N*T);
v=zeros(sizepop,M*N*T);
fitness=zeros(sizepop,1);
for i = 1:sizepop
    %��ʼ��ȺΪ��С����b
    for j=1:M*N*T
        pop(i,j)=b;
        v(i,j)=2*rand(1);%��ʼ���ٶ�
    end
    % ������Ӧ��
    fitness(i,1) = funcm(pop(i,:),F,T,M,N)+func1(pop(i,:),R,V,T,N,M);   %Ⱦɫ�����Ӧ��
   
end


%% V. ���弫ֵ��Ⱥ�弫ֵ
[bestfitness, bestindex] = min(fitness);
zbest = pop(bestindex,:);   %ȫ�����
gbest = pop;    %�������
fitnessgbest = fitness;   %���������Ӧ��ֵ
fitnesszbest = bestfitness;   %ȫ�������Ӧ��ֵ

%% VI. ����Ѱ��
w=winit;
Et=0;
Fg=zeros(1,maxgen);
Cm=zeros(sizepop,1);
C1=zeros(sizepop,1);
for i = 1:maxgen
     Fg(i) = fitnesszbest; 
         % ��̬����Ȩ��
         w=funcw(gbest,zbest,M*N*T,sizepop,Fg,i,winit);
%         w=winit-(winit-wend)*(i-1)/(maxgen-1);
    for j = 1:sizepop
        % �ٶȸ���
        v(j,:) = round(w*v(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + v(j,:); 
        

         
        G=funG(pop(j,:),T,N,M,1./A);
        uu=G-H;

        
        % �޸���Ⱥ��ʹ����Լ��1
%         if(max(max(G-H))>0)
            for t=1:T
             for j2=1:N
               for i2=1:M
                   if(uu(i2,t)>0)
                  pop(j,(t-1)*M*N+(j2-1)*M+i2)=min(pop(j,(t-1)*M*N+(j2-1)*M+i2)-v(j,(t-1)*M*N+(j2-1)*M+i2),H(i2,t)*A(i2,j2)/M);   
                   else
                        if(pop(j,(t-1)*M*N+(j2-1)*M+i2)<b)
                        pop(j,(t-1)*M*N+(j2-1)*M+i2)=0;
                        end 
                   end
               end
            end
            end
 

        % ��Ӧ��ֵ����
        fitness(j,1) = funcm(pop(j,:),F,T,M,N)+func1(pop(j,:),R,V,T,N,M);
        % �������Ÿ���
        if fitness(j,1) < fitnessgbest(j,1)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j,1) = fitness(j,1);
        end
             % Ⱥ�����Ÿ���
        if fitness(j,1) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j,1);
            

        end
    end
   
           
end

%% VII.������
display(fitnesszbest);
z(r)=fitnesszbest;
plot(1,z(r),'o');
hold on
        G1=funG(zbest,T,N,M,1./A);
        uu1=G1-H;
end
% figure
% plot(Fg)
% grid on;
% title('���Ÿ�����Ӧ��','fontsize',12);
% xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);
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
