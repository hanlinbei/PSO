%% �Ľ�����Ⱥ�㷨�������ģ��
%% I. ��ջ���
clc
clear
%% II. ��������
load F1.mat
load A1.mat
load V1.mat
load dk1.mat
load D1.mat
load H1.mat
%% III. ������ʼ��
% D1(1,1)=1500;
c1 = 2.0;
c2 = 2.0;
winit =1.2;
wend = 0.4;
M = 3;  % װ��������
N = 3;    % ������Ʒ����
T = 3;     % �ƻ�����
L = 3;       %������
b= 10;     %��С��������

maxgen = 500;   % ��������  
sizepop = 100;   %��Ⱥ��ģ
popmax = 480;
popmin = 0;

a = zeros(L,T) ;%����Ok�ڵ�t���Ƿ񽻻�����
for k=1:L
    for t=1:T
        if(t==dk1(k,:))
            a(k,t)=1;
        end
    end
end

R = zeros(N,T) ;%����t��ʱ��Ҫ������j�ֲ�Ʒ������ 
for j=1:N
    s=0;
 for q=1:t
    for k=1:L
        s =s+D1(k,j)*a(k,q);
    end
        R(j,q)=s;
 end
end

%% IV. ������ʼ���Ӻ��ٶ�
for r=1:10
pop=zeros(sizepop,M*N*T);
v=zeros(sizepop,M*N*T);
fitness=zeros(sizepop,1);
for i = 1:sizepop
    %��ʼ��ȺΪ��С����b
    for j=1:M*N*T
        pop(i,j)=10;
        v(i,j)=10*rand(1);%��ʼ���ٶ�
    end
    % ������Ӧ��
    fitness(i,1) = funcm(pop(i,:),F1,T,M,N)+func1(pop(i,:),R,V1,T,N,M);   %Ⱦɫ�����Ӧ��
   
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
        v(j,:) = round(1*v(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + v(j,:);        
         
        G=funG(pop(j,:),T,N,M,1./A1);
        uu=G-H1;

        
        % �޸���Ⱥ��ʹ����Լ��1
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

        % ��Ӧ��ֵ����
        fitness(j,1) = funcm(pop(j,:),F1,T,M,N)+func1(pop(j,:),R,V1,T,N,M);
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
% title('���Ÿ�����Ӧ��','fontsize',12);
% xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);
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

