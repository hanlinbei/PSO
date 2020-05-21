%% �Ľ�����Ⱥ�㷨�����Ӳ�Ʒģ��
tic
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
load FgL.mat
%% III. ������ʼ��
% D1(1,1)=1500;
c1 = 2.0;
c2 = 2.0;
winit =1.0;
wend = 0.4;
M = 4;  % װ��������
N = 6;    % ������Ʒ����
T = 30;     % �ƻ�����
L = 8;       %������
b= 100;     %��С��������

maxgen = 500;   % ��������  
sizepop = 30;   %��Ⱥ��ģ
popmax = 480;
popmin = 0;

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
%v=zeros(sizepop,M*N*T);
v=round(2*rand(sizepop,M*N*T));
fitness=zeros(sizepop,1);
for i = 1:sizepop
    %��ʼ��ȺΪ��С����b
    for j=1:M*N*T
        %v(i,j)=2*rand(1);%��ʼ���ٶ�
        pop(i,j)=b+v(i,j);
        
    end

    G=funG(pop(i,:),T,N,M,1./A);
        uu=G-H;
        
        for t=1:T
             for j2=1:N
               for i2=1:M
                   if(uu(i2,t)>0)
                   pop(i,(t-1)*M*N+(j2-1)*M+i2)=min(pop(i,(t-1)*M*N+(j2-1)*M+i2)-v(i,(t-1)*M*N+(j2-1)*M+i2),round(H(i2,t)*A(i2,j2)/N)); 
                   end
%                    else
                    if(pop(i,(t-1)*M*N+(j2-1)*M+i2)<b)
                    pop(i,(t-1)*M*N+(j2-1)*M+i2)=0;
                    end                               
                end
              end
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

%  G1=funG(zbest,T,N,M,1./A);
%         uu1=G1-H;

%% VI. ����Ѱ��
w=winit;
Et=0;
Fg=zeros(1,maxgen);
Cm=zeros(sizepop,1);
C1=zeros(sizepop,1);
for i = 1:maxgen
     Fg(i) = fitnesszbest; 
%          ��̬����Ȩ��
          w=funcw(gbest,zbest,M*N*T,sizepop,Fg,i,winit);
          
%          w=winit-(winit-wend)*(i-1)/(maxgen-1);
w1(i)=w;
    for j = 1:sizepop
        % �ٶȸ���
        v(j,:) = round(w*v(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:)));
        % ��Ⱥ����
        pop(j,:) = pop(j,:) + v(j,:);        
         
        G=funG(pop(j,:),T,N,M,1./A);
        uu=G-H;

        
        % �޸���Ⱥ��ʹ����Լ��1
%          if(max(max(G-H))>0)

                
            for t=1:T
             for j2=1:N
               for i2=1:M
                   if(uu(i2,t)>0)
                   pop(j,(t-1)*M*N+(j2-1)*M+i2)=min(pop(j,(t-1)*M*N+(j2-1)*M+i2)-v(j,(t-1)*M*N+(j2-1)*M+i2),round(H(i2,t)*A(i2,j2)/N)); 
                   end
%                    else
                    if(pop(j,(t-1)*M*N+(j2-1)*M+i2)<b)
                    pop(j,(t-1)*M*N+(j2-1)*M+i2)=0;
                    end 
                               
                end
              end
            end
     
         %pop(j,find(pop(j,:)>popmax)) = popmax;
       % pop(j,find(pop(j,:)<popmin)) = popmin;

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
zbest;

% o=[77 78 78 77 78 78 133 133 134 77 78 78 77 78 78 133 133 134 78 78 78 78 78 78 133 133 134];
%         G1=funG(o,T,N,M,1./A1);
%         uu1=G1-H1;
cmm=funcm(zbest,F,T,M,N);
c11=func1(zbest,R,V,T,N,M);
figure
plot(Fg,'r');
% hold on
% plot(FgL,'g');
grid on;
title('���Ÿ�����Ӧ��','fontsize',12);
xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);

        G1=funG(zbest,T,N,M,1./A);
        uu1=G1-H;


% B=zeros(N,T);
% 
%     for t=1:T
%         for j=1:N
%            B(j,t)=zbest(1,(t-1)*M*N+(j-1)*M+1); 
%         end
%         
%     end

out=zeros(T,N);
for t=1:T
    for j=1:N
        for i=1:M
           out(t,j)=out(t,j)+zbest(1,(t-1)*M*N+(j-1)*M+i); 
        end
    end
end

out4=zeros(T,N);
for t=1:T
    for j=1:N
        if(t==1)
           out4(t,j)=out(t,j);
        else           
          out4(t,j)=out(t,j)+out4(t-1,j); 
        end
    end
end

out2=zeros(N);
 for j=1:N
    for t=1:T
        for i=1:M
        out2(j)=out2(j)+zbest(1,(t-1)*M*N+(j-1)*M+i);
        end
    end
 end
 
 out3=zeros(M,T);
 for j=1:N
    for t=1:T
        for i=1:M
        out3(i,t)=out3(i,t)+zbest(1,(t-1)*M*N+(j-1)*M+i);
        end
    end
 end
 
% A=[3 5 6;7 5 9;8 4 2];
% 
% B=max(0,A);
% for j=1:N
% han(j,:)=zbest(1,(6-1)*M*N+(j-1)*M+1);
% end

toc
