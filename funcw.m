function  w= funcw(x,p,n,m,Fg,t,winit)
%�������ڼ������Ȩ��
%x           input           ��������
%p           input           ����Ⱥ������ֵ
%n           input           ����Ⱥά��
%m           input           ����Ⱥ��ģ
%Fg          input           ������������Ӧ��
%t           input           ��������
%winit       input           w��ʼ��ֵ
%w           output          ��̬����Ȩ��

Dtt=zeros(m,1);
s=0;
for i=1:m
    for d=1:n
        s=s+power(p(1,d)-x(i,d),2);
    end
    Dtt(i,1)=sqrt(s);
end
if(max(Dtt)==0)
    Dt=0;
else
    Dt=(sum(Dtt)/m)/(max(Dtt));
    
end

if(t==1)
    Et=0;
else
    Et=Fg(1,t)/Fg(1,t-1);
end
%         plot(t, Dt,'+');
%         plot(t, Et,'*');
%         hold on;
 w=winit*(1-Dt)*(1-Et);
    

 




