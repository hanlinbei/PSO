function  w= funcw(x,p,n,m,Fg,t,winit)
%函数用于计算惯性权重
%x           input           输入粒子
%p           input           粒子群体最优值
%n           input           粒子群维度
%m           input           粒子群规模
%Fg          input           迭代过程中适应度
%t           input           迭代次数
%winit       input           w初始化值
%w           output          动态惯性权重

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
    

 




