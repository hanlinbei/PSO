function C1 = func1(x,R,V,T,N,M)
%函数用于计算惩罚费用
%x           input           输入粒子
%R           input           第t天需要交货j类产品数量
%T           input           计划周期
%M           input           装配线条数
%N           input           产品种数
%V           input           惩罚系数
%C1           output         惩罚费用 
%% 到第t天时j种产品的生产量
Q=zeros(N,T);
for j=1:N
    s=0;
  for q=1:T
    for i=1:M
        s=s+x(:,(q-1)*M*N+(j-1)*M+i);
    end
        Q(j,q)=s;
  end
end
%% 拖期费用
y=0;
for t=1:T
    for j=1:N
            y=y+V(1,j)*max(R(j,t)-Q(j,t),0);
            C1=y;
    end
end
