function G = funG(x,T,N,M,g)
%函数用于计算惩罚费用
%x           input           输入粒子
%A           input           生产能力矩阵
%T           input           计划周期
%M           input           装配线条数
%N           input           产品种数
%C1           output         惩罚费用 
G=zeros(M,T);
s=0;
for i=1:M
    for t=1:T
        for j=1:N
           s=s+g(i,j)*x(1,(t-1)*M*N+(j-1)*M+i);
        end
        G(i,t)=s;
        s=0;
    end
end