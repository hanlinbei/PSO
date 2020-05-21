function Cm = funcm(x,f,T,M,N)
%函数用于计算加工费用
%x           input           输入粒子
%f           input           装配线装配单位产品所需的费用
%T           input           计划周期
%M           input           装配线条数
%N           input           产品种数
%Cm           output         加工费用 
y=0;
for t=1:T
    for j=1:N
        for i=1:M
            y=y+f(i,j)*x(1,(t-1)*M*N+(j-1)*M+i);
            Cm=y;
        end
    end
end


