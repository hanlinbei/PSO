function G = funG(x,T,N,M,g)
%�������ڼ���ͷ�����
%x           input           ��������
%A           input           ������������
%T           input           �ƻ�����
%M           input           װ��������
%N           input           ��Ʒ����
%C1           output         �ͷ����� 
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