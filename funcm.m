function Cm = funcm(x,f,T,M,N)
%�������ڼ���ӹ�����
%x           input           ��������
%f           input           װ����װ�䵥λ��Ʒ����ķ���
%T           input           �ƻ�����
%M           input           װ��������
%N           input           ��Ʒ����
%Cm           output         �ӹ����� 
y=0;
for t=1:T
    for j=1:N
        for i=1:M
            y=y+f(i,j)*x(1,(t-1)*M*N+(j-1)*M+i);
            Cm=y;
        end
    end
end


