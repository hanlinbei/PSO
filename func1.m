function C1 = func1(x,R,V,T,N,M)
%�������ڼ���ͷ�����
%x           input           ��������
%R           input           ��t����Ҫ����j���Ʒ����
%T           input           �ƻ�����
%M           input           װ��������
%N           input           ��Ʒ����
%V           input           �ͷ�ϵ��
%C1           output         �ͷ����� 
%% ����t��ʱj�ֲ�Ʒ��������
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
%% ���ڷ���
y=0;
for t=1:T
    for j=1:N
            y=y+V(1,j)*max(R(j,t)-Q(j,t),0);
            C1=y;
    end
end
