function [K,Adj,v]=szabo(N)

x=linspace(-0.5*pi,0.5*pi,N);
y=linspace(-0.5*pi,0.5*pi,N);
e=exp(1);
for i = 1:length(x)
    for j = 1:length(y)
        v(i,j)=-10*(e.^(-2*(x(i)+1).^2-2.*(y(j)-1).^2)+e.^(-2.*(x(i)+0.8).^2-2.*(y(j)+1).^2)+e.^(-2.*(x(i)-1).^2-2.*(y(j)+0.5).^2));
    end
end
v=v-min(min(v));
v=v;
% want to create a 1D rate matrix from this potential energy
% have to create a mapping from 2D states in to a single index
K=zeros(N^2);
A=1;
kbT=0.596;
for j=1:length(y)
    for i=1:length(x)
        new_index=(j-1)*N+i;
        up=(j-1)*N+i-1;
        down=(j-1)*N+i+1;
        left=(j-2)*N+i;
        right=(j)*N+i;
        if (i-1)>0
            K(new_index,up)=A*exp((v(i,j)-v(i-1,j))/(2*kbT));
            K(up,new_index)=A*exp(-(v(i,j)-v(i-1,j))/(2*kbT));
            %keyboard
        end
        if (i+1)<(N+1)
            K(new_index,down)=A*exp((v(i,j)-v(i+1,j))/(2*kbT));
            K(down,new_index)=A*exp(-(v(i,j)-v(i+1,j))/(2*kbT));
            %keyboard
        end
        if (j-1)>0
            K(new_index,left)=A*exp((v(i,j)-v(i,j-1))/(2*kbT));
            K(left,new_index)=A*exp(-(v(i,j)-v(i,j-1))/(2*kbT));
            %keyboard
        end
        if (j+1)<(N+1)
            K(new_index,right)=A*exp((v(i,j)-v(i,j+1))/(2*kbT));
            K(right,new_index)=A*exp(-(v(i,j)-v(i,j+1))/(2*kbT));
            %keyboard
        end
    end
end
K=K';
for i=1:size(K,1)
    K(i,i)=-sum(K(:,i));
end
K=K';

Adj=K./K;
Adj(isnan(Adj))=0;

figure()
contour(v)