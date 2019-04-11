function [A]=erdosrenyi_N(N,Wvec)
%
% nodes: number of nodes in the 2 cluster graph
% N1: number of nodes in each cluster 
% W: is the probability of inter/intracluster connections
%

nodes=sum(N);
for i=1:length(N)
    p(i)=N(i)/nodes;
end
x(1:nodes)=2;
x(1:N(1))=1;
for i=2:length(N)
   x((sum(N(1:i-1))+1):sum(N(1:i)))=i; 
end

W=ones(length(N))*Wvec(2);
for i=1:length(N)
   W(i,i)=Wvec(1);
end
A=zeros(nodes);


for i=1:nodes
    for j=i:nodes
        P(i,j) = (4)*W(x(i),x(j))/(nodes*(p(x(i))*p(x(j))));
        r = rand(1);
        if r<P(i,j)
            A(i,j)=1;
        end
    end
end
for i=1:nodes
    for j=1:i
        A(i,j)=A(j,i);
    end
end

% this bit adds symmetry
%A=(A+A')/2;
%A=ceil(A);


rem=[];
for i=1:nodes
    if max(A(i,:))==0 || max(A(:,i))==0
        rem=[rem i];
    end
end
A(rem,:)=[];
A(:,rem)=[];
rem

end