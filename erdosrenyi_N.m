function [K,Adj]=erdosrenyi_N(N,Wvec)
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
Adj=zeros(nodes);


for i=1:nodes
    for j=i:nodes
        P(i,j) = (4)*W(x(i),x(j))/(nodes*(p(x(i))*p(x(j))));
        r = rand(1);
        if r<P(i,j)
            Adj(i,j)=1;
        end
    end
end
for i=1:nodes
    for j=1:i
        Adj(i,j)=Adj(j,i);
    end
end

% this bit adds symmetry
%A=(A+A')/2;
%A=ceil(A);


rem=[];
for i=1:nodes
    if max(Adj(i,:))==0 || max(Adj(:,i))==0
        rem=[rem i];
    end
end
Adj(rem,:)=[];
Adj(:,rem)=[];

% make a figure of the graph we have generated
%G = digraph(Adj); % directed graph from adjacency matrix generated
%h=plot(G);

% create a matrix which describes the rate of transition between states
% evenly divide rate among number of states linked to

for i = 1:nodes
    for j = 1:nodes
        K(i,j)=Adj(i,j)/sum(Adj(:,j));
    end
end
for i = 1:nodes
    K(i,i)=0;
    K(i,i)=-sum(K(i,:));
end

end