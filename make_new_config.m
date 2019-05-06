function [A_new]=make_new_config(i,A,Adj)

N=size(Adj,1);

single_node=0;
A_tmp=squeeze(A(:,:,i));
count=0;
for ii=1:size(A,2)
    if sum(A(:,ii,i))<2
        count=count+1;
        single_node(count)=find(A_tmp(:,ii));
    end
end

% check which nodes are connected to other clusters
edge=(Adj-eye(length(Adj)))*squeeze(A_tmp);
edge=(1-A_tmp).*edge;
edge_nodes=mod(find(edge),N);
edge_nodes(edge_nodes==0)=N;
if count~=0
    for jj=single_node
        edge_nodes(edge_nodes==jj)=[];
    end
end

% from this list of nodes, choose one at random and propose a flip
switch_node=datasample(edge_nodes,1);
oldC=find(A(switch_node,:,i));
edge_cluster=find(edge(switch_node,:));
newC=datasample(edge_cluster,1);
A_new(:,:,i)=A(:,:,i);
A_new(switch_node,oldC,i)=0;
A_new(switch_node,newC,i)=1;

end