% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
close all

%%%% I've got code that does the diffusive clustering but it's spending a
%%%% lot of its time computing kemeny for clusters with 0 microstates. Need
%%%% to make sure that if a cluster has only 1 microstate then that cannot
%%%% be removed.

% Create graph
% choose a number of edges
% Draw a random pair of states and connect them (if they are not already
% connected) until theres sufficient edges

% M = 8; % number of edges within a cluster
% 
% [A]=make_graph(M); % make random graph
% %[A]=make_lin_graph(M); % make linear chain graph
% 
% %if i want a multi state-esque network
% multi_state = 2; % number of clusters
% ncon = 1; % number of intercluster links
% for i = 1:multi_state-1
%     [A1]=make_graph(M);
%     A2 = blkdiag(A,A1);
%     for j=1:ncon
%         rand_a = randi(size(A,1),1);
%         rand_b = randi([size(A,1)+1,size(A2,1)],1);
%         if A2(rand_a,rand_b)==0
%             A2(rand_a,rand_b)=1;
%             A2(rand_b,rand_a)=1;
%         end
%     end
%     A = A2;
% end
% N=size(A,1); % total number of nodes in the network


nodes=[10,10,10,10,10,10];
A=erdosrenyi_N(nodes,[0.8,0.025]);
N=length(A);
% make a figure of the graph we have generated
G = digraph(A); % directed graph from adjacency matrix generated
h=plot(G);
keyboard
% create a matrix which describes the rate of transition between states
% evenly divide rate among number of states linked to
for i = 1:N
    for j = 1:N
        K(i,j)=A(i,j)/sum(A(:,j));
    end
end

for i = 1:N
    K(i,i)=0;
    K(i,i)=-sum(K(:,i));%.*A(i,:));
end
Adj=A;
% do spectral decomposition
[Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K');
kemeny = sum(-1./Keigs(2:end));

% what are the endpoints suggested by the second eigenvector
[a,b1] = min(K_eig_R(:,2));
[a,b2] = max(K_eig_R(:,2));

one_vec=ones(1,length(K));
INV_K=(inv(eq*one_vec-K'));


%%
% in this piece of code I want to compare:
% choosing end states from committor clustering
% followed by diffusive boundaries
tic
kem_max=0;
secvec=1;
if secvec==0
    MM=jjhunter(expm(K));
    maximum = max(max(MM));
    [x,y]=find(MM==maximum);
    end_points=[x,y];
elseif secvec==1
    end_points=[b1,b2];
end


% now for each other state I want to compute the commitor probability to
% reach one state or the other first, details on what i'm doing are here:
% www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
[committor]=compute_commit(K',end_points);

[~,tmp2]=sort(committor);

% Now that I've identified the optimal wells, I will do a diffusive
% boundary search for the optimal clustering
T=N*100;
t=0;
% let's set up n_sim simulateneous searches at different temperatures
n_sim=50;

% My way of choosing the initial grouping is pretty arbitrary and maybe not
% the best way. Something to come back to.
NCLUS=6; %number of clusters to use
boundary(1,:)=randi([2,floor(N/(NCLUS-1))],[1,n_sim]);
for i=1:n_sim
    for j=2:(NCLUS-1)
        boundary(j,i)=randi([boundary((j-1),i)+2,(floor(j*N/(NCLUS-1))-1)],1);
    end
end

A=zeros([N,3,n_sim]);
for i=1:n_sim
    A(tmp2(1:boundary(1,i)),1,i)=1;
    for j=2:(NCLUS-1)
        A(tmp2(boundary((j-1),i)+1:boundary(j,i)),j,i)=1;
    end
    A(tmp2(boundary(end,i)+1:end),NCLUS,i)=1;
end
% temperatures for sims
for i=1:n_sim
    i
    [kemeny_latest(i)]=kemeny_boundary(K,INV_K,eq,A(:,:,i));
end

switchcount=0;
prop_count=0;
temp=linspace(0.025,1,n_sim);
% I could calculate the acceptance ration on the fly and adjust the
% temperature as I gok

best_kem_yet=0;
%keyboard
while t<T
    t=t+1;
    tic
    % for each simulation
    for i=1:n_sim
        % pick a connection between differing clusters
        % find two nodes from different clusters and propose flipping one
        % of the nodes
        A_tmp=squeeze(A(:,:,i));
        count=0;
        for ii=1:size(A,2)
            if sum(A(:,ii,i))<2
                count=count+1;
                single_node(count)=find(A_tmp(:,ii));
            end
        end
        edge=(Adj-eye(length(Adj)))*squeeze(A_tmp);
        edge=(1-A_tmp).*edge;
        edge_nodes=mod(find(edge),N);
        edge_nodes(edge_nodes==0)=N;
        if count~=0
            for jj=single_node
                edge_nodes(edge_nodes==jj)=[];
            end
        end
        switch_node=datasample(edge_nodes,1); % choose a node that is
        oldC=find(A(switch_node,:,i));
        % connected to a different cluster
        edge_cluster=find(edge(switch_node,:)); % check which clusters it connected to
        newC=datasample(edge_cluster,1); % pick a cluster
        %switch_node2=Adj(switch_node,:).*A(:,switch_cluster,i); 
        %switch_node2=datasample(switch_node2,1);
        % second node to do the switch with
        %keyboard
        % propose a flip to the boundary
        A_new(:,:,i)=A(:,:,i);
        A_new(switch_node,oldC,i)=0;
        A_new(switch_node,newC,i)=1;
        try
            %A_new(:,:,i);
            [kemenyR_new]=kemeny_boundary(K,INV_K,eq,A_new(:,:,i));
            if kemenyR_new>kemeny_latest(i)
                switchcount=switchcount+1;
                kemeny_latest(i)=kemenyR_new;
                A(:,:,i)=A_new(:,:,i);
            elseif kemenyR_new<=kemeny_latest(i)
                val=rand(1);
                condition = exp((kemenyR_new-kemeny_latest(i))*(1/temp(i)));
                %keyboard
                if val<condition
                    switchcount=switchcount+1;
                    kemeny_latest(i)=kemenyR_new;
                    A(:,:,i)=A_new(:,:,i);
                end
            end
        catch
            kemenyR_newX = 0;
        end
    end
    % after a certain number of steps attempt some changes
    if mod(t,10)==0
        %prop_count=prop_count+1;
        for j=1:n_sim
            propose=randi([1,n_sim],2);
            ii=propose(1);
            jj=propose(2);
            if kemeny_latest(ii)>kemeny_latest(jj)
                %switchcount=switchcount+1;
                kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                A(:,:,[ii jj])=A(:,:,[jj, ii]);
            else
                val=rand(1);
                condition = exp((kemeny_latest(jj)-kemeny_latest(ii))*((1/temp(jj))-(1/temp(ii))));
                if condition<val
                    %switchcount=switchcount+1;
                    kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                    A(:,:,[ii jj])=A(:,:,[jj, ii]);
                end
            end
        end
    end
    
    % I'm going to update the temperature during the simulation to try and
    % keep the 50% acceptance
    if mod(t,10)==0
       t
       px=switchcount/(n_sim*10); % currect accepted fraction since last update
       switchcount=0;
       correction=log(px)/log(0.5);
       temp=temp*correction;
    end
    
    % after each timestep, sweep all the states and look for the optimum
    [a,b]=max(kemeny_latest);
    if a>best_kem_yet
        opt_bound=A(:,:,b);
        best_kem_yet=a;
    end
    timetaken=toc;
end
best_split=opt_bound;
diffusion=toc;
keyboard
%%
% colour based on 'best' splitting
figure()
subplot(1,2,1)
h = plot(G);
for i=1:NCLUS
    highlight(h,find(best_split(:,i)),'NodeColor',[rand(1,3)])
end
title('Best splitting','FontSize', 18)
%
% colour based on 'best' splitting
% for i=1:10
%     figure(i+5)
%     best_split=A(:,:,i);
%     h = plot(G);
%     highlight(h,find(best_split(:,1)),'NodeColor','g')
%     highlight(h,find(best_split(:,2)),'NodeColor','r')
%     highlight(h,find(best_split(:,3)),'NodeColor','k')
%     title('Best splitting','FontSize', 18)
% end
%
% Here I want to do the exhaustive version where I try every combination of
% states to cluster in to.
%tic

% kem_max2=0;
% count=0;
% for i=1:(N-2)
%     clus1=nchoosek([1:N],i);
%     [i]
%     tic
%     for k=1:size(clus1,1)
%         clus1_tmp=clus1(k,:);
%         for j=1:N-i-1
%             clus2=nchoosek(setdiff([1:N],clus1_tmp),j);
%             for kk=1:size(clus2,1)
%                 count=count+1;
%                 A=zeros(N,3);
%                 clus2_tmp=clus2(kk,:);
%                 clus3_tmp=setdiff([1:N],[clus1_tmp,clus2_tmp]);
%                 A(clus1_tmp,1)=1;
%                 A(clus2_tmp,2)=1;
%                 A(clus3_tmp,3)=1;
%                 [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K',INV_K, eq, A);
%                 [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
%                 kemenyR=sum(-1./Reigs(2));%:end));
%                 if kemenyR>kem_max2
%                     kem_max2 = kemenyR;
%                     best_split_exhaus = Aclus;
%                 end
%             end
%         end
%     end
%     time = toc
% end
%exhaust_time=toc;

% third of all, I want to do a semi-exhaustive search where I take the 1-D
% ordering and search every 1D clustering

% This bit needs to be adapted to do general number of clusters
kem_max=0;
for i1=2:N-2
    for j1=i+1:N-1
        [i1,j1]
        A=zeros(N,3);
        A(tmp2(1:i1),1)=1;
        A(tmp2(i1+1:j1),2)=1;
        A(tmp2(j1+1:end),3)=1;
        try
            [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K',INV_K, eq, A);
            [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
            kemenyR=sum(-1./Reigs(2:end));
            kemenyR=sum(-1./Reigs(2));
        catch
            kemenyR=0;
        end
        % analyse eigenvalues vectors of clustered matrix R
        
        % compute reduced kemeny (for 2 state clustering this is just the
        % same as the slowest timescale
        %:end));
        
        %  save the info for the iteration which maximises kemeny (this may
        %  not be unique), will only save the first choice found

        if kemenyR>kem_max
            kem_max = kemenyR;
            best_wells = end_points;
            best_split_exhaus = Aclus;
            best_bound=[i1,j1];
        end
    end
end

% colour based on exhaustive
%figure()
subplot(1,2,2)
h = plot(G);
for i=1:NCLUS
    highlight(h,find(best_split_exhaus(:,i)),'NodeColor',[rand(1,3)])
end
title('Exhaustive splitting','FontSize', 18)

% % colour based on second eigenvector end points
% figure()
% h = plot(G);
% highlight(h,find(sec_eig_split(:,1)),'NodeColor','g')
% highlight(h,find(sec_eig_split(:,2)),'NodeColor','r')
% title('Second eigenvector end point splitting','FontSize', 18)
