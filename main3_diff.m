% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
close all
% Create graph
% choose a number of edges
% Draw a random pair of states and connect them (if they are not already
% connected) until theres sufficient edges

M = 3; % number of edges within a cluster

[A]=make_graph(M); % make random graph
%[A]=make_lin_graph(M); % make linear chain graph

%if i want a multi state-esque network
multi_state = 3; % number of clusters
ncon = 1; % number of intercluster links
for i = 1:multi_state-1
    [A1]=make_graph(M);
    A2 = blkdiag(A,A1);
    for j=1:ncon
        rand_a = randi(size(A,1),1);
        rand_b = randi([size(A,1)+1,size(A2,1)],1);
        if A2(rand_a,rand_b)==0
            A2(rand_a,rand_b)=1;
            A2(rand_b,rand_a)=1;
        end
    end
    A = A2;
end
N=size(A,1); % total number of nodes in the network

% make a figure of the graph we have generated
G = digraph(A); % directed graph from adjacency matrix generated
h=plot(G);

% create a matrix which describes the rate of transition between states
% evenly divide rate among number of states linked to
for i = 1:N
    for j = 1:N
        K(i,j)=A(i,j)/sum(A(:,j));
    end
end

for i = 1:N
    K(i,i)=-sum(K(i,:).*A(i,:));
end

% do spectral decomposition
[Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K');
kemeny = sum(-1./Keigs(2:end));

% what are the endpoints suggested by the second eigenvector
[a,b1] = min(K_eig_R(:,2));
[a,b2] = max(K_eig_R(:,2));


%%
% in this piece of code I want to compare:
% choosing end states from committor clustering
% followed by diffusive boundaries
tic
kem_max=0;
MM=jjhunter(expm(K));
maximum = max(max(MM));
[x,y]=find(MM==maximum);
end_points=[x,y];

% now for each other state I want to compute the commitor probability to
% reach one state or the other first, details on what i'm doing are here:
% www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
[committor]=compute_commit(K,end_points);

[~,tmp2]=sort(committor);

% Now that I've identified the optimal wells, I will do a diffusive
% boundary search for the optimal clustering
T=1000;
t=0;
% let's set up n_sim simulateneous searches at different temperatures
n_sim=10;
boundary(1,:)=randi([2,N-2],[1,n_sim]);
for i=1:n_sim
    boundary(2,i)=randi([boundary(1,i)+1,N-1],1);
end
% temperatures for sims
for i=1:n_sim
    boundary(:,i)
    [kemeny_latest(i)]=kemeny_boundary(K,eq,boundary(:,i),tmp2);
end
switchcount=0;
prop_count=0;
temp=1:10;
best_kem_yet=0;
while t<T
    t=t+1;
    % for each simulation
    for i=1:n_sim
        % pick a boundary
        bound = randi([1,2]);
        % propose a change to the boundary
        bound_new(:,i) = boundary(:,i);
        bound_new(bound,i) = bound_new(bound,i)+randi([0,1])*2-1;
        try
            [kemenyR_new]=kemeny_boundary(K,eq,bound_new(:,i),tmp2);
        catch
            break
        end
        if kemenyR_new>kemeny_latest(i)
            switchcount=switchcount+1;
            kemeny_latest(i)=kemenyR_new;
            boundary(:,i)=bound_new(:,i);
        elseif kemenyR_new<=kemeny_latest(i)
            val=rand(1);
            condition = exp((kemenyR_new-kemeny_latest(i))*(1/temp(i)));
            if val<condition
                switchcount=switchcount+1;
                kemeny_latest(i)=kemenyR_new;
                boundary(:,i)=bound_new(:,i);
            end
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
                boundary(:,[ii jj])=boundary(:,[jj, ii]);
            else
                val=rand(1);
                condition = exp((kemeny_latest(jj)-kemeny_latest(ii))*((1/temp(ii))-(1/temp(jj))));
                if condition<val
                    %switchcount=switchcount+1;
                    kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                    boundary(:,[ii jj])=boundary(:,[jj, ii]);
                end
            end
        end
    end
    
    % after each timestep, sweep all the states and look for the optimum
    [a,b]=max(kemeny_latest);
    if a>best_kem_yet
        opt_bound=boundary(:,b);
        best_kem_yet=a;
    end
    
end
A=zeros(N,3);
A(tmp2(1:opt_bound(1)),1)=1;
A(tmp2(opt_bound(1)+1:opt_bound(2)),2)=1;
A(tmp2(opt_bound(2)+1:end),3)=1;
best_split=A;
diffusion=toc;
keyboard
% Here I want to do the exhaustive version where I try every combination of
% states to cluster in to.
tic

kem_max2=0;
count=0;
for i=1:(N-2)
    clus1=nchoosek([1:N],i);
    [i]
    for k=1:size(clus1,1)
        clus1_tmp=clus1(k,:);
        for j=1:N-i-1
            clus2=nchoosek(setdiff([1:N],clus1_tmp),j);
            for kk=1:size(clus2,1)
                count=count+1;
                A=zeros(N,3);
                clus2_tmp=clus2(kk,:);
                clus3_tmp=setdiff([1:N],[clus1_tmp,clus2_tmp]);
                A(clus1_tmp,1)=1;
                A(clus2_tmp,2)=1;
                A(clus3_tmp,3)=1;
                [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K', eq, A);
                [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
                kemenyR=sum(-1./Reigs(2:end));
                if kemenyR>kem_max2
                    kem_max2 = kemenyR;
                    best_split_exhaus = Aclus;
                end
            end
        end
    end
end
exhaust_time=toc;

% third of all, I want to do a semi-exhaustive search where I take the 1-D
% ordering and search every 1D clustering
for i1=2:N-2
    for j1=i+1:N-1

        [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K', eq, A);
        
        % analyse eigenvalues vectors of clustered matrix R
        [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
        
        % compute reduced kemeny (for 2 state clustering this is just the
        % same as the slowest timescale
        kemenyR=sum(-1./Reigs(2:end));
        
        %  save the info for the iteration which maximises kemeny (this may
        %  not be unique), will only save the first choice found
        if kemenyR>kem_max
            kem_max = kemenyR;
            best_wells = end_points;
            best_split = Aclus;
        end
    end
end


% colour based on 'best' splitting
figure()
subplot(1,2,1)
h = plot(G);
highlight(h,find(best_split(:,1)),'NodeColor','g')
highlight(h,find(best_split(:,2)),'NodeColor','r')
title('Best splitting','FontSize', 18)

% colour based on exhaustive
%figure()
subplot(1,2,2)
h = plot(G);
highlight(h,find(best_split_exhaus(:,3)),'NodeColor','g')
highlight(h,find(best_split_exhaus(:,1)),'NodeColor','r')
title('Exhaustive splitting','FontSize', 18)

% % colour based on second eigenvector end points
% figure()
% h = plot(G);
% highlight(h,find(sec_eig_split(:,1)),'NodeColor','g')
% highlight(h,find(sec_eig_split(:,2)),'NodeColor','r')
% title('Second eigenvector end point splitting','FontSize', 18)