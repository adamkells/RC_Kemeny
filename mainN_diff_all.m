% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
close all
potential_type=1;
% each element of the vector is the number of nodes within each cluster
if potential_type==0
    nodes=[5,5]; % length of vector defines number of clusters
    [K,Adj]=erdosrenyi_N(nodes,[0.7,0.1]); % this creates a random erdos renyi graph
elseif potential_type==1
    [K,Adj,v]=szabo(15);
end
K=K';
% the first input is the vector of nodes from before, the second input is
% a 2 element vector with i) probability of node connection within cluster
% and ii) probability of node connection between clusters
% the output is the adjency matrix of the graph
N=size(K,1); % N is the total number of nodes
N2=sqrt(N);
% do spectral decomposition of the rate matrix of the system
[Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K);
kemeny = sum(-1./Keigs(2:end)); % kemeny constant of system
slow_rels = -1./Keigs(2:end); % relaxation processes

G = digraph(Adj); % directed graph from adjacency matrix generated
h=plot(G);
%K=K';
%% CLUSTERING %%%%%%%%%%%%%%%

% In this section, the aim is to combine parallel tempering with
% variational dynamics to find the optimal clustering of the nodes in the
% network.

% these are just some quantities I'll need to use within a big for loop
% later so I'll do them here to save on computation
one_vec=ones(1,length(K));
INV_K=(inv(eq*one_vec-K));

NCLUS=4; %number of clusters to look for

for iii=1:NCLUS
    color_scheme(iii,:)=rand(1,3);
end

T=N*100; % the total amount of time that I will simulate for

% next set up n_sim simultaneous searches at different temperatures
n_sim=50;

% choice of reduction method, 0 for Hummer-Szabo, 1 for local equilibrium
counter=0;
px=0;
bloop=0;
for red_method=0:1
    display(num2str(red_method))
    % choice of variational parameter, 0 for kemeny, 1 for tau_2, 2 for kemeny-1
    for param=0:2
        counter=counter+1;
        % The first step to try and find a 1D ordering of my nodes. To do this I
        % find the two most kinetically distinct states and then order the
        % remaining states based on their committor probability to reach one state
        % over the other.
        
        % There are two ways to choose my end points: either the pair of states
        % with the largest mean first passage time between them (MFPT) or the most
        % negative and positive elements of the second eigenvector
        [end_points]=ep_choice(K,K_eig_R);
        
        % now for each other state I want to compute the commitor probability to
        % reach one state or the other first. Details on what i'm doing are here:
        % www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
        [committor]=compute_commit(K',end_points); % find committor for each state
        [~,tmp2]=sort(committor); % order all nodes from the committor

        % My way of choosing the initial grouping of states is pretty arbitrary and
        % maybe not the best way. Something to come back to. I just pick some
        % random boundaries along the initial 1D ordering. I let each simulation
        % start at a different configuration
        % This code picks some random boundaries
        boundary(1,:)=randi([2,floor(N/(NCLUS-1))-1],[1,n_sim]);
        for i=1:n_sim
            for j=2:(NCLUS-1)
                boundary(j,i)=randi([boundary((j-1),i)+2,(floor(j*N/(NCLUS-1))-1)],1);
            end
        end
        
        % This code creates a clustering matrix A from the boundaries which
        % specifies which cluster each node belongs to
        A=zeros([N,NCLUS,n_sim]);
        for i=1:n_sim
            A(tmp2(1:boundary(1,i)),1,i)=1;
            for j=2:(NCLUS-1)
                A(tmp2(boundary((j-1),i)+1:boundary(j,i)),j,i)=1;
            end
            A(tmp2(boundary(end,i)+1:end),NCLUS,i)=1;
        end
        %keyboard
        % calculate the kemeny value of each starting configuration
        for i=1:n_sim
            [kemeny_latest(i)]=kemeny_boundary(K,INV_K,eq,A(:,:,i),red_method,param);
        end

        % for the parallel tempering, I need to choose temperatures for sims
        % as temperature has no physical meaning here, the temperatures will be
        % changed over the course of the simulation to maintain an acceptance ratio
        % of 50%. So I can pretty much choose whatever numbers here.
        temp=linspace(1,10,n_sim);
        temp=temp*0.01;
        % count of how many proposed moves are accepted, used to calculate the
        % neccesary change to the temperature
        switchcount=0;
        attswitch=0;
        switchcount2(1:n_sim)=0;
        switchcount3(1:n_sim)=0;
        bad_suggest(1:n_sim)=0;
        % optimal value of the variational parameter so far (largest of starting values)
        best_kem_yet=max(kemeny_latest);
        % preset the modularity value
        Q=0;
        % I might in a later version, make a movie showing the convergence of the
        % clusters and this will be used to index the frames of the movie
        framecount=0;
        % starting time
        t=0;
        while t<T
            % counter to display progress
            if mod(t,100)==0
                %keyboard
                display(['Completion is:', num2str(100*t/T), '%']);
                display(['Best value is:', num2str(best_kem_yet)]);
                display(['Modularity is:', num2str(Q)]);
                display(['Acceptance is:', num2str(px)]);
                display(['Bloop is:', num2str(bloop)]);
            end
            t=t+1;
            for i=1:n_sim
                % this bit of code is a bit dense (and probably suboptimally
                % written) basically it looks at which nodes are connected to nodes
                % of a different cluster. It then randomly chooses one of these
                % nodes to propose a flip. It then creates a new configuration
                % A_new.
                [A_new]=make_new_config(i,A,Adj);

                % here we calculate the value of the variational parameter for the
                % new configuration and if it improves the value we accept it, if
                % it makes the value worse we accept with some probability
                % this approach is inspired by parallel tempering
                try
                    [kemenyR_new]=kemeny_boundary(K,INV_K,eq,A_new(:,:,i),red_method,param);
                    if kemenyR_new>kemeny_latest(i)
                        switchcount=switchcount+1;
                        switchcount2(i)=switchcount2(i)+1;
                        kemeny_latest(i)=kemenyR_new;
                        A(:,:,i)=A_new(:,:,i);
                    elseif kemenyR_new<=kemeny_latest(i)
                        attswitch=attswitch+1;
                        val=rand(1);
                        condition = exp((kemenyR_new-kemeny_latest(i))*(1/temp(i)));
                        if condition==0
                            bad_suggest(i)=bad_suggest(i)+1;
                        end
                        if val<condition
                            switchcount=switchcount+1;
                            switchcount3(i)=switchcount3(i)+1;
                            kemeny_latest(i)=kemenyR_new;
                            A(:,:,i)=A_new(:,:,i);
                        end
                    end
                catch
                    bloop=bloop+1;
                    keyboard
                end
            end
            
            % every 10 steps attempt some interchanges between the configurations
            % at different temperatures
            % again, using the parallel tempering approach, we accept these swaps
            % with some probability
            if mod(t,10)==0
                for j=1:(n_sim-1)
                    ii=j+1;
                    jj=j;
                    %keyboard
                    if kemeny_latest(ii)>kemeny_latest(jj)
                        kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                        A(:,:,[ii jj])=A(:,:,[jj, ii]);
                    else
                        val=rand(1);
                        condition = exp((kemeny_latest(jj)-kemeny_latest(ii))*((1/temp(jj))-(1/temp(ii))));
                        if condition<val
                            kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                            A(:,:,[ii jj])=A(:,:,[jj, ii]);
                        end
                    end
                end
            end
            %keyboard
            % next I update temperature during the simulation to keep 50% acceptance
            if mod(t,200)==0
                px=sum(switchcount3)/(attswitch); % currect accepted fraction since last update
                %keyboard
                attswitch=0;
                switchcount3(1:n_sim)=0; % reset the counter of accepted configurations
                correction=log(px)/log(0.5);
                temp=temp*correction;
            end
            
%             if t==100
%                 % This code creates a clustering matrix A from the boundaries which
%                 % specifies which cluster each node belongs to
%                 A=zeros([N,NCLUS,n_sim]);
%                 for iii=1:n_sim
%                     A(tmp2(1:boundary(1,iii)),1,iii)=1;
%                     for jjj=2:(NCLUS-1)
%                         A(tmp2(boundary((jjj-1),iii)+1:boundary(jjj,iii)),jjj,iii)=1;
%                     end
%                     A(tmp2(boundary(end,iii)+1:end),NCLUS,iii)=1;
%                 end
%                 
%                 % calculate the kemeny value of each starting configuration
%                 for iii=1:n_sim
%                     [kemeny_latest(iii)]=kemeny_boundary(K,INV_K,eq,A(:,:,iii),red_method,param);
%                 end
%             end
            
            % after each timestep, sweep all the states and look for the optimum
            [a,b]=max(kemeny_latest);
            % if the optimum is better than the previous best, update
            if a>best_kem_yet
                best_split=A(:,:,b);
                best_kem_yet=a;
                %framecount=framecount+1;
                %frames(:,:,framecount) = opt_bound;
                
                % compute x from A
                for iii=1:size(A,1)
                    x(iii)=find(best_split(iii,:));
                end
                
                % Computation of the modularity of the best clustering
                m = sum(sum(Adj));
                Q = 0;
                COMu = unique(x);
                for jj=1:length(COMu)
                    Cj = find(x==COMu(jj));
                    Ec = sum(sum(Adj(Cj,Cj)));
                    Et = sum(sum(Adj(Cj,:)));
                    if Et>0
                        Q = Q + Ec/m-(Et/m)^2;
                    end
                end
            end
        end
        
        
        if potential_type==0
            figure(counter+1)
            h = plot(G);
            title('Best splitting','FontSize', 18)
            txt = ['Parameter: ' num2str(best_kem_yet)];
            text(-4,4,txt)
            txt = ['Modularity: ' num2str(Q)];
            text(-4,3,txt)
            saveas(gcf,[num2str(counter) '_figure.png'],'png')
            cc(:,:,counter)=best_split;
        elseif potential_type==1
            cluster_V=zeros(N2);
            for jjj=1:NCLUS
                [a,b]=find(best_split(:,jjj));
                new_ind1=mod(a,N2);
                new_ind1(new_ind1==0)=N2;
                new_ind2=1+(a-new_ind1)/N2;
                for iii=1:length(a)
                    cluster_V(new_ind1(iii),new_ind2(iii))=jjj;
                end
            end
            cc(:,:,counter)=best_split;
            dd(:,:,counter)=cluster_V;
        end
    end
end

%%
conv_check()