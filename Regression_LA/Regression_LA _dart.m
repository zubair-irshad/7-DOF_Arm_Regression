
% Load data (q,qdot,qdotdot,torque) available online and save it
% A=load('sarcos_inv.mat');
% 
% B=A.sarcos_inv;
% 
% all_q=B(1:1000,1:7);
% 
% all_dq=B(1:1000,8:14);
% 
% all_dqq=B(1:1000,15:21);
% 
% all_torque=B(1:1000,22:28);

%W=[];

%T=[];

phi_mat=load('phi.txt');

all_torque=load('dataTorque_RHS.txt');

[m,n]=size(phi_mat);

[o,p]=size(all_torque);


T=[];

for i = 1:o
    ta=all_torque(i,:)';
    tau=[T;ta];
    T=tau;
end


% Load data from simulation

% all_q=load('dataQ_train.txt');
% all_dq=load('dataQdot_train.txt');
% all_dqq=load('dataQdotdot_train.txt');
% all_torque=load('dataTorque_train.txt');
% 
% 
% % Unwrap q,dq,dqq,all_torque into understandable parameters
% 
% for i = 1:2000
%    q=all_q(i,:);
%    dq=all_dq(i,:);
%    dqq=all_dqq(i,:);
%    
%    ph=PHI(q,dq,dqq);
%    phi_mat=[W;ph];
%    
%    W=phi_mat;
%    
%    
%    ta=all_torque(i,:)';
%    
%    tau=[T;ta];
%    
%    T=tau;
%    
% end

[m,n]=size(phi_mat);


% std_param is the numerical value of parameters with '1' being M1 '2'
% being 'MX1' and so on...

std_param=linspace(1,n,n);

r=rank(phi_mat);
%%%

% [Q,R]=qr(phi_mat);
% 
% k=diag(R);
% 
% indices=zeros(n,1);
% 
% % count number of zeros in indices (independent parameters)
% 
% nindependant=sum(indices==0);
% 
% ndependant=size(k)-size(indices);
% 
% for i =1:size(k)
%     tao=1000*eps*abs(max(k));
%     if abs(k(i))<=tao
%         indices(i)=1;
%     end    
% end

%%%

ic=zeros(n,1);

for i =1:n
    if norm(phi_mat(:,i)) == 0;
        ic(i)=1;
    end
end

un_ident=find(ic==1);

phi_mat1=phi_mat;

phi_mat(:,un_ident)=[];

std_param(:,un_ident)=[];



ident=[];
ind_dep=[];

[o,p]=size(phi_mat);

L=zeros(o,p);

phi_matprior=phi_mat;

for i = 1:p
    
    fci=phi_mat(:,i);
    
    phi_mat(:,i)=[];
    
    rwi=rank(phi_mat);
    
    if rwi <r
        ident(end+1)=std_param(i);
        L(:,i)=fci;
        
    else
        ind_dep(end+1)=i;
        
    end
        
    if i==1
        
        B=[fci phi_mat];
        
        phi_mat=B;
        
    end
    
    if i>1 && i<n
        B=[phi_mat(:,1:i-1) fci phi_mat(:,i:p-1)];
        
        phi_mat=B;
        
    end
    
    if i==n
        B=[phi_mat fci];
        phi_mat=B;
    end
    
end

%%%%

% ic2=zeros(p,1);
% 
% for i =1:p
%     if norm(L(:,i)) == 0;
%         ic2(i)=1;
%     end
% end        
% 
% ind=find(ic2==1);


%%%%

L(:,ind_dep)=[];

DBS=ident;

rank_L=rank(L);

[~,q]=size(ind_dep);

ident_LC=[];


for i =1:q
    
    rank_L=rank(L);
    
    fc=phi_mat(:,ind_dep(i));
    
    L1=[L fc];
    
    rank_L1=rank(L1);
    
    if rank_L1>rank_L
        
        DBS(end+1)= std_param(ind_dep(i));
   
        L(:,end+1)=fc;
    
    else
        
        ident_LC(end+1)=std_param(ind_dep(i));
        
    end
    
    
end

s_LC=size(ident_LC);
LC_theta=zeros(s_LC);

ident_theta=pinv(L'*L)*L'*tau;

for i = 1:s_LC
    
    fc_LC=phi_mat1(:,ident_LC(i));
    
    alpha=pinv(L'*L)*L'*fc_LC;
    
    LC_theta(i)=ident_theta'*alpha; 
end

i_theta=[ident_theta DBS'];
li_theta=[LC_theta' ident_LC']; 



% alpha = 0.01;
% num_iters = 400;
% 
% theta = zeros(77, 1);
% [theta, J_history] = gradientDescentMulti(phi_base, tau, theta, alpha, num_iters);
% 
