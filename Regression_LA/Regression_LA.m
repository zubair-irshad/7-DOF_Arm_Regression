clear
clc

% Load data
A=load('sarcos2.mat');
B=A.sar;
[row, col] = size(B);

phi_mat = load('C:/Users/Areeb Mehmood/Desktop/phifile_comp');
all_curr=B(1:row,22:28);
% Create holder for all_torques
all_torque = all_curr*0;

W=[];
T=[];

km(1) = 31.4e-3;
km(2) = 31.4e-3;
km(3) = 38e-3;
km(4) = 38e-3;
km(5) = 16e-3;
km(6) = 16e-3;
km(7) = 16e-3;

G_R(1) = 596;
G_R(2) = 596;
G_R(3) = 625;
G_R(4) = 625;
G_R(5) = 552;
G_R(6) = 552;
G_R(7) = 552;

% Convert currents to torques
for i = 1:row
        all_torque(i,:) = all_curr(i,:).*km.*G_R;
end

for i = 1:row
%    q=all_q(i,:);
%    dq=all_dq(i,:);
%    dqq=all_dqq(i,:);
%    
%    ph=PHI(q,dq,dqq);
%    phi_mat=[W;ph];
%    
%    W=phi_mat;
   
   
   ta=all_torque(i,:)';
   
   tau=[T;ta];
   
   T=tau;
   
end

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
