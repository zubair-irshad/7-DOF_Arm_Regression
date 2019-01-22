clear;
clc;

all_curr = load('C:/Users/Areeb Mehmood/Desktop/dataCur_NF.txt');
phi_mat = load('C:/Users/Areeb Mehmood/Desktop/phifile_comp2');
% beta_k contains Mass and COM parameters that we estimated from COM experiments
% Treating these values as known (beta_k)
beta_k = load('C:/Users/Areeb Mehmood/Desktop/betaK.txt'); 
[row, col] = size(all_curr);

% Create holder for all_torques
all_torque = all_curr*0;

W=[];
T=[];

% km = zeros(7,1);
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

% Stack <#datapoints> rows of <7> torques into one <7>*<#datapoints> column
for i = 1:row
   ta=all_torque(i,:)';
   tau=[T;ta];
   T=tau;
end

phi_k = [];
phi_u = [];
% Separate phi_k and phi_u out of phi_mat
for i=1:13:91
   phi_k = [phi_k  phi_mat(:,i) phi_mat(:,i+1) phi_mat(:,i+2) phi_mat(:,i+3)];
end
for i=5:13:91
   phi_u = [phi_u  phi_mat(:,i) phi_mat(:,i+1) phi_mat(:,i+2) phi_mat(:,i+3) phi_mat(:,i+4) phi_mat(:,i+5) phi_mat(:,i+6) phi_mat(:,i+7) phi_mat(:,i+8)];
end

% Subtracting known phi*beta from existing tau to get tau_prime
tau_prime = tau - phi_k*beta_k';

% Running regression on unknown betas
ridge_u = (phi_u'*phi_u);
beta_u = pinv(ridge_u)*phi_u'*tau_prime;

% Putting known and calculated beta parameters in single vector
betaVec = zeros(91,1);
for i=0:6
   betaVec(13*i+1:13*i+13) = [beta_k(4*i+1:4*i+4)'; beta_u(9*i+1:9*i+9)];
end

betaVec = betaVec';
save('C:/Users/Areeb Mehmood/Desktop/betaVec.txt','betaVec','-ascii');

% ridge=(phi_mat'*phi_mat);%+(1e-5*eye(n));
% param=pinv(ridge)*phi_mat'*tau;
% % Matrix format for params 
% paramMat = [param(1:13) param(14:26) param(27:39) param(40:52) param(53:65) param(66:78) param(79:91)]

% std_param is the numerical value of parameters with '1' being M1 '2'
% being 'MX1' and so on...

% std_param=linspace(1,n,n);

% r=rank(phi_mat);

% rr=phi_mat'*phi_mat;
% e=eig(rr);
% min_e=min(e);
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
% 
% ic=zeros(n,1);
% 
% for i =1:n
%     if norm(phi_mat(:,i)) == 0;
%         ic(i)=1;
%     end
% end
% 
% un_ident=find(ic==1);
% 
% phi_mat1=phi_mat;
% 
% phi_mat(:,un_ident)=[];
% 
% std_param(:,un_ident)=[];
% 
% 
% 
% ident=[];
% ind_dep=[];
% 
% [o,p]=size(phi_mat);
% 
% L=zeros(o,p);
% 
% phi_matprior=phi_mat;
% 
% for i = 1:p
%     
%     fci=phi_mat(:,i);
%     
%     phi_mat(:,i)=[];
%     
%     rwi=rank(phi_mat);
%     
%     if rwi <r
%         ident(end+1)=std_param(i);
%         L(:,i)=fci;
%         
%     else
%         ind_dep(end+1)=i;
%         
%     end
%         
%     if i==1
%         
%         B=[fci phi_mat];
%         
%         phi_mat=B;
%         
%     end
%     
%     if i>1 && i<n
%         B=[phi_mat(:,1:i-1) fci phi_mat(:,i:p-1)];
%         
%         phi_mat=B;
%         
%     end
%     
%     if i==n
%         B=[phi_mat fci];
%         phi_mat=B;
%     end
%     
% end
% 
% %%%%
% 
% % ic2=zeros(p,1);
% % 
% % for i =1:p
% %     if norm(L(:,i)) == 0;
% %         ic2(i)=1;
% %     end
% % end        
% % 
% % ind=find(ic2==1);
% 
% 
% %%%%
% 
% L(:,ind_dep)=[];
% 
% DBS=ident;
% 
% rank_L=rank(L);
% 
% [~,q]=size(ind_dep);
% 
% ident_LC=[];
% 
% 
% for i =1:q
%     
%     rank_L=rank(L);
%     
%     fc=phi_mat(:,ind_dep(i));
%     
%     L1=[L fc];
%     
%     rank_L1=rank(L1);
%     
%     if rank_L1>rank_L
%         
%         DBS(end+1)= std_param(ind_dep(i));
%    
%         L(:,end+1)=fc;
%     
%     else
%         
%         ident_LC(end+1)=std_param(ind_dep(i));
%         
%     end
%     
%     
% end
% 
% s_LC=size(ident_LC);
% LC_theta=zeros(s_LC);
% 
% ident_theta=pinv(L'*L)*L'*tau;
% 
% for i = 1:s_LC
%     
%     fc_LC=phi_mat1(:,ident_LC(i));
%     
%     alpha=pinv(L'*L)*L'*fc_LC;
%     
%     LC_theta(i)=ident_theta'*alpha; 
% end
% 
% i_theta=[ident_theta DBS'];
% li_theta=[LC_theta' ident_LC']; 
% 


% alpha = 0.01;
% num_iters = 400;
% 
% theta = zeros(77, 1);
% [theta, J_history] = gradientDescentMulti(phi_base, tau, theta, alpha, num_iters);
% 
