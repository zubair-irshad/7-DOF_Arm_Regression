clear
clc
% Load data (q,qdot,qdotdot,torque) available online and save it
A=load('sarcos2.mat');
B=A.sar;
[row, col] = size(B);

% Read in all data points, but only want to consider first 1000
row = 1000;
all_curr=B(1:row,22:28);

% Create holder for all_torques
all_torque = all_curr*0;

% W=[];
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

phi_mat = load('C:/Users/Areeb Mehmood/Desktop/phifile_comp');

% Reduce amount of data points from 5429 to 1000
% all_torque(row+1:5429,:) = [];
%Reduce size of phi_mat
phi_mat((row+1)*7:5429*7,:) = [];
[m,n]=size(phi_mat);

% Stack <#datapoints> rows of <7> torques into one <7>*<#datapoints> column
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


% std_param is the numerical value of parameters with '1' being M1 '2'
% being 'MX1' and so on...

std_param=linspace(1,n,n);

r=rank(phi_mat);

ic=zeros(n,1);

for i =1:n
    if norm(phi_mat(:,i)) == 0
        ic(i)=1;
    end
end

un_ident=find(ic==1);

phi_mat(:,un_ident)=[];

std_param(:,un_ident)=[];

%%%

[Q,R]=qr(phi_mat);
k=diag(R);

[o,p]=size(phi_mat);

indices=zeros(p,1);

% count number of zeros in indices (independent parameters)

theta_ind=[];
theta_oth=[];

for i =1:size(k)
    tao=1000*eps*abs(max(k));
    if abs(k(i))<=tao
        indices(i)=1;
        theta_oth(end+1)=i;
    else
        theta_ind(end+1)=i;
    end    
end

n_ind=find(indices==0);

n_dep=find(indices==1);

R1=zeros(m,size(n_ind,1));

R2=zeros(m,size(n_dep,1));

for i=1:size(n_ind)
    R1(:,i)=phi_mat(:,n_ind(i));
end

for i=1:size(n_dep)
    R2(:,i)=phi_mat(:,n_dep(i));
end


WP=[R1 R2];

[Q_new,R_new]=qr(WP);

R1_new=R_new(1:r,1:r);

R2_new=R_new(1:r,r+1:p);

ABC=inv(R1_new)*R2_new;