% Neural network for Wiener-Hopf kernel factorization
% Author: Liang Sicong
% Date: November 13, 2023
% Reference: @article{hurd1981scattering, title={Scattering by hard and soft parallel half-planes}, 
% author={Hurd, RA and L{"u}neburg, E}, journal={Canadian Journal of Physics}, 
% volume={59}, number={12}, pages={1879--1885}, year={1981}, publisher={NRC Research Press Ottawa, Canada} }

%% Parameters
clear all;close all;
N1=10000; %% Integral path sampling size
alpha_c=zeros(N1,1); %% Integral path
P=zeros(2,N1); %% P(\alpha)
U=zeros(2,N1); %% U(\alpha), analytic in the upper half-plane
L=zeros(2,N1); %% L(\alpha), analytic in the lower half-plane
A1=zeros(1,N1); %% A1(\alpha),
A2=zeros(1,N1); %% A2(\alpha),
A3=zeros(1,N1); %% A3(\alpha),
B2=zeros(1,N1); %% B2(\alpha),
k=5+0.01*1i;    %% normalized wavenumber
h=2.0;d=h/2;    %% distance between parallel plates

theta0=pi/3;   %% Incident angle
k0=k*cos(theta0); %% k0
beta0=k*sin(theta0);  %% beta_0
psi0=exp(1i*beta0*d);  %% psi_0
c1=beta0/(2*pi*1i*psi0) %% constant of P(\alpha)
c2=-psi0/(2*pi*1i)      %% constant of P(\alpha)
%% Integral path
cc3=imag(k);
x_circle_plus=zeros(1,1);
y_circle_plus=zeros(1,1);
x_circle_neg=zeros(1,1);
y_circle_neg=zeros(1,1);
rr=imag(k);
%%% Caucy integral path near the poles
for ii=1:2000
    rr1=rr/2;
    tk1=ii*2*pi/2000;
    x_circle_plus(ii)=rr1*cos(tk1)+real(k);
    y_circle_plus(ii)=rr1*sin(tk1)+imag(k);
    x_circle_neg(ii)=rr1*cos(tk1)+real(-k);
    y_circle_neg(ii)=rr1*sin(tk1)+imag(-k);
end
%%% Integral path on the real axis
dh_zone1=0.05;
dh_zone2=0.01;
dh_zone3=0.0001;
alpha_2_zone1=-50:dh_zone1:-20-dh_zone1;
alpha_2_zone2=-20:dh_zone2:-real(k0)-0.1213-dh_zone2;
alpha_2_zone3=-real(k0)-0.1213:dh_zone3:-real(k0)+0.1213-dh_zone3;
alpha_2_zone4=-real(k0)+0.1213:dh_zone2:0-dh_zone2;
alpha2_real_left=[alpha_2_zone1,alpha_2_zone2,alpha_2_zone3,alpha_2_zone4];
alpha2_real_right=-flip(alpha2_real_left);
alpha_2=[alpha2_real_left,alpha2_real_right];
rr=abs(imag(k));
left_index=rr/(dh_zone3)-1;% -50:50
circle_index=size(alpha_2_zone1,2)+size(alpha_2_zone2,2)+size(alpha_2_zone3,2)/2+1;
for ii=circle_index-left_index:circle_index+left_index
    len=abs(alpha_2(ii)-real(-k0));

    imagx=sqrt(rr^2-len^2);
    alpha_2(ii)=alpha_2(ii)-1i*imagx;

end
%%% Put the Caucy integral and Wiener-Hopf integral path together
x_circle_plus=[x_circle_plus,real(alpha_2)];
y_circle_plus=[y_circle_plus,imag(alpha_2)];
x_circle_neg=[x_circle_neg,real(alpha_2)];
y_circle_neg=[y_circle_neg,imag(alpha_2)];
x_circle_plus=transpose(x_circle_plus);
y_circle_plus=transpose(y_circle_plus);
x_circle_neg=transpose(x_circle_neg);
y_circle_neg=transpose(y_circle_neg);
%%% 1:2000: Caucy integral path, 2001:12000 Wiener-Hopf integral path
save(['PQVR_V8.mat'],'x_circle_plus','y_circle_plus','x_circle_neg','y_circle_neg'); 
%%% Integral path coordinates
x_pre=real(alpha_2.');
y_pre=imag(alpha_2.');
 save(['xy_pre.mat'],'x_pre','y_pre'); 
%%% k0 coordinates
 x_k0=real(-k0);y_k0=imag(-k0);
  save(['k0_pre.mat'],'x_k0','y_k0'); 
%%% \delta \alpha of the integral path
alpha_c=alpha_2;
dalpha=0;
for ii=1:N1-1
    dalpha(ii)=alpha_c(ii+1)-alpha_c(ii);
end
dalpha(N1)=0;

%% Load factorization result of MINN.ipynb on the integral path 
G_PLUS=zeros(2,2,N1);G_NEG=zeros(2,2,N1);G_PLUS2=zeros(2,2,N1);
load('G_mat w_5.mat');% G_{+}^{-1} and G_{-}
G_PLUS(1,1,:)=pr+1i*pi;clear pi;
G_PLUS(1,2,:)=qr+1i*qi;
G_NEG(1,1,:)=pr2+1i*pi2;
G_NEG(1,2,:)=qr2+1i*qi2;
G_PLUS(2,1,:)=rr+1i*ri;
G_PLUS(2,2,:)=sr+1i*si;
G_NEG(2,1,:)=rr2+1i*ri2;
G_NEG(2,2,:)=sr2+1i*si2;
G_PLUS0=G_PLUS;
for ii=1:N1
    G_PLUS2(:,:,ii)=inv(squeeze(G_PLUS(:,:,ii)));
end
G_PLUS=G_PLUS2;  % G_{+}
%% Load factorization result of MINN.ipynb at k0
GK=zeros(2,2);
load('GK_mat w_5.mat');% G_{+}^{-1}(k0)
GK(1,1)=prk+1i*pik;
GK(1,2)=qrk+1i*qik;
GK(2,1)=rrk+1i*rik;
GK(2,2)=srk+1i*sik;
GK=inv(GK);  % G_{+}(k0)
%% Solutions of Hurd (1981) parallel plates scattering problem
PK=getp(-k0,k,d,theta0);
for ii=1:N1
    P(:,ii)=getp(alpha_c(ii),k,d,theta0);
    U(:,ii)=getu(alpha_c(ii),squeeze(G_PLUS(:,:,ii)),squeeze(G_NEG(:,:,ii)),P(:,ii),k,d,GK,PK,theta0);
    L(:,ii)=getL(alpha_c(ii),squeeze(G_PLUS(:,:,ii)),squeeze(G_NEG(:,:,ii)),P(:,ii),k,d,GK,PK,theta0);
end
beta_mat=zeros(N1,N1);
for ii=1:N1
    alpha=alpha_c(ii);
    U1=U(1,ii);U2=U(2,ii);L1=L(1,ii);L2=L(2,ii);
    b_plus=beta_plus(alpha_c(ii),k);
    b_neg=beta_neg(alpha_c(ii),k);
    beta=b_plus*b_neg;
    beta2=(alpha^2-k^2)^0.5;
    beta0=k*sin(theta0);
    psik=exp(-beta2*d);
    psi=exp(1i*beta*d);
    psi0=exp(1i*beta0*d);
    A1(1,ii)=beta^(-1)*(psi)^(-1)*U1+beta^(-1)*(psi)^(-1)*beta0*psi0^(-1)/(2*pi*1i*(alpha+k0));
    A3(1,ii)=(psi)^(-1)*U2-(psi)^(-1)*psi0/(2*pi*1i*(alpha+k0));
    B2(1,ii)=(psi^2)/(1+psi^4)*(-A1(1,ii)+psi^2*A3(1,ii)); 
    A2(1,ii)=(psi^2)/(1+psi^4)*(psi^2*A1(1,ii)+A3(1,ii));
    PSIT(1,ii)=psi;
    PSI0(1,ii)=psi0;
    ALPHA(1,ii)=alpha;
end

%% Compute near-field result of zones 1,2,3
N2=100;
N3=100;
LL=4;
v1=zeros(N3,N2);
v1i=zeros(N3,N2);
vin1=zeros(N3,N2);
y_mat1=zeros(N2,1);
begin1=1;end1=N1;
parfor ii=1:N2+1
    z_mat1=zeros(N3,1);
    for jj=1:N3
        yy=d+(LL-d)/N2*(ii-1);
        zz=-LL+(LL*2)/N3*jj;
        jifen=0;
        y_mat1(ii)=yy;z_mat1(jj)=zz;
        for jk=begin1:end1
            alpha=alpha_c(jk);
            b_plus=beta_plus(alpha_c(jk),k);
            b_neg=beta_neg(alpha_c(jk),k);
            beta=b_plus*b_neg;
            beta2=(alpha^2-k^2)^0.5;
            jifen=jifen+A1(1,jk)*(exp(1i*alpha*zz+1i*beta*(yy)))*dalpha(jk);
        end
        v1i(jj,ii)=(exp(-1i*k*(yy*sin(theta0)+zz*cos(theta0))));
        v1(jj,ii)=jifen;
    end
end
v2=zeros(N3,N2);
v2i=zeros(N3,N2);
y_mat2=zeros(N2,1);

parfor ii=1:N2+1
    z_mat2=zeros(N3,1);
    for jj=1:N3
        yy=-d+(2*d)/N2*(ii-1);
        zz=-LL+(LL*2)/N3*jj;
        jifen=0;
        y_mat2(ii)=yy;z_mat2(jj)=zz;
        for jk=begin1:end1
            alpha=alpha_c(jk);
            b_plus=beta_plus(alpha_c(jk),k);
            b_neg=beta_neg(alpha_c(jk),k);
            beta2=(alpha^2-k^2)^0.5;
            beta=b_plus*b_neg;
            jifen=jifen+A2(1,jk)*(exp(1i*alpha*zz+1i*(beta)*(yy)))*dalpha(jk)+B2(1,jk)*(exp(1i*alpha*zz-1i*(beta)*(yy)))*dalpha(jk);

        end
        v2i(jj,ii)=(exp(-1i*k*(yy*sin(theta0)+zz*cos(theta0))));
        v2(jj,ii)=jifen;
    end
end
v3=zeros(N3,N2);
v3i=zeros(N3,N2);
y_mat3=zeros(N2,1);

parfor ii=1:N2+1
    z_mat3=zeros(N3,1);
    for jj=1:N3
        yy=-LL+(-d+LL)/N2*(ii-1);
        zz=-LL+(LL*2)/N3*jj;
        jifen=0;
        y_mat3(ii)=yy;z_mat3(jj)=zz;
        for jk=begin1:end1
            alpha=alpha_c(jk);
            b_plus=beta_plus(alpha_c(jk),k);
            b_neg=beta_neg(alpha_c(jk),k);
            beta2=(alpha^2-k^2)^0.5;
            beta=b_plus*b_neg;
            jifen=jifen+(   A3(1,jk)*(exp(-1i*(beta)*(yy))*exp(1i*alpha*zz)) )*dalpha(jk);
        end
        v3i(jj,ii)=(exp(-1i*k*(yy*sin(theta0)+zz*cos(theta0))));
        v3(jj,ii)=jifen;
    end
end
for jj=1:N3
    zz=-LL+(LL*2)/N3*jj;
    z_mat1(jj)=zz;
    zz=-LL+(LL*2)/N3*jj;
    z_mat3(jj)=zz;
    zz=-LL+(LL*2)/N3*jj;
    z_mat2(jj)=zz;
end



%% View the near-field results
[Y,Z] = meshgrid(y_mat1,z_mat1);
s=surf(Z,Y,real(v1));colormap(gca,jet);view(2)
s.EdgeColor = 'none';
hold on;
[Y,Z] = meshgrid(y_mat2,z_mat2);
s=surf(Z,Y,real(v2));colormap(gca,jet);view(2)
s.EdgeColor = 'none';
hold on;
[Y,Z] = meshgrid(y_mat3,z_mat3);
s=surf(Z,Y,real(v3));colormap(gca,jet);view(2)
s.EdgeColor = 'none';
hold on;
y_dd=0.95:0.1/100:1.05;
z_dd=0.0:4/100:4.0;
[Y,Z] = meshgrid(y_dd,z_dd);
s=surf(Z,Y,abs(Z+1),'FaceColor','k');view(2)
s.EdgeColor = 'none';
hold on;
y_dd=-1.05:0.1/100:-0.95;
z_dd=0.0:4/100:4.0;
[Y,Z] = meshgrid(y_dd,z_dd);
s=surf(Z,Y,abs(Z+1),'FaceColor','k');view(2)
s.EdgeColor = 'none';
axis equal;

% colorbar;
xlim([-4,4])
ylim([-4,4])
caxis([-1.2 1.2])
set(gcf,'Position',[100 100 300 300]);
% title('v(\theta)')
  
set(gca,'FontSize',12,'Fontwei','bold','Fontname', 'Times New Roman')
%%
%%% Compute P(\alpha)
function pp=getp(alpha,k,d,theta0)
    beta0=k*sin(theta0);
    psi0=exp(1i*beta0*d);
    b_plus=beta_plus(alpha,k);
    pp1=beta0/psi0/b_plus/(2*pi*1i);
    pp2=-b_plus*psi0/(2*pi*1i);
    pp=[pp1;pp2];
end
%%% Compute U(\alpha)
function uu=getU(alpha,G_PLUS,G_NEG,P,k,d,GK,PK,theta0)

    beta0=k*sin(theta0);
    k0=k*cos(theta0);
    psi0=exp(1i*beta0*d);
    b_plus=beta_plus(alpha,k);
    A=[b_plus 0;0 b_plus^(-1)];
    B=G_PLUS*inv(GK)*PK-P;
    uu=(alpha+k0)^(-1)*A*B;
end
%%% Compute L(\alpha)
function LL=getL(alpha,G_PLUS,G_NEG,P,k,d,GK,PK,theta0)


    beta0=k*sin(theta0);
    k0=k*cos(theta0);
    psi0=exp(1i*beta0*d);
    b_plus=beta_plus(alpha,k);
    b_neg=beta_neg(alpha,k);
    A=[b_neg^(-1) 0;0 b_neg];
    B=inv(G_NEG)*inv(GK)*PK;
    LL=2*(alpha+k0)^(-1)*A*B;
end
%%% Compute beta_{+}(\alpha)
function pp=beta_plus(alpha,k)

    b_plus=(k+alpha).^0.5;
    b_neg=(k-alpha).^0.5;
             beta=b_plus*b_neg;
            if imag(beta)<0
                pp=-b_plus;
            else
                pp=b_plus;
                
            end
end
%%% Compute beta_{-}(\alpha)
function pp=beta_neg(alpha,k)

    pp=(k-alpha).^0.5;
end