%% Function to calculate strengthening from a maraging steel
% _See [1] T.H. Simm, L. Sun, D.R. Galvin, E.P. Gilbert, D. Alba Venero, Y. Li, T.L. Martin, P.A.J. Bagot, M.P. Moody, P. Hill, H.K.D.H. Bhadeshia, S. Birosca, M.J. Rawson, K.M. Perkins, D.A. Venero, Y. Li, T.L. Martin, P.A.J. Bagot, M.P. Moody, P. Hill, H.K.D.H. Bhadeshia, S. Birosca, M.J. Rawson, K.M. Perkins, A SANS and APT study of precipitate evolution and strengthening in a maraging steel, Mater. Sci. Eng. A. 702 (2017) 414–424. doi:10.1016/j.msea.2017.07.013.__
function PrecipitateStrengthEquations(x0, betax, lavesx, Dbeta, Dlaves,atyp)
%% Function main
% fitting variables 1=laves, 2=beta, 3=ratio2, 4=lath, 5=SSS, 6=0 for
% both, -1 for NiAl and 1 for laves
if nargin==0 %loads data
    atyp=2;%for alloy 9922
    x0 = [1.9267    0.9565    1.4241  839.7736];%the parameters to describe strength contributions
    load('9922dataOuts')%the size, volume fraction and ageing time of the precipitates
end
% the strength contributions constants 
lavesOro=x0(1);
betaST=x0(2);
lath_iron=x0(4)*(1/3.2);
ratio2=x0(3);
% if we want to extract individual components
chois = 'alls';       
xx=1;

%% strength functions

%For Laves phase 
HVL=Orowan2(Dlaves,lavesx);HVL(1,:)=0;
HVL=lavesOro*HVL;%laves orowan

%For Beta phase
[HVB_ , Fsh, Foro]=ShearOrowan(Dbeta,betax,ratio2); 
HVB_=HVB_*betaST;
%Shear and Orowan parts to Beta
Fsh =Fsh*betaST;     Foro = Foro*betaST;

% solid solution strengthening
SSS=solidsolution(betax, lavesx, atyp);%1= NiAl 2=Laves

% select what to get out
switch chois
    case 'beta'
        F_=(  HVB_);
    case 'besh'
        F_=Fsh;
    case 'boro'
        F_=Foro;
    case 'lave'
        F_=(  HVL);
    case 'sols'
        F_=(   SSS );
    case 'alls'
        F_=(  lath_iron + HVL + HVB_ + SSS );
end
    
Fout(:,:,xx)=F_;
if size(Fout,2)==3
    Fout(:,4)=0*Fout(:,3);
end

%% Plots

h=figure(10);

plot(time,Dbeta(:,1),'bv-','linewidth',2)
hold on
plot(time,Dbeta(:,2),'go-','linewidth',2)
plot(time,Dbeta(:,3),'r^-','linewidth',2)

plot(time,Dlaves(:,1),'bv--','linewidth',2)
plot(time,Dlaves(:,2),'go--','linewidth',2)
plot(time,Dlaves(:,3),'r^--','linewidth',2)
xlim([0 24])
xlabel('Ageing time (h)')
ylabel('Precipitate size ECD (nm)')
legend('NiAl 520 \circC','NiAl 540 \circC','NiAl 560 \circC',...
    'Laves 520 \circC','Laves 540 \circC','Laves 560 \circC')
%
h=figure(11);

plot(time,betax(:,1),'bv-','linewidth',2)
hold on
plot(time,betax(:,2),'go-','linewidth',2)
plot(time,betax(:,3),'r^-','linewidth',2)

plot(time,lavesx(:,1),'bv--','linewidth',2)
plot(time,lavesx(:,2),'go--','linewidth',2)
plot(time,lavesx(:,3),'r^--','linewidth',2)
xlim([0 24])
xlabel('Ageing time (h)')
ylabel('Precipitate volume fraction')
legend('NiAl 520 \circC','NiAl 540 \circC','NiAl 560 \circC',...
    'Laves 520 \circC','Laves 540 \circC','Laves 560 \circC')
%
h=figure(1);

plot(time,Fout(:,1),'bv-','linewidth',2)
hold on
plot(time,Fout(:,2),'go-','linewidth',2)
plot(time,Fout(:,3),'r^-','linewidth',2)
xlim([0 24])
xlabel('Ageing time (h)')
ylabel('Strength combined (HV)')
legend('520 \circC','540 \circC','560 \circC')
%
h=figure(2);

plot(time,SSS(:,1),'bv-','linewidth',2)
hold on
plot(time,SSS(:,2),'go-','linewidth',2)
plot(time,SSS(:,3),'r^-','linewidth',2)
xlim([0 24])
xlabel('Ageing time (h)')
ylabel('Solid Solution Strength (HV)')
legend('520 \circC','540 \circC','560 \circC')
%
h=figure(3);

plot(time,HVL(:,1),'bv-','linewidth',2)
hold on
plot(time,HVL(:,2),'go-','linewidth',2)
plot(time,HVL(:,3),'r^-','linewidth',2)
xlim([0 24])
xlabel('Ageing time (h)')
ylabel('Laves Strength (HV)')
legend('520 \circC','540 \circC','560 \circC')
%
h=figure(4);

plot(time,HVB_(:,1),'bv-','linewidth',2)
hold on
plot(time,HVB_(:,2),'go-','linewidth',2)
plot(time,HVB_(:,3),'r^-','linewidth',2)
xlim([0 24])
xlabel('Ageing time (h)')
ylabel('NiAl Strength (HV)')
legend('520 \circC','540 \circC','560 \circC')
end
%% Shear and Orowan functions for Beta phase
function [HVB_, HVBSi__,HVBOi__ ]=ShearOrowan(Dbeta,betax,ratio2)

if nargin==2
    ratio2=1.15;
end

for ii=1:size(Dbeta,2)%temp
    for ij=1:size(Dbeta,1)%time ignoring 0h
%mean values        
DbetaM=Dbeta(ij,ii);betaxM=betax(ij,ii);
%create a distribution of sizes
xsz=DbetaM-1:.1:DbetaM+1;
fDbetaD=normpdf(xsz , DbetaM, 1);
fDbetaD=fDbetaD/sum(fDbetaD);
for jj=1:length(xsz)
    HVBOi_(jj)=Oro(xsz(jj),betaxM);
    HVBSi_(jj)=ratio2*Shh(xsz(jj),betaxM);
    % find beta strengthening- from lower of shear or orowan
    [HVBi_(jj), posmn(jj)]=min([HVBSi_(jj),HVBOi_(jj)]);

end
HVB_(ij,ii)=sum(HVBi_.*fDbetaD);
HVBSi__(ij,ii)=sum(HVBSi_.*fDbetaD);
HVBOi__(ij,ii)=sum(HVBOi_.*fDbetaD);
    end
end

HVB_((isnan(HVB_)==1))=0;

if isnan(HVBSi_(1))==1
    HVBSi_(1)=0;
end
if isnan(HVBOi_(1))==1
    HVBOi_(1)=0;
end

end
%% Orowan equation for Beta
function FO=Oro(ECD,vf)
D=ECD*1e-9;        
M=3.1;
b=0.254e-9;G=70e9;mu=0.31;
rs=(2/3)^0.5*D/2;
A=0.81*M*G*b;B=log(2*rs/b);
C=4*pi*rs;
D= (pi./(4)).^.5 - .2 ; E=(1-mu)^0.5;
FO= (vf).^.5*1e-6*A.*B ./ (C.*D.*E);
FO=(1/3.2)*FO/4;
end

%% Shear Equation for Beta
function FSh=Shh(ECD,vf)

% Friedel from Schitzner
b=0.254e-9;G=70e9;%mu=0.31;
r=ECD*0.5e-9;
gamma=0.256;
S=G*b^2/2;
wr=.82;wq=.75;
A=1/b;
B=2*S/(pi*wq)^.5;
C=(gamma*wr*r/S).^1.5;
D=(vf.^.5)./r;
FSh=1e-6*A.*B.*C.*D;
FSh=FSh*4*(1/3.2);
end
%% Orowan Equation for Laves phase
function F=Orowan2(D,vf)
%         METALLURGICAL TRANSACTIONS A VOLUME 7A, FEBRUARY 1976-293
D=D*1e-9;        
M=3.1;
b=0.254e-9;G=70e9;mu=0.31;
rs=(2/3)^0.5*D/2;
A=0.81*M*G*b;B=log(2*rs/b);
C=4*pi*rs;
D= (pi./(4*vf)).^.5 - 1 ; E=(1-mu)^0.5;
F= 1e-6*A.*B ./ (C.*D.*E);
%%
F=(1/3.2)*F/4;
end
%% Solid Solution Strength Equations
function [SSstrength, mWTpc]=solidsolution(vf1, vf2, alloy)%1= NiAl 2=Laves
eleme={'Al'; 'Co'; 'Cr';'Fe'; 'Mo'; 'Ni';'W'};
[p1ATMpc ~]=precip_ATMpc(1);%1= NiAl 2=Laves
[p2ATMpc ~]=precip_ATMpc(2);%1= NiAl 2=Laves

for ij=1:size(vf1,1)%ageing time
for ii=1:size(vf1,2)%ageing temp.
    [mATMpc ,~]=matrix_ATMpc(vf1(ij,ii),vf2(ij,ii),p1ATMpc,p2ATMpc,alloy);
    mWTpc(:,ii,ij)=convert_atm2weight(mATMpc)*100;
    SSstrength(ij,ii)=70+mWTpc(1,ii,ij)*46.5 + mWTpc(5,ii,ij)*46.5 - mWTpc(6,ii,ij)*14 + mWTpc(7,ii,ij)*46.5 ;
end
end

SSstrength = SSstrength*(1/3.2);
end
    function mWTpc=convert_atm2weight(mATMpc)
        atm_mass=[26.981539 58.933 51.996 55.845 95.94 58.6934 183.84];
        mWTpc=mATMpc.*atm_mass;
        mWTpc=mWTpc/sum(mWTpc);
    end

%%
%%
function [ATMpc_ ,cell_det]=matrix_ATMpc(vf1,vf2,p1ATMpc,p2ATMpc,alloy)
switch alloy
    case 1%lowAl
ATMpc(1)=2.56/100;ATMpc(2)=7.96/100;ATMpc(3)=10.91/100;
ATMpc(4)=69.35/100;ATMpc(5)=1.62/100;ATMpc(6)=6.85/100;ATMpc(7)=0.75/100;
    case 2%9922
ATMpc(1)=3.58/100;ATMpc(2)=7.91/100;ATMpc(3)=9.72/100;
ATMpc(4)=68.31/100;ATMpc(5)=1.17/100;ATMpc(6)=8.71/100;ATMpc(7)=0.60/100;
    case 3%F1E
ATMpc(1)=3.66/100;ATMpc(2)=7.95/100;ATMpc(3)=10.74/100;
ATMpc(4)=68.49/100;ATMpc(5)=1.62/100;ATMpc(6)=6.80/100;ATMpc(7)=0.75/100;
end
cell_det(1)=2.875^3*1e-24;%vol unit cell
cell_det(2)=2;%atoms/cell
%%
vfM=1-vf1-vf2;
ATMpc_ = (ATMpc - p1ATMpc*vf1 - p2ATMpc*vf2)/vfM ;


for n=1:length(ATMpc)
    if ATMpc_(n)<0
%         eleme={'Al'; 'Co'; 'Cr';'Fe'; 'Mo'; 'Ni';'W'};
        
%         disp(['laves error <0% ', eleme{n},' vf: ',num2str(vf2),'  ',num2str(100*ATMpc_)])
        ATMpc_(n)=0;
    end
end
ATMpc_=ATMpc_/sum(ATMpc_);
end
function [ATMpc ,cell_det]=precip_ATMpc(Ni_lav)
%% given precipitate type gets atomic % and unit cell
if Ni_lav==1%element order eleme={'1)Al, 2)Co 3)Cr 4)Fe 5)Mo 6)Ni 7)W'}
%NiAl
% ATMpc(1)=16.10/100;ATMpc(2)=5.62/100;ATMpc(3)=6.00/100;
% ATMpc(4)=43.16/100;ATMpc(5)=0.66/100;ATMpc(6)=27.60/100;ATMpc(7)=0.03/100;
%Nial v2
ATMpc(1)=35.99/100;ATMpc(2)=4.07/100;ATMpc(3)=4.28/100;
ATMpc(4)=28.04/100;ATMpc(5)=0.49/100;ATMpc(6)=26.34/100;ATMpc(7)=0.16/100;

cell_det(1)=2.875^3*1e-24;%vol unit cell
cell_det(2)=2;%atoms/cell
else
%laves
ATMpc(1)=0.44/100;ATMpc(2)=5.96/100;ATMpc(3)=22.31/100;
ATMpc(4)=42.14/100;ATMpc(5)=17.09/100;ATMpc(6)=3.46/100;ATMpc(7)=6.37/100;
cell_det(1)=1.589e-22;%vol unit cell
cell_det(2)=12;%atoms/cell
end
end