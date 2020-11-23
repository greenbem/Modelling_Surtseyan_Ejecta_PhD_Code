function SursteyEjectaModel

% Considering a spherical ejecta with a centred inclusion
% and the model determined in the notes with a coordinate
% transformation to fix the moving boundary and using 
% method of lines on the spacial terms to solve.

% The space is split into two sections the inclusion and
% the magma and each section has a different coordinate
% spacing that meet at the boundary between the two. Also
% in each section has an arctan coordinate transform applied
% to it in order to retain the uniformity of the steps but
% allow for small step sizes at the boundary.This allow for
% smooth movement from the smallest to largest mesh sizes 

%% Typical Values
format long
M = 18*10^(-3);% Molecular mass of water
hsl = 2300000 ;% the latent heat of vaporisation
Rg= 8.314 ;% Gas Constant
Tm= 1275 ;% Original magma temperature
T0=373 ;% Original inclusion temperature
R1=0.01 ;% Radius of the inclusion
R2=0.1; % Radius of the magma ball
k=10^(-14); % permeability of magma
muv=3*10^(-5); % viscosity
p0=10^(5); % atmospheric pressure
rhom=2750; % Basalt density Range: 2650-2800
cpm=840; % specific heat Basalt
cps=2000; % specific heat of water vapour
rhol=1000; % Density of water
cpl=4200; % specific heat water
phi=0.4 ; % Porosity
rhol0=rhol; % Water Density ND
rhos0=(p0*M)/(Rg*Tm); % Vapour Density ND
Ke= 2; %thermal conductivity magma and vapour
Kel=3; % Effective thermal conductivity in the inclusion
pc=((1-phi)*rhom*cpm+phi*rhos0*cps); % specific heat and
% density comb magma vapour
pci=((1-phi)*rhom*cpm+phi*rhol*cpl);% specific heat and
% density comb inclusion 
t0=(phi*rhol0*hsl*R1^2)/(3*Ke*(Tm-T0));% time scaling
vs0=(R2*rhol0)/(t0*rhos0); % velocity of vapour
epsilon= R1/R2;% Original boundary position(does not move)
p0m=10^(5); % Initial magma pressure 
pre=10^(-5); %Size of the pore space;


%% Creating the mesh
% Using an arctan transformation mesh set up
% Setting m and n to give correct spacing
% In the inclusion and magma
n=200;
m=800;

smallestMesh=pre/R2; 
% the smallest desired mesh size is set to twice the
% typical diameter of pores, taken to be 10 microns = 1E-05
% now set aa, given n=m, so that the smallest mesh size
% is roughly smallestMesh:

aa = 2*(1+smallestMesh*n)/( n*(n+1)*smallestMesh*pi);

%constants needed later
Bm = atan(aa*(m+1)); 
Bn = atan(aa*(n+1)); 

chi=linspace(0,m+1,m+2); % 0,1,2,3, ... n+1 in value
psi=linspace(0,n+1,n+2); % 0,1,2,3, ... m+1 in value
zeta=epsilon*atan(aa*chi)/atan(aa*(m+1));  
% zeta in slurry goes from 0 to epsilon, tightens on epsilon
xi= 1 + (1-epsilon)*atan(aa*(psi-n-1))/atan(aa*(n+1)); 
% xi in hot magma goes from epsilon to 1, 
% and tightens mesh on epsilon
Cn = (Bn./(aa.*(1-epsilon))).*(1+(aa.*(psi-n-1)).^2); 
% vector term that arises from chain rule

Cm = (Bm./(aa.*epsilon))*(1+(aa.*chi).^2); 
% vector term that arises from chain rule                                          

if (xi(2)-xi(1)) < smallestMesh/2   
% possibly n is too large, smallest spacing is too small
    disp('bailing out, smallest mesh size is too small maybe n is too large?')
    display((xi(2)-xi(1)))
    display((zeta(m+2)-zeta(m+1)))
    disp('Magma Step Size Error')
end


 if (zeta(m+2)-zeta(m+1)) < smallestMesh/2   
 % possibly m is too large, smallest spacing is too small
    disp('bailing out, smallest mesh size is too small,maybe m is too large?')
    display((xi(2)-xi(1)))
    display((zeta(m+2)-zeta(m+1)))
    disp('Slurry Step Size Error')
 end

%% Parameters from the ND Model
H=(M*hsl)/(Rg*T0);
f9=(Kel*t0)/(R2^2*pci);
f11=(t0*k*p0)/(phi*muv*R2^2);
f15=(Ke*t0)/(R2^2*pc);
f3=(vs0*R2*phi*rhos0*hsl)/(Ke*Tm);
f10=(k*p0)/(muv*R2*phi*vs0);

%% Setting up a time mesh

ts=0 ;% start time
te=0.1; % end time
tsteps=500;% Number of time steps
tspan= [0 logspace(-7, -1,tsteps-1)];

%% Initial Values of y set up

y0= zeros(1,2*n+m+1);
y0(1:m)=375/1275; % Initial Temp in Inclusion
y0(m+1:n+m)=1; % Initial Temp in Magma
y0(n+m+1:2*n+m)=1; % Initial density in magma
y0(2*n+m+1)=R1/R2;% position of flash front


% Setting up the different temp, density and
% position of the boundary matrices
Tinc=zeros(1,m+2);
Tmagma= zeros(1,n+2);
rhomagma=zeros(1,n+2);
Pmagma=zeros(1,n+2);

%% boundary condition at the surface of the 
% magma from fixed end point

Tmagma(n+2)=375/1275;
Pmagma(n+2)=p0m/p0;
rhomagma(n+2)= Pmagma(n+2)/Tmagma(n+2);



%% Setting up to solve
options=odeset('RelTol', 10^(-8), 'AbsTol', 10^(-8), ...
'events', @functionevents, 'OutputFcn', @odeprint);

sol=ode15s(@SurtseyModelODE,[ts te],y0,options);


tspread=tspan;
%provides the solutions for the wanted times
Y=deval(sol,tspread); 

% %Splitting Y into different parts
 sval= Y(end,:);
 % svalR=transpose(sval);
%display(sval)
 Tincvals= Y(1:m,:);
%display(Tincvals)
 Tmagmavals= Y(m+1:n+m, :);
%display(Tmagmavals)
 %rhomagmavals=Y(j+m+1:2*j+m, :); 
%display(rhomagmavals)
 Pmagmavals=Y(n+m+1:2*n+m, :).*Y(m+1:n+m, :);
%display(Pmagmavals)


 
%% Boundary  
TBvals=zeros(1,tsteps);
PBvals=zeros(1,tsteps);
for i=1:tsteps
    TBvals(i)=functionboundary(Tmagmavals(1,i), ...
    Tincvals(m,i),Pmagmavals(1,i),sval(i));
    PBvals(i)=exp(H*((TBvals(i)-(T0/Tm))/TBvals(i)));
end

plot(log10(tspan),(PBvals),'-o')
 title('Pressure vs Time at boundary')
 xlabel('log time')
 ylabel('pressure')
 %% Function
function f=SurtseyModelODE(~,y)
 
% This function generates the right-hand sides of the 
% ODEs using the method of lines. This is discretised
% into two sections the inclusion and the magma sections. 
% The input is the column vector y and t with is a scalar

% y is a input column vector of solutions. Note f should
% therefore also be a column vector. There are m 
% temperature values from inside the inclusion not 
% including a boundary condition at either end (as they 
% are not calculated with the right hand sides they
% are imputed and n temperature values from the magma
% which also do not include the boundaries. Then there
% are a further n density values from inside the magma,
% boundary not included. So y is a column vector with 
% the size m+2n+1

% Input values
%Values of Tinc not including boundary at ends
Tinc(2:m+1)=y(1:m); 
%Values of Tmagma not including boundary at ends
Tmagma(2:n+1)=y((m+1):(n+m)); 
%Values of density not including boundary at ends
rhomagma(2:n+1)=y(((n+m)+1):(2*n+m));
% position of the flash front 
s=y(2*n+m+1); 
 %Pressure Ideal Gas Law
Pmagma(2:(n+1))=rhomagma(2:(n+1)).*Tmagma(2:(n+1));

%set up the boundary conditions 
% the interior boundary at the centre of the 
% inclusion Tinc(1)=y(1) has a temperature flux of zero. 

Tinc(1)= Tinc(2);

%boundary condition between inclusion and magma
tau = functionboundary(Tmagma(2),Tinc(m+1),Pmagma(2),s);

Tinc(m+2)=tau;
Tmagma(1)=tau;

% The boundary condition of the density between
% the inclusion and magma

Pmagma(1)= exp(H*(1-T0/(Tm*Tmagma(1))));
rhomagma(1)= Pmagma(1)/Tmagma(1);

%setting up diff
% from 1 to m+1 is a a upstream difference at Tinc(2) 
% with the first term of Tinc(2)-Tinc(1) taking 2 to m+2
% this gives a downstream difference of Tinc(3)-Tinc(2)
% same idea with the others

DTinc=diff(Tinc); 
DPmagma= diff(Pmagma); 
DTmagma=diff(Tmagma);
Drhomagma=diff(rhomagma);

%derivatives used in right hand sides of equations note
%the differences in chi and psi are all 1

% pressure in magma 1st derivative with i-1 &i
DPN=DPmagma(2:n+1); 
DPP=DPmagma(1:n);
% Temp in Inc 1st derivative with i & i-1
DTiP=DTinc(1:m); 
DTiN=DTinc(2:m+1);
% Temp in Magma 1st derivative with i & i-1
DTmP=DTmagma(1:n); 
DTmN=DTmagma(2:n+1);
% Density in Magma 1st derivative with i &i-1
DrhoP=Drhomagma(1:n);
DrhoN=Drhomagma(2:n+1);


%second derivatives
%Tinc second derivative
DDTi=zeros(1,m);
for i1=2:m+1
  DDTi(i1-1)=(Tinc(i1+1)-2.*Tinc(i1)+Tinc(i1-1));
end
%Tmag second derivative
DDTm=zeros(1,n);
for i2=2:n+1
  DDTm(i2-1)=(Tmagma(i2+1)-2.*Tmagma(i2)+Tmagma(i2-1));
end

%Pmag second derivative
DDP=zeros(1,n);
for i3=2:n+1
  DDP(i3-1)=(Pmagma(i3+1)-2.*Pmagma(i3)+Pmagma(i3-1));
end


%% Calculating sdot

sdot=functionvelocity(Tmagma(2),Tmagma(1),s, ...
    Tinc(m+2),Tinc(m+1));

%Now to calculate the f function for the Tinc

fTinc=functioninctemp(s,DTiP,DTiN,sdot,DDTi);

%Now calculate the f function for the Tmagma

fTmagma=functionmagma(s,DTmP,DTmN,sdot,DDTm,xi(2:n+1));

%Now Calculate the f function of rhomagma

frhomagma=functionrhomagma(s,rhomagma(1:n+1),DPN, ...
        DrhoN,DrhoP,sdot,DPP,xi(2:n+1),n, DPmagma);


% Now to form the right hand side f
fr=[fTinc,fTmagma,frhomagma,sdot];
f=transpose(fr);
end

%% stop function
function [value,isterminal,direction]= functionevents(~,y)
szero=y(end)-10^(-4)*R1; %value that equals zero 
stop=1; %halt the integration
side=0; % when coming from either direction, +1 
        % increasing function -1 decreasing
value=szero;
isterminal=stop;
direction=side;
end


%% boundary function
function [tau]=functionboundary(Th,Ti,Pm,sb)
    function [result]=bfun(tau)
    epsfrac= (1-epsilon)/(1-sb); 
    esfrac=epsilon/sb; 
    Pb=exp(H*((tau-(T0/Tm))/tau));
    LHS=f10*f3*Pb*Cn(1)*epsfrac*(Pm-Pb);
    RHS=-tau*(epsfrac*Cn(1)*(Th-tau)-Cm(m+1)*esfrac*(tau-Ti));
    result= LHS-RHS; 
    end

Min=0.1;
r1=bfun(Min); r2=bfun(0.5); 
%T=0.5 is approximately critical temperature
            if r1*r2 > 0   % same sign, not straddling!!
                disp('not straddling zero')
                disp(sb)
                disp(Th)
                disp(Pm)
                disp(Ti)
                % LookAtFlash
                Tflash= fzero(@bfun,Ti); 
                % look for a solution near the slurry T
            else % we do have straddle:
            Tflash= fzero(@bfun,[Min 0.5]); 
            % solve BC at flash for temperature there
            end  
tau=Tflash;
end 
%% s dot function
function [sdot]=functionvelocity(Tm2,Tm1,sp,Tie,Tin)
%Parameters that are needed in the function 
epsfrac1= (1-epsilon)./(1-sp);
esfrac1=epsilon./sp; 

%Calculate sdot
fluxmagma=(Tm2-Tm1);
fluxinc=(Tie-Tin);
sdot=-(1/f3).*(Cn(1)*epsfrac1.*fluxmagma-Cm(m+1)*esfrac1.*fluxinc); 
end

%% T Inc Function
function [fTinc]=functioninctemp(s,DTiP,DTiN,sdot,DDTi)
esfrac=epsilon/s; 
VT=(((esfrac.^2)*.2.*f9)./zeta(2:m+1))+((zeta(2:m+1).*sdot)./s);
a43=((esfrac^2)*f9);
A1=VT.*Cm(2:m+1)+a43.*((2.*Cm(2:m+1).*Bm.*aa.*chi(2:m+1))./epsilon);
a44P=min(A1,0);
a44N=max(A1,0);
WDTi=a44P.*DTiP+a44N.*DTiN; %Winding
fTinc= (a43.*Cm(2:m+1).^2).*DDTi+WDTi;        
    
end
%% T magma Function
function [fTmagma]=functionmagma(s,DTmP,DTmN,sdot,DDTm,xi)
epsfrac= (1-epsilon)/(1-s);
gfrac=((1-s).*xi+(s-epsilon))./(1-epsilon);
WT=(((epsfrac)*2*f15)./(gfrac))-(((xi-1).*sdot)./(1-s));
a53=((epsfrac^2)*f15);
A2=WT.*Cn(2:n+1)+a53.*((2.*aa.*Cn(2:n+1).*Bn.*(psi(2:n+1)-n-1)) ... 
        ./(1-epsilon));
a54P=min(A2,0);
a54N=max(A2,0);
WDTm=a54P.*DTmP+a54N.*DTmN; %Winding

fTmagma= (Cn(2:n+1).^2).*a53.*DDTm+WDTm;   
end

%% rho magma function
function [frhomagma]=functionrhomagma(s,rhomagma,DPN,DrhoN, ...
                        DrhoP,sdot,DPP,xi,n,DPmagma)
epsfrac= (1-epsilon)/(1-s); 
gfrac=((1-s).*xi+(s-epsilon))./(1-epsilon);
Xp=(f11.*2.*rhomagma(2:n+1).*epsfrac./gfrac);
Xrho=-((xi-1)./(1-s))*sdot;
a61=((f11.*(epsfrac).^2));
A3=Xrho.*Cn(2:n+1);
A4=Xp.*Cn(2:n+1)+((a61.*rhomagma(2:n+1).*Cn(2:n+1).*2.*aa.* ...
    Bn.*(psi(2:n+1)-n-1))/(1-epsilon));

a62P=min(A3,0);
a62N=max(A3,0);
WDrhoZ=a62P.*DrhoP+a62N.*DrhoN; 

a64P=min(A4,0);
a64N=max(A4,0);
A4WD=a64P.*DPP+a64N.*DPN; 

DD=diff(rhomagma.*DPmagma);
frhomagma= (Cn(2:n+1).^2).*a61.*DD+WDrhoZ+A4WD;
end
end