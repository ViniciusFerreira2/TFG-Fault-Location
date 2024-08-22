close all
clear all

% dados gerais do relé e sistema

f0=60;
w=2*pi*f0;
T=1/f0;
fs=16*f0; % freq de amostragem
nac=fs/f0;
dt=1/fs; % intervalo
ph=0;
h=3; %nac-1;

ncpf=6;
Iatual=4;
Imax=sqrt(2)*Iatual;

% ajustes do relé de sobrecorrente

curve=1;
dial=0.1;
Ipk=5; % partida

% constantes de ajuste 

if curve == 1
   kk=dial*0.14*fs;  % curva NI
elseif curve == 2
   kk=dial*13.5*fs;  % curva MI
elseif curve == 3
   kk=dial*80*fs;  % curva EI
else
   kk=dial*120*fs;  % curva TL
end;
   
% dados do curto-circuito   
   
   pucc=15; % pu de curto-circuito
   %pucc=8; %23.0769
   
   T1=0.1; % constante de tempo até ponto de falta
   
   teta=atan(w*T1)*180/pi   % ANGULO DO SISTEMA
   
   psi=0; % INSTANTE DA FALTA
   %psi=teta;    % INSTANTE DA FALTA
   
   delta=psi-teta;
   delta=delta*pi/180; % ângulo de caracterização
   

% Cnstantes de programação

k=1;
mm=1;
fl=0;

% b=(1/w)^2;
c=1/sqrt(2);

for k=1:nac-1
    I1(k)=Imax*sin(w*dt*(k-1));
    Ih(k)=ph*Imax*sin(h*w*dt*(k-1));
    I(k)=I1(k)+Ih(k);
end;
%k=2;
%I(k)=Imax*sin(w*dt*(k-1));
n=k;
nk=0;
s=0;
p=0;
ni=0;

% ciclo de processamento


while fl==0
   
   n=n+1;
   if n < (ncpf*nac)+1 
      mm=1;
      del=0;
      nk=0;
   else
      mm=pucc;
      del=delta;
      nk=nk+1;
   end;
   
   Ialt(n) = Imax*sin(w*dt*(n-1)+del)*mm;
   Idc(n) = -Imax*mm*exp(-dt*nk/T1)*sin(del);
   Ih(n) = ph*Imax*sin(h*(w*dt*(n-1)+del))*mm;
   I(n)=Ialt(n)+Idc(n)+Ih(n);
   
   %pause
   
% algoritmo da DFT ciclo completo
somaycI=0;
somaysI=0;     
for i=ni:nac-1+ni
   somaycI=somaycI+I(i+1)*cos(2*pi*i/nac);
   somaysI=somaysI+I(i+1)*sin(2*pi*i/nac);
end;
ni=ni+1;   
YcI(n)=(2/nac)*(somaycI);
YsI(n)=(2/nac)*(somaysI);

Ief(n)=c*sqrt(YcI(n)^2 + YsI(n)^2);

% algoritmo do relé  
   M(n-1)=Ief(n-1)/Ipk;
   
   if M(n-1) <= 1
      s=0;
   else
      if curve == 1
         s=s+(M(n-1))^0.02-1;
      elseif curve == 2
         s=s+M(n-1)-1;
      elseif curve == 3
         s=s+M(n-1)^2-1;
      else
         s=s+M(n-1)-1;
      end
      
      p=p+1;
      if s>=kk
         fl=1;
      end;
   end;
end;

% TRIP - resultados


t=dt*[0:n-1];

% apenas para calcular com a formula, considerando o ultimo M ou M
% constante
Mt=(pucc*Iatual)/Ipk;
top=p*dt

if curve == 1
   temp=dial*0.14/(Mt^0.02-1)
elseif curve == 2
   temp=dial*13.5/(Mt-1)   
elseif curve == 3
   temp=dial*80/(Mt^2-1)   
else
   temp=dial*120/(Mt-1)   
end

% ********************
fs=16*60;

M=[0 M];
figure;
plot(t,Ief,'g');
ylabel("Valor RMS");
xlabel("tempo s")
figure;
plot(t,M,'g');
ylabel("Múltiplo M");
xlabel("tempo s")
figure;
plot(t,Ialt,'r',t,Idc,'k',t,I,'b'); %,t,v,'k');
ylabel("Icc Amp");
xlabel("tempo s");
legend("Iac","Idc","Icc");
figure;
plot(t,Ief,'m',t,I,'b');
ylabel("Icc Amp");
xlabel("tempo s")
legend("Irms","Icc");



