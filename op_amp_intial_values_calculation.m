%Vellios Georgios-Serafeim
%AEM: 9471

x=71;
CL=(2+0.01*x)*10^-12;
SR=(18+0.01*x)*10^6;
Vdd = (1.8 + 0.003*x);
Vss = -Vdd;
GB = (7 + 0.01*x)*10^6;
GBrad=GB*2*pi;
AdB = (20 + 0.01*x);
P = (50 + 0.01*x)*10^-3;
Vinmax=0.1;
Vinmin=-0.1;

kntonos=175*10^-6;
VTOn=0.7860;
lamdan=0.05;

kptonos=60*10^-6;
VTOp=-0.9056;
lamdap=0.15;

Cox=4.6*10^-3;


%step1 thewrw oti einai 1u gia kales pra3eis. arkei na nai panw apo 0.35
L=1 * 10^-6;

%step2 to pernw pio megalo apo 0.22CL
Cc=0.22*CL +0.01*CL;

%step3
I5=SR*Cc;

I3=I5/2;
I2=I3;

%step4
S3=I5/((kptonos)*((Vdd-Vinmax-abs(VTOp-0.15)+(VTOn-0.15))^2));

if S3<1
    S3=1;
end

S4=S3;
W3=ceal(S3) *L;
W4=W3;


%step 5 
p3=(sqrt(2*kptonos*I3*S3)/(2*0.667*W3*L*Cox));
p3hz=p3*0.1591549;
if p3hz>10*GB
    disp('p3>10GB');
else
    disp('p3 not geater than 10GB');
end

%step6 
gm1=GBrad*Cc;
gm2=gm1;
S2=(((gm2)^2)/(kntonos*I5));
S1=S2;
W2=ceal(S2) *L;
W1=W2;

%step 7
Vds5=Vinmin-Vss-sqrt(I5/(kntonos*S1))-(VTOn+0.15);

if Vds5>0.1
    disp('Vds5 is greater than 100mV');
else
    disp('Vds5 is lesser than 100mV');
end

S5=2*I5/(kntonos*((Vds5)^2));
W5=ceal(S5) *L;
W8=W5;

%step 8
gm6=2.2*gm2*(CL/Cc);
gm4=sqrt(2*kptonos*S4*I2);
S6=S4*(gm6/gm4);
I6=((gm6)^2)/(2*kptonos*S6);

W4=ceal(S4) *L;
W6=ceal(S6) *L;

%step 10
S7=(I6/I5)*S5;
W7=ceal(S7) *L;

%step 11
A=(2*gm2*gm6)/(I5*(lamdan+lamdap)*I6*(lamdap+lamdan));
Pdiss=(I5+I6)*(Vdd+ abs(Vss));

if 20*log10(A) >AdB
    disp('Calculated A is greater than A initial');
else
     disp('Calculated A is lesser than A initial');
end

if Pdiss<P
    disp('Pdiss is lesser than P');
else
    disp('Pdiss is greater than P');
end
