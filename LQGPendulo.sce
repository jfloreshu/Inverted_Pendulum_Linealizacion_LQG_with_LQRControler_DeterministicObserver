
// load the data
load("pendulojLTILab1.sod","X","U","sys")

Ap=sys.A;
Bp=sys.B;
Cp=sys.C;
Dp=sys.D;

//checking controllability and observability
[i,j] = size(Ap);
// e=[B, AB, A^2 B,..., A^(n-1) B]
e = cont_mat(sys.A,sys.B);
rankC=rank(e);
if i == rankC then
    disp('Continuous System is Controllable');
end

// o=[C; CA; CA^2;...; CA^(n-1) ]
o = obsv_mat(sys.A, sys.C);
rankO=rank(o);
if j == rankO then
    disp('Continuous System is Observable');
end

tranM=ss2tf(sys); // Matriz de transferencia
disp('Matriz de Transferencia',tranM);

//tfc11 = tranM(1,1);
//tfc22 = tranM(2,2); 

/* Plot singular values of LTI the model */
tr = trzeros(sys)
w = logspace(-3,3);
sv = svplot(sys,w);
scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot of the system","Frequency (rad/s)", "Amplitude (dB)");

//Obtencion de los polors  zeros del modelo de software
scf(2);
plzr(sys);
xtitle("Poles and zeros plot of system","Real", "Imaginarie");
///////////////////////////////////---------------------------

//Augment Plant with Integrators at Plant Input
[ns,nc]=size(Bp); //ns= number of inputs; nc=number of control
Ai=[Ap             Bp;
    0*ones(nc,ns) 0*ones(nc,nc)];

Bi=[0*ones(ns,nc); eye(nc)];
    
Ci=[Cp 0*ones(2,nc)];

Di=0*ones(2,nc);

sysi=syslin('c',Ai,Bi,Ci,Di);

I=eye(2,2);

/* Plot singular values  */
tri = trzeros(sysi)
w = logspace(-3,3);
svi = svplot(sysi,w);
scf(3);
plot2d("ln", w, 20*log(svi')/log(10))
xgrid(12)
xtitle("Design Plant with integrator:Singular Values","Frequency (rad/s)", "Amplitude (dB)");

//Obtenciion de los polors  zeros del modelo de software
scf(4);
plzr(sysi);
xtitle("Design Plant with integrator:poles and zeros","Real", "Imaginarie");

//lqr controller calculation
//We use the ricatti equation for calculate de gain of the lqr controller
//for this we have  A'*X+X*A-X*B*X+C=0 for function X=riccati(A,B,C,'c','eigen')
C=0.9*Ci'*Ci;        //State Weighting Matrix
rho=2e-0;       //Cheap control recovery parameter 
                //The smaller the parameter, the better the recovery.
R = rho*eye(nc);//Control Weigthing Matrix


//now we calculate B
B=Bi*inv(R)*Bi';

A=Ai;

//Solv the ricatti equation
X=riccati(A,B,C,'c','eigen');

//the value of the gain G
G=inv(R)*Bi'*X; //<--this value is important mtfk

//[G1, X1]=lqr(sysi,C,R);

//---------------------------------------------------------------------
//Design of Target Loop Singular Values Using Kalman Filter

//computing H

//H=(ppol(Ai',Ci',[-2.5,-1.5,-1,-0.5,-1]))';  
H=(ppol(Ai',Ci',[-2.5,-1.5,-1,-0.5,-1]))';  

//-----------------------------------------------------------------------
//Compensador
Ak = [ Ai-Bi*G-H*Ci  0*ones(ns+nc,nc)
       G          0*ones(nc,nc) ]

Bk = [ H
       0*ones(1,2) ]

Ck = [0*ones(nc, ns+nc) eye(nc,nc) ]

Dk = [0*ones(1,2)]

sysk=syslin('c',Ak,Bk,Ck,Dk);

//Discretization------------------------------------
T=1/9600;
sysdk=cls2dls(sysk,T);

sysdktf=ss2tf(sysdk);

disp('Compensador Discretizado',sysdktf);
//compensador Discretizado
Adk=sysdk.A;
Bdk=sysdk.B;
Cdk=sysdk.C;
Ddk=sysdk.D;
//Planta discretizada

sysdp=cls2dls(sys,1/9600);
Adp=sysdp.A;
Bdp=sysdp.B;
Cdp=sysdp.C;
Ddp=sysdp.D;

//----------------------------------------
//Open loop 
Aol = [ Adp                     Bdp*Cdk
       0*ones(ns+nc+nc,ns)    Adk    ]

Bol = [ 0*ones(ns,2)
       Bdk ]
    
Col = [ Cdp  0*ones(2,ns+nc+nc) ]

Dol = 0*ones(2,2);

sysol = syslin('d',Aol,Bol,Col,Dol);


//----------------------------------------
//Response in closed loop
syscl = syslin('d',Aol-Bol*Col, Bol, Col, 0*eye(2,2));

//Obtenciion de los polors  zeros del modelo de software
scf(5);
plzr(syscl);

//----------Response to step
n=[0:T:30];
//input defined by a time function
tranclM=ss2tf(syscl);
tfcl11 = tranclM(1,1);
tfcl21 = tranclM(2,1); 

deff(['[u1,u2]=timefun(t)'],['u1=0','u2=0']) //*%pi/180
scf(6);
//plot2d(t',(csim(zeros(t),syscl,eye(10,1)))')
y1 = flts([zeros(n);zeros(n)],syscl,ones(10,1));
plot2d(n',y1')
e=gce();e.children.polyline_style=2;
xtitle("Response of the model for 1st reference","Amplitud", "t(s)");
legend('Response for Position','Response for Angle')

