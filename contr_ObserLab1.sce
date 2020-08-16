// load the data
load("pendulojLTILab1.sod","X","U","sys")
// load the data
load("pendulojLTIanaliticoLab1.sod","sysa")

// obtaingin the matrices A,B,C,D
A=sys.A
B=sys.B
C=sys.C
D=sys.D

Aa=sysa.A
Ba=sysa.B
Ca=sysa.C
Da=sysa.D

// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(A,B)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(A, C)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a
/* Plot singular values of LTI the model */
G = syslin('c', A, B, C, D);
tr = trzeros(G)
w = logspace(-3,3);
sv = svplot(G,w);

//valores analiticos
Ga = syslin('c', Aa, Ba, Ca, Da);
tra = trzeros(Ga)
wa = logspace(-3,3);
sva = svplot(Ga,wa);
scf(1);
plot2d("ln", w, [20*log(sv')/log(10) 20*log(sva')/log(10)],leg="Caso Soft@Caso Analitico")
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
/* Scaling */

//Obtenciion de los polors  zeros del modelo de software
scf(2);
plzr(sys);

//Obtencion de los polos  zeros del modelo analitico
scf(3);
plzr(sysa);

//Obtencion de las funciones de transferencia
[h]=ss2tf(sys)

[ha]=ss2tf(sysa)
