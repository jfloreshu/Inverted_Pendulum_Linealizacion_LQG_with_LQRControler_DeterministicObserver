//Parameters of Pendulum to Labo1

M_0=0.696 //Mass of the cart(Kg)
M_1=0.017  //Pedulum mass(Kg)
b=0.001 //Friction coefficient(Ns/m)
I=0.0011 //Inertia of the pendulum(Kgm2)
g=9.8 //Grvitational acceleration(m/s2)
l_s=0.3 //Pendulum length(m)
C=0.00 //Friccion seca o de Coulomb

//asumciones del laboratorio1 para reducir 
p=I*( M_0+M_1)+M_0*M_1*l_s^2;

THETA= I+ M_1*l_s^2;
alpha=M_1*l_s;
M=M_0+M_1;

bet=THETA.*M-alpha.^2;
