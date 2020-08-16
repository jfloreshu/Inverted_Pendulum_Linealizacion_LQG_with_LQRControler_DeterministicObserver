function xdot=pendulojLab1(u1,u2,u3,u4,u5)
// This is nonlinear PenduluModel

// Load the parameters
exec('pendulojParametersLab1.sce', -1);

// state variables
x=u1;		
thet=u2;
v=u3;
vartheta=u4;

// control variables
F=u5;	// fuerza


// Modelo MoDiCA-X
//Estado 1-posicion
e1dot= v;
//Estado 2-Angulo theta
e2dot= vartheta;
//Estado 3-Velocidad lineal respecto a 1
e3dot= ((alpha^2*g*(sin(thet))*(cos(thet)))/(THETA*M-alpha^2*(cos(thet))^2))-...
            ((b*THETA*v)/(THETA*M-alpha^2*(cos(thet))^2))-...
        ((alpha*C*cos(thet)*vartheta)/(THETA*M-alpha^2*(cos(thet))^2))-...
            ((alpha*sin(thet)*vartheta^2)/(THETA*M-alpha^2*(cos(thet))^2))+...
        ((THETA*F)/(THETA*M-alpha^2*(cos(thet))^2));

//Estado 4-Velocidad Angular respecto de 3
e4dot= ((alpha*M*g*sin(thet))/(THETA*M-alpha^2*(cos(thet))^2))-...
            ((alpha*b*cos(thet)*v)/(THETA*M-alpha^2*(cos(thet))^2))-...
        ((M*C*vartheta)/(THETA*M-alpha^2*(cos(thet))^2))-...
 ((alpha^2*sin(thet)*cos(thet)*vartheta^2)/(THETA*M-alpha^2*(cos(thet))^2))+...
        ((alpha*cos(thet)*F)/(THETA*M-alpha^2*(cos(thet))^2));

//Salida xdot
xdot =[e1dot;e2dot;e3dot;e4dot];

endfunction
