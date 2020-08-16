//valores Linealizados analiticamente caso del 1er lab

// Load the parameters
exec('pendulojParametersLab1.sce', -1);


Aa=[0 0 1 0;...
    0 0 0 1;...
    0 ((alpha^2*g)/(bet)) -((THETA*b)/(bet)) 0;...
    0 ((alpha*M*g)/(bet)) -((alpha*b)/(bet)) 0];
   
Ba=[0;...
   0;...
   THETA/bet;...
   alpha/bet];
   
Ca=[1 0 0 0;
    0 1 0 0];
  
Da=[0;
   0];
   
sysa=syslin('c', Aa, Ba, Ca, Da);
// save the data
save("pendulojLTIanaliticoLab1.sod","sysa")
