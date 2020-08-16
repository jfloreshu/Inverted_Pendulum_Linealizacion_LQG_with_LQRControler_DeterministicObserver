// linearizing.sce

// steps before
// 1. Open edsonj.m and execute from editor "Execute > file with echo".
// 2. Open edsonjXcosLincos.zcos run the simulation using "srtat" button.
// 3. Then run this file

// Search the SUPERBLOCK in Xcos
for i=1:length(scs_m.objs)
    if typeof(scs_m.objs(i))=="Block" & scs_m.objs(i).gui=="SUPER_f" then
        scs_m = scs_m.objs(i).model.rpar;
        break;
    end
end

// Set the equilibrium point, in this parto we use x=0 and theta=pi
X=[0.001;0.001;0.001;0.001];
U=[0.001];

// linearize the model
sys = lincos(scs_m,X,U);

// obtaingin the matrices A,B,C,D
A=sys.A
B=sys.B
C=sys.C
D=sys.D

// save the data
save("pendulojLTILab1.sod","X","U","sys")

// load the data
//load("edsonjLTI.sod","X","U","sys")
