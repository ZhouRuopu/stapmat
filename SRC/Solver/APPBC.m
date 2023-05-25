function APPBC()

global sdata;
global cdata;

BC = sdata.BC; NDOF = sdata.NDOF;
STIFF = sdata.STIFF; R = sdata.R;
NUMNP = cdata.NUMNP;

inf = 1e16;
for N = 1:NUMNP
    for dof = 1:NDOF
        if (BC(dof,N) == 1)
           STIFF(6*N-6 + dof, 6*N-6 + dof) = inf * STIFF(6*N-6 + dof, 6*N-6 + dof);
           R(6*N-6 + dof) = 0;
        end
    end
end

sdata.STIFF = STIFF; sdata.R = R;