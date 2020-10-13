function [S1r,S2r] = smpropag_fw(S1,S2,p1,p2,A1,A2,f1,f2)

% Forward propagation of scattering matrix. 

X = p1\p2; Y = A1\A2;
I = (X+Y)/2;
J = (X-Y)/2;

S1r = (I - diag(f1)*S2*J)\(diag(f1)*S1);
S2r = (I - diag(f1)*S2*J)\(diag(f1)*S2*I - J)*diag(f2);

end