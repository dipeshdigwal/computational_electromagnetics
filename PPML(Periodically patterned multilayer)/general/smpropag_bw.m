function [S3r,S4r] = smpropag_bw(S3,S4,p1,p2,A1,A2,f1,f2)

% Backward propagation of scattering matrix. 

X = p1\p2; Y = A1\A2;
L = (X+Y)/2;
M = (X-Y)/2;

S3r = (L - diag(f1)*S3*M)\(diag(f1)*S3*L - M)*diag(f2);
S4r = (L - diag(f1)*S3*M)\(diag(f1)*S4);

end