function n = refractive_index_BK7(lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
B1=1.03961212;
B2=0.231792344;
B3=1.01046945;
C1=0.00600069867; %in micrometer-square
C2=0.0200179144;  %in micrometer-square
C3=103.560653;    %in micrometer-square

P1=(B1*(lambda^2))/((lambda^2)-C1);
P2=(B2*(lambda^2))/((lambda^2)-C2);
P3=(B3*(lambda^2))/((lambda^2)-C3);
n=sqrt(1+P1+P2+P3);
end

