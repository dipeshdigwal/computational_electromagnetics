lambda=linspace(400e-9,700e-9,100)
y=LD(lambda,'Ag','LD')
plot(lambda,real(y))
hold on
plot(lambda,imag(y))
hold off