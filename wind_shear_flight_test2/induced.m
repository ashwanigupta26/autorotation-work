function f = induced(lambda,nu,alpha,Ct)
f=lambda - nu.*tan(alpha) - Ct./(2*sqrt(lambda^2 + nu^2));
end