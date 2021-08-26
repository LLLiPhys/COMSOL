function f=fermi(x)
f=1.0.*(x<=-200)+1./(exp(x)+1).*(abs(x)<200)+0.0.*(x>=200);
end