x0 = [5, 0.25, 3];
fun = @(x) mainF(1,x);
options = optimoptions('fmincon','UseParallel',true);
[x, fval] = fmincon(fun, x0, [], [], [], [], [], [], [], options)