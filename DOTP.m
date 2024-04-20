clear all
clc


disp("============================================================================================")
disp("Design Optimization of Telecommunication Pole")
disp("Liam R. Spring | University of Connecticut <3 | Department of Mechanical Engineering")
disp("Advisor: Prof. Matthew D. Stuber")
disp("============================================================================================")
disp(" ")



gamma = 78500;

f = @(x)(gamma * fun1(x));      %objective function (weight function)

x0 = [0.4,0.005,0.02];          %x = [d_t,t,tau]


J = zeros(2,1);

J(1) = 463.1;          %alpha
J(2) = 0.25;           %beta


nonlcon = @(x)con_nl(x,J);

lb = [0.3;0.0032;0];
ub = [1;0.0254;0.05];


options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
[x,fval,exitflag,output,lambda] = fmincon(f,x0,[],[],[],[],lb,ub,nonlcon,options);



disp("============================================================================================")
disp(" ")

disp("Optimum solution in the form x* = [d_t[m],t[m],tau[m/m]]:")

format shortG

x
disp("-----------------------------------------------")
disp(" ")
disp("Cost function value at x* in N:")

fval

disp("-----------------------------------------------")
disp(" ")
disp("Lagrange multipliers:")
mu = lambda.ineqnonlin

disp("============================================================================================")
disp(" ")
SD(x,J)



function [c,ceq] = con_nl(x,J)    %nonlinear constraints: home to goofy equations

d_t = x(1);
t = x(2);
tau = x(3);

A_lc = 0.3;
A_t = 10;
C_lc = 1;
C_p = 0.75;
C_t = 1;
E = 210e9;
F_v = 10400;
H = 30;
p_lc = 400;
gamma = 78500;
sig_a = 150e6;
v_a = 0.6;

alpha = J(1);
beta = J(2);

q = @(z)(alpha*z.^beta);
F_h = C_t.*A_t.*q(H);
A = @(y)(pi/4).*((d_t + 2.*tau*(H-y).^2)-(d_t + 2.*tau.*(H-y)-2.*t).^2);
p_v = @(y)(p_lc + gamma*A(y));                  
N = @(z)F_v + integral(p_v,0,H-z);
p_h = @(y)((A_lc.*C_lc + (d_t + 2.*tau.*(y)).*C_p)).*q(H-y);
core1 = @(y)(p_h(y).*y);
core2 = @(y)(p_h(y));
M = @(z)((F_h.*(H-z)) - integral(core1,0,H-z) + H*integral(core2,0,H-z) - z*integral(core2,0,H-z));
I = @(z)((pi/64).*((d_t + 2.*tau*(H-z)).^4 - (d_t + 2.*tau.*(H-z)-2.*t).^4));
S = @(z)((2.*I(z))/(d_t + 2.*tau*(H-z)));
sig = @(z)((N(z)/A(z))+(M(z)/S(z)));

def  = @(z,v)[v(2);(-M(z)/(E*I(z)))];    %Solving ODE:  EIv"(z) = -M(x)  w/ I.C.s: v(0) = v'(0) = 0        
[t,v] = ode45(def,[0 H],[0; 0]);

%nonlinear equality constriants
ceq = [];

%nonlinear inequality constriants
c = zeros(2,1);

c(1) = sig(0)/sig_a - 1;
c(2) = abs(v(end,1))/v_a - 1;

%plot(t,v(:,1))

end



function SD = SD(x,J)

d_t = x(1);
t = x(2);
tau = x(3);

A_lc = 0.3;
A_t = 10;
C_lc = 1;
C_p = 0.75;
C_t = 1;
E = 210e9;
F_v = 10400;
H = 30;
p_lc = 400;
gamma = 78500;

alpha = J(1);
beta = J(2);

q = @(z)(1*(alpha*z.^beta));
F_h = C_t.*A_t.*q(H);
A = @(y)(pi/4).*((d_t + 2.*tau*(H-y).^2)-(d_t + 2.*tau.*(H-y)-2.*t).^2);
p_v = @(y)(p_lc + gamma*A(y));                  
N = @(z)F_v + integral(p_v,0,H-z);
p_h = @(y)((A_lc.*C_lc + (d_t + 2.*tau.*(y)).*C_p)).*q(H-y);
core1 = @(y)(p_h(y).*y);
core2 = @(y)(p_h(y));
M = @(z)((F_h.*(H-z)) - integral(core1,0,H-z) + H*integral(core2,0,H-z) - z*integral(core2,0,H-z));
I = @(z)((pi/64).*((d_t + 2.*tau*(H-z)).^4 - (d_t + 2.*tau.*(H-z)-2.*t).^4));
S = @(z)((2.*I(z))/(d_t + 2.*tau*(H-z)));
sig = @(z)((N(z)/A(z))+(M(z)/S(z)));

def  = @(z,v)[v(2);(-M(z)/(E*I(z)))];        
[t,v] = ode45(def,[0 H],[0; 0]);



format shortG

disp("Stress at base in Pa:")
disp(" ")
disp(sig(0))
disp("-----------------------------------------------")
disp(" ")
disp("Tip deflection, v(30) in m:")
disp(" ")
disp(abs(v(end,1)));
disp("============================================================================================")

end



function calc1 = fun1(x)      %A(y) function

H = 30;

d_t = x(1);
t = x(2);
tau = x(3);

A = @(y)(pi/4*((d_t + 2.*tau*(H-y)).^2 - ((d_t + 2.*tau.*(H-y))-2*t).^2));

calc1 = integral(A,0,H);

end


