set encoding utf8

#x = k
set yrange [-1.1:3]
set xrange [-300:300]
#set ytics 0 font ",13"
#set xtics 0 font ",13" 
set xlabel "kD" font ",12"# x軸のラベルを設定
set ylabel "C/U" font ",12" # y軸のラベルを設定
R_a=0.7
R_b=1.0
R_c=1.3



#density ratio
den1 = 0.1

#d2/D
d2 = 0.6
d1 = 1-d2

#U2/U
U2_a = (1 + d2*(den1 - 1))/(1 + d2*(R_a*den1 - 1))
U2_b = (1 + d2*(den1 - 1))/(1 + d2*(R_b*den1 - 1))
U2_c = (1 + d2*(den1 - 1))/(1 + d2*(R_c*den1 - 1))

#U1/U
U1_a = R_a*U2_a
U1_b = R_b*U2_b
U1_c = R_c*U2_c


#delta_rho_g/U^2
rho_g_a=(U1_a**2)*(den1/(1-d2))+(U2_a**2)*(1/d2)
rho_g_b=(U1_b**2)*(den1/(1-d2))+(U2_b**2)*(1/d2)
rho_g_c=(U1_c**2)*(den1/(1-d2))+(U2_c**2)*(1/d2)

#T1(x)=T1/ρ2
T1(x)=den1*x/tanh(x*(1-d2))

#T2(x)=T2/ρ2
T2(x)=x/tanh(x*d2)

#d(T1)/dk
d_T1(x)=den1*(coth(x*d1) - x*d1*((csch(x*d1))**2))/d1
#d(T2)/dk
d_T2(x)=(coth(x*d2) - x*d2*((csch(x*d2))**2))/d2

#A
A(x)=T1(x)+T2(x)

d_A(x)=d_T1(x)+d_T2(x)

#coth(x)
coth(x)=cosh(x)/sinh(x)

#csch(x)=1/sinh(x)
csch(x)=1/sinh(x)


#Cp(x)=Cp(k)/U
Cp_a(x) = (U1_a*T1(x)+U2_a*T2(x),\
    -sqrt((T1(x)+T2(x))*rho_g_a-T1(x)*T2(x)*(U1_a-U2_a)**2))/(T1(x)+T2(x))
Cp_b(x) = (U1_b*T1(x)+U2_b*T2(x),\
    -sqrt((T1(x)+T2(x))*rho_g_b-T1(x)*T2(x)*(U1_b-U2_b)**2))/(T1(x)+T2(x))
Cp_c(x) = (U1_c*T1(x)+U2_c*T2(x),\
    -sqrt((T1(x)+T2(x))*rho_g_c-T1(x)*T2(x)*(U1_c-U2_c)**2))/(T1(x)+T2(x))

Cp_0(x) = (U1_b*(den1/(1-d2))+U2_b*(1/d2),\
    -sqrt((den1/(1-d2)+(1/d2))*((U1_b**2)*(den1/(1-d2))+(U2_b**2)*(1/d2))-(den1/(1-d2))*(1/d2)*(U1_b-U2_b)**2))/(den1/(1-d2)+(1/d2))

#Cg(x)=Cg(k)/U
Cg_a(x) = Cp_a(x) + x*((U1_a*d_T1(x)+U2_a*d_T2(x))/A(x)\
            -(U1_a*T1(x)+U2_a*T2(x))*d_A(x)/A(x)**2\
            -(d_A(x)*rho_g_a - (d_T1(x)*T2(x) + T1(x)*d_T2(x))*((U1_a-U2_a)**2))\
            /(2*A(x)*sqrt(A(x)*rho_g_a - T1(x)*T2(x)*((U1_a-U2_a)**2)))\
            +d_A(x)*(A(x)*rho_g_a - T1(x)*T2(x)*((U1_a-U2_a)**2))/A(x)**2)
Cg_b(x) = Cp_b(x) + x*((U1_b*d_T1(x)+U2_b*d_T2(x))/A(x)\
            -(U1_b*T1(x)+U2_b*T2(x))*d_A(x)/A(x)**2\
            -(d_A(x)*rho_g_b - (d_T1(x)*T2(x) + T1(x)*d_T2(x))*((U1_b-U2_b)**2))\
            /(2*A(x)*sqrt(A(x)*rho_g_b - T1(x)*T2(x)*((U1_b-U2_b)**2)))\
            +d_A(x)*(A(x)*rho_g_b - T1(x)*T2(x)*((U1_b-U2_b)**2))/A(x)**2)
Cg_c(x) = Cp_c(x) + x*((U1_c*d_T1(x)+U2_c*d_T2(x))/A(x)\
            -(U1_c*T1(x)+U2_c*T2(x))*d_A(x)/A(x)**2\
            -(d_A(x)*rho_g_c - (d_T1(x)*T2(x) + T1(x)*d_T2(x))*((U1_c-U2_c)**2))\
            /(2*A(x)*sqrt(A(x)*rho_g_c - T1(x)*T2(x)*((U1_c-U2_c)**2)))\
            +d_A(x)*(A(x)*rho_g_c - T1(x)*T2(x)*((U1_c-U2_c)**2))/A(x)**2)
    
        

set xrange [-10:10] # k=0付近を広げる
set samples 1000 
set key outside
set key spacing 1.5
plot Cp_0(x) title "Cp(k=0,R=1)" w l lc 'black' lw 2 smooth csplines,\
    Cg_a(x) title "Cg(k,R=0.7)" w l lc 'blue'lw 2 smooth csplines,\
    Cg_b(x) title "Cg(k,R=1.0)" w l lc 'green' lw 2 smooth csplines,\
    Cg_c(x) title "Cg(k,R=1.3)" w l lc 'red' lw 2 smooth csplines,\
    Cp_a(x) title "Cp(k,R=0.7)" w l lc 'blue' lw 2 dt 2 smooth csplines,\
    Cp_b(x) title "Cp(k,R=1.0)" w l lc 'green' lw 2 dt 2 smooth csplines,\
    Cp_c(x) title "Cp(k,R=1.3)" w l lc 'red' lw 2 dt 2 smooth csplines

