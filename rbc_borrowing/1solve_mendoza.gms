*--------------------------------------------------------
* Simple RBC
*--------------------------------------------------------
scalar starttime; starttime = jnow;

Parameters
beta_true,gamma_true,alphak_true,delta_true,rhoz_true,sigmaz_true, a_true, rhor_true, sigmar_true
;
beta_true=0.98;
gamma_true=0.75;
alphak_true=0.36;
delta_true=0.025;
rhoz_true=0;
rhor_true=0;
sigmaz_true=0.04;
a_true=0;
*rhor_true=0.99;
sigmar_true=0.003;


Parameters
beta         'discount rate'
gamma        'risk preference'
alphak       'capital share'
delta        'depreciation'
rhoz         'TFP persistence'
sigmaz       'std TFP'
a            'capital adjustment cost'
rhor         'interest rate persistence'
sigmar       'interest rate std'
kbar         'steady state capital'
ybar         'steady state production'
cbar         'steady state consumption'
ibar         'steady state investment'
bbar         'steady state bonds'
sbar
rbar         'average interest rate'
alphar1      'stochastic constant of the interest rate paid by the agent'
alphar2      'slope of the interest rate function (polynomial of degree 2)'
;

beta=beta_true;
gamma=gamma_true;
alphak=alphak_true;
delta=delta_true;
rhoz=rhoz_true;
rhor=rhor_true;
sigmaz=sigmaz_true;
a=a_true;
*rhor = rhor_true;
sigmar = sigmar_true;
kbar=1;
ybar=(1/beta+delta-1)/alphak*kbar**alphak ;
ibar = kbar*delta;
bbar = 1 ;
rbar = 1.015 ;
alphar2 = (1/beta-rbar)/2 ;
alphar1 = rbar-alphar2*bbar*bbar ;
cbar=bbar*(1/rbar-1) - ibar+ybar;
sbar = ybar - ibar - cbar ;
display ybar, ibar, cbar, sbar, bbar, alphar1, alphar2;


*---------------------------------------------------------------------------
* Solution Step 1: Construct a tensor product grid for computing a solution
*---------------------------------------------------------------------------
sets nk  number of nodes for capital /1*10/ ;
sets nz  number of nodes for productivity shock / 1*10/ ;
sets nb  number of nodes for bonds /1*10/;
sets nr  number of nodes for interest rate shock /1*10/ ;
alias (nk,n_gridk) ;
alias (nz,n_gridz)  ;

* Gaussion Herminte quadrature std=1
Set Qn 'Number of integration nodes in Gauss Hermite quadrature' /1*5/;
alias(Qn,nzb,nrb);
Parameters
epsi_nodes(Qn)
weight_nodes(Qn)
;

epsi_nodes("1")= -2.85697001387281; epsi_nodes("2")= -1.35562617997427; epsi_nodes("3")= 0; epsi_nodes("4")= 1.35562617997427; epsi_nodes("5")= 2.85697001387281;
weight_nodes("1")=0.0112574113277207;weight_nodes("2")=0.222075922005613; weight_nodes("3")=0.533333333333333;weight_nodes("4")=0.222075922005613 ;weight_nodes("5")=0.0112574113277207 ;


Parameters
kmin
kmax
bmin
bmax
zexpmin
zexpmax
rshockmin
rshockmax
k0(nk)
z0exp(nz)
z0(nz)
rshock0(nr)
r0(nr)
b0(nb)
;

* Unidimensional grid
kmin = kbar*0.8;
kmax = kbar*1.2;
k0(nk)  = kmin + (kmax-kmin)/(card(nk)-1) * (ord(nk)-1) ;

bmin = bbar*0.8;
bmax = bbar*1.2;
b0(nb)  = bmin + (bmax-bmin)/(card(nb)-1) * (ord(nb)-1) ;

zexpmin = 0.8;
zexpmax = 1.2;
z0exp(nz) = zexpmin + (zexpmax-zexpmin)/(card(nz)-1) * (ord(nz)-1) ;
z0(nz) = log(z0exp(nz));


*rmin = rbar*0.99 ;
*rmax = rbar*1.01 ;
*r0(nr)  = rmin + (rmax-rmin)/(card(nr)-1) * (ord(nr)-1) ;

rshockmin = 0.99;
rshockmax = 1.01;
rshock0(nr)  = rshockmin + (rshockmax-rshockmin)/(card(nr)-1) * (ord(nr)-1) ;
rshock0(nr)=log(rshock0(nr)) ;
r0(nr)=alphar1*exp(rshock0(nr)) ;

Parameter
z1(nz,nzb);
z1(nz,nzb) = rhoz*z0(nz)+sigmaz*epsi_nodes(nzb);


Parameter
r1(nr,nrb);
r1(nr,nrb) = alphar1*exp(rhor*rshock0(nr)+sigmar*epsi_nodes(nrb));

display r0, r1;




*if 3rd order, ncoeff=35
Sets ncoeff number of coefficients   /1*35/;

Parameters
y0ini(nk,nb,nz,nr)
c0ini(nk,nb,nz,nr)
i0ini(nk,nb,nz,nr)
vf(nk,nb,nz,nr)  ;

y0ini(nk,nb,nz,nr) = (1/beta+delta-1)/alphak * exp(z0(nz)) * k0(nk)**alphak ;
c0ini(nk,nb,nz,nr) = y0ini(nk,nb,nz,nr) * cbar/ybar ;
i0ini(nk,nb,nz,nr) = y0ini(nk,nb,nz,nr) * ibar/ybar ;

Variables
coeff(ncoeff)
noise(nk,nb,nz,nr)
obj
;


equation
eqvf(nk,nb,nz,nr)
eqobj
;

eqvf(nk,nb,nz,nr)..           vf(nk,nb,nz,nr)  =E=   coeff("1")
                                  + coeff("2")*k0(nk) + coeff("3")*b0(nb)+ coeff("4")*exp(z0(nz))+ coeff("5")*r0(nr)
                                  + coeff("6")*k0(nk)*k0(nk)+ coeff("7")*k0(nk)*b0(nb)+ coeff("8")*k0(nk)*exp(z0(nz))+ coeff("9")*k0(nk)*r0(nr)
                                  + coeff("10")*b0(nb)*b0(nb)+ coeff("11")*b0(nb)*exp(z0(nz))+ coeff("12")*b0(nb)*r0(nr)
                                  + coeff("13")*exp(z0(nz))*exp(z0(nz))+ coeff("14")*exp(z0(nz))*r0(nr)
                                  + coeff("15")*r0(nr)*r0(nr)
*                                  + coeff("16")*k0(nk)*k0(nk)*k0(nk)+ coeff("17")*k0(nk)*k0(nk)*b0(nb)+ coeff("18")*k0(nk)*k0(nk)*exp(z0(nz))+ coeff("19")*k0(nk)*k0(nk)*r0(nr)
*                                  + coeff("20")*k0(nk)*b0(nb)*b0(nb)+ coeff("21")*k0(nk)*b0(nb)*exp(z0(nz))+ coeff("22")*k0(nk)*b0(nb)*r0(nr)
*                                  + coeff("23")*k0(nk)*exp(z0(nz))*exp(z0(nz))+ coeff("24")*k0(nk)*exp(z0(nz))*r0(nr)
*                                  + coeff("25")*k0(nk)*r0(nr)*r0(nr)
*                                  + coeff("26")*b0(nb)*b0(nb)*b0(nb)+ coeff("27")*b0(nb)*b0(nb)*exp(z0(nz))+ coeff("28")*b0(nb)*b0(nb)*r0(nr)
*                                  + coeff("29")*b0(nb)*exp(z0(nz))*exp(z0(nz))+ coeff("30")*b0(nb)*exp(z0(nz))*r0(nr)
*                                  + coeff("31")*b0(nb)*r0(nr)*r0(nr)
*                                  + coeff("32")*exp(z0(nz))*exp(z0(nz))*exp(z0(nz))+ coeff("33")*exp(z0(nz))*exp(z0(nz))*r0(nr)
*                                  + coeff("34")*exp(z0(nz))*r0(nr)*r0(nr)
*                                  + coeff("35")*r0(nr)*r0(nr)*r0(nr)
                                  + noise(nk,nb,nz,nr)

                                  ;

EQobj..                 obj =E= sum((nk,nb,nz,nr), noise(nk,nb,nz,nr)*noise(nk,nb,nz,nr));

Model leastsquare/
eqvf
eqobj
/
;

Parameters
mu              the severity of the borrowing constraint
coeffini(ncoeff)  ;
* we assume that the capital is valuable, not debt
coeffini("2") = 1 ;
coeffini("3") = -1 ;
mu =10 ;

Positive variable
k1in(nk,nb,nz,nr)
b1in(nk,nb,nz,nr)
c0(nk,nb,nz,nr)
i0(nk,nb,nz,nr)
;
Variables
obje(nk,nb,nz,nr)
obj1
;

equations
eqobj1
eqk1in(nk,nb,nz,nr)
eqb1in(nk,nb,nz,nr)
eqobje(nk,nb,nz,nr)
eqconstraint(nk,nb,nz,nr)
;



eqk1in(nk,nb,nz,nr)..       k1in(nk,nb,nz,nr) =E=  (1-delta)*k0(nk)+i0(nk,nb,nz,nr) ;
eqb1in(nk,nb,nz,nr)..       b1in(nk,nb,nz,nr) =E= ( -(1/beta+delta-1)/alphak * exp(z0(nz)) * k0(nk)**alphak
                                                      +c0(nk,nb,nz,nr)+i0(nk,nb,nz,nr)    + a/2*(k1in(nk,nb,nz,nr)-k0(nk))*(k1in(nk,nb,nz,nr)-k0(nk))/k0(nk)
                                                  )+b0(nb)
                                                  *(r0(nr) + alphar2*b0(nb)*b0(nb) )  ;

eqconstraint(nk,nb,nz,nr)..   b1in(nk,nb,nz,nr) =L= mu*k1in(nk,nb,nz,nr)*(r0(nr) + alphar2*b0(nb)*b0(nb) )  ;

eqobje(nk,nb,nz,nr)..         obje(nk,nb,nz,nr) =E=
                                         c0(nk,nb,nz,nr)**(1-gamma) / (1-gamma)
                                        + beta*sum( (nzb,nrb),weight_nodes(nzb)*weight_nodes(nrb)*
                                                              (  coeffini("1")
                                  + coeffini("2")*k1in(nk,nb,nz,nr)  + coeffini("3")*b1in(nk,nb,nz,nr) + coeffini("4")*exp(z1(nz,nzb))+ coeffini("5")*r1(nr,nrb)
                                  + coeffini("6")*k1in(nk,nb,nz,nr) *k1in(nk,nb,nz,nr) + coeffini("7")*k1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) + coeffini("8")*k1in(nk,nb,nz,nr) *exp(z1(nz,nzb))+ coeffini("9")*k1in(nk,nb,nz,nr) *r1(nr,nrb)
                                  + coeffini("10")*b1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) + coeffini("11")*b1in(nk,nb,nz,nr) *exp(z1(nz,nzb))+ coeffini("12")*b1in(nk,nb,nz,nr) *r1(nr,nrb)
                                  + coeffini("13")*exp(z1(nz,nzb))*exp(z1(nz,nzb))+ coeffini("14")*exp(z1(nz,nzb))*r1(nr,nrb)
                                  + coeffini("15")*r1(nr,nrb)*r1(nr,nrb)
*                                  + coeffini("16")*k1in(nk,nb,nz,nr) *k1in(nk,nb,nz,nr) *k1in(nk,nb,nz,nr) + coeffini("17")*k1in(nk,nb,nz,nr) *k1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) + coeffini("18")*k1in(nk,nb,nz,nr) *k1in(nk,nb,nz,nr) *exp(z1(nz,nzb))+ coeffini("19")*k1in(nk,nb,nz,nr) *k1in(nk,nb,nz,nr) *r1(nr,nrb)
*                                  + coeffini("20")*k1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) + coeffini("21")*k1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) *exp(z1(nz,nzb))+ coeffini("22")*k1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) *r1(nr,nrb)
*                                  + coeffini("23")*k1in(nk,nb,nz,nr) *exp(z1(nz,nzb))*exp(z1(nz,nzb))+ coeffini("24")*k1in(nk,nb,nz,nr) *exp(z1(nz,nzb))*r1(nr,nrb)
*                                  + coeffini("25")*k1in(nk,nb,nz,nr) *r1(nr,nrb)*r1(nr,nrb)
*                                  + coeffini("26")*b1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) + coeffini("27")*b1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) *exp(z1(nz,nzb))+ coeffini("28")*b1in(nk,nb,nz,nr) *b1in(nk,nb,nz,nr) *r1(nr,nrb)
*                                  + coeffini("29")*b1in(nk,nb,nz,nr) *exp(z1(nz,nzb))*exp(z1(nz,nzb))+ coeffini("30")*b1in(nk,nb,nz,nr) *exp(z1(nz,nzb))*r1(nr,nrb)
*                                  + coeffini("31")*b1in(nk,nb,nz,nr) *r1(nr,nrb)*r1(nr,nrb)
*                                  + coeffini("32")*exp(z1(nz,nzb))*exp(z1(nz,nzb))*exp(z1(nz,nzb))+ coeffini("33")*exp(z1(nz,nzb))*exp(z1(nz,nzb))*r1(nr,nrb)
*                                  + coeffini("34")*exp(z1(nz,nzb))*r1(nr,nrb)*r1(nr,nrb)
*                                  + coeffini("35")*r1(nr,nrb)*r1(nr,nrb)*r1(nr,nrb)
                                  )
                                                    )
                                       ;

eqobj1..                     obj1 =E= sum((nk,nb,nz,nr), obje(nk,nb,nz,nr) ) ;

Model findc0i0/
eqk1in
eqb1in
eqobje
eqobj1
eqconstraint
/;

c0.l(nk,nb,nz,nr)=c0ini(nk,nb,nz,nr);
c0.lo(nk,nb,nz,nr) = 0.001 ;
i0.fx(nk,nb,nz,nr)=i0ini(nk,nb,nz,nr);
k1in.up(nk,nb,nz,nr) = kmax*1.2 ;
k1in.lo(nk,nb,nz,nr) = kmin*0.8 ;
b1in.up(nk,nb,nz,nr) = bmax*1.2 ;
b1in.lo(nk,nb,nz,nr) = bmin*0.8 ;


solve findc0i0 using nlp maximizing obj1;

i0.lo(nk,nb,nz,nr) = -inf ;
i0.l(nk,nb,nz,nr) = i0ini(nk,nb,nz,nr) ;
i0.up(nk,nb,nz,nr) = +inf ;

sets iterr/1*2000 / ;
Parameters
stopiter
diffvf(iterr)
vfold(nk,nb,nz,nr)
vfnew(nk,nb,nz,nr)
;

stopiter = 0 ;
loop(iterr$(stopiter EQ 0),

vf(nk,nb,nz,nr) = obje.l(nk,nb,nz,nr) ;
vfold(nk,nb,nz,nr) = obje.l(nk,nb,nz,nr) ;
solve leastsquare minimizing obj using nlp;

coeffini(ncoeff) = coeff.l(ncoeff) ;
coeffini(ncoeff)=0.5*coeff.l(ncoeff)+0.5*coeffini(ncoeff);
solve findc0i0 using nlp maximizing obj1;
vfnew(nk,nb,nz,nr) = obje.l(nk,nb,nz,nr) ;

diffvf(iterr)= (sum((nk,nb,nz,nr),abs(1-vfnew(nk,nb,nz,nr)/ vfold(nk,nb,nz,nr))));
stopiter = 1$(abs(diffvf(iterr)) LE 1E-8) ;

) ;

display diffvf, stopiter ;

