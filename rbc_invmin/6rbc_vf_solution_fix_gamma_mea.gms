*--------------------------------------------------------
* one sector model
*--------------------------------------------------------
scalar starttime; starttime = jnow;

Parameters
beta_true,gamma_true,alphak_true,delta_true,rho_true,sigma_true;
beta_true=0.99;
gamma_true=2;
alphak_true=0.36;
delta_true=0.025;
rho_true=0.85;
sigma_true=0.04;

Set t    / 1*100/;
Set boot  /1*100/;

Parameters
beta,gamma,alphak,delta,rho,sigma,kbar,cbar;
beta=beta_true;
gamma=gamma_true;
alphak=alphak_true;
delta=delta_true;
rho=rho_true;
sigma=sigma_true;
kbar=1;
cbar=(1/beta+delta-1)/alphak - delta;

*--------------------------------------------------
* parameters for projection
*--------------------------------------------------
sets nk  number of nodes for capital /1*5/ ;
sets nz  number of nodes for productivity shock / 1*10/ ;
alias (nk,nka,nkb) ;
alias (nz,nza,nzb,nz)  ;

Parameters
err0(nz)
prob(nz)
;

Parameters
k0(nk)      initial capital stock
z0exp(nz)      initial market price level
z1(nz,nzb)  next period  price level
z0expmax
z0expmin
k0max
k0min
prob(nz)
;


*gaussion quadrature std=1
*err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;
*prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;

** Gaussian quadrature
* std=1
*gaussion quadrature std=1
*err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;
*prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;

*err0("1")= -2.85697001387281; err0("2")= -1.35562617997427; err0("3")= 0; err0("4")= 1.35562617997427; err0("5")= 2.85697001387281;
*prob("1")=0.0112574113277207;prob("2")=0.222075922005613; prob("3")=0.533333333333333;prob("4")=0.222075922005613 ;prob("5")=0.0112574113277207 ;

* std=1
err0("1")=-4.85946282833231; err0("2")=-3.58182348355193; err0("3")=-2.48432584163896; err0("4")=-1.46598909439116;err0("5")=-0.484935707515498;
err0("6")=0.484935707515498;err0("7")=1.46598909439116;err0("8")=2.48432584163896;err0("9")=3.58182348355193;err0("10")=4.85946282833231;

prob("1") = 4.31065263071829e-06;prob("2") = 0.000758070934312218;prob("3") = 0.0191115805007703;prob("4") = 0.135483702980268;prob("5") = 0.344642334932019;
prob("6") = 0.344642334932019;prob("7") = 0.135483702980268;prob("8") = 0.0191115805007703;prob("9") = 0.000758070934312218;prob("10") = 4.31065263071829e-06;


k0max = kbar*1.1;
k0min = kbar*0.9 ;
k0(nk)  = k0min + (k0max-k0min)/(card(nk)-1) * (ord(nk)-1) ;

z0expmax = 1.1;
z0expmin = 0.9;
z0exp(nz)  = z0expmin + (z0expmax-z0expmin)/(card(nz)-1) * (ord(nz)-1) ;


*--------------------------------------------------
* value function coefficients from MaliarJudd2016
*--------------------------------------------------
Parameters
mus0
mus1
mus2
mus11
mus12
mus22
mus111
mus112
mus122
mus222
mus1111
mus1112
mus1122
mus1222
mus2222
vfs(t)
;

* no binding constraint
*mus0=-1934.34628076383;
*mus1=718.707415294864;
*mus2=441.254481673093;
*mus11=-306.909205844206;
*mus12=-198.494283185508;
*mus22=-151.579808698686;
*mus111=58.3862770368426;
*mus112=40.5711908272261;
*mus122=29.4726952446437;
*mus222=25.9841273503470;

mus0=-2064.17757419394;
mus1=1040.49735048140;
mus2=640.142526058677 ;
mus11=-673.650781405834 ;
mus12=-433.831657068556 ;
mus22=-333.567211847312;
mus111=257.758002870836 ;
mus112=178.136199175614 ;
mus122=128.837555331704 ;
mus222=114.890261486110 ;
mus1111=-42.2210833814525;
mus1112=-30.8388594819872 ;
mus1122=-23.0860115802875 ;
mus1222=-18.0724899668120  ;
mus2222=-17.8573192761438 ;

$ontext
* with binding constraint
mus0=-1837.53758882935 ;
mus1=455.467852904367 ;
mus2 =344.048425870813 ;
mus11=-97.6881594180974 ;
mus12=-67.6087438111959 ;
mus22=-78.1951270498086 ;
$offtext

*---------------------------
* data simulation
*---------------------------
Parameters
vfs_deriv(t,boot)
zss(t,boot)
kss(t,boot)
css(t,boot)
yss(t,boot)
uss(t,boot)
z1ss(t,nzb,boot)
k1ss(t,nzb,boot)
vfs_deriv1
u1ss(t,nzb,boot)
c1ss(t,nzb,boot)
Residuals (t,boot)
Mean_Residuals(boot)
Residuals1 (t,boot)
Mean_Residuals1(boot)


;

zss(t,boot)=0;

loop(boot,
loop(t$(ord(t) LT card(t)),
zss(t+1,boot) = rho*zss(t,boot) + normal(0,0.04) ;
) ;
*execseed = 1 + gmillisec(jnow);
) ;

parameters
meanzss(boot)
stdzss(boot) ;
meanzss(boot) = sum(t, zss(t,boot))/card(t) ;
stdzss(boot) = (sum(t, (zss(t,boot)-meanzss(boot))*(zss(t,boot)-meanzss(boot)))/card(t) )**0.5 ;
display meanzss, stdzss ;
zss(t,boot) = zss(t,boot) - meanzss(boot) ;
meanzss(boot) = sum(t, zss(t,boot))/card(t) ;
stdzss(boot) = (sum(t, (zss(t,boot)-meanzss(boot))*(zss(t,boot)-meanzss(boot)))/card(t) )**0.5 ;
display meanzss, stdzss ;

kss(t,boot) = kbar;

loop(t$(ord(t) LE card(t)),
vfs_deriv(t,boot) = mus1
              + mus11*2*kss(t,boot)
              + mus12*exp(zss(t,boot))
              + mus111*3*kss(t,boot)*kss(t,boot)
              + mus112*2*kss(t,boot)*exp(zss(t,boot))
              + mus122  *exp(zss(t,boot))*exp(zss(t,boot))
              + mus1111*4*kss(t,boot)*kss(t,boot)*kss(t,boot)
              + mus1112*3*kss(t,boot)*kss(t,boot)*exp(zss(t,boot))
              + mus1122*2*kss(t,boot)*exp(zss(t,boot))*exp(zss(t,boot))
              + mus1222  *exp(zss(t,boot))*exp(zss(t,boot))*exp(zss(t,boot))
;

yss(t,boot) = ((1/beta+delta-1)/alphak)*exp(zss(t,boot))*kss(t,boot)**alphak  ;
uss(t,boot) = vfs_deriv(t,boot) / (1-delta + ((1/beta+delta-1)/alphak)*exp(zss(t,boot))*alphak*kss(t,boot)**(alphak-1)) ;
css(t,boot) = uss(t,boot)**(-1/gamma) ;
kss(t+1,boot) = (1-delta)*kss(t,boot) + yss(t,boot) - css(t,boot);

*CHECK BELOW THE APPROXIMATION ERROR
* Future period quantities in n_nodes integration nodes (k1,a1)
z1ss(t,nzb,boot) = rho*zss(t,boot)+sigma*err0(nzb);
k1ss(t,nzb,boot) =  (1-delta)*kss(t,boot) + yss(t,boot) - css(t,boot);
vfs_deriv1(t,nzb,boot) = mus1
              + mus11*2*k1ss(t,nzb,boot)
              + mus12*exp(z1ss(t,nzb,boot))
              + mus111*3*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)
              + mus112*2*k1ss(t,nzb,boot)*exp(z1ss(t,nzb,boot))
              + mus122  *exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))
              + mus1111*4*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)
              + mus1112*3*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)*exp(z1ss(t,nzb,boot))
              + mus1122*2*k1ss(t,nzb,boot)*exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))
              + mus1222  *exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))
;
u1ss(t,nzb,boot)=vfs_deriv1(t,nzb,boot) / (1-delta + ((1/beta+delta-1)/alphak)*exp(z1ss(t,nzb,boot))*alphak*k1ss(t,nzb,boot)**(alphak-1)) ;
c1ss(t,nzb,boot) = u1ss(t,nzb,boot)**(-1/gamma) ;


Residuals (t,boot)=sum(nzb, prob(nzb)*beta*c1ss(t,nzb,boot)**(-gamma)
                            *(1-delta+ ((1/beta+delta-1)/alphak)*exp(z1ss(t,nzb,boot))*alphak*k1ss(t,nzb,boot)**(alphak-1)
                             )
                      )
                  - css(t,boot)**(-gamma);

*Residuals1 (t,boot) = sum(nzb, prob(nzb)*beta*vfs_deriv1(t,nzb,boot)) - css(t,boot)**(-gamma);
)
;

Mean_Residuals(boot)=sum(t, Residuals(t,boot))/card(t);
*Mean_Residuals1(boot)=sum(t, Residuals1(t,boot))/card(t);
*Mean_Residuals(boot)=log10(sum(t, Residuals(t,boot))/card(t));
*Max_Residuals(boot)=max(Residuals(t,boot));

display vfs_deriv,css,kss, yss,zss,Residuals,Mean_Residuals ;


Parameters
cs(t)
ys(t)
;

cs(t) = css(t,"1");
ys(t) = yss(t,"1");

Parameters
cmean
ymean ;
cmean = sum(t, cs(t))/card(t) ;
ymean = sum(t, ys(t))/card(t) ;
Parameters
cs0(t)
;
cs0(t) = cs(t) ;
display cs ;
cs(t) = cs(t) + 0.02*normal(0,1)*cmean ;
ys(t) = ys(t) + 0.02*normal(0,1)*ymean ;
display cs;
parameters
diffc(t) ;
diffc(t) =cs(t)-cs0(t) ;
display diffc ;

$ontext
*----------------------
* selection of "smolyak" points at the second order of approximation for d=2 and mu=2 ;
*----------------------
Sets jk /1*3 / ;
alias(jk,jp) ;
Parameters
sk(jk,nk) ;
sk(jk,nk) = 1$(ord(nk) LE (2**(ord(jk)-1) + 1) ) ;
sk("1",nk) = 1$(ord(nk) EQ 1) ;
Parameters
sp(jp,np) ;
sp(jp,np) = 1$(ord(np) LE (2**(ord(jp)-1) + 1) ) ;
sp("1",np) = 1$(ord(np) EQ 1) ;


Parameters
smoltemp(nk,np)
smol(nk,np)
smolcount ;
smoltemp(nk,np) = sum((jk,jp)$(         ( (ord(jk)+ord(jp)) GE 2)
                                    and ( (ord(jk)+ord(jp)) LE 4)
                                         ), sk(jk,nk)*sp(jp,np) )  ;
smol(nk,np) = 1$smoltemp(nk,np) ;
*smol(nk,np)=1;
smolcount = sum((nk,np), smol(nk,np) ) ;
display smol, smolcount ;
$offtext


*----------------------
* estimation GME
*----------------------

parameters
err0(nz)
z0(nz)
prob(nz)
kmax
kmin
zmax
zmin
predict
predict1
;


sets k number of support values / 1*3 / ;

Parameters
betasv(k)
alphaksv(k)
gammasv(k)
deltasv(k)
rhozsv(k)
sigmazsv(k)
errzsv(k)           TFP shock
erreulersv(k)
erreulernesv(k)
errysv(k)
invmin
;

Parameters
betae
gammae
;

Variables
entropie
alphake
deltae
rhoze
sigmaze

ke(t)
kfin
ze(t)


zne(t,nz)
errz(t)
eulerer(t)
eulererne(t,nz)
erry(t)

mu0
mu1
mu2
mu11
mu12
mu22
mu111
mu112
mu122
mu222
mu1111
mu1112
mu1122
mu1222
mu2222
mu11111
mu11112
mu11122
mu11222
mu12222
mu22222
;

Positive variables
palphaksv(k)
pdeltasv(k)
prhozsv(k)
psigmazsv(k)

perrzsv(t,k)
perreulerer(t,k)
perreulerneer(t,k,nz)
perrysv(t,k)

vf(t)
vf_deriv(t)
vfne(t,nz)
vfne_deriv(t,nz)

mut(t)
;

equations
eqalphake
eqdeltae
eqrhoze
eqsigmaze
eqpalphaksv
eqpdeltasv
eqprhozsv
eqpsigmazsv

eqke(t)
eqzne(t,nz)
eqze(t)

eqcs(t)
eqys(t)
eqentropie

eqmeanze
eqstdze

eqerrz(t)
eqperrzsv(t)

eqvf(t)
eqvf_deriv(t)
eqvfne(t,nz)
eqvfne_deriv(t,nz)
eqvf_FOC(t)

eqeulerer(t)
eqpeulerersv(t)
eqeulerneer(t,nz)
eqpeulerneersv(t,nz)

eqmut(t)
eqcons(t)

eqerry(t)
eqperrysv(t)
;


eqentropie..     entropie =e=     - predict*sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - predict*sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - predict*sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - predict*sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - predict*sum((k,t), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - sum((k,t), perreulerer(t,k)*log(1.e-5+perreulerer(t,k)))
                                  - sum((k,t,nz), perreulerneer(t,k,nz)*log(1.e-5+perreulerneer(t,k,nz)))
                                  - sum((k,t), perrysv(t,k)*log(1.e-5+perrysv(t,k)))

;

$ontext
eqvf(t)..
                                   vf(t)                   =E=   mu0
                                                            + mu1*ke(t)
                                                            + mu2*exp(ze(t))
                                                            + mu11*ke(t)*ke(t)
                                                            + mu22*exp(ze(t))*exp(ze(t))
                                                            + mu12*ke(t)*exp(ze(t))
                                                            + mu111*ke(t)*ke(t)*ke(t)
                                                            + mu112*ke(t)*ke(t)*exp(ze(t))
                                                            + mu122*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + mu222*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                            + mu1111*ke(t)*ke(t)*ke(t)*ke(t)
                                                            + mu1112*ke(t)*ke(t)*ke(t)*exp(ze(t))
                                                            + mu1122*ke(t)*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + mu1222*ke(t)*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                            + mu2222*exp(ze(t))*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                            ;
$offtext

eqvf_deriv(t)..
                                    vf_deriv(t) + eulerer(t)   =E= mu1
                                                            + mu11*2*ke(t)
                                                            + mu12  *exp(ze(t))
                                                            + mu111*3*ke(t)*ke(t)
                                                            + mu112*2*ke(t)*exp(ze(t))
                                                            + mu122  *exp(ze(t))*exp(ze(t))
                                                            + mu1111*4*ke(t)*ke(t)*ke(t)
                                                            + mu1112*3*ke(t)*ke(t)*exp(ze(t))
                                                            + mu1122*2*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + mu1222  *exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                             ;

eqcs(t)..                          cs(t)**(-gammae) =e= (vf_deriv(t) ) /((1-deltae + ((1/betae+deltae-1)/alphake)*exp(ze(t))*alphake*ke(t)**(alphake-1))) ;

*eqcs(t)..                          vf(t) =e= ((cs(t))**(1-gammae)-1)/(1-gammae) + betae* sum(nz,prob(nz)*vfne(t,nz)) ;

$ontext
eqvfne(t,nz)..
                                    vfne(t,nz)   =E=  (mu0
                                                            + mu1*ke(t+1)
                                                            + mu2*exp(zne(t,nz))
                                                            + mu11*ke(t+1)*ke(t+1)
                                                            + mu22*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu12*ke(t+1)*exp(zne(t,nz))
                                                            + mu111*ke(t+1)*ke(t+1)*ke(t+1)
                                                            + mu112*ke(t+1)*ke(t+1)*exp(zne(t,nz))
                                                            + mu122*ke(t+1)*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu222*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu1111*ke(t+1)*ke(t+1)*ke(t+1)*ke(t+1)
                                                            + mu1112*ke(t+1)*ke(t+1)*ke(t+1)*ze(t+1)
                                                            + mu1122*ke(t+1)*ke(t+1)*ze(t+1)*ze(t+1)
                                                            + mu1222*ke(t+1)*ze(t+1)*ze(t+1)*ze(t+1)
                                                            + mu2222*ze(t+1)*ze(t+1)*ze(t+1)*ze(t+1)
                                                             )$(ord(t) LT card(t))
                                                             +
                                                             (mu0
                                                            + mu1*kfin
                                                            + mu2*exp(zne(t,nz))
                                                            + mu11*kfin*kfin
                                                            + mu22*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu12*kfin*exp(zne(t,nz))
                                                            + mu111*kfin*kfin*kfin
                                                            + mu112*kfin*kfin*exp(zne(t,nz))
                                                            + mu122*kfin*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu222*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu1111*kfin*kfin*kfin*kfin
                                                            + mu1112*kfin*kfin*kfin*ze(t+1)
                                                            + mu1122*kfin*kfin*ze(t+1)*ze(t+1)
                                                            + mu1222*kfin*ze(t+1)*ze(t+1)*ze(t+1)
                                                            + mu2222*ze(t+1)*ze(t+1)*ze(t+1)*ze(t+1)
                                                             )$(ord(t) EQ card(t))
                                                            ;
$offtext

eqvfne_deriv(t,nz)..
                                     vfne_deriv(t,nz) + eulererne(t,nz) =E= (mu1
                                                            + mu11*2*ke(t+1)
                                                            + mu12*exp(zne(t,nz))
                                                            + mu111*3*ke(t+1)*ke(t+1)
                                                            + mu112*2*ke(t+1)*exp(zne(t,nz))
                                                            + mu122*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu1111*4*ke(t+1)*ke(t+1)*ke(t+1)
                                                            + mu1112*3*ke(t+1)*ke(t+1)*exp(zne(t,nz))
                                                            + mu1122*2*ke(t+1)*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu1222*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            )$(ord(t) LT card(t))
                                                            +
                                                            (mu1
                                                            + mu11*2*kfin
                                                            + mu12*exp(zne(t,nz))
                                                            + mu111*3*kfin*kfin
                                                            + mu112*2*kfin*exp(zne(t,nz))
                                                            + mu122*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu1111*4*kfin*kfin*kfin
                                                            + mu1112*3*kfin*kfin*exp(zne(t,nz))
                                                            + mu1122*2*kfin*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + mu1222*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            )$(ord(t) EQ card(t))
                                                             ;


eqvf_FOC(t)..         cs(t)**(-gammae)  =E= betae*sum(nz,prob(nz)*(vfne_deriv(t,nz)+ 0*eulererne(t,nz)));


* Modeling next period TFP expectation
eqzne(t,nz)..        zne(t,nz) =E=  rhoze * ze(t) + sigmaze*err0(nz);

* Modeling TFP evolution process
eqze(t)..            ze(t+1)   =E=  rhoze * ze(t) + errz(t);


eqke(t)..            (ke(t+1) - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                                               )$(ord(t) LT card(t))
                                                               + (kfin    - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                                               )$(ord(t) EQ card(t)) =E= 0  ;

eqmeanze..            sum(t, errz(t))/card(t) =E= 0;

eqstdze..             (sum(t, errz(t)*errz(t))/card(t))**0.5 =E= sigmaze  ;

eqys(t)..             ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake  + erry(t);


eqalphake..       alphake   =E= sum(k, palphaksv(k)*alphaksv(k)) ;
eqpalphaksv..     1         =E= sum(k, palphaksv(k)) ;
eqdeltae..        deltae    =E= sum(k, pdeltasv(k)*deltasv(k)) ;
eqpdeltasv..      1         =E= sum(k, pdeltasv(k)) ;
eqrhoze..         rhoze   =E= sum(k, prhozsv(k)*rhozsv(k)) ;
eqprhozsv..       1         =E= sum(k, prhozsv(k)) ;
eqsigmaze..       sigmaze   =E= sum(k, psigmazsv(k)*sigmazsv(k)) ;
eqpsigmazsv..     1         =E= sum(k, psigmazsv(k)) ;

eqerrz(t)..           errz(t)    =E= sum(k, perrzsv(t,k)*errzsv(k)) ;
eqperrzsv(t)..        1          =E= sum(k, perrzsv(t,k)) ;

eqeulerer(t)..          eulerer(t) =E= sum(k, perreulerer(t,k)*erreulersv(k)) ;
eqpeulerersv(t)..       1          =E= sum(k, perreulerer(t,k)) ;

eqeulerneer(t,nz)..          eulererne(t,nz) =E= sum(k, perreulerneer(t,k,nz)*erreulernesv(k)) ;
eqpeulerneersv(t,nz)..       1          =E= sum(k, perreulerneer(t,k,nz)) ;

eqerry(t)..           erry(t)    =E= sum(k, perrysv(t,k)*errysv(k)) ;
eqperrysv(t)..        1          =E= sum(k, perrysv(t,k)) ;


model estimation /
eqalphake
eqdeltae
eqpalphaksv
eqpdeltasv
eqrhoze
eqprhozsv
eqsigmaze
eqpsigmazsv
eqke
eqze
eqcs
eqys
eqzne
eqentropie
eqmeanze
eqstdze
eqerrz
eqperrzsv
*eqvf
eqvf_deriv
*eqvfne
eqvfne_deriv
eqvf_FOC
eqeulerer
eqpeulerersv
eqeulerneer
eqpeulerneersv
eqerry
eqperrysv
/ ;


*initiate chebyshev coefficient for estimation
betae = 0.99;
gammae = 2;

*$ontext
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.36 ;
alphaksv("3")   = 0.5 ;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.025 ;
deltasv("3")    = 0.05 ;
rhozsv("1")    = 0.7 ;
rhozsv("2")    = 0.85 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.04 ;
sigmazsv("3")   = 0.08;
*$offtext

$ontext
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.5 ;
alphaksv("3")   = 0.8 ;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.05 ;
deltasv("3")    = 0.10 ;
rhozsv("1")    = 0.6 ;
rhozsv("2")    = 0.8 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2;
$offtext


z0(nz) = err0(nz)*sigmazsv("3")  ;

Parameters
kbar
cbar
;


zmax = z0("10")*1.5 ;
zmin = -zmax ;

*zmax = log(1.5);
*zmin = -zmax ;

errzsv("1") = zmin ;
errzsv("3") = zmax ;

parameters
deltaeboot(boot)
alphakeboot(boot)
rhozeboot(boot)
sigmazeboot(boot)
modelboot(boot)
;

parameter
vfderivss ;

loop(boot$(ord(boot) ),
cs(t) = css(t,boot);
ys(t) = yss(t,boot);

* we reinitialize endogenous variables
deltae.l  = deltasv("2") ;
alphake.l = alphaksv("2") ;
rhoze.l    = rhozsv("2") ;
sigmaze.l  = sigmazsv("2") ;

deltae.lo  = deltasv("1") ;
alphake.lo = alphaksv("1") ;
rhoze.lo    = rhozsv("1") ;
sigmaze.lo  = sigmazsv("1") ;

deltae.up  = deltasv("3") ;
alphake.up = alphaksv("3") ;
rhoze.up    = rhozsv("3") ;
sigmaze.up  = sigmazsv("3") ;

kbar   = 1;
cbar   = (1/betae+deltae.l-1)/alphake.l - deltae.l ;
kmax     = kbar*1.5;
kmin     = kbar*0.5;
ke.l(t)           = kbar  ;
kfin.l            = kbar ;
ke.lo(t)          = kmin;
ke.up(t)          = kmax ;
ze.l(t)           = 0;
ze.up(t)          = zmax ;
ze.lo(t)          = zmin ;
zne.l(t,nz)       = 0 ;
zne.up(t,nz)       = zmax ;
zne.lo(t,nz)       = zmin ;
errz.l(t)         = 0 ;
eulerer.l(t)      = 0 ;
eulererne.l(t,nz) = 0 ;
palphaksv.l(k)          = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
psigmazsv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perreulerer.l(t,k)      = 1/card(k) ;
perreulerneer.l(t,k,nz) = 1/card(k) ;
perrysv.l(t,k)          = 1/card(k) ;

vf_deriv.l(t) =  cs(t)**(-gammae)*
            (1-deltae.l + ((1/betae+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vf_deriv.lo(t) = 0;
vfne_deriv.l(t,nz) =  cs(t)**(-gammae)*
            (1-deltae.l + ((1/betae+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vfne_deriv.lo(t,nz) = 0;

vfderivss = (sum(t,cs(t))/card(t) )**(-gammae)*
            (1-deltae.l + ((1/betae+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));


mu1.l = vfderivss;
*display mu1.l;
$ontext
mu0.l=mus0;
mu1.l=mus1;
mu2.l=mus2;
mu11.l=mus11;
mu12.l=mus12;
mu22.l=mus22;
mu111.l=mus111;
mu112.l=mus112;
mu122.l=mus122;
mu222.l=mus222;
mu1111.l=mus1111;
mu1112.l=mus1112 ;
mu1122.l=mus1122 ;
mu1222.l=mus1222 ;
mu2222.l=mus2222;
$offtext

predict = 0.01;

erreulersv("1")         = -1e-1*vfderivss;
erreulersv("3")         = 1e-1*vfderivss;
erreulernesv("1")       = -1e-1*vfderivss;
erreulernesv("3")       = 1e-1*vfderivss ;
errysv("1")             = -0.05*sum(t, ys(t))/card(t);
errysv("3")             = 0.05*sum(t, ys(t))/card(t);
solve estimation using nlp maximising entropie ;

deltaeboot(boot) = deltae.l ;
alphakeboot(boot) = alphake.l ;
rhozeboot(boot) = rhoze.l   ;
sigmazeboot(boot) = sigmaze.l  ;
modelboot(boot) = estimation.solvestat ;

) ;

parameters
alphakem
deltaem
rhozem
sigmazem
alphakestd
gammaestd
deltaestd
betaestd
rhozestd
sigmazestd
;

alphakem = sum(boot, alphakeboot(boot))/card(boot) ;
alphakestd = (sum(boot, (alphakeboot(boot)-alphakem)*(alphakeboot(boot)-alphakem) )/card(boot) )**0.5 ;
deltaem = sum(boot, deltaeboot(boot))/card(boot) ;
deltaestd = (sum(boot, (deltaeboot(boot)-deltaem)*(deltaeboot(boot)-deltaem) )/card(boot) )**0.5 ;
rhozem = sum(boot, rhozeboot(boot))/card(boot) ;
rhozestd = (sum(boot, (rhozeboot(boot)-rhozem)*(rhozeboot(boot)-rhozem) )/card(boot) )**0.5 ;
sigmazem = sum(boot, sigmazeboot(boot))/card(boot) ;
sigmazestd = (sum(boot, (sigmazeboot(boot)-sigmazem)*(sigmazeboot(boot)-sigmazem) )/card(boot) )**0.5 ;

parameters
alphakemse
deltaemse
rhozemse
sigmazemse
;

alphakemse = (power((alphakem-alphak_true),2)+alphakestd**2)**0.5;
deltaemse = (power((deltaem-delta_true),2)+deltaestd**2)**0.5;
rhozemse = (power((rhozem-rho_true),2)+rhozestd**2)**0.5;
sigmazemse = (power((sigmazem-sigma_true),2)+sigmazestd**2)**0.5;

display alphakem, alphakestd, alphakemse ;
display deltaem, deltaestd, deltaemse ;
display rhozem, rhozestd, rhozemse ;
display sigmazem, sigmazestd, sigmazemse ;

scalar elapsed; elapsed = (jnow - starttime)*24*3600;

display deltaeboot, alphakeboot, rhozeboot, sigmazeboot, modelboot,elapsed ;

Parameters
res(boot,*);

res(boot,"alphae") = alphakeboot(boot) ;
res(boot,"deltae") = deltaeboot(boot) ;
res(boot,"rhoe") = rhozeboot(boot);
res(boot,"sigmae") = sigmazeboot(boot);
res(boot,"solvestat") = modelboot(boot);

execute_unload 'res_vf_solution_fix_gamma_mea.gdx';
$libinclude xlexport res res4.xlsx res!a1:h101




