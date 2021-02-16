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
beta,gamma,alphak,delta,rho,sigma,kbar,cbar,psib;
beta=beta_true;
gamma=gamma_true;
alphak=alphak_true;
delta=delta_true;
rho=rho_true;
sigma=sigma_true;
kbar=1;
cbar=(1/beta+delta-1)/alphak - delta;
psib=0.95;

*-----------------------------------------------------------
* Simulate data from the vf coefficients from MaliarJudd2016
*-----------------------------------------------------------
* parameters for projection
sets nk  number of nodes for capital /1*5/ ;
sets nz  number of nodes for productivity shock / 1*5/ ;
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

err0("1")= -2.85697001387281; err0("2")= -1.35562617997427; err0("3")= 0; err0("4")= 1.35562617997427; err0("5")= 2.85697001387281;
prob("1")=0.0112574113277207;prob("2")=0.222075922005613; prob("3")=0.533333333333333;prob("4")=0.222075922005613 ;prob("5")=0.0112574113277207 ;

* std=1
*err0("1")=-4.85946282833231; err0("2")=-3.58182348355193; err0("3")=-2.48432584163896; err0("4")=-1.46598909439116;err0("5")=-0.484935707515498;
*err0("6")=0.484935707515498;err0("7")=1.46598909439116;err0("8")=2.48432584163896;err0("9")=3.58182348355193;err0("10")=4.85946282833231;

*prob("1") = 4.31065263071829e-06;prob("2") = 0.000758070934312218;prob("3") = 0.0191115805007703;prob("4") = 0.135483702980268;prob("5") = 0.344642334932019;
*prob("6") = 0.344642334932019;prob("7") = 0.135483702980268;prob("8") = 0.0191115805007703;prob("9") = 0.000758070934312218;prob("10") = 4.31065263071829e-06;

k0max = kbar*1.1;
k0min = kbar*0.9 ;
k0(nk)  = k0min + (k0max-k0min)/(card(nk)-1) * (ord(nk)-1) ;

z0expmax = 1.1;
z0expmin = 0.9;
z0exp(nz)  = z0expmin + (z0expmax-z0expmin)/(card(nz)-1) * (ord(nz)-1) ;


*--------------------------------------------------
* value function coefficients from MaliarJudd2016
*--------------------------------------------------
Sets ncoeff number of coefficients   /1*21/;
Parameters
coeffs(ncoeff)
vfs(t)
;


* with binding constraint
coeffs("1") = -3231.35215387506 ;
coeffs("2") =1383.96121955887 ;
coeffs("3") =4645.27101925873 ;
coeffs("4") =-795.746871155866 ;
coeffs("5") =-1168.75953038294 ;
coeffs("6") =-5601.64342010823 ;
coeffs("7") =288.265892493659 ;
coeffs("8") =319.891012520627 ;
coeffs("9") =690.960008538866 ;
coeffs("10") =3229.02122436838 ;
coeffs("11") =-46.4089018676276 ;
coeffs("12") =-43.6216367357155 ;
coeffs("13") =-72.0711597124692 ;
coeffs("14") =-165.780577855195 ;
coeffs("15") =-713.916083993840 ;

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
invini_check(t,boot)


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
vfs_deriv(t,boot) = coeffs("2")
              + coeffs("4")*2*kss(t,boot)
              + coeffs("5")*exp(zss(t,boot))
              + coeffs("7")*3*kss(t,boot)*kss(t,boot)
              + coeffs("8")*2*kss(t,boot)*exp(zss(t,boot))
              + coeffs("9")  *exp(zss(t,boot))*exp(zss(t,boot))
              + coeffs("11")*4*kss(t,boot)*kss(t,boot)*kss(t,boot)
              + coeffs("12")*3*kss(t,boot)*kss(t,boot)*exp(zss(t,boot))
              + coeffs("13")*2*kss(t,boot)*exp(zss(t,boot))*exp(zss(t,boot))
              + coeffs("14")  *exp(zss(t,boot))*exp(zss(t,boot))*exp(zss(t,boot))
;

yss(t,boot) = ((1/beta+delta-1)/alphak)*exp(zss(t,boot))*kss(t,boot)**alphak  ;
uss(t,boot) = vfs_deriv(t,boot) / (1-delta + ((1/beta+delta-1)/alphak)*exp(zss(t,boot))*alphak*kss(t,boot)**(alphak-1)) ;
css(t,boot) = uss(t,boot)**(-1/gamma) ;
css(t,boot) = min(css(t,boot), yss(t,boot) - psib*delta);
kss(t+1,boot) = (1-delta)*kss(t,boot) + yss(t,boot) - css(t,boot);

*CHECK BELOW THAT INVESTMENT IS POSITIVE
invini_check(t,boot) = kss(t+1,boot)-(1-delta)*kss(t,boot)- psib*delta ;

$ontext
*CHECK BELOW THE APPROXIMATION ERROR
* Future period quantities in n_nodes integration nodes (k1,a1)
z1ss(t,nzb,boot) = rho*zss(t,boot)+sigma*err0(nzb);
k1ss(t,nzb,boot) =  (1-delta)*kss(t,boot) + yss(t,boot) - css(t,boot);
vfs_deriv1(t,nzb,boot) = coeffs("2")
              + coeffs("4")*2*k1ss(t,nzb,boot)
              + coeffs("5")*exp(z1ss(t,nzb,boot))
              + coeffs("7")*3*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)
              + coeffs("8")*2*k1ss(t,nzb,boot)*exp(z1ss(t,nzb,boot))
              + coeffs("9")  *exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))
              + coeffs("11")*4*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)
              + coeffs("12")*3*k1ss(t,nzb,boot)*k1ss(t,nzb,boot)*exp(z1ss(t,nzb,boot))
              + coeffs("13")*2*k1ss(t,nzb,boot)*exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))
              + coeffs("14")  *exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))*exp(z1ss(t,nzb,boot))
;
u1ss(t,nzb,boot)=vfs_deriv1(t,nzb,boot) / (1-delta + ((1/beta+delta-1)/alphak)*exp(z1ss(t,nzb,boot))*alphak*k1ss(t,nzb,boot)**(alphak-1)) ;
c1ss(t,nzb,boot) = u1ss(t,nzb,boot)**(-1/gamma) ;


Residuals (t,boot)=sum(nzb, prob(nzb)*beta*c1ss(t,nzb,boot)**(-gamma)
                            *(1-delta+ ((1/beta+delta-1)/alphak)*exp(z1ss(t,nzb,boot))*alphak*k1ss(t,nzb,boot)**(alphak-1)
                             )
                      )
                  - css(t,boot)**(-gamma);

*Residuals1 (t,boot) = sum(nzb, prob(nzb)*beta*vfs_deriv1(t,nzb,boot)) - css(t,boot)**(-gamma);
$offtext
)
;

*Mean_Residuals(boot)=sum(t, Residuals(t,boot))/card(t);


display vfs_deriv,css,kss, yss,zss,invini_check ;



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
cs(t) = cs(t) + 0.005*normal(0,1)*cmean ;
ys(t) = ys(t) + 0.005*normal(0,1)*ymean ;
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
eulererne(t)
erry(t)

coeffe(ncoeff)






;

Positive variables
palphaksv(k)
pdeltasv(k)
prhozsv(k)
psigmazsv(k)

perrzsv(t,k)
perreulerer(t,k)
perreulerneer(t,k)
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
eqeulerneer(t)
eqpeulerneersv(t)

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
                                  - sum((k,t), perreulerneer(t,k)*log(1.e-5+perreulerneer(t,k)))
                                  - sum((k,t), perrysv(t,k)*log(1.e-5+perrysv(t,k)))

;

eqvf_deriv(t)..
                                    vf_deriv(t)    =E= coeffe("2")
                                                            + coeffe("4")*2*ke(t)
                                                            + coeffe("5")  *exp(ze(t))
                                                            + coeffe("7")*3*ke(t)*ke(t)
                                                            + coeffe("8")*2*ke(t)*exp(ze(t))
                                                            + coeffe("9")  *exp(ze(t))*exp(ze(t))
                                                            + coeffe("11")*4*ke(t)*ke(t)*ke(t)
                                                            + coeffe("12")*3*ke(t)*ke(t)*exp(ze(t))
                                                            + coeffe("13")*2*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + coeffe("14")  *exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                             ;

eqcs(t)..                            (cs(t)- eulerer(t))**(-gammae) =e= (vf_deriv(t) + (1-deltae)*mut(t))/((1-deltae + ((1/betae+deltae-1)/alphake)*exp(ze(t))*alphake*ke(t)**(alphake-1))) ;

*eqcs(t)..                          vf(t) =e= ((cs(t))**(1-gammae)-1)/(1-gammae) + betae* sum(nz,prob(nz)*vfne(t,nz)) ;


eqvfne_deriv(t,nz)..
                                     vfne_deriv(t,nz)  =E= (coeffe("2")
                                                            + coeffe("4")*2*ke(t+1)
                                                            + coeffe("5")*exp(zne(t,nz))
                                                            + coeffe("7")*3*ke(t+1)*ke(t+1)
                                                            + coeffe("8")*2*ke(t+1)*exp(zne(t,nz))
                                                            + coeffe("9")*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("11")*4*ke(t+1)*ke(t+1)*ke(t+1)
                                                            + coeffe("12")*3*ke(t+1)*ke(t+1)*exp(zne(t,nz))
                                                            + coeffe("13")*2*ke(t+1)*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("14")*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            )$(ord(t) LT card(t))
                                                            +
                                                            (coeffe("2")
                                                            + coeffe("4")*2*kfin
                                                            + coeffe("5")*exp(zne(t,nz))
                                                            + coeffe("7")*3*kfin*kfin
                                                            + coeffe("8")*2*kfin*exp(zne(t,nz))
                                                            + coeffe("9")*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("11")*4*kfin*kfin*kfin
                                                            + coeffe("12")*3*kfin*kfin*exp(zne(t,nz))
                                                            + coeffe("13")*2*kfin*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("14")*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            )$(ord(t) EQ card(t))
                                                             ;


eqvf_FOC(t)..                         (cs(t)- eulererne(t))**(-gammae)  =E= betae*sum(nz,prob(nz)*(vfne_deriv(t,nz))) + mut(t) ;

*to make sure the data satify the constraint
eqcons(t)..          ke(t+1)-(1-deltae)*ke(t) =G= invmin;

eqmut(t)..           mut(t)*(ys(t)-cs(t)-invmin) =E= 0;

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

eqeulerneer(t)..          eulererne(t) =E= sum(k, perreulerneer(t,k)*erreulernesv(k)) ;
eqpeulerneersv(t)..       1          =E= sum(k, perreulerneer(t,k)) ;

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
eqmut
*eqcons
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


z0(nz) = err0(nz)*sigmazsv("3");

Parameters
kbar
cbar
;


zmax = z0("5")*1.5 ;
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
eulererne.l(t) = 0 ;
palphaksv.l(k)          = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
psigmazsv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perreulerer.l(t,k)      = 1/card(k) ;
perreulerneer.l(t,k) = 1/card(k) ;
perrysv.l(t,k)          = 1/card(k) ;

vf_deriv.l(t) =  cs(t)**(-gammae)*
            (1-deltae.l + ((1/betae+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vf_deriv.lo(t) = 0;
vfne_deriv.l(t,nz) =  cs(t)**(-gammae)*
            (1-deltae.l + ((1/betae+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vfne_deriv.lo(t,nz) = 0;

vfderivss = (sum(t,cs(t))/card(t) )**(-gammae)*
            (1-deltae.l + ((1/betae+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));



*coeffe("1").l = 100;
coeffe.l("2") = vfderivss;
*display coeffe("2").l;
$ontext
coeffe.l(ncoeff)=coeffs(ncoeff);
$offtext

mut.l(t) = 0;
mut.lo(t) = 0;
mut.up(t) = +inf;
invmin = 0.025*0.95;

predict = 0.01;

erreulersv("1")         = -0.2*sum(t,cs(t))/card(t);
erreulersv("3")         = 0.2*sum(t,cs(t))/card(t);
erreulernesv("1")       = -0.2*sum(t,cs(t))/card(t);
erreulernesv("3")       = 0.2*sum(t,cs(t))/card(t);
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

parameters
alphakebias
deltaebias
rhozebias
sigmazebias
;

alphakebias = (alphakem-alphak_true)/alphak_true;
deltaebias = (deltaem-delta_true)/delta_true;
rhozebias = (rhozem-rho_true)/rho_true;
sigmazebias = (sigmazem-sigma_true)/sigma_true;


display alphakem, alphakestd, alphakemse ;
display deltaem, deltaestd, deltaemse ;
display rhozem, rhozestd, rhozemse ;
display sigmazem, sigmazestd, sigmazemse ;

display deltaeboot, alphakeboot, rhozeboot, sigmazeboot, modelboot, elapsed   ;

scalar elapsed; elapsed = (jnow - starttime)*24*3600;

Parameters
res(boot,*);

res(boot,"alphae") = alphakeboot(boot) ;
res(boot,"deltae") = deltaeboot(boot) ;
res(boot,"rhoe") = rhozeboot(boot);
res(boot,"sigmae") = sigmazeboot(boot);
res(boot,"solvestat") = modelboot(boot);

Parameters
res_table(*,*);

res_table("alpha","True")=alphak_true;
res_table("alpha","S.D.")=alphakestd;
res_table("alpha","Bias")=alphakebias;
res_table("alpha","MSE")=alphakemse;

res_table("delta","True")=delta_true;
res_table("delta","Mean")=deltaem;
res_table("delta","S.D.")=deltaestd;
res_table("delta","Bias")=deltaebias;
res_table("delta","MSE")=deltaemse;

res_table("rho","True")=rho_true;
res_table("rho","Mean")=rhozem;
res_table("rho","S.D.")=rhozestd;
res_table("rho","Bias")=rhozebias;
res_table("rho","MSE")=rhozemse;

res_table("sigma","True")=sigma_true;
res_table("sigma","Mean")=sigmazem;
res_table("sigma","S.D.")=sigmazestd;
res_table("sigma","Bias")=sigmazebias;
res_table("sigma","MSE")=sigmazemse;

*execute_unload 'res-invmin-vf-solution-fix-gamma_mea.gdx';
*$libinclude xlexport res res3.xlsx res!a1:h101



