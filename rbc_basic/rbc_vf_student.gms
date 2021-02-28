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
z0exp(nz)   initial market price level
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

* no binding constraint
coeffs("1")=-2064.17757419394;
coeffs("2")=1040.49735048140;
coeffs("3")=640.142526058677 ;
coeffs("4")=-673.650781405834 ;
coeffs("5")=-433.831657068556 ;
coeffs("6")=-333.567211847312;
coeffs("7")=257.758002870836 ;
coeffs("8")=178.136199175614 ;
coeffs("9")=128.837555331704 ;
coeffs("10")=114.890261486110 ;
coeffs("11")=-42.2210833814525;
coeffs("12")=-30.8388594819872 ;
coeffs("13")=-23.0860115802875 ;
coeffs("14")=-18.0724899668120  ;
coeffs("15")=-17.8573192761438 ;

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

* generate tfp/price shock from student distribution
$funclibin stolib stodclib
function randstudentt  /stolib.dstudentt /;
set i /i1*i20/
set j /randstudentt/
parameter randx(i,j)    numbers from distributions;
randx(i,"randstudentt")=randstudentt(5);
display randx;

loop(boot,
loop(t$(ord(t) LT card(t)),
zss(t+1,boot) = rho*zss(t,boot) + 0.04*randstudentt(5) ;
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
kss(t+1,boot) = (1-delta)*kss(t,boot) + yss(t,boot) - css(t,boot);

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
invmin
;

Variables
entropie
betae
alphake
gammae
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

coeffe(ncoeff)
;

Positive variables
pbetasv(k)
palphaksv(k)
pgammasv(k)
pdeltasv(k)
prhozsv(k)
psigmazsv(k)

perrzsv(t,k)
perreulerer(t,k)
perreulerneer(t,k)

vf(t)
vf_deriv(t)
vfne(t,nz)
vfne_deriv(t,nz)

mut(t)
;

equations
eqbetae
eqalphake
eqgammae
eqdeltae
eqrhoze
eqsigmaze
eqpalphaksv
eqpgammasv
eqpdeltasv
eqpbetasv
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

eqeulererm
eqeulerneerm
;


eqentropie..     entropie =e=     - predict*sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - predict*sum(k, pgammasv(k)*LOG(1.e-5+pgammasv(k)))
                                  - predict*sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - predict*sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - predict*sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - predict*sum(k, pbetasv(k)*LOG(1.e-5+pbetasv(k)))
                                  - predict*sum((k,t), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - 100*sum((k,t), perreulerer(t,k)*log(1.e-5+perreulerer(t,k)))
                                  - 100*sum((k,t), perreulerneer(t,k)*log(1.e-5+perreulerneer(t,k)))
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

eqcs(t)..                          (cs(t) - eulerer(t))**(-gammae) =e= (vf_deriv(t) ) /((1-deltae + ((1/betae+deltae-1)/alphake)*exp(ze(t))*alphake*ke(t)**(alphake-1))) ;

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


eqvf_FOC(t)..         (cs(t)-eulererne(t))**(-gammae)  =E= betae*sum(nz,prob(nz)*(vfne_deriv(t,nz)+ 0*eulererne(t)));


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

eqys(t)..             ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake  ;


eqbetae..         betae    =E= sum(k, pbetasv(k)*betasv(k)) ;
eqpbetasv..       1         =E= sum(k, pbetasv(k)) ;
eqgammae..        gammae    =E= sum(k, pgammasv(k)*gammasv(k)) ;
eqpgammasv..      1         =E= sum(k, pgammasv(k)) ;
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
eqeulererm..            sum(t,eulerer(t)) =E= 0 ;


eqeulerneer(t)..          eulererne(t) =E= sum(k, perreulerneer(t,k)*erreulernesv(k)) ;
eqpeulerneersv(t)..       1          =E= sum(k, perreulerneer(t,k)) ;
eqeulerneerm ..           sum(t,eulererne(t)) =E= 0 ;


model estimation /
eqalphake
eqgammae
eqdeltae
eqpalphaksv
eqpgammasv
eqpdeltasv
eqrhoze
eqprhozsv
eqsigmaze
eqpsigmazsv
eqbetae
eqpbetasv
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
eqeulererm
eqeulerneerm
/ ;


*initiate chebyshev coefficient for estimation
$ontext
betasv("1")   = 0.98 ;
betasv("2")   = 0.99 ;
betasv("3")   = 0.999 ;
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.36 ;
alphaksv("3")   = 0.5 ;
gammasv("1")    = 0.1 ;
gammasv("2")    = 2 ;
gammasv("3")    = 4;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.025 ;
deltasv("3")    = 0.05 ;
rhozsv("1")    = 0.7 ;
rhozsv("2")    = 0.85 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.04 ;
sigmazsv("3")   = 0.08;
$offtext

*$ontext
betasv("1")   = 0.98 ;
betasv("2")   = 0.99 ;
betasv("3")   = 0.999 ;
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.5 ;
alphaksv("3")   = 0.8 ;
gammasv("1")    = 0.1 ;
gammasv("2")    = 1.5 ;
gammasv("3")    = 3;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.05 ;
deltasv("3")    = 0.10 ;
rhozsv("1")    = 0.6 ;
rhozsv("2")    = 0.8 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2;
*$offtext


z0(nz) = err0(nz)*sigmazsv("3")   ;

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
gammaeboot(boot)
deltaeboot(boot)
alphakeboot(boot)
betaeboot(boot)
rhozeboot(boot)
sigmazeboot(boot)
modelboot(boot)
;

parameter
vfderivss ;

loop(boot$(ord(boot)),
cs(t) = css(t,boot);
ys(t) = yss(t,boot);

* we reinitialize endogenous variables
gammae.l  = gammasv("2") ;
deltae.l  = deltasv("2") ;
alphake.l = alphaksv("2") ;
betae.l   = betasv("2") ;
rhoze.l    = rhozsv("2") ;
sigmaze.l  = sigmazsv("2") ;

gammae.lo  = gammasv("1") ;
deltae.lo  = deltasv("1") ;
alphake.lo = alphaksv("1") ;
betae.lo   = betasv("1") ;
rhoze.lo    = rhozsv("1") ;
sigmaze.lo  = sigmazsv("1") ;

gammae.up  = gammasv("3") ;
deltae.up  = deltasv("3") ;
alphake.up = alphaksv("3") ;
betae.up   = betasv("3") ;
rhoze.up    = rhozsv("3") ;
sigmaze.up  = sigmazsv("3") ;

kbar   = 1;
cbar   = (1/betae.l+deltae.l-1)/alphake.l - deltae.l ;
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
pbetasv.l(k)            = 1/card(k) ;
palphaksv.l(k)          = 1/card(k) ;
pgammasv.l(k)           = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
psigmazsv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perreulerer.l(t,k)      = 1/card(k) ;
perreulerneer.l(t,k) = 1/card(k) ;

vf_deriv.l(t) =  cs(t)**(-gammae.l)*
            (1-deltae.l + ((1/betae.l+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vf_deriv.lo(t) = 0;
vfne_deriv.l(t,nz) =  cs(t)**(-gammae.l)*
            (1-deltae.l + ((1/betae.l+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vfne_deriv.lo(t,nz) = 0;

vfderivss = (sum(t,cs(t))/card(t) )**(-gammae.l)*
            (1-deltae.l + ((1/betae.l+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));



*coeffe("1").l = 100;
coeffe.l("2") = vfderivss;
*display coeffe("2").l;
$ontext
coeffe.l(ncoeff)=coeffs(ncoeff);
$offtext

predict = 0.01;

erreulersv("1")         = -0.2*sum(t,cs(t))/card(t);
erreulersv("3")         = 0.2*sum(t,cs(t))/card(t);
erreulernesv("1")       = -0.2*sum(t,cs(t))/card(t);
erreulernesv("3")       = 0.2*sum(t,cs(t))/card(t);
solve estimation using nlp maximising entropie ;

gammaeboot(boot) = gammae.l ;
deltaeboot(boot) = deltae.l ;
alphakeboot(boot) = alphake.l ;
betaeboot(boot) = betae.l   ;
rhozeboot(boot) = rhoze.l   ;
sigmazeboot(boot) = sigmaze.l  ;
modelboot(boot) = estimation.solvestat ;

) ;

parameters
alphakem
gammaem
deltaem
betaem
rhozem
sigmazem
alphakestd
gammaestd
deltaestd
betaestd
rhozestd
sigmazestd
;

alphakem = sum(boot$((modelboot(boot) LE 2) ), alphakeboot(boot))/sum(boot$((modelboot(boot) LE 2) ), 1) ;
alphakestd = (sum(boot$((modelboot(boot) LE 2) ), (alphakeboot(boot)-alphakem)*(alphakeboot(boot)-alphakem) )/sum(boot$((modelboot(boot) LE 2) ),1) )**0.5 ;
gammaem = sum(boot$((modelboot(boot) LE 2) ), gammaeboot(boot))/sum(boot$((modelboot(boot) LE 2) ),1) ;
gammaestd = (sum(boot$((modelboot(boot) LE 2) ), (gammaeboot(boot)-gammaem)*(gammaeboot(boot)-gammaem) )/sum(boot$((modelboot(boot) LE 2) ),1) )**0.5 ;
deltaem = sum(boot$((modelboot(boot) LE 2) ), deltaeboot(boot))/sum(boot$((modelboot(boot) LE 2) ),1) ;
deltaestd = (sum(boot$((modelboot(boot) LE 2) ), (deltaeboot(boot)-deltaem)*(deltaeboot(boot)-deltaem) )/sum(boot$((modelboot(boot) LE 2) ),1) )**0.5 ;
betaem = sum(boot$((modelboot(boot) LE 2) ), betaeboot(boot))/sum(boot$((modelboot(boot) LE 2) ),1) ;
betaestd = (sum(boot$((modelboot(boot) LE 2) ), (betaeboot(boot)-betaem)*(betaeboot(boot)-betaem) )/sum(boot$((modelboot(boot) LE 2) ),1) )**0.5 ;
rhozem = sum(boot$((modelboot(boot) LE 2) ), rhozeboot(boot))/sum(boot$((modelboot(boot) LE 2) ),1) ;
rhozestd = (sum(boot$((modelboot(boot) LE 2) ), (rhozeboot(boot)-rhozem)*(rhozeboot(boot)-rhozem) )/sum(boot$((modelboot(boot) LE 2) ),1) )**0.5 ;
sigmazem = sum(boot$((modelboot(boot) LE 2) ), sigmazeboot(boot))/sum(boot$((modelboot(boot) LE 2) ),1) ;
sigmazestd = (sum(boot$((modelboot(boot) LE 2) ), (sigmazeboot(boot)-sigmazem)*(sigmazeboot(boot)-sigmazem) )/sum(boot$((modelboot(boot) LE 2) ),1) )**0.5 ;

parameters
alphakemse
gammaemse
betaemse
deltaemse
rhozemse
sigmazemse
;

alphakemse = (power((alphakem-alphak_true),2)+alphakestd**2)**0.5;
gammaemse = (power((gammaem-gamma_true),2)+gammaestd**2)**0.5;
deltaemse = (power((deltaem-delta_true),2)+deltaestd**2)**0.5;
betaemse = (power((betaem-beta_true),2)+betaestd**2)**0.5;
rhozemse = (power((rhozem-rho_true),2)+rhozestd**2)**0.5;
sigmazemse = (power((sigmazem-sigma_true),2)+sigmazestd**2)**0.5;

parameters
alphakebias
gammaebias
betaebias
deltaebias
rhozebias
sigmazebias
;

alphakebias = (alphakem-alphak_true)/alphak_true;
gammaebias = (gammaem-gamma_true)/gamma_true;
betaebias = (betaem-beta_true)/beta_true;
deltaebias = (deltaem-delta_true)/delta_true;
rhozebias = (rhozem-rho_true)/rho_true;
sigmazebias = (sigmazem-sigma_true)/sigma_true;


display alphakem, alphakestd, alphakemse ;
display gammaem, gammaestd, gammaemse ;
display deltaem, deltaestd, deltaemse ;
display betaem, betaestd, betaemse ;
display rhozem, rhozestd, rhozemse ;
display sigmazem, sigmazestd, sigmazemse ;


scalar elapsed; elapsed = (jnow - starttime)*24*3600;

display gammaeboot, deltaeboot, alphakeboot, betaeboot, rhozeboot, sigmazeboot, modelboot, elapsed   ;


Parameters
res(boot,*);

res(boot,"betae") = betaeboot(boot) ;
res(boot,"gammae") = gammaeboot(boot) ;
res(boot,"alphae") = alphakeboot(boot) ;
res(boot,"deltae") = deltaeboot(boot) ;
res(boot,"rhoe") = rhozeboot(boot);
res(boot,"sigmae") = sigmazeboot(boot);
res(boot,"solvestat") = modelboot(boot);

Parameters
res_table(*,*);
res_table("beta","True")=beta_true;
res_table("beta","Mean")=betaem;
res_table("beta","S.D.")=betaestd;
res_table("beta","Bias")=betaebias;
res_table("beta","MSE")=betaemse;

res_table("gamma","True")=gamma_true;
res_table("gamma","Mean")=gammaem;
res_table("gamma","S.D.")=gammaestd;
res_table("gamma","Bias")=gammaebias;
res_table("gamma","MSE")=gammaemse;

res_table("alpha","True")=alphak_true;
res_table("alpha","Mean")=alphakem;
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

execute_unload 'rbc-vf-student-other.gdx',res,res_table,elapsed;




