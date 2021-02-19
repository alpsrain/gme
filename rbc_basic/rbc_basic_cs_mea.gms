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

Parameters
beta,gamma,alphak,delta,rho,sigma,kbar,cbar;
beta=beta_true;
gamma=gamma_true;
alphak=alphak_true;
delta=delta_true;
rho=rho_true;
sigma=sigma_true;

Set t    / 1*10000/
    name /css, yss,kss,zss/
;

*Parameters
*series(t,name) ;
*$libinclude xlimport series SimGrowthMaliar.xlsx k20nodes!a1:e10001
*$libinclude xlimport series SimGrowthMaliar.xlsx cons_invmin80!a1:e10001


Parameters
series(t,name)
$gdxIn simdata.gdx
$LOAD series
$GDXIN

Parameters
cs(t)
ys(t)
zs(t)
invs(t)
ks(t)
;

cs(t)=series(t,"css");
ys(t)=series(t,"yss");
zs(t)=series(t,"zss");
ks(t)=series(t,"kss");
invs(t)= ys(t)-cs(t);

zs(t) = log(zs(t));



Parameters
cmean
ymean
invmean
;
cmean = sum(t, cs(t))/card(t) ;
ymean = sum(t, ys(t))/card(t) ;
invmean = sum(t,invs(t))/card(t);
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

Parameters
invs2(t)
diffinv
diffinv_mean
diffinv_std;
;
invs2(t) = ys(t)-cs(t);
diffinv(t) = invs2(t) - invs(t);
diffinv_mean = sum(t,diffinv(t))/card(t);
diffinv_std = (sum(t,(diffinv(t)-diffinv_mean)*(diffinv(t)-diffinv_mean))/card(t))**0.5;

display diffinv_std, invmean,cmean, ymean;
sets nk  number of nodes for capital /1*5/ ;
sets nz  number of nodes for productivity shock / 1*5/ ;
alias (nk,nka,nkb) ;
alias (np,npa,npb,nz,nza,nzb)  ;

Parameters
err0(nz)
prob(nz)
;

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

*------------------------------
* Compute Approximation errors
*------------------------------
*We need to take VF coefficients from Maliar because data is generating with VFI
Sets ncoeff number of coefficients   /1*21/;
Parameters
coeffs(ncoeff)
vfs(t)
;

coeffs("1") = -2180.15171495657 ;
coeffs("2")=1400.01608514611 ;
coeffs("3")=862.029760805889 ;
coeffs("4")=-1220.14528059860 ;
coeffs("5")=-783.588984561966 ;
coeffs("6")=-604.033513698975 ;
coeffs("7")=703.633878591297  ;
coeffs("8")=485.034417022351 ;
coeffs("9")=349.366106273512 ;
coeffs("10")=313.066641734630 ;
coeffs("11")=-231.209152546818 ;
coeffs("12")=-168.697609655073 ;
coeffs("13")=-125.177845924741 ;
coeffs("14")=-98.1276469190940 ;
coeffs("15")=-97.4290009177044 ;
coeffs("16")=32.8440027073006 ;
coeffs("17")=24.9985797392119 ;
coeffs("18")=19.3862071531376 ;
coeffs("19")=15.0111656816389 ;
coeffs("20")=12.6763031927360 ;
coeffs("21")=13.4660064430677 ;

Parameters
z1s
k1s
vfs_deriv1
u1s(t,nzb)
c1s(t,nzb)
Residuals (t)
Mean_Residuals;

* Future period quantities in n_nodes integration nodes (k1,a1)
z1s(t,nzb) = rho*zs(t)+sigma*err0(nzb);
k1s(t,nzb) =  (1-delta)*ks(t) + ys(t) - cs(t);
vfs_deriv1(t,nzb) = coeffs("2")
              + coeffs("4")*2*k1s(t,nzb)
              + coeffs("5")*exp(z1s(t,nzb))
              + coeffs("7")*3*k1s(t,nzb)*k1s(t,nzb)
              + coeffs("8")*2*k1s(t,nzb)*exp(z1s(t,nzb))
              + coeffs("9")  *exp(z1s(t,nzb))*exp(z1s(t,nzb))
              + coeffs("11")*4*k1s(t,nzb)*k1s(t,nzb)*k1s(t,nzb)
              + coeffs("12")*3*k1s(t,nzb)*k1s(t,nzb)*exp(z1s(t,nzb))
              + coeffs("13")*2*k1s(t,nzb)*exp(z1s(t,nzb))*exp(z1s(t,nzb))
              + coeffs("14")  *exp(z1s(t,nzb))*exp(z1s(t,nzb))*exp(z1s(t,nzb))
              + coeffs("16")*5*k1s(t,nzb)*k1s(t,nzb)*k1s(t,nzb)*k1s(t,nzb)
              + coeffs("17")*4*k1s(t,nzb)*k1s(t,nzb)*k1s(t,nzb)*exp(z1s(t,nzb))
              + coeffs("18")*3*k1s(t,nzb)*k1s(t,nzb)*exp(z1s(t,nzb))*exp(z1s(t,nzb))
              + coeffs("19")*2*k1s(t,nzb)*exp(z1s(t,nzb))*exp(z1s(t,nzb))*exp(z1s(t,nzb))
              + coeffs("20")*  exp(z1s(t,nzb))*exp(z1s(t,nzb))*exp(z1s(t,nzb))*exp(z1s(t,nzb))
;
u1s(t,nzb) = vfs_deriv1(t,nzb) / (1-delta + ((1/beta+delta-1)/alphak)*exp(z1s(t,nzb))*alphak*k1s(t,nzb)**(alphak-1)) ;
c1s(t,nzb) = u1s(t,nzb)**(-1/gamma) ;


Residuals (t)=sum(nzb, prob(nzb)*beta*c1s(t,nzb)**(-gamma)
                            *(1-delta+ ((1/beta+delta-1)/alphak)*exp(z1s(t,nzb))*alphak*k1s(t,nzb)**(alphak-1)
                             )
                      )
                  - cs(t)**(-gamma);
Mean_Residuals = sum(t,Residuals(t))/card(t);

display Residuals, Mean_Residuals;


*----------------------
* estimation GME
*----------------------
* selection of "smolyak" points at the second order of approximation for d=2 and mu=2 ;
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



parameters
z0(nz)
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

cne(t,nz)
zne(t,nz)
errz(t)
eulerer(t)
eulererne(t,nz)
erry(t)

phiz(nza,t)
phizne(nza,t,nz)
phik(nka,t)
coeffe(nka,nza)
phikfin(nka)
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
perreulerneer(t,k,nz)
perrysv(t,k)
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
eqcne(t,nz)
eqze(t)

eqcs(t)
eqys(t)
eqentropie

eqmeanze
eqstdze

eqerrz(t)
eqperrzsv(t)

eqphiz(nza,t)
eqphizne(nza,t,nz)
eqphik(nka,t)
eqpolicy(t)
eqphikfin(nka)
eqmeaneulerer

eqeulerer(t)
eqpeulerersv(t)
eqeulerneer(t,nz)
eqpeulerneersv(t,nz)

eqerry(t)
eqperrysv(t)

;

Parameters
tmin
tmax ;


eqentropie..     entropie =e=     - predict*sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - predict*sum(k, pgammasv(k)*LOG(1.e-5+pgammasv(k)))
                                  - predict*sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - predict*sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - predict*sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - predict*sum(k, pbetasv(k)*LOG(1.e-5+pbetasv(k)))
                                  - predict*sum((k,t)$((ord(t) GE tmin) and (ord(t) LE Tmax) ), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - sum((k,t)$((ord(t) GE tmin) and (ord(t) LE Tmax) ), perreulerer(t,k)*log(1.e-5+perreulerer(t,k)))
                                  - sum((k,t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) ), perreulerneer(t,k,nz)*log(1.e-5+perreulerneer(t,k,nz)))
                                  - sum((k,t)$((ord(t) GE tmin) and (ord(t) LE Tmax) ), perrysv(t,k)*log(1.e-5+perrysv(t,k)))
;

eqcs(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..         (   - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/betae+deltae-1)/alphake)*exp(zne(t,nz))*alphake*ke(t+1)**(alphake-1)
                                                                                                            ) )
                                                            )$(ord(t) LT tmax)
                                                          +(  - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/betae+deltae-1)/alphake)*exp(zne(t,nz))*alphake*kfin**(alphake-1)
                                                                                                             ))
                                                            )$(ord(t) EQ tmax)
                                                            =E= 0 ;

* Modeling next period TFP expectation
eqzne(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..        zne(t,nz) =E=  rhoze * ze(t) + sigmaze*err0(nz);

* Modeling TFP evolution process
eqze(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..            ze(t+1)   =E=  rhoze * ze(t) + errz(t);


eqke(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..            (ke(t+1) - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                      )$(ord(t) LT tmax)
                                    + (kfin    - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                      )$(ord(t) EQ tmax) =E= 0  ;

eqmeanze..                          sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), errz(t))/(tmax-tmin) =E= 0;

eqstdze..                           (sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), errz(t)*errz(t))/(tmax-tmin))**0.5 =E= sigmaze  ;

eqys(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..             ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake
                                                               + erry(t)
                                                               ;


eqphiz(nza,t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                     phiz(nza,t)     =E=                                                1$(ord(nza) eq 1)
                                                        + (2*ze(t)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax))$(ord(nza) eq 2)
                                                           + (2*phiz("2",t)*phiz(nza-1,t)-phiz(nza-2,t))$(ord(nza) GE 3) ;

eqphik(nka,t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                     phik(nka,t)     =E=                                                1$(ord(nka) eq 1)
                                                        + (2*ke(t)/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax))$(ord(nka) eq 2)
                                                           + (2*phik("2",t)*phik(nka-1,t)-phik(nka-2,t))$(ord(nka) GE 3) ;

eqphikfin(nka)..                    phikfin(nka)   =E=                                                 1$(ord(nka) eq 1)
                                                         + (2*kfin/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax))$(ord(nka) eq 2)
                                                        + (2*phikfin("2")*phikfin(nka-1)-phikfin(nka-2))$(ord(nka) GE 3) ;

eqphizne(nza,t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
                                    phizne(nza,t,nz)     =E=                                           1$(ord(nza) eq 1)
                                                    + (2*zne(t,nz)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax))$(ord(nza) eq 2)
                                            + (2*phizne("2",t,nz)*phizne(nza-1,t,nz)-phizne(nza-2,t,nz))$(ord(nza) GE 3) ;

eqpolicy(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                        (cs(t)- eulerer(t))**(-gammae) =E=  sum((nka,nza)$smol(nka,nza), coeffe(nka,nza)*phik(nka,t) * phiz(nza,t))  ;

eqcne(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                        (cne(t,nz)- eulererne(t,nz))**(-gammae) =E=  ( sum((nka,nza)$smol(nka,nza), coeffe(nka,nza)*phik(nka,t+1) * phizne(nza,t,nz))
                                                               )$(ord(t) LT tmax)
                                                             + ( sum((nka,nza)$smol(nka,nza), coeffe(nka,nza)*phikfin(nka) * phizne(nza,t,nz))
                                                               )$(ord(t) EQ tmax)
                                    ;

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

eqerrz(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..           errz(t)    =E= sum(k, perrzsv(t,k)*errzsv(k)) ;
eqperrzsv(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..        1          =E= sum(k, perrzsv(t,k)) ;

eqeulerer(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..          eulerer(t) =E= sum(k, perreulerer(t,k)*erreulersv(k)) ;
eqpeulerersv(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..       1          =E= sum(k, perreulerer(t,k)) ;

eqeulerneer(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..          eulererne(t,nz) =E= sum(k, perreulerneer(t,k,nz)*erreulernesv(k)) ;
eqpeulerneersv(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..       1          =E= sum(k, perreulerneer(t,k,nz)) ;

eqerry(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..           erry(t)    =E= sum(k, perrysv(t,k)*errysv(k)) ;
eqperrysv(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..        1          =E= sum(k, perrysv(t,k)) ;


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
eqphik
eqphizne
eqphiz
eqcne
eqpolicy
eqphikfin
eqeulerer
eqpeulerersv
eqeulerneer
eqpeulerneersv
eqerry
eqperrysv
/ ;


*initiate chebyshev coefficient for estimation
*$ontext
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
*$offtext

$ontext
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
rhozsv("1")    = 0.5 ;
rhozsv("2")    = 0.75 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2;
$offtext


Parameters
kbar
cbar
;

z0(nz) = err0(nz)*sigmazsv("3");

zmax = z0("5")*1.5 ;
zmin = -zmax ;

errzsv("1")             = zmin ;
errzsv("3")             = zmax ;

kbar   = 1;
kmax     = kbar*1.5 ;
kmin     = kbar*0.5 ;


sets boot / 1*100/ ;
parameters
gammaeboot(boot)
deltaeboot(boot)
alphakeboot(boot)
betaeboot(boot)
rhozeboot(boot)
sigmazeboot(boot)
modelboot(boot)
;


scalar tstep /100/;
loop(boot$(ord(boot)  ),
tmin = 1+(ord(boot)-1)*tstep ;
tmax = tstep+(ord(boot)-1)*tstep ;

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
kmax     = kbar*1.5 ;
kmin     = kbar*0.5 ;
ke.l(t)           = kbar  ;
kfin.l            = kbar ;
ke.lo(t)          = kmin;
ke.up(t)          = kmax ;
ze.l(t)           = 0;
ze.up(t)          = zmax ;
ze.lo(t)          = zmin ;
cne.l(t,npa)      = cbar ;
cne.lo(t,npa)     = 0.5*cbar ;
cne.up(t,npa)     = 1.5*cbar ;
zne.l(t,nz)       = 0 ;
errz.l(t)         = 0 ;
eulerer.l(t)      = 0 ;
eulererne.l(t,nz) = 0 ;
erry.l(t)         = 0 ;

pbetasv.l(k)            = 1/card(k) ;
palphaksv.l(k)          = 1/card(k) ;
pgammasv.l(k)           = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
psigmazsv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perreulerer.l(t,k)      = 1/card(k) ;
perreulerneer.l(t,k,nz) = 1/card(k) ;
perrysv.l(t,k)          = 1/card(k) ;

phik.l("1",t) = 1;
phik.l("2",t) = 2*ke.l(t)/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax);
loop(nka$(ord(nka) GE 3),
phik.l(nka,t) = 2*phik.l("2",t)*phik.l(nka-1,t)-phik.l(nka-2,t) ;
     );
phikfin.l("1") = 1 ;
phikfin.l("2") = 2*kbar/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax);
loop(nka$(ord(nka) GE 3),
phikfin.l(nka) = 2*phikfin.l("2")*phikfin.l(nka-1)-phikfin.l(nka-2) ;
     );
phiz.l("1",t) = 1;
phiz.l("2",t) = 2*ze.l(t)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax);
loop(nza$(ord(nza) GE 3),
phiz.l(nza,t) = 2*phiz.l("2",t)*phiz.l(nza-1,t)-phiz.l(nza-2,t) ;
     );
phizne.l("1",t,nz) = 1;
phizne.l("2",t,nz) = 2*ze.l(t)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax);
loop(nza$(ord(nza) GE 3),
phizne.l(nza,t,nz) = 2*phizne.l("2",t,nz)*phizne.l(nza-1,t,nz)-phizne.l(nza-2,t,nz) ;
     );

coeffe.l(nka,nza) =  0 ;
coeffe.l("1","1") = (sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin) )**(-gammae.l);

predict = 0.1 ;

erreulersv("1")         = -0.3*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulersv("3")         = 0.3*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulernesv("1")       = -0.3*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulernesv("3")       = 0.3*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
errysv("1")             = -0.1*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), ys(t))/(tmax-tmin);
errysv("3")             = 0.1*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), ys(t))/(tmax-tmin);
*solve estimation using nlp maximising entropie;

smol(nk,np)=1;
solve estimation using nlp maximising entropie;


gammaeboot(boot) = gammae.l ;
deltaeboot(boot) = deltae.l ;
alphakeboot(boot) = alphake.l ;
betaeboot(boot) = betae.l   ;
rhozeboot(boot) = rhoze.l   ;
sigmazeboot(boot) = sigmaze.l  ;
modelboot(boot) = estimation.solvestat ;

) ;

scalar elapsed; elapsed = (jnow - starttime)*24*3600;


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

alphakem = sum(boot$(modelboot(boot) LE 2), alphakeboot(boot))/card(boot$(modelboot(boot) LE 2)) ;
alphakestd = (sum(boot$(modelboot(boot) LE 2), (alphakeboot(boot)-alphakem)*(alphakeboot(boot)-alphakem) )/card(boot$(modelboot(boot) LE 2)) )**0.5 ;
gammaem = sum(boot$(modelboot(boot) LE 2), gammaeboot(boot))/card(boot$(modelboot(boot) LE 2)) ;
gammaestd = (sum(boot$(modelboot(boot) LE 2), (gammaeboot(boot)-gammaem)*(gammaeboot(boot)-gammaem) )/card(boot$(modelboot(boot) LE 2)) )**0.5 ;
deltaem = sum(boot$(modelboot(boot) LE 2), deltaeboot(boot))/card(boot$(modelboot(boot) LE 2)) ;
deltaestd = (sum(boot$(modelboot(boot) LE 2), (deltaeboot(boot)-deltaem)*(deltaeboot(boot)-deltaem) )/card(boot$(modelboot(boot) LE 2)) )**0.5 ;
betaem = sum(boot$(modelboot(boot) LE 2), betaeboot(boot))/card(boot$(modelboot(boot) LE 2)) ;
betaestd = (sum(boot$(modelboot(boot) LE 2), (betaeboot(boot)-betaem)*(betaeboot(boot)-betaem) )/card(boot$(modelboot(boot) LE 2)) )**0.5 ;
rhozem = sum(boot$(modelboot(boot) LE 2), rhozeboot(boot))/card(boot$(modelboot(boot) LE 2)) ;
rhozestd = (sum(boot$(modelboot(boot) LE 2), (rhozeboot(boot)-rhozem)*(rhozeboot(boot)-rhozem) )/card(boot$(modelboot(boot) LE 2)) )**0.5 ;
sigmazem = sum(boot$(modelboot(boot) LE 2), sigmazeboot(boot))/card(boot$(modelboot(boot) LE 2)) ;
sigmazestd = (sum(boot$(modelboot(boot) LE 2), (sigmazeboot(boot)-sigmazem)*(sigmazeboot(boot)-sigmazem) )/card(boot$(modelboot(boot) LE 2)) )**0.5 ;

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

display gammaeboot, deltaeboot, alphakeboot, betaeboot, rhozeboot, sigmazeboot, modelboot, elapsed   ;

scalar elapsed; elapsed = (jnow - starttime)*24*3600;

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

execute_unload 'rbc-cs-mea.gdx';
*$libinclude xlexport res res.xlsx res!a1:h101





