*--------------------------------------------------------
* one sector model
*--------------------------------------------------------
scalar starttime; starttime = jnow;

Set t    / 1*10000/
    name /css, yss,kss,zss/
;

Parameters
series(t,name) ;
$libinclude xlimport series SimGrowthMaliar.xlsx k20nodes!a1:e10001
*$libinclude xlimport series SimGrowthMaliar.xlsx cons_invmin80!a1:e10001

*Parameters
*series(t,name)

*$gdxIn SimGrowthMaliar_macro_low.gdx
*$LOAD series
*$GDXIN


Parameters
cs(t)
ys(t)
;

cs(t)=series(t,"css");
ys(t)=series(t,"yss");
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
*cs(t) = cs(t) + 0.025*normal(0,1)*cmean ;
*ys(t) = ys(t) + 0.025*normal(0,1)*ymean ;
display cs;
parameters
diffc(t) ;
diffc(t) =cs(t)-cs0(t) ;
display diffc ;

*----------------------
* estimation GME
*----------------------
sets nk  number of nodes for capital /1*5 / ;
sets nz  number of nodes for productivity shock / 1*5 / ;
alias (nk,nka,nkb) ;
alias (np,npa,npb,nz,nza,nzb)  ;

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

err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;
prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;

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


mu0
mu1
mu2
mu11
mu12
mu22
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

eqpolicy(t)
eqeulerer(t)
eqpeulerersv(t)
eqeulerneer(t,nz)
eqpeulerneersv(t,nz)
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
;

eqcs(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                         (   - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/betae+deltae-1)/alphake)*exp(zne(t,nz))*alphake*ke(t+1)**(alphake-1)
                                                                                                            )
                                                                      )
                                   )$(ord(t) LT tmax)
                                  +(  - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/betae+deltae-1)/alphake)*exp(zne(t,nz))*alphake*kfin**(alphake-1)
                                                                                                             )
                                                                      )
                                   )$(ord(t) EQ tmax)
                                    =E= 0 ;

* Modeling next period TFP expectation
eqzne(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                        zne(t,nz) =E=  rhoze * ze(t) + sigmaze*err0(nz);

* Modeling TFP evolution process
eqze(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                            ze(t+1)   =E=  rhoze * ze(t) + errz(t);


eqke(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                             (ke(t+1) - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                      )$(ord(t) LT tmax)
                                    + (kfin    - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                      )$(ord(t) EQ tmax) =E= 0  ;

eqmeanze..                          sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), errz(t))/(tmax-tmin) =E= 0;

eqstdze..                           (sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), errz(t)*errz(t))/(tmax-tmin))**0.5 =E= sigmaze  ;

eqys(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                           ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake;                                                                         ;


eqpolicy(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                        (cs(t)+ eulerer(t)) **(-gammae)  =E=  mu0
                                                            + mu1*ke(t)
                                                            + mu2*ze(t)
                                                            + 0.5*mu11*ke(t)*ke(t)
                                                            + 0.5*mu22*ze(t)*ze(t)
                                                            + mu12*ke(t)*ze(t)
                                                            ;

eqcne(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..                        (-(cne(t,nz)+ eulererne(t,nz))**(-gammae) + mu0
                                                            + mu1*ke(t+1)
                                                            + mu2*zne(t,nz)
                                                            + 0.5*mu11*ke(t+1)*ke(t+1)
                                                            + 0.5*mu22*zne(t,nz)*zne(t,nz)
                                                            + mu12*ke(t+1)*zne(t,nz)
                                                  )$(ord(t) LT tmax)
                                   + (-(cne(t,nz)+ eulererne(t,nz))**(-gammae) + mu0
                                                            + mu1*kfin
                                                            + mu2*zne(t,nz)
                                                            + 0.5*mu11*kfin*ke(t+1)
                                                            + 0.5*mu22*zne(t,nz)*zne(t,nz)
                                                            + mu12*kfin*zne(t,nz)
                                                  )$(ord(t) EQ tmax)
                                   =E= 0 ;

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
eqcne
eqpolicy
eqeulerer
eqpeulerersv
eqeulerneer
eqpeulerneersv
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
rhozsv("1")    = 0.5 ;
rhozsv("2")    = 0.75 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2;
*$offtext


*gaussion quadrature std=1
err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;
prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;

z0(nz) = err0(nz)*sigmazsv("3") ;


Parameters
kbar
cbar
;


zmax = z0("5")*1.5 ;
zmin = -zmax ;

errzsv("1")             = zmin ;
errzsv("3")             = zmax ;

sets boot / 1*20/ ;
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
loop(boot$(ord(boot) eq 2  ),
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
pbetasv.l(k)            = 1/card(k) ;
palphaksv.l(k)          = 1/card(k) ;
pgammasv.l(k)           = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
psigmazsv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perreulerer.l(t,k)      = 1/card(k) ;
perreulerneer.l(t,k,nz) = 1/card(k) ;

mu0.l = (sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin) )**(-gammae.l) ;
*mu11.l = 0.1;
*mu22.l = 0.1;
*mu12.l = 0.1;


predict = 0.01 ;

erreulersv("1")         = -0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulersv("3")         = 0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulernesv("1")       = -0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulernesv("3")       = 0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
solve estimation using nlp maximising entropie;

gammaeboot(boot) = gammae.l ;
deltaeboot(boot) = deltae.l ;
alphakeboot(boot) = alphake.l ;
betaeboot(boot) = betae.l   ;
rhozeboot(boot) = rhoze.l   ;
sigmazeboot(boot) = sigmaze.l  ;
modelboot(boot) = estimation.solvestat ;

) ;

Parameters
beta_true,gamma_true,alphak_true,delta_true,rho_true,sigma_true;
beta_true=0.99;
gamma_true=2;
alphak_true=0.36;
delta_true=0.025;
rho_true=0.85;
sigma_true=0.04;

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

alphakem = sum(boot, alphakeboot(boot))/card(boot) ;
alphakestd = (sum(boot, (alphakeboot(boot)-alphakem)*(alphakeboot(boot)-alphakem) )/card(boot) )**0.5 ;
gammaem = sum(boot, gammaeboot(boot))/card(boot) ;
gammaestd = (sum(boot, (gammaeboot(boot)-gammaem)*(gammaeboot(boot)-gammaem) )/card(boot) )**0.5 ;
deltaem = sum(boot, deltaeboot(boot))/card(boot) ;
deltaestd = (sum(boot, (deltaeboot(boot)-deltaem)*(deltaeboot(boot)-deltaem) )/card(boot) )**0.5 ;
betaem = sum(boot, betaeboot(boot))/card(boot) ;
betaestd = (sum(boot, (betaeboot(boot)-betaem)*(betaeboot(boot)-betaem) )/card(boot) )**0.5 ;
rhozem = sum(boot, rhozeboot(boot))/card(boot) ;
rhozestd = (sum(boot, (rhozeboot(boot)-rhozem)*(rhozeboot(boot)-rhozem) )/card(boot) )**0.5 ;
sigmazem = sum(boot, sigmazeboot(boot))/card(boot) ;
sigmazestd = (sum(boot, (sigmazeboot(boot)-sigmazem)*(sigmazeboot(boot)-sigmazem) )/card(boot) )**0.5 ;

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

display alphakem, alphakestd, alphakemse ;
display gammaem, gammaestd, gammaemse ;
display deltaem, deltaestd, deltaemse ;
display betaem, betaestd, betaemse ;
display rhozem, rhozestd, rhozemse ;
display sigmazem, sigmazestd, sigmazemse ;

display gammaeboot, deltaeboot, alphakeboot, betaeboot, rhozeboot, sigmazeboot, modelboot ;

*execute_unload 'macro_low.gdx';

*execute_unload '1sector_monomial_low.gdx';




