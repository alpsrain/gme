*--------------------------------------------------------
* one sector model
*--------------------------------------------------------
scalar starttime; starttime = jnow;
sets nk  number of nodes for capital /1*5 / ;
sets nz  number of nodes for productivity shock / 1*5 / ;
alias (nk,nka,nkb) ;
alias (np,npa,npb,nz,nza,nzb)  ;

Set
    ttot simulated years / 1*500/
    t(ttot)              / 251*300/
    name /css, yss/
;

*Parameters
*seri(ttot,*)          ;
*$libinclude xlimport seri SimGrowthMaliar.xlsx onesector!a1:c501

Parameters
series(ttot,name)

$gdxIn SimGrowthMaliar_low.gdx
$LOAD series
$GDXIN


Parameters
cs(t)
ys(t)
;

cs(t)=series(t,"css");
ys(t)=series(t,"yss");

display cs, ys;

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

coeffsv(k,nka,nza)
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
pcoeffsv(k,nka,nza)
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
eqt(t,nz)
eqt1(t)
eqmeaneulerer
eqeulerer(t)
eqpeulerersv(t)
eqeulerneer(t,nz)
eqpeulerneersv(t,nz)
;


eqentropie..     entropie =e=     - sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - sum(k, pgammasv(k)*LOG(1.e-5+pgammasv(k)))
                                  - sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - sum(k, pbetasv(k)*LOG(1.e-5+pbetasv(k)))
                                  - sum((k,t), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - predict*sum((k,t), perreulerer(t,k)*log(1.e-5+perreulerer(t,k)))
                                  - predict1*sum((k,t,nz), perreulerneer(t,k,nz)*log(1.e-5+perreulerneer(t,k,nz)))
;

eqcs(t)..                         (   - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/betae+deltae-1)/alphake)*exp(zne(t,nz))*alphake*ke(t+1)**(alphake-1)
                                                                                                            )
                                                                      )
                                   )$(ord(t) LT card(t))
                                  +(  - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/betae+deltae-1)/alphake)*exp(zne(t,nz))*alphake*kfin**(alphake-1)
                                                                                                             )
                                                                      )
                                   )$(ord(t) EQ card(t))
                                    =E= 0 ;

* Modeling next period TFP expectation
eqzne(t,nz)..                        zne(t,nz) =E=  rhoze * ze(t) + sigmaze*err0(nz);

* Modeling TFP evolution process
eqze(t)..                            ze(t+1)   =E=  rhoze * ze(t) + errz(t);


eqke(t)..                             (ke(t+1) - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                      )$(ord(t) LT card(t))
                                    + (kfin    - ( (1-deltae)*ke(t) + ys(t) - cs(t))
                                      )$(ord(t) EQ card(t)) =E= 0  ;

eqmeanze..                          sum(t, errz(t))/card(t) =E= 0;

eqstdze..                           (sum(t, errz(t)*errz(t))/card(t))**0.5 =E= sigmaze  ;

eqys(t)..                           ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake;                                                                         ;


eqpolicy(t)..                        (cs(t)+ eulerer(t)) **(-gammae)  =E=  mu0
                                                            + mu1*ke(t)
                                                            + mu2*ze(t)
                                                            + 0.5*mu11*ke(t)*ke(t)
                                                            + 0.5*mu22*ze(t)*ze(t)
                                                            + mu12*ke(t)*ze(t)
                                                            ;

eqcne(t,nz)..                        (-(cne(t,nz)+ eulererne(t,nz))**(-gammae) + mu0
                                                            + mu1*ke(t+1)
                                                            + mu2*zne(t,nz)
                                                            + 0.5*mu11*ke(t+1)*ke(t+1)
                                                            + 0.5*mu22*zne(t,nz)*zne(t,nz)
                                                            + mu12*ke(t+1)*zne(t,nz)
                                                  )$(ord(t) LT card(t))
                                   + (-(cne(t,nz)+ eulererne(t,nz))**(-gammae) + mu0
                                                            + mu1*kfin
                                                            + mu2*zne(t,nz)
                                                            + 0.5*mu11*kfin*ke(t+1)
                                                            + 0.5*mu22*zne(t,nz)*zne(t,nz)
                                                            + mu12*kfin*zne(t,nz)
                                                  )$(ord(t) EQ card(t))
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

eqerrz(t)..           errz(t)    =E= sum(k, perrzsv(t,k)*errzsv(k)) ;
eqperrzsv(t)..        1          =E= sum(k, perrzsv(t,k)) ;

eqeulerer(t)..          eulerer(t) =E= sum(k, perreulerer(t,k)*erreulersv(k)) ;
eqpeulerersv(t)..       1          =E= sum(k, perreulerer(t,k)) ;

eqeulerneer(t,nz)..          eulererne(t,nz) =E= sum(k, perreulerneer(t,k,nz)*erreulernesv(k)) ;
eqpeulerneersv(t,nz)..       1          =E= sum(k, perreulerneer(t,k,nz)) ;


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
betasv("1")   = 0.92 ;
betasv("2")   = 0.95 ;
betasv("3")   = 0.99 ;
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.4 ;
alphaksv("3")   = 0.6 ;
gammasv("1")    = 0.1 ;
gammasv("2")    = 1.8;
gammasv("3")    = 4 ;
deltasv("1")    = 0.01 ;
deltasv("2")    = 0.05 ;
deltasv("3")    = 0.10 ;
rhozsv("1")    = 0.7 ;
rhozsv("2")    = 0.8 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2 ;

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

*gaussion quadrature std=1
err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;
prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;

z0(nz) = err0(nz)*sigmazsv("2") ;


Parameters
kbar
cbar
;

kbar   = 1;
cbar   = (1/betae.l+deltae.l-1)/alphake.l - deltae.l ;

zmax = z0("5")*1.5 ;
zmin = -zmax ;


errzsv("1")             = zmin ;
errzsv("3")             = zmax ;
erreulersv("1")         = -0.5*sum(t, cs(t))/card(t);
erreulersv("3")         = 0.5*sum(t, cs(t))/card(t);
erreulernesv("1")       = -0.5*sum(t, cs(t))/card(t);
erreulernesv("3")       = 0.5*sum(t, cs(t))/card(t);


ke.l(t)           = kbar  ;
kfin.l            = kbar ;
ke.lo(t)          = kbar*0.5 ;
ke.up(t)          = kbar*1.5 ;
ze.l(t)           = 0;
ze.up(t)          = zmax ;
ze.lo(t)          = zmin ;
cne.l(t,npa)      = cbar ;
cne.lo(t,npa)     = 0.5*cbar ;
cne.up(t,npa)     = 1.5*cbar ;

mu0.l = (sum(t, cs(t))/card(t) )**(-gammae.l) ;
*mu11.l = 0.1;
*mu22.l = 0.1;
*mu12.l = 0.1;

predict = 1 ;
predict1 = 1 ;

solve estimation using nlp maximising entropie;

erreulersv("1")         = -0.05*(sum(t, cs(t))/card(t) );
erreulersv("3")         = 0.05*(sum(t, cs(t))/card(t) );
erreulernesv("1")       = -0.05*(sum(t, cs(t))/card(t) );
erreulernesv("3")       = 0.05*(sum(t, cs(t))/card(t) );

solve estimation using nlp maximising entropie;

*execute_unload '1sector_compare.gdx';




