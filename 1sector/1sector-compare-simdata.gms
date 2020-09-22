*--------------------------------------------------------
* one sector model
*--------------------------------------------------------
scalar starttime; starttime = jnow;

sets nk  number of nodes for capital /1*5 / ;
sets nz  number of nodes for productivity shock / 1*5 / ;
alias (nk,nka,nkb) ;
alias (np,npa,npb,nz,nza,nzb)  ;

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
smoltemp(nk,np) = sum((jk,jp)$(     ( (ord(jk)+ord(jp)) GE 2)
                                    and ( (ord(jk)+ord(jp)) LE 4)
                                         ), sk(jk,nk)*sp(jp,np) )  ;
smol(nk,np) = 1$smoltemp(nk,np) ;
*smol(nk,np)=1;
smolcount = sum((nk,np), smol(nk,np) ) ;
display smol, smolcount ;



Set
    ttot simulated years / 1*500/
    t(ttot)              / 1*50/
    name /css, yss/
;

*Parameters
*seri(ttot,*)          ;
*$libinclude xlimport seri SimGrowthMaliar.xlsx onesector!a1:c501

Parameters
series(ttot,name)

$gdxIn SimGrowthMaliar_Method6_D2.gdx
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


eqphiz(nza,t)..                     phiz(nza,t)     =E=                                                1$(ord(nza) eq 1)
                                                        + (2*ze(t)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax))$(ord(nza) eq 2)
                                                           + (2*phiz("2",t)*phiz(nza-1,t)-phiz(nza-2,t))$(ord(nza) GE 3) ;

eqphik(nka,t)..                     phik(nka,t)     =E=                                                1$(ord(nka) eq 1)
                                                        + (2*ke(t)/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax))$(ord(nka) eq 2)
                                                           + (2*phik("2",t)*phik(nka-1,t)-phik(nka-2,t))$(ord(nka) GE 3) ;

eqphikfin(nka)..                    phikfin(nka)   =E=                                                 1$(ord(nka) eq 1)
                                                         + (2*kfin/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax))$(ord(nka) eq 2)
                                                        + (2*phikfin("2")*phikfin(nka-1)-phikfin(nka-2))$(ord(nka) GE 3) ;

eqphizne(nza,t,nz)..
                                    phizne(nza,t,nz)     =E=                                           1$(ord(nza) eq 1)
                                                    + (2*zne(t,nz)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax))$(ord(nza) eq 2)
                                            + (2*phizne("2",t,nz)*phizne(nza-1,t,nz)-phizne(nza-2,t,nz))$(ord(nza) GE 3) ;

eqpolicy(t)..                        (cs(t)+ eulerer(t))**(-gammae)  =E=  sum((nka,nza)$smol(nka,nza), coeffe(nka,nza)*phik(nka,t) * phiz(nza,t))  ;

eqcne(t,nz)..                        (cne(t,nz)+ eulererne(t,nz))**(-gammae) =E=  ( sum((nka,nza)$smol(nka,nza), coeffe(nka,nza)*phik(nka,t+1) * phizne(nza,t,nz))
                                                               )$(ord(t) LT card(t))
                                                             + ( sum((nka,nza)$smol(nka,nza), coeffe(nka,nza)*phikfin(nka) * phizne(nza,t,nz))
                                                               )$(ord(t) EQ card(t))
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

eqerrz(t)..           errz(t)    =E= sum(k, perrzsv(t,k)*errzsv(k)) ;
eqperrzsv(t)..        1          =E= sum(k, perrzsv(t,k)) ;
eqmeaneulerer..         sum(t, eulerer(t))/card(t) =E= 0;

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
eqphik
eqphizne
eqphiz
eqcne
eqpolicy
eqphikfin
eqmeaneulerer
eqeulerer
eqpeulerersv
eqeulerneer
eqpeulerneersv
/ ;


*initiate chebyshev coefficient for estimation
*$ontext
betasv("1")   = 0.92 ;
betasv("2")   = 0.95 ;
betasv("3")   = 0.99 ;
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.4 ;
alphaksv("3")   = 0.6 ;
gammasv("1")    = 0.1 ;
gammasv("2")    = 2;
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
*$offtext


gammae.l  = gammasv("2") ;
deltae.l  = deltasv("2") ;
alphake.l = alphaksv("2") ;
betae.l   = betasv("2") ;
rhoze.l    = rhozsv("2") ;
sigmaze.l  = sigmazsv("2") ;

*gammae.l  = 2 ;
*deltae.l  = 0.025 ;
*alphake.l = 0.36 ;
*betae.l   = 0.95 ;
*rhoze.l    = 0.85 ;
*sigmaze.l  = 0.04 ;

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
*z0(nz) = gq(nz)*0.1 ;

Parameters
kbar
cbar
;

kbar   = 1;
cbar   = (1/betae.l+deltae.l-1)/alphake.l - deltae.l ;

kmax     = kbar*1.5 ;
kmin     = kbar*0.5 ;

zmax = z0("5")*1.5 ;
zmin = -zmax ;

*$ontext
errzsv("1")             = zmin ;
errzsv("3")             = zmax ;
erreulersv("1")         = -0.5*sum(t, cs(t))/card(t);
erreulersv("3")         = 0.5*sum(t, cs(t))/card(t);
erreulernesv("1")       = -0.5*sum(t, cs(t))/card(t);
erreulernesv("3")       = 0.5*sum(t, cs(t))/card(t);
*$offtext
display errzsv,erreulersv,erreulernesv;

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

$ontext
*-----------------------------------------------
* solve the model to find approximated coeffe.l
*-----------------------------------------------
parameters
gamma, delta, alphak, beta, rho, sigma, A, cbar0
;
gamma  = gammasv("2") ;
delta  = deltasv("2") ;
alphak = alphaksv("2") ;
beta   = betasv("2") ;
rho    = rhozsv("2") ;
sigma  = sigmazsv("2") ;
A      = ((1/beta+delta-1)/alphak);
cbar0   = A - delta ;

Parameters
k0(nk), phi0k(nka,nk),phi0z(nza,nz),z1(nz,nzb),phi1z(nza,nz,nzb), z0max, z0min, k0max, k0min, z0cs(nz)
;

k0max = kmax ;
k0min = kmin ;
k0(nk)  = (k0max+k0min)/2 + (k0max-k0min)/2*  cos(  pi * (1+ (0.5 - ord(nk))/card(nk) ) ) ;
Parameters
k0temp(nk) ;
k0temp("1") = k0("3") ;
k0temp("2") = k0("1") ;
k0temp("3") = k0("5") ;
k0temp("4") = k0("2") ;
k0temp("5") = k0("4") ;
k0(nk) = k0temp(nk) ;

phi0k("1",nk) = 1;
phi0k("2",nk) = 2*k0(nk)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max);
loop(nka$(ord(nka) GE 3),
phi0k(nka,nk) = 2*phi0k("2",nk)*phi0k(nka-1,nk)-phi0k(nka-2,nk);
     );

z0max = zmax ;
z0min = zmin ;
phi0z("1",nz) = 1;
phi0z("2",nz) = cos((1+ (0.5 - ord(nz))/card(nz))*pi);
Parameters
phi0ztemp(npa,np) ;
phi0ztemp("2","1") = phi0z("2","3") ;
phi0ztemp("2","2") = phi0z("2","1") ;
phi0ztemp("2","3") = phi0z("2","5") ;
phi0ztemp("2","4") = phi0z("2","2") ;
phi0ztemp("2","5") = phi0z("2","4") ;
phi0z("2",nz) = phi0ztemp("2",nz) ;
loop(nza$(ord(nza) GE 3),
phi0z(nza,nz) = 2*phi0z("2",nz)*phi0z(nza-1,nz)-phi0z(nza-2,nz) ;
     );

z0cs(nz)    = (z0max+z0min+(z0max-z0min)*cos((1+ (0.5 - ord(nz))/card(nz))*pi))/2;
Parameters
z0cstemp(np) ;
z0cstemp("1") = z0cs("3") ;
z0cstemp("2") = z0cs("1") ;
z0cstemp("3") = z0cs("5") ;
z0cstemp("4") = z0cs("2") ;
z0cstemp("5") = z0cs("4") ;
z0cs(nz) = z0cstemp(nz) ;

z1(nz,nzb) = rho*z0cs(nz) + z0(nzb);
phi1z("1",nz,nzb) = 1 ;
phi1z("2",nz,nzb) = 2*z1(nz,nzb)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax);
loop(nza$(ord(nza) GE 3),
phi1z(nza,nz,nzb) = 2*phi1z("2",nz,nzb)*phi1z(nza-1,nz,nzb)-phi1z(nza-2,nz,nzb) ;
     );


*----------------------
* projection definition
*----------------------
Variables
c(nk,nz)
k1(nk,nz)
c1(nk,nz,nzb)
phi1k(nka,nk,nz)
coeff(nk,nz)
;

Equations
eqc(nk,nz)
eqk1(nk,nz)
eqcoeff(nk,nz)
eqc1(nk,nz,nzb)
eqphi1k(nka,nk,nz)
  ;


eqk1(nk,nz)$smol(nk,nz)..            k1(nk,nz) =E= (1-delta)*k0(nk) + ((1/beta+delta-1)/alphak)*exp(z0cs(nz))*k0(nk)**(alphak) - c(nk,nz) ;

eqphi1k(nka,nk,nz)$smol(nk,nz)..     phi1k(nka,nk,nz) =E=                                                                  1$(ord(nka) eq 1)
                                                     + ( 2*k1(nk,nz)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max))$(ord(nka) eq 2)
                                                    + (2*phi1k("2",nk,nz)*phi1k(nka-1,nk,nz)-phi1k(nka-2,nk,nz))$(ord(nka) GE 3) ;


eqc(nk,nz)$smol(nk,nz)..             c(nk,nz)**(-gamma)  =E= sum((nka,nza)$smol(nka,nza), coeff(nka,nza)*phi0k(nka,nk) * phi0z(nza,nz) ) ;

eqc1(nk,nz,nzb)$smol(nk,nz)..        c1(nk,nz,nzb)**(-gamma) =E= sum((nka,nza)$smol(nka,nza), coeff(nka,nza)*phi1k(nka,nk,nz) * phi1z(nza,nz,nzb) ) ;

eqcoeff(nk,nz)$smol(nk,nz)..         c(nk,nz)**(-gamma)/ beta =E= sum(nzb, prob(nzb) * (
                                                                        c1(nk,nz,nzb)**(-gamma)*( 1 - delta + ((1/beta+delta-1)/alphak)*exp(z1(nz,nzb))*alphak*k1(nk,nz)**(alphak-1) )

                                                                        )
                                                           ) ;

model projection /
eqc.c
eqk1.k1
eqcoeff.coeff
eqc1.c1
eqphi1k.phi1k
/ ;

coeff.l(nka,nza) = 0.1 ;
coeff.l("1","1") = cbar0**(-gamma) ;
c.l(nk,nz) = cbar0 ;
k1.l(nk,nz) = kbar ;

phi1k.l("1",nk,nz) =  1 ;
phi1k.l("2",nk,nz) =  2*k1.l(nk,nz)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max)  ;
loop(nka$(ord(nka) GE 3),
phi1k.l(nka,nk,nz) = 2*phi1k.l("2",nk,nz)*phi1k.l(nka-1,nk,nz)-phi1k.l(nka-2,nk,nz) ;
     );
c1.l(nk,nz,nzb) = cbar0 ;

solve projection using mcp ;

rho    = 0.5 ;
z1(nz,nzb) = rho*z0cs(nz) + z0(nzb);
phi1z("1",nz,nzb) = 1 ;
phi1z("2",nz,nzb) = 2*z1(nz,nzb)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax);
loop(nza$(ord(nza) GE 3),
phi1z(nza,nz,nzb) = 2*phi1z("2",nz,nzb)*phi1z(nza-1,nz,nzb)-phi1z(nza-2,nz,nzb) ;
     );
solve projection using mcp ;

rho    = 0.85 ;
z1(nz,nzb) = rho*z0cs(nz) + z0(nzb);
phi1z("1",nz,nzb) = 1 ;
phi1z("2",nz,nzb) = 2*z1(nz,nzb)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax);
loop(nza$(ord(nza) GE 3),
phi1z(nza,nz,nzb) = 2*phi1z("2",nz,nzb)*phi1z(nza-1,nz,nzb)-phi1z(nza-2,nz,nzb) ;
     );
solve projection using mcp ;

rho    = rhozsv("2");
z1(nz,nzb) = rho*z0cs(nz) + z0(nzb);
phi1z("1",nz,nzb) = 1 ;
phi1z("2",nz,nzb) = 2*z1(nz,nzb)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax);
loop(nza$(ord(nza) GE 3),
phi1z(nza,nz,nzb) = 2*phi1z("2",nz,nzb)*phi1z(nza-1,nz,nzb)-phi1z(nza-2,nz,nzb) ;
     );
solve projection using mcp ;

$offtext

coeffe.l("1","1") = (sum(t, cs(t))/card(t) )**(-gammae.l) ;
*coeffe.l(nka,nza)$smol(nka,nza) = coeff.l(nka,nza)$smol(nka,nza);

predict = 1 ;
predict1 = 1 ;

solve estimation using nlp maximising entropie;
*smol(nka,nza) = 1 ;
*solve estimation using nlp maximising entropie;

erreulersv("1")         = -0.1*(sum(t, cs(t))/card(t) );
erreulersv("3")         = 0.1*(sum(t, cs(t))/card(t) );
erreulernesv("1")       = -0.1*(sum(t, cs(t))/card(t) );
erreulernesv("3")       = 0.1*(sum(t, cs(t))/card(t) );

smol(nka,nza) = 1 ;
solve estimation using nlp maximising entropie;

*execute_unload '1sector_compare.gdx';




