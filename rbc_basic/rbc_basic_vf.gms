*--------------------------------------------------------
* one sector model
*--------------------------------------------------------
scalar starttime; starttime = jnow;


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
;

cs(t)=series(t,"css");
ys(t)=series(t,"yss");

display cs, ys;


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

Sets ncoeff number of coefficients   /1*15/;
Parameters
coeffs(ncoeff)
vfs(t)
;

* no binding constraint
*coeffs("1")=-1934.34628076383;
*coeffs("2")=718.707415294864;
*coeffs("3")=441.254481673093;
*coeffs("4")=-306.909205844206;
*coeffs("5")=-198.494283185508;
*coeffs("6")=-151.579808698686;
*coeffs("7")=58.3862770368426;
*coeffs("8")=40.5711908272261;
*coeffs("9")=29.4726952446437;
*coeffs("10")=25.9841273503470;

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


*----------------------
* estimation GME
*----------------------
sets nk  number of nodes for capital /1*5 / ;
sets nz  number of nodes for productivity shock / 1*5/ ;
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

vf(t)
vf_deriv(t)
vfne(t,nz)
vfne_deriv(t,nz)

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

eqcs1(t)

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
                                  - sum((k,t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) ), perreulerneer(t,k)*log(1.e-5+perreulerneer(t,k)))
;

$ontext
eqvf(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
                                vf(t)     =E=               coeffe("1")
                                                            + coeffe("2")*ke(t)
                                                            + coeffe("3")*exp(ze(t))
                                                            + coeffe("4")*ke(t)*ke(t)
                                                            + coeffe("6")*exp(ze(t))*exp(ze(t))
                                                            + coeffe("5")*ke(t)*exp(ze(t))
                                                            + coeffe("7")*ke(t)*ke(t)*ke(t)
                                                            + coeffe("8")*ke(t)*ke(t)*exp(ze(t))
                                                            + coeffe("9")*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + coeffe("10")*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                            + coeffe("11")*ke(t)*ke(t)*ke(t)*ke(t)
                                                            + coeffe("12")*ke(t)*ke(t)*ke(t)*exp(ze(t))
                                                            + coeffe("13")*ke(t)*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + coeffe("14")*ke(t)*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                            + coeffe("15")*exp(ze(t))*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                            ;
$offtext

eqvf_deriv(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
                                 vf_deriv(t)             =E= coeffe("2")
                                                            + coeffe("4")*2*ke(t)
                                                            + coeffe("5")  *exp(ze(t))
                                                            + coeffe("7")*3*ke(t)*ke(t)
                                                            + coeffe("8")*2*ke(t)*exp(ze(t))
                                                            + coeffe("9")  *exp(ze(t))*exp(ze(t))
                                                            + coeffe("11")*4*ke(t)*ke(t)*ke(t)
                                                            + coeffe("12")*3*ke(t)*ke(t)*exp(ze(t))
                                                            + coeffe("13")*2*ke(t)*exp(ze(t))*exp(ze(t))
                                                            + coeffe("14")  *exp(ze(t))*exp(ze(t))*exp(ze(t))
*                                                            + coeffe("16")*5*ke(t)*ke(t)*ke(t)*ke(t)
*                                                            + coeffe("17")*4*ke(t)*ke(t)*ke(t)*exp(ze(t))
*                                                            + coeffe("18")*3*ke(t)*ke(t)*exp(ze(t))*exp(ze(t))
*                                                            + coeffe("19")*2*ke(t)*exp(ze(t))*exp(ze(t))*exp(ze(t))
*                                                            + coeffe("20")  *exp(ze(t))*exp(ze(t))*exp(ze(t))*exp(ze(t))
                                                                             ;

eqcs(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
                                    vf_deriv(t) =e= (cs(t)- eulerer(t))**(-gammae)*(1-deltae + ((1/betae+deltae-1)/alphake)*exp(ze(t))*alphake*ke(t)**(alphake-1)) ;


*eqcs(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
*                                    vf(t) =e= ((cs(t)+ eulerer(t))**(1-gammae)-1)/(1-gammae) + betae* sum(nz,prob(nz)*vfne(t,nz)) ;

$ontext
eqvfne(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
                                    vfne(t,nz)   =E=  (coeffe("1")
                                                            + coeffe("2")*ke(t+1)
                                                            + coeffe("3")*exp(zne(t,nz))
                                                            + coeffe("4")*ke(t+1)*ke(t+1)
                                                            + coeffe("6")*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("5")*ke(t+1)*exp(zne(t,nz))
                                                            + coeffe("7")*ke(t+1)*ke(t+1)*ke(t+1)
                                                            + coeffe("8")*ke(t+1)*ke(t+1)*exp(zne(t,nz))
                                                            + coeffe("9")*ke(t+1)*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("10")*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("11")*ke(t+1)*ke(t+1)*ke(t+1)*ke(t+1)
                                                            + coeffe("12")*ke(t+1)*ke(t+1)*ke(t+1)*exp(zne(t,nz))
                                                            + coeffe("13")*ke(t+1)*ke(t+1)*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("14")*ke(t+1)*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("15")*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                             )$(ord(t) LT tmax)
                                                             +
                                                             (coeffe("1")
                                                            + coeffe("2")*kfin
                                                            + coeffe("3")*exp(zne(t,nz))
                                                            + coeffe("4")*kfin*kfin
                                                            + coeffe("6")*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("5")*kfin*exp(zne(t,nz))
                                                            + coeffe("7")*kfin*kfin*kfin
                                                            + coeffe("8")*kfin*kfin*exp(zne(t,nz))
                                                            + coeffe("9")*kfin*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("10")*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("11")*kfin*kfin*kfin*kfin
                                                            + coeffe("12")*kfin*kfin*kfin*exp(zne(t,nz))
                                                            + coeffe("13")*kfin*kfin*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("14")*kfin*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                            + coeffe("15")*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))*exp(zne(t,nz))
                                                             )$(ord(t) EQ tmax)
                                                            ;
$offtext

eqvfne_deriv(t,nz)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
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
                                                            )$(ord(t) LT tmax)
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
                                                            )$(ord(t) EQ tmax)
                                                                             ;


eqvf_FOC(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..
                                   -(cs(t)-eulererne(t))**(-gammae) + betae*sum(nz,prob(nz)*(vfne_deriv(t,nz)))  =E=0 ;


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

eqys(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..            ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake ;                                                                         ;




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

eqeulerneer(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..          eulererne(t) =E= sum(k, perreulerneer(t,k)*erreulernesv(k)) ;
eqpeulerneersv(t)$((ord(t) GE tmin) and (ord(t) LE Tmax) )..       1          =E= sum(k, perreulerneer(t,k)) ;


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
gammasv("2")    = 3 ;
gammasv("3")    = 6;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.05 ;
deltasv("3")    = 0.10 ;
rhozsv("1")    = 0.5 ;
rhozsv("2")    = 0.75 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.15;
$offtext


*gaussion quadrature std=1
*err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;
*prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;

err0("1")= -2.85697001387281; err0("2")= -1.35562617997427; err0("3")= 0; err0("4")= 1.35562617997427; err0("5")= 2.85697001387281;
prob("1")=0.0112574113277207;prob("2")=0.222075922005613; prob("3")=0.533333333333333;prob("4")=0.222075922005613 ;prob("5")=0.0112574113277207 ;

*err0("1") = -4.145 ; err0("2") = -2.802 ; err0("3") = -1.637 ; err0("4") = -0.539 ;
*err0("8") = 4.145 ; err0("7") = 2.802 ; err0("6") = 1.637 ; err0("5") = 0.539 ;
*prob("1") = 0.0001 ; prob("2") = 0.0096 ; prob("3") = 0.1172 ; prob("4") = 0.3730 ;
*prob("8") = 0.0001 ; prob("7") = 0.0096 ; prob("6") = 0.1172 ; prob("5") = 0.3730 ;

z0(nz) = err0(nz)*sigmazsv("3");

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

parameter
vfderivss ;

scalar tstep /100/;
loop(boot$(ord(boot) ge 10 and  ord(boot) le 20 ),
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
zne.l(t,nz)       = 0 ;
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


vfderivss = (sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin) )**(-gammae.l)*
            (1-deltae.l + ((1/betae.l+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
display vfderivss  ;
*coeffe("1").l = 100;
*coeffe("2").l = vfderivss;
*display coeffe("2").l;

$ontext
coeffe.l(ncoeff)=coeffs(ncoeff);
$offtext

coeffe.l(ncoeff)=0;
*coeffe.l("1")=(sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin) )**(-gammae.l) ;
coeffe.l("2")=vfderivss;

predict = 0.01;

$ontext
erreulersv("1")         = -1e-1*vfderivss;
erreulersv("3")         = 1e-1*vfderivss ;
erreulernesv("1")       = -1e-1*vfderivss ;
erreulernesv("3")       = 1e-1*vfderivss;
$offtext

*$ontext
erreulersv("1")         = -0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulersv("3")         = 0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulernesv("1")       = -0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
erreulernesv("3")       = 0.2*sum(t$((ord(t) GE tmin) and (ord(t) LE Tmax) ), cs(t))/(tmax-tmin);
*$offtext

solve estimation using nlp maximising entropie;

*$ontext
erreulersv("1")         = -0.1;
erreulersv("3")         = 0.1 ;
erreulernesv("1")       = -0.1 ;
erreulernesv("3")       = 0.1;
*$offtext

*solve estimation using nlp maximising entropie;


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

*execute_unload 'res_smolyak_100_true.gdx';





