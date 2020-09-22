*--------------------------------------------------------
* one sector model
*--------------------------------------------------------
*------------------------------
* definition of deep parameters
*------------------------------
scalar starttime; starttime = jnow;

Parameters
gamma_true    risk aversion
delta_true    capital depreciation
alphak_true   marginal capital productivity
beta_true     time preference
rho_true      shock autocorrelation
sigma_true    standard error of the innovation
A_true        technology level to normalize capital (to be consistant with the Judd2017)
;

gamma_true  = 2;
delta_true  = 0.025 ;
alphak_true = 0.36 ;
beta_true   = 0.95 ;
rho_true    = 0;
sigma_true  = 1 ;
A_true = (1/beta_true+delta_true-1)/alphak_true;


Parameters
gamma    risk aversion
delta    capital depreciation
alphak   marginal capital productivity
beta     time preference
rho      shock autocorrelation
sigma    standard error of the innovation
kbar     determinist capital steady state
cbar     determinist consumption steady state
ybar     output steady state
A        technology level to normalize capital (to be consistant with the Judd2017)
;

gamma  = gamma_true ;
delta  = delta_true ;
alphak = alphak_true ;
beta   = beta_true ;

*gamma  = 0.9;
*delta  = 0.05 ;
*alphak = 0.5 ;
*beta   = 0.95 ;

rho    = rho_true;
sigma  = sigma_true ;
A      = (1/beta+delta-1)/alphak;
kbar   = 1 ;
cbar   = A - delta ;
ybar   = A ;
display kbar, cbar, ybar ;


*---------------------------
* parameters for projection
*---------------------------
sets nk  number of nodes for capital /1*5 / ;
sets np  number of nodes for price shock / 1*5 / ;
Set boot /1*1/  ;
Sets t simulated years / 51*250/ ;

alias (nk,nka,nkb) ;
alias (np,npa,npb,nz,nza,nzb)  ;

Parameters
k0(nk)      initial capital stock
phi0k(nka,nk)
p0(np)      initial market price level
phi0p(npa,np)
p0n(np,npb)  next period  price level
phip0n(npa,np,npb)
pmax
pmin
k0max
k0min
prob(np)
p0cs(np)    initial market level by cs
;

k0max = kbar*1.2 ;
k0min = kbar*0.8 ;
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
phi0k(nka,nk) = 2*phi0k("2",nk)*phi0k(nka-1,nk)-phi0k(nka-2,nk) ;
     );


** Gaussian quadrature
* std = 0.1
p0("1")=-0.2857; p0("2")=-0.1356; p0("3")=0; p0("4")=0.1356;p0("5")=0.2857;
prob("1")=0.0113; prob("2")=0.2221; prob("3")=0.5333; prob("4")=0.2221; prob("5")=0.0113;
p0(np)=p0(np)*0.8;
pmax = p0("5")*1.2 ;
pmin = -pmax ;
phi0p("1",np) = 1;
phi0p("2",np) = cos((1+ (0.5 - ord(np))/card(np))*pi);
Parameters
phi0ptemp(npa,np) ;
phi0ptemp("2","1") = phi0p("2","3") ;
phi0ptemp("2","2") = phi0p("2","1") ;
phi0ptemp("2","3") = phi0p("2","5") ;
phi0ptemp("2","4") = phi0p("2","2") ;
phi0ptemp("2","5") = phi0p("2","4") ;
phi0p("2",np) = phi0ptemp("2",np) ;
loop(npa$(ord(npa) GE 3),
phi0p(npa,np) = 2*phi0p("2",np)*phi0p(npa-1,np)-phi0p(npa-2,np) ;
     );
p0cs(np)    = (pmax+pmin+(pmax-pmin)*cos((1+ (0.5 - ord(np))/card(np))*pi))/2;
Parameters
p0cstemp(np) ;
p0cstemp("1") = p0cs("3") ;
p0cstemp("2") = p0cs("1") ;
p0cstemp("3") = p0cs("5") ;
p0cstemp("4") = p0cs("2") ;
p0cstemp("5") = p0cs("4") ;
p0cs(np) = p0cstemp(np) ;

*next period price p0n
p0n(np,npb) = rho*p0cs(np) + sigma*p0(npb);
phip0n("1",np,npb) = 1 ;
phip0n("2",np,npb) = 2*p0n(np,npb)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phip0n(npa,np,npb) = 2*phip0n("2",np,npb)*phip0n(npa-1,np,npb)-phip0n(npa-2,np,npb) ;
     );

display p0n ;


$ontext
*----------------------
* smolyak determination
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
smoltemp(nk,np) = sum((jk,jp)$(     ( (ord(jk)+ord(jp)) GE 2)
                                    and ( (ord(jk)+ord(jp)) LE 4)
                                                   ), sk(jk,nk)*sp(jp,np) )  ;
smol(nk,np) = 1$smoltemp(nk,np) ;
smol(nk,np)=1;
smolcount = sum((nk,np), smol(nk,np) ) ;
display smol, smolcount ;


*------------------
* model definition
*------------------
Variables
c(nk,np)
k0n(nk,np)
c0n(nk,np,npb)
phik0n(nka,nk,np)
coeff(nk,np)
;

Equations
eqc(nk,np)
eqk0n(nk,np)
eqcoeff(nk,np)
eqc0n(nk,np,npb)
eqphik0n(nka,nk,np)
  ;


* Complete polynomail definition

*eqk0n(nk,np)..            k0n(nk,np) =E= (1-delta)*k0(nk) + exp(p0cs(np))*k0(nk)**(alphak) - c(nk,np) ;

*eqphik0n(nka,nk,np)..     phik0n(nka,nk,np) =E=                                                                  1$(ord(nka) eq 1)
*                                                     + ( 2*k0n(nk,np)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max))$(ord(nka) eq 2)
*                                                    + (2*phik0n("2",nk,np)*phik0n(nka-1,nk,np)-phik0n(nka-2,nk,np))$(ord(nka) GE 3) ;

*eqc(nk,np)..             c(nk,np)**(-gamma) =E= sum((nka,npa), coeff(nka,npa)*phi0k(nka,nk) * phi0p(npa,np) )  ;

*eqc0n(nk,np,npb)..       c0n(nk,np,npb)**(-gamma) =E= sum((nka,npa), coeff(nka,npa)*phik0n(nka,nk,np) * phip0n(npa,np,npb) ) ;;

*eqcoeff(nk,np)..         c(nk,np)**(-gamma) / beta =E= sum(npb, prob(npb) * (
*                                                                        c0n(nk,np,npb)**(-gamma)*( 1 - delta + exp(p0n(np,npb))*alphak*k0n(nk,np)**(alphak-1) )
*
*                                                                        )
*                                                           ) ;

*Smolyak definition

eqk0n(nk,np)$smol(nk,np) ..            k0n(nk,np) =E= (1-delta)*k0(nk) + ((1/beta+delta-1)/alphak)*exp(p0cs(np))*k0(nk)**(alphak) - c(nk,np) ;

eqphik0n(nka,nk,np)$smol(nk,np) ..     phik0n(nka,nk,np) =E=                                                                  1$(ord(nka) eq 1)
                                                     + ( 2*k0n(nk,np)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max))$(ord(nka) eq 2)
                                                    + (2*phik0n("2",nk,np)*phik0n(nka-1,nk,np)-phik0n(nka-2,nk,np))$(ord(nka) GE 3) ;

eqc(nk,np)$smol(nk,np) ..             c(nk,np)**(-gamma) =E= sum((nka,npa)$smol(nka,npa), coeff(nka,npa)*phi0k(nka,nk) * phi0p(npa,np) )  ;

eqc0n(nk,np,npb)$smol(nk,np) ..       c0n(nk,np,npb)**(-gamma) =E= sum((nka,npa)$smol(nka,npa), coeff(nka,npa)*phik0n(nka,nk,np) * phip0n(npa,np,npb) ) ;;

eqcoeff(nk,np)$smol(nk,np) ..         c(nk,np)**(-gamma) / beta =E= sum(npb, prob(npb) * (
                                                                        c0n(nk,np,npb)**(-gamma)*( 1 - delta + ((1/beta+delta-1)/alphak)*exp(p0n(np,npb))*alphak*k0n(nk,np)**(alphak-1) )

                                                                        )
                                                           ) ;


model projection /
eqc.c
eqk0n.k0n
eqcoeff.coeff
eqc0n.c0n
eqphik0n.phik0n
/ ;

coeff.l(nka,npa)$smol(nka,npa) = 0.1 ;
coeff.l("1","1") = cbar**(-gamma) ;
c.l(nk,np) = cbar ;
k0n.l(nk,np) = kbar ;

phik0n.l("1",nk,np) =  1 ;
phik0n.l("2",nk,np) =  2*k0n.l(nk,np)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max)  ;
loop(nka$(ord(nka) GE 3),
phik0n.l(nka,nk,np) = 2*phik0n.l("2",nk,np)*phik0n.l(nka-1,nk,np)-phik0n.l(nka-2,nk,np) ;
     );
c0n.l(nk,np,npb) = cbar ;

solve projection using mcp ;

rho    = 0.5 ;
p0n(np,npb) = rho*p0cs(np) + p0(npb);
phip0n("1",np,npb) = 1 ;
phip0n("2",np,npb) = 2*p0n(np,npb)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phip0n(npa,np,npb) = 2*phip0n("2",np,npb)*phip0n(npa-1,np,npb)-phip0n(npa-2,np,npb) ;
     );
solve projection using mcp ;

rho    = 0.7 ;
p0n(np,npb) = rho*p0cs(np) + p0(npb);
phip0n("1",np,npb) = 1 ;
phip0n("2",np,npb) = 2*p0n(np,npb)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phip0n(npa,np,npb) = 2*phip0n("2",np,npb)*phip0n(npa-1,np,npb)-phip0n(npa-2,np,npb) ;
     );
solve projection using mcp ;

rho    = 0.8 ;
p0n(np,npb) = rho*p0cs(np) + p0(npb);
phip0n("1",np,npb) = 1 ;
phip0n("2",np,npb) = 2*p0n(np,npb)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phip0n(npa,np,npb) = 2*phip0n("2",np,npb)*phip0n(npa-1,np,npb)-phip0n(npa-2,np,npb) ;
     );
solve projection using mcp ;

rho    = 0.85 ;
p0n(np,npb) = rho*p0cs(np) + p0(npb);
phip0n("1",np,npb) = 1 ;
phip0n("2",np,npb) = 2*p0n(np,npb)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phip0n(npa,np,npb) = 2*phip0n("2",np,npb)*phip0n(npa-1,np,npb)-phip0n(npa-2,np,npb) ;
     );
solve projection using mcp ;

Parameters
coeffs(nka,npa) ;
coeffs(nka,npa) = coeff.l(nka,npa)$(abs(coeff.l(nka,npa)) GE 1.0E-3) ;

display coeffs
$offtext



$ontext

*----------------
* data simulation
*----------------
Parameters
pss(t,boot)
*lambdass(t,boot)
css(t,boot)
yss(t,boot)
kss(t,boot)
invss(t,boot)

*invss(t)
phi0kss(nka,t,boot)
phik0nss(nka,t,boot)
*lambda1ss(t,npb,boot)
c0nss(t,boot)
verif(t,boot)
phi0pss(npa,t,boot)
phip0nss(npa,t,boot);

pss(t,boot) = 0 ;
*execseed = 1e8*(frac(jnow));

loop(boot,
loop(t$(ord(t) LT card(t)),
pss(t+1,boot) = rho*pss(t,boot) + normal(0,0.04) ;
) ;
*execseed = 1 + gmillisec(jnow);
) ;

*display pss

parameters
meanpss(boot)
stdpss(boot) ;
meanpss(boot) = sum(t, pss(t,boot))/card(t) ;
stdpss(boot) = (sum(t, (pss(t,boot)-meanpss(boot))*(pss(t,boot)-meanpss(boot)))/card(t) )**0.5 ;
display meanpss, stdpss ;
pss(t,boot) = pss(t,boot) - meanpss(boot) ;
meanpss(boot) = sum(t, pss(t,boot))/card(t) ;
stdpss(boot) = (sum(t, (pss(t,boot)-meanpss(boot))*(pss(t,boot)-meanpss(boot)))/card(t) )**0.5 ;
display meanpss, stdpss ;

kss(t,boot) = kbar ;

loop(t$(ord(t) LE card(t)),

phi0kss("1",t,boot) = 1;
phi0kss("2",t,boot) = 2*kss(t,boot)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max);
loop(nka$(ord(nka) GE 3),
phi0kss(nka,t,boot) = 2*phi0kss("2",t,boot)*phi0kss(nka-1,t,boot)-phi0kss(nka-2,t,boot) ;
     );

phi0pss("1",t,boot) = 1;
phi0pss("2",t,boot) = 2*pss(t,boot)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phi0pss(npa,t,boot) = 2*phi0pss("2",t,boot)*phi0pss(npa-1,t,boot)-phi0pss(npa-2,t,boot) ;
     );

*lambdass(t,boot) = sum((nka,npa), coeffs(nka,npa)*phi0kss(nka,t,boot)*phi0pss(npa,t,boot) ) ;
css(t,boot)      = (sum((nka,npa)$smol(nka,npa), coeffs(nka,npa)*phi0kss(nka,t,boot)*phi0pss(npa,t,boot) ))**(-1/gamma) ;
yss(t,boot)      = exp(pss(t,boot))*kss(t,boot)**alphak ;
invss(t,boot)    = exp(pss(t,boot))*kss(t,boot)**(alphak) - css(t,boot) ;
kss(t+1,boot)    = (1-delta)*kss(t,boot) + exp(pss(t,boot))*kss(t,boot)**(alphak) - css(t,boot) ;

phik0nss("1",t,boot) = 1;
phik0nss("2",t,boot) = 2*kss(t+1,boot)/(k0max-k0min)+ (k0max+k0min)/(k0min-k0max);
loop(nka$(ord(nka) GE 3),
phik0nss(nka,t,boot) = 2*phik0nss("2",t,boot)*phik0nss(nka-1,t,boot)-phik0nss(nka-2,t,boot) ;
     );

phip0nss("1",t,boot) = 1;
phip0nss("2",t,boot) = 2*pss(t+1,boot)/(pmax-pmin)+ (pmax+pmin)/(pmin-pmax);
loop(npa$(ord(npa) GE 3),
phip0nss(npa,t,boot) = 2*phip0nss("2",t,boot)*phip0nss(npa-1,t,boot)-phip0nss(npa-2,t,boot) ;
     );

*c0nss(t,npb,boot) = (sum((nka,npa)$smol(nka,npa) , coeffs(nka,npa)*phik0nss(nka,t,boot)*phip0n(npa,"1",npb) ))**(-1/gamma) ;

*verif(t,boot)$(ord(t) LT card(t)) = -css(t,boot)**(-gamma)/ beta + sum(npb, prob(npb) * (
*                                                                        c0nss(t,npb,boot)**(-gamma)*( 1 - delta + exp(p0n("1",npb))*alphak*kss(t+1,boot)**(alphak-1) )
*
*                                                                        )

c0nss(t,boot) = (sum((nka,npa)$smol(nka,npa) , coeffs(nka,npa)*phik0nss(nka,t,boot)*phip0nss(npa,t,boot)))**(-1/gamma) ;

verif(t,boot)$(ord(t) LT card(t)) = -css(t,boot)**(-gamma)/ beta + (c0nss(t,boot)**(-gamma)*( 1 - delta + exp(pss(t+1,boot))*alphak*kss(t+1,boot)**(alphak-1) )
                                                           ) ;

) ;


display css, c0nss, kss, invss, yss, pss, verif,smol,smolcount,coeff.l ;


*$libinclude xlexport css results.xlsx datacs!a1:k0n01
*$libinclude xlexport yss results.xlsx datacs!m1:w101

parameters
cs(t)
ys(t)
ps(t)
;

cs(t) = css(t,"1");
ys(t) = yss(t,"1");
ps(t) = pss(t,"1");

display cs, ys

$offtext

*----------------------------------------------------------------------
* Import data (Data generated using Judd et al (2017)'s Matlab program)
*----------------------------------------------------------------------
*Sets
*t          /1*100/
*name       /css, yss, kss, zss/
*;

*Parameters
*series(t,name)
*$call GDXXRW SimGrowthMaliar.xls trace=3 par=series rng=onesector!a1:e101

*$GDXIN SimGrowthMaliar.gdx

Set
*   t    /1*100/
   name /css, yss/
;

Parameters
series(t,name)

$gdxIn SimGrowthMaliar_Method6.gdx
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
err0(np)   ;

*sigmap = 1
*err0("1") = -4.1445 ; err0("2") = -2.8025 ; err0("3") = -1.6365 ; err0("4") = -0.5391 ;
*err0("8") = 4.1445  ; err0("7") = 2.8025  ; err0("6") = 1.6365  ; err0("5") = 0.5391 ;

err0("1")=-2.8570; err0("2")=-1.3556; err0("3")=0;  err0("4")=1.3556; err0("5")= 2.8570;

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
*$ontext
betasv("1")   = 0.9 ;
betasv("2")   = 0.95 ;
betasv("3")   = 0.98 ;
alphaksv("1")   = 0.1 ;
alphaksv("2")   = 0.5 ;
alphaksv("3")   = 0.7 ;
gammasv("1")    = 0.9 ;
gammasv("2")    = 2 ;
gammasv("3")    = 3 ;
deltasv("1")    = 0.01 ;
deltasv("2")    = 0.05 ;
deltasv("3")    = 0.10 ;
rhozsv("1")    = 0.7 ;
rhozsv("2")    = 0.9 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2 ;
*$offtext
$ontext
betasv("1")   = 0.9 ;
betasv("2")   = 0.95 ;
betasv("3")   = 0.999 ;
alphaksv("1")   = 0 ;
alphaksv("2")   = 0.5 ;
alphaksv("3")   = 1 ;
gammasv("1")    = 0.01 ;
gammasv("2")    = 1 ;
gammasv("3")    = 3 ;
deltasv("1")    = 0.01 ;
deltasv("2")    = 0.10 ;
deltasv("3")    = 0.15 ;
rhozsv("1")    = 0.01 ;
rhozsv("2")    = 0.5 ;
rhozsv("3")    = 0.99 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.1 ;
sigmazsv("3")   = 0.2 ;
$offtext
coeffsv("1",nka,nza)    = -100 ;
coeffsv("2",nka,nza)    = 0 ;
coeffsv("3",nka,nza)    = 100 ;
errzsv("1")             = -1 ;
errzsv("3")             = 1;
erreulersv("1")         = -1 ;
erreulersv("3")         = 1;
erreulernesv("1")       = -1;
erreulernesv("3")       = 1;

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

*phiz(nza,t)
*phizne(nza,t,nz)
*phik(nka,t)
*coeffe(nka,nza)
*phikfin(nka)

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

*eqcoeffe
*eqpcoeffsv
*eqphiz(nza,t)
*eqphizne(nza,t,nz)
*eqphik(nka,t)
*eqphikfin(nka)
eqpolicy(t)
eqt(t,nz)
eqt1(t)
eqmeaneulerer
eqeulerer(t)
eqpeulerersv(t)
eqeulerneer(t,nz)
eqpeulerneersv(t,nz)
;

parameter
kmax
kmin
zmax
zmin;

kmax     = kbar*1.2 ;
kmin     = kbar*0.8 ;
zmax = p0("5")*1.2 ;
zmin = -pmax ;

parameters
predict ;
predict = 100;

eqentropie..     entropie =e=     - sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - sum(k, pgammasv(k)*LOG(1.e-5+pgammasv(k)))
                                  - sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - sum(k, pbetasv(k)*LOG(1.e-5+pbetasv(k)))
                                  - predict*sum((k,t), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - predict*sum((k,t), perreulerer(t,k)*log(1.e-5+perreulerer(t,k)))
                                  - predict*sum((k,t,nz), perreulerneer(t,k,nz)*log(1.e-5+perreulerneer(t,k,nz)))
;

eqcs(t)..                         (   - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/beta+delta-1)/alphak)*exp(zne(t,nz))*alphake*ke(t+1)**(alphake-1)
                                                                                                            )
                                                                      )
                                   )$(ord(t) LT card(t))
                                  +(  - cs(t)**(-gammae) + betae * sum(nz, prob(nz)* cne(t,nz)**(-gammae) * (1-deltae
                                                         + ((1/beta+delta-1)/alphak)*exp(zne(t,nz))*alphake*kfin**(alphake-1)
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

eqys(t)..                           ys(t) =E= ((1/beta+delta-1)/alphak)*exp(ze(t))*ke(t)**alphake;                                                                         ;


*eqphiz(nza,t)..                     phiz(nza,t)     =E=                                                1$(ord(nza) eq 1)
*                                                        + (2*ze(t)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax))$(ord(nza) eq 2)
*                                                           + (2*phiz("2",t)*phiz(nza-1,t)-phiz(nza-2,t))$(ord(nza) GE 3) ;
*
*eqphik(nka,t)..                     phik(nka,t)     =E=                                                1$(ord(nka) eq 1)
*                                                        + (2*ke(t)/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax))$(ord(nka) eq 2)
*                                                           + (2*phik("2",t)*phik(nka-1,t)-phik(nka-2,t))$(ord(nka) GE 3) ;
*
*eqphikfin(nka)..                    phikfin(nka)   =E=                                                 1$(ord(nka) eq 1)
*                                                         + (2*kfin/(kmax-kmin)+ (kmax+kmin)/(kmin-kmax))$(ord(nka) eq 2)
*                                                        + (2*phikfin("2")*phikfin(nka-1)-phikfin(nka-2))$(ord(nka) GE 3) ;
*
*eqphizne(nza,t,nz)..
*                                    phizne(nza,t,nz)     =E=                                           1$(ord(nza) eq 1)
*                                                    + (2*zne(t,nz)/(zmax-zmin)+ (zmax+zmin)/(zmin-zmax))$(ord(nza) eq 2)
*                                            + (2*phizne("2",t,nz)*phizne(nza-1,t,nz)-phizne(nza-2,t,nz))$(ord(nza) GE 3) ;

eqpolicy(t)..                        cs(t)**(-gammae)  =E=  mu0
                                                            + mu1*ke(t)
                                                            + mu2*ze(t)
                                                            + 0.5*mu11*ke(t)*ke(t)
                                                            + 0.5*mu22*ze(t)*ze(t)
                                                            + mu12*ke(t)*ze(t)
                                                            + eulerer(t) ;

eqcne(t,nz)..                        (-cne(t,nz)**(-gammae) + mu0
                                                            + mu1*ke(t+1)
                                                            + mu2*zne(t,nz)
                                                            + 0.5*mu11*ke(t+1)*ke(t+1)
                                                            + 0.5*mu22*zne(t,nz)*zne(t,nz)
                                                            + mu12*ke(t+1)*zne(t,nz)
                                                  )$(ord(t) LT card(t))
                                   + (-cne(t,nz)**(-gammae) + mu0
                                                            + mu1*kfin
                                                            + mu2*zne(t,nz)
                                                            + 0.5*mu11*kfin*ke(t+1)
                                                            + 0.5*mu22*zne(t,nz)*zne(t,nz)
                                                            + mu12*kfin*zne(t,nz)
                                                  )$(ord(t) EQ card(t))   + eulererne(t,nz)
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

*eqcoeffe(nka,nza)..     coeffe(nka,nza)=E= sum(k, pcoeffsv(k,nka,nza)*coeffsv(k,nka,nza)) ;
*eqpcoeffsv(nka,nza)..   1         =E= sum(k, pcoeffsv(k,nka,nza)) ;

eqmeaneulerer..         sum(t, eulerer(t))/card(t) =E= 0;

eqeulerer(t)..          eulerer(t) =E= sum(k, perreulerer(t,k)*erreulersv(k)) ;
eqpeulerersv(t)..       1          =E= sum(k, perreulerer(t,k)) ;

eqeulerneer(t,nz)..          eulererne(t,nz) =E= sum(k, perreulerneer(t,k,nz)*erreulernesv(k)) ;
eqpeulerneersv(t,nz)..       1          =E= sum(k, perreulerneer(t,k,nz)) ;

eqt(t,nz)..                 cne(t,nz)-eulererne(t,nz) =G= 0.001 ;
eqt1(t)..                   cs(t) - eulerer(t) =G= 0.001 ;

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
*eqcoeffe
*eqpcoeffsv
*eqphik
*eqphizne
*eqphiz
*eqphikfin
eqcne
eqpolicy
eqt
eqt1
eqmeaneulerer
*eqeulerer
*eqpeulerersv
*eqeulerneer
*eqpeulerneersv
/ ;


model estimationt /
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
*eqcoeffe
*eqpcoeffsv
*eqphik
*eqphizne
*eqphikfin
*eqphiz
eqcne
eqpolicy
eqt
eqt1
eqmeaneulerer
eqeulerer
eqpeulerersv
eqeulerneer
eqpeulerneersv
/ ;


*initiate chebyshev coefficient for estimation
gammae.l  = gamma_true ;
deltae.l  = delta_true ;
alphake.l = alphak_true ;
betae.l   = beta_true ;
rhoze.l    = 0.85;
sigmaze.l  = 0.04 ;
*gammae.l  = 0.9;
*deltae.l  = 0.05 ;
*alphake.l = 0.5 ;
*betae.l   = 0.95 ;
*rhoze.l    = 0.7;
*sigmaze.l  = 1 ;
kbar   = 1;
cbar   = (1/beta+delta-1)/alphak - delta ;

k0max = kbar*1.2 ;
k0min = kbar*0.8 ;
k0(nk)  = (k0max+k0min)/2 + (k0max-k0min)/2*  cos(  pi * (1+ (0.5 - ord(nk))/card(nk) ) ) ;

ke.l(t)           = kbar  ;
kfin.l            = kbar ;
ke.lo(t)          = kbar*0.5 ;
ke.up(t)          = kbar*1.5 ;
ze.l(t)           = 0;
ze.up(t)          = zmax ;
ze.lo(t)          = zmin ;
cne.l(t,npa)      = cbar ;
cne.lo(t,npa)     = 0.1*cbar ;
cne.up(t,npa)     = 5*cbar ;

mu0.l = cs("51");
mu11.l = 0.05;
mu22.l = 0.05;
mu12.l = 0.05;

erreulersv("1")         = -1;
erreulersv("3")         = 1;
erreulernesv("1")       = -1;
erreulernesv("3")       = 1;
solve estimation using nlp maximising entropie;
solve estimationt using nlp maximising entropie;

erreulersv("1")         = -0.01;
erreulersv("3")         = 0.01;
erreulernesv("1")       = -0.01;
erreulernesv("3")       = 0.01;
*solve estimation using nlp maximising entropie;
*solve estimationt using nlp maximising entropie;

*execute_unload '1sector_complete.gdx';



