
*----------------
* data simulation
*----------------
Set ttot    / 1*150/
t(ttot)     / 51*150 / ;
Set boot  /1*100 /;


Parameters
err0(nzb)
errr0(nrb)
z0e(nzb)
probz(nzb)
r0e(nrb)
probr(nrb)
;
*gaussion quadrature std=1
err0("1")= -2.85697001387281; err0("2")= -1.35562617997427; err0("3")= 0; err0("4")= 1.35562617997427; err0("5")= 2.85697001387281;
probz("1")=0.0112574113277207;probz("2")=0.222075922005613; probz("3")=0.533333333333333;probz("4")=0.222075922005613 ;probz("5")=0.0112574113277207 ;

errr0("1")= -2.85697001387281; errr0("2")= -1.35562617997427; errr0("3")= 0; errr0("4")= 1.35562617997427; errr0("5")= 2.85697001387281;
probr("1")=0.0112574113277207; probr("2")=0.222075922005613; probr("3")=0.533333333333333;probr("4")=0.222075922005613 ;probr("5")=0.0112574113277207 ;


Parameters
zss(ttot,boot)
rss(ttot,boot)
rshockss(ttot,boot)
kss(ttot,boot)
iss(ttot,boot)
bss(ttot,boot)
css(ttot,boot)
yss(ttot,boot)
verifkss(ttot,boot)
verifbss(ttot,boot)

;

*rbar = 1.0125;
zss(ttot,boot) = 0 ;
rshockss(ttot,boot) = 0 ;
*execseed = 1e8*(frac(jnow));

loop(boot,
loop(ttot,
zss(ttot+1,boot) = rhoz*zss(ttot,boot) + sigmaz*normal(0,1) ;
rshockss(ttot,boot) = sigmar*normal(0,1) ;
) ;
*execseed = 1 + gmillisec(jnow);
) ;

*display zss
parameters
meanzss(boot)
stdzss(boot) ;
meanzss(boot) = sum(ttot, zss(ttot,boot))/card(ttot) ;
stdzss(boot) = (sum(ttot, (zss(ttot,boot)-meanzss(boot))*(zss(ttot,boot)-meanzss(boot)))/card(ttot) )**0.5 ;
display meanzss, stdzss ;
zss(ttot,boot) = zss(ttot,boot) - meanzss(boot) ;
meanzss(boot) = sum(ttot, zss(ttot,boot))/card(ttot) ;
stdzss(boot) = (sum(ttot, (zss(ttot,boot)-meanzss(boot))*(zss(ttot,boot)-meanzss(boot)))/card(ttot) )**0.5 ;
display meanzss, stdzss ;
display zss ;

parameters
meanrshockss(boot)
stdrshockss(boot) ;
meanrshockss(boot) = sum(ttot, rshockss(ttot,boot))/card(t) ;
stdrshockss(boot) = (sum(ttot, (rshockss(ttot,boot)-meanrshockss(boot))*(rshockss(ttot,boot)-meanrshockss(boot)))/card(ttot) )**0.5 ;
display meanrshockss, stdrshockss ;
rshockss(ttot,boot) = rshockss(ttot,boot) - meanrshockss(boot) ;
meanrshockss(boot) = sum(ttot, rshockss(ttot,boot))/card(ttot) ;
stdrshockss(boot) = (sum(ttot, (rshockss(ttot,boot)-meanrshockss(boot))*(rshockss(ttot,boot)-meanrshockss(boot)))/card(ttot) )**0.5 ;
display meanrshockss, stdrshockss ;


display sigmar ;

rss(ttot,boot) = alphar1*exp(rshockss(ttot,boot)) ;
kss(ttot,boot) = kbar ;
bss(ttot,boot) = bbar ;


parameters
z0sim
r0sim
k0sim
b0sim
;

variables
i0sim
c0sim
k1sim
b1sim
objsim
;

equations
eqk1sim
eqb1sim
eqobjsim
;


eqk1sim..       k1sim =E=  (1-delta)*k0sim+i0sim ;
eqb1sim..       b1sim =E= ( -(1/beta+delta-1)/alphak * exp(z0sim) * k0sim**alphak
                                                      +c0sim+i0sim   + a/2*(k1sim-k0sim)*(k1sim-k0sim)/k0sim
                                                  )
                                                 +b0sim  *(r0sim + alphar2*b0sim*b0sim )  ;

eqobjsim..         objsim =E=
                                         c0sim**(1-gamma) / (1-gamma)
                                        + beta*sum( (nzb,nrb),weight_nodes(nzb)*weight_nodes(nrb)*
                                                              (  coeffini("1")
                                  + coeffini("2")*k1sim  + coeffini("3")*b1sim + coeffini("4")*exp(z1("1",nzb))+ coeffini("5")*r1("1",nrb)
                                  + coeffini("6")*k1sim *k1sim + coeffini("7")*k1sim *b1sim + coeffini("8")*k1sim *exp(z1("1",nzb))+ coeffini("9")*k1sim *r1("1",nrb)
                                  + coeffini("10")*b1sim *b1sim + coeffini("11")*b1sim *exp(z1("1",nzb))+ coeffini("12")*b1sim *r1("1",nrb)
                                  + coeffini("13")*exp(z1("1",nzb))*exp(z1("1",nzb))+ coeffini("14")*exp(z1("1",nzb))*r1("1",nrb)
                                  + coeffini("15")*r1("1",nrb)*r1("1",nrb)
                                  )
                                                    )
                                       ;



model simulation /
eqk1sim
eqb1sim
eqobjsim / ;


k0sim = kbar ;
b0sim = bbar ;
z0sim = zss("1","1");
r0sim = rss("1","1");
z0sim = 0 ;
r0sim = alphar1 ;


k1sim.l = kbar ;
b1sim.l = bbar ;
c0sim.l = cbar ;
i0sim.l = ibar ;
c0sim.lo = cbar*0.1 ;

solve simulation using nlp maximising objsim ;

*zss("1",boot) = -0.1 ;
*rss("1",boot) = 1.04 ;

loop(boot,
k0sim = kbar ;
b0sim = bbar ;
loop(ttot,
z0sim = zss(ttot,boot);
r0sim = rss(ttot,boot);
solve simulation using nlp maximising objsim ;

css(ttot,boot)=c0sim.l ;
yss(ttot,boot)=((1/beta+delta-1)/alphak)*exp(z0sim)*k0sim**alphak;
iss(ttot,boot) = k1sim.l-(1-delta)*k0sim ;
kss(ttot,boot) = k0sim ;
bss(ttot,boot) = b0sim ;
k0sim = k1sim.l ;
b0sim = b1sim.l ;

);
);

display css, yss, iss, kss, bss ;



Parameters
boundmax
boundmin
titi(boot) ;
boundmax = 1.1  ;
boundmin  = 0.9 ;
titi(boot) = 1$ ( (smax(t, kss(t,boot)) GE kmax*boundmax)  or (smax(t, bss(t,boot)) GE bmax*boundmax)  or (smin(t, kss(t,boot)) LE kmin*boundmin)  or (smin(t, bss(t,boot)) LE bmin*boundmin)) ;
*titi(boot) = 0 ;


$exit

parameters
cs(t)
ys(t)
zs(t)
is(t)
;

cs(t) = css(t,"1");
ys(t) = yss(t,"1");
zs(t) = zss(t,"1");
is(t) = iss(t,"1");


*----------------------
* estimation GME
*----------------------
parameters
kemax
kemin
zemax
zemin
bemax
bemin
remax
remin
predict
;


sets k number of support values / 1*3 / ;

Parameters
betasv(k)
alphaksv(k)
gammasv(k)
deltasv(k)
rhozsv(k)
sigmazsv(k)
sigmarsv(k)
alphar1sv(k)
alphar2sv(k)
asv(k)

errzsv(k)
errrsv(k)
errisv(k)
erreulerksv(k)
erreulerbsv(k)
;

Variables
entropie
betae
alphake
gammae
deltae
rhoze
sigmaze
sigmare
alphar1e
alphar2e
ae

ke(t)
kfin
be(t)
bfin
ze(t)
re(t)


zne(t,nzb)
rne(t,nrb)
errz(t)
errr(t)
erri(t)
eulererk(t)
eulererb(t)
eulererkb(t)
eulererbb(t)

coeffe(ncoeff)

vf_derivk(t)
vf_derivb(t)
vfne_derivk(t,nzb,nrb)
vfne_derivb(t,nzb,nrb)

;

Positive variables
pbetasv(k)
palphaksv(k)
pgammasv(k)
pdeltasv(k)
prhozsv(k)
psigmazsv(k)
psigmarsv(k)
palphar1sv(k)
palphar2sv(k)
pasv(k)

perrzsv(t,k)
perrisv(t,k)
perrrsv(t,k)
perreulererk(t,k)
perreulererb(t,k)
perreulererkb(t,k)
perreulererbb(t,k)

;

equations
eqbetae
eqalphake
eqgammae
eqdeltae
eqrhoze
eqsigmaze
eqsigmare
eqalphare1
eqalphare2
eqpalphaksv
eqpgammasv
eqpdeltasv
eqpbetasv
eqprhozsv
eqpsigmazsv
eqpsigmarsv
eqpalphar1sv
eqpalphar2sv
eqae
eqpasv

eqke(t)
eqbe(t)
eqzne(t,nzb)
eqze(t)
eqrne(t,nrb)
eqre(t)

eqcs(t)
eqys(t)
eqis(t)
eqentropie

eqmeanze
eqstdze
eqmeanre
eqstdre

eqerrz(t)
eqperrzsv(t)
eqerrr(t)
eqperrrsv(t)
eqerri(t)
eqperrisv(t)

eqvfne_derivk(t,nzb,nrb)
eqvfne_derivb(t,nzb,nrb)
eqvf_FOCk(t)
eqvf_FOCb(t)

eqeulererk(t)
eqpeulererksv(t)
eqeulererb(t)
eqpeulererbsv(t)
eqeulererkb(t)
eqpeulererksvb(t)
eqeulererbb(t)
eqpeulererbsvb(t)

eqvf_derivk(t)
eqvf_derivb(t)
eqvf_derivkb(t)
eqvf_derivbb(t)

eqeulererkm
eqeulererbm
eqeulererkbm
eqeulererbbm
;


eqentropie..     entropie =E=     - predict*sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - predict*sum(k, pgammasv(k)*LOG(1.e-5+pgammasv(k)))
                                  - predict*sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - predict*sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - predict*sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - predict*sum(k, psigmarsv(k)*LOG(1.e-5+psigmarsv(k)))
                                  - predict*sum(k, palphar1sv(k)*LOG(1.e-5+palphar1sv(k)))
                                  - predict*sum(k, palphar2sv(k)*LOG(1.e-5+palphar2sv(k)))
                                  - predict*sum(k, pbetasv(k)*LOG(1.e-5+pbetasv(k)))
                                  - predict*sum(k, pasv(k)*LOG(1.e-5+pasv(k)))
                                  - predict*sum((k,t), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - predict*sum((k,t), perrrsv(t,k)*log(1.e-5+perrrsv(t,k)))

                                  - 0*sum((k,t), perrisv(t,k)*log(1.e-5+perrisv(t,k)))
                                  - 100*sum((k,t), perreulererk(t,k)*log(1.e-5+perreulererk(t,k)))
                                  - 100*sum((k,t), perreulererb(t,k)*log(1.e-5+perreulererb(t,k)))
                                  - 100*sum((k,t), perreulererkb(t,k)*log(1.e-5+perreulererkb(t,k)))
                                  - 100*sum((k,t), perreulererbb(t,k)*log(1.e-5+perreulererbb(t,k)))
                                  ;


eqvf_derivk(t).. vf_derivk(t)  =E= coeffe("2")
                  + coeffe("6")*2*ke(t)+ coeffe("7")*be(t) + coeffe("8")*exp(ze(t)) +coeffe("9")*re(t)
*                  + coeffe("16")*3*ke(t)*ke(t)+ coeffe("17")*2*ke(t)*be(t)+ coeffe("18")*2*ke(t)*exp(ze(t))+ coeffe("19")*2*ke(t)*re(t)
*                  + coeffe("20")*be(t)*be(t)+ coeffe("21")*be(t)*exp(ze(t))+ coeffe("22")*be(t)*re(t)
*                  + coeffe("23")*exp(ze(t))*exp(ze(t))+ coeffe("24")*exp(ze(t))*re(t)
*                  + coeffe("25")*re(t)*re(t)
;
eqvf_derivkb(t).. vf_derivk(t) =E= betae*sum((nzb,nrb),probz(nzb)*probr(nrb)*
                                        ( vfne_derivk(t,nzb,nrb)*(1-deltae)
                                          - vfne_derivb(t,nzb,nrb)*
                                             ( ((1/betae+deltae-1)/alphake)*exp(ze(t))*alphake*ke(t)**(alphake-1)
                                               + ae/2/ke(t)/ke(t)*(is(t) - deltae*ke(t))*(is(t) + deltae*ke(t))
                                             )
                                        ) ) + eulererkb(t);



eqvf_derivb(t)..  vf_derivb(t) =E= coeffe("3")
                  + coeffe("7")*ke(t) + coeffe("10")*2*be(t) + coeffe("11")*exp(ze(t)) + coeffe("12")*re(t)
*                  + coeffe("17")*ke(t)*ke(t)
*                  + coeffe("20")*ke(t)*2*be(t)+ coeffe("21")*ke(t)*exp(ze(t))+ coeffe("22")*ke(t)*re(t)
*                  + coeffe("26")*3*be(t)*be(t)+ coeffe("27")*2*be(t)*exp(ze(t))+ coeffe("28")*2*be(t)*re(t)
*                  + coeffe("29")*exp(ze(t))*exp(ze(t))+ coeffe("30")*exp(ze(t))*re(t)
*                  + coeffe("31")*re(t)*re(t)
;

eqvf_derivbb(t).. vf_derivb(t) =E= betae*sum((nzb,nrb),probz(nzb)*probr(nrb)*
                                          vfne_derivb(t,nzb,nrb)*( 2*alphar2e*be(t)*be(t)
                                                                  +re(t)+alphar2e*be(t)*be(t)
                                                                 )

                                        ) + eulererbb(t);




eqvf_FOCk(t)..       cs(t)**(-gammae)  + betae*sum((nzb,nrb),probz(nzb)*probr(nrb)*(vfne_derivb(t,nzb,nrb)) ) =E= eulererk(t) ;

eqvf_FOCb(t)..       sum((nzb,nrb),probz(nzb)*probr(nrb)*(  vfne_derivb(t,nzb,nrb)*(1+ae/ke(t)*(is(t)-deltae*ke(t)))  +
                                                        vfne_derivk(t,nzb,nrb)
                                                      )
                       ) =E= eulererb(t)  ;



eqvfne_derivk(t,nzb,nrb)..  vfne_derivk(t,nzb,nrb)   =E= (coeffe("2")
                                                              + coeffe("6")*2*ke(t+1)+ coeffe("7")*be(t+1)+ coeffe("8")*exp(zne(t,nzb))+ coeffe("9")*rne(t,nrb)
*                                                              + coeffe("16")*3*ke(t+1)*ke(t+1)+ coeffe("17")*2*ke(t+1)*be(t+1)+ coeffe("18")*2*ke(t+1)*exp(zne(t,nzb))+ coeffe("19")*2*ke(t+1)*rne(t,nrb)
*                                                              + coeffe("20")*be(t+1)*be(t+1)+ coeffe("21")*be(t+1)*exp(zne(t,nzb))+ coeffe("22")*be(t+1)*rne(t,nrb)
*                                                              + coeffe("23")*exp(zne(t,nzb))*exp(zne(t,nzb))+ coeffe("24")*exp(zne(t,nzb))*rne(t,nrb)
*                                                              + coeffe("25")*rne(t,nrb)*rne(t,nrb)
                                                              )$(ord(t) LT card(t))
                                                              +
                                                               (coeffe("2")
                                                              + coeffe("6")*2*kfin + coeffe("7")*bfin + coeffe("8")*exp(zne(t,nzb))+ coeffe("9")*rne(t,nrb)
*                                                              + coeffe("16")*3*kfin*kfin+ coeffe("17")*2*kfin*bfin+ coeffe("18")*2*kfin*exp(zne(t,nzb))+ coeffe("19")*2*kfin*rne(t,nrb)
*                                                              + coeffe("20")*bfin*bfin+ coeffe("21")*bfin*exp(zne(t,nzb))+ coeffe("22")*bfin*rne(t,nrb)
*                                                              + coeffe("23")*exp(zne(t,nzb))*exp(zne(t,nzb))+ coeffe("24")*exp(zne(t,nzb))*rne(t,nrb)
*                                                              + coeffe("25")*rne(t,nrb)*rne(t,nrb)
                                                              )$(ord(t) EQ card(t))
                                                             ;

eqvfne_derivb(t,nzb,nrb)..  vfne_derivb(t,nzb,nrb)  =E= (coeffe("3")
                                                               + coeffe("7")*ke(t+1)+ coeffe("10")*2*be(t+1)+ coeffe("11")*exp(zne(t,nzb))+ coeffe("12")*rne(t,nrb)
*                                                               + coeffe("17")*ke(t+1)*ke(t+1)
*                                                               + coeffe("20")*ke(t+1)*2*be(t+1)+ coeffe("21")*ke(t+1)*exp(zne(t,nzb))+ coeffe("22")*ke(t+1)*rne(t,nrb)
*                                                               + coeffe("26")*3*be(t+1)*be(t+1)+ coeffe("27")*2*be(t+1)*exp(zne(t,nzb))+ coeffe("28")*2*be(t+1)*rne(t,nrb)
*                                                               + coeffe("29")*exp(zne(t,nzb))*exp(zne(t,nzb))+ coeffe("30")*exp(zne(t,nzb))*rne(t,nrb)
*                                                               + coeffe("31")*rne(t,nrb)*rne(t,nrb)
                                                               )$(ord(t) LT card(t))
                                                            +
                                                                (coeffe("3")
                                                               + coeffe("7")*kfin + coeffe("10")*2*bfin + coeffe("11")*exp(zne(t,nzb))+ coeffe("12")*rne(t,nrb)
*                                                               + coeffe("17")*kfin*kfin
*                                                               + coeffe("20")*kfin*2*bfin+ coeffe("21")*kfin*exp(zne(t,nzb))+ coeffe("22")*kfin*rne(t,nrb)
*                                                               + coeffe("26")*3*bfin*bfin+ coeffe("27")*2*bfin*exp(zne(t,nzb))+ coeffe("28")*2*bfin*rne(t,nrb)
*                                                               + coeffe("29")*exp(zne(t,nzb))*exp(zne(t,nzb))+ coeffe("30")*exp(zne(t,nzb))*rne(t,nrb)
*                                                               + coeffe("31")*rne(t,nrb)*rne(t,nrb)
                                                                )$(ord(t) EQ card(t))
                                                             ;

* Modeling next period TFP expectation
eqzne(t,nzb)..        zne(t,nzb) =E=  rhoze * ze(t) + sigmaze*err0(nzb);

* Modeling TFP evolution process
eqze(t)..            ze(t)   =E=  rhoze * ze(t) + errz(t);

* Modeling next period r expectation
eqrne(t,nrb)..        rne(t,nrb) =E=  alphar1e*exp(sigmare*errr0(nrb));

* Modeling r evolution process
eqre(t)..            re(t)   =E=  alphar1e*exp(errr(t));

eqke(t)..            (is(t) + ae/2*(is(t)-deltae*ke(t))*(is(t)-deltae*ke(t))/ke(t)
                           + cs(t) - ys(t) -be(t+1)+(re(t)+alphar2e*be(t)*be(t))*be(t) )$(ord(t) LT card(t))
                                                               +
                     (is(t) + ae/2*(is(t)-deltae*ke(t))*(is(t)-deltae*ke(t))/ke(t)
                              + cs(t) - ys(t) -bfin+(re(t)+alphar2e*be(t)*be(t))*be(t)   )$(ord(t) EQ card(t))  =E= 0
                     ;

eqis(t)..            is(t)=E= (ke(t+1) - (1-deltae)*ke(t))$(ord(t) LT card(t))
                              +  (kfin - (1-deltae)*ke(t))$(ord(t) EQ card(t)) + erri(t);

eqmeanze..            sum(t, errz(t))/card(t) =E= 0;

eqstdze..             (sum(t, errz(t)*errz(t))/card(t))**0.5 =E= sigmaze  ;

eqmeanre..            sum(t, errr(t))/card(t) =E= 0;

eqstdre..             (sum(t, errr(t)*errr(t))/card(t))**0.5 =E= sigmare  ;

eqys(t)..             ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake ;

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
eqsigmare..       sigmare   =E= sum(k, psigmarsv(k)*sigmarsv(k)) ;
eqpsigmarsv..     1         =E= sum(k, psigmarsv(k)) ;
eqalphare1..       alphar1e   =E= sum(k, palphar1sv(k)*alphar1sv(k)) ;
eqpalphar1sv..     1         =E= sum(k, palphar1sv(k)) ;
eqalphare2..       alphar2e   =E= sum(k, palphar2sv(k)*alphar2sv(k)) ;
eqpalphar2sv..     1         =E= sum(k, palphar2sv(k)) ;
eqae..             ae    =E= sum(k, pasv(k)*asv(k)) ;
eqpasv..           1         =E= sum(k, pasv(k)) ;

eqerrz(t)..           errz(t)    =E= sum(k, perrzsv(t,k)*errzsv(k)) ;
eqperrzsv(t)..        1          =E= sum(k, perrzsv(t,k)) ;

eqerrr(t)..           errr(t)    =E= sum(k, perrrsv(t,k)*errrsv(k)) ;
eqperrrsv(t)..        1          =E= sum(k, perrrsv(t,k)) ;

eqerri(t)..           erri(t)    =E= sum(k, perrisv(t,k)*errisv(k)) ;
eqperrisv(t)..        1          =E= sum(k, perrisv(t,k)) ;

eqeulererk(t)..          eulererk(t) =E= sum(k, perreulererk(t,k)*erreulerksv(k)) ;
eqpeulererksv(t)..       1          =E= sum(k, perreulererk(t,k)) ;
eqeulererkm..            sum(t,eulererk(t)) =E= 0 ;

eqeulererb(t)..          eulererb(t) =E= sum(k, perreulererb(t,k)*erreulerbsv(k)) ;
eqpeulererbsv(t)..       1          =E= sum(k, perreulererb(t,k)) ;
eqeulererbm..            sum(t,eulererb(t)) =E= 0 ;

eqeulererkb(t)..          eulererkb(t) =E= sum(k, perreulererkb(t,k)*erreulerksv(k)) ;
eqpeulererksvb(t)..       1          =E= sum(k, perreulererkb(t,k)) ;
eqeulererkbm..            sum(t,eulererkb(t)) =E= 0 ;

eqeulererbb(t)..          eulererbb(t) =E= sum(k, perreulererbb(t,k)*erreulerbsv(k)) ;
eqpeulererbsvb(t)..       1          =E= sum(k, perreulererbb(t,k)) ;
eqeulererbbm..            sum(t,eulererbb(t)) =E= 0 ;


model estimation /
eqalphake
eqpalphaksv
eqgammae
eqpgammasv
eqdeltae
eqpdeltasv
eqrhoze
eqprhozsv
eqsigmaze
eqpsigmazsv
eqsigmare
eqpsigmarsv
eqbetae
eqpbetasv
eqalphare1
eqpalphar1sv
eqalphare2
eqpalphar2sv

eqke
eqre
eqze
eqis
eqys
eqzne
eqrne
eqentropie
eqmeanze
eqstdze
eqmeanre
eqstdre
eqerrz
eqperrzsv
eqerrr
eqperrrsv
eqerri
eqperrisv
eqvfne_derivk
eqvfne_derivb
eqvf_FOCk
eqvf_FOCb
eqeulererk
eqpeulererksv
eqeulererb
eqpeulererbsv
eqeulererkb
eqpeulererksvb
eqeulererbb
eqpeulererbsvb

eqeulererkm
eqeulererbm
eqeulererkbm
eqeulererbbm

eqvf_derivk
eqvf_derivb
eqvf_derivkb
eqvf_derivbb
eqae
eqpasv

/
;

*initiate chebyshev coefficient for estimation
*$ontext
betasv("1")   = 0.97 ;
betasv("2")   = 0.98 ;
betasv("3")   = 0.99 ;
alphaksv("1")   = 0.1 ;
alphaksv("2")   = 0.25 ;
alphaksv("3")   = 0.5 ;
gammasv("1")    = 0.1 ;
gammasv("2")    = 1.5 ;
gammasv("3")    = 3 ;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.05 ;
deltasv("3")    = 0.10 ;
rhozsv("1")    = 0 ;
rhozsv("2")    = 0 ;
rhozsv("3")    = 0 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.075 ;
sigmazsv("3")   = 0.15 ;
sigmarsv("1")   = 0.001 ;
sigmarsv("2")   = 0.005 ;
sigmarsv("3")   = 0.010;
alphar1sv("1")   = 1.001 ;
alphar1sv("2")   = 1.025 ;
alphar1sv("3")   = 1.05;
alphar2sv("1")   = 0.001 ;
alphar2sv("2")   = 0.0075;
alphar2sv("3")   = 0.015;
asv("1")   = 0 ;
asv("2")   = 0.5 ;
asv("3")   = 1 ;

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
sigmazsv("3")   = 0.2;
$offtext

z0e(nzb) = err0(nzb) * sigmazsv("2");
r0e(nrb) = errr0(nrb)* sigmarsv("2");

Parameters
kbar
cbar
;

zemax = z0e("5")*1.5 ;
zemin = -zemax ;

errzsv("1") = zemin ;
errzsv("3") = zemax ;

display errzsv;

*remax = rbar*exp(rshockmax)*1.5 ;
*remin = rbar*exp(-rshockmax)*1.5 ;

errrsv("1") = -r0e("5")*1.5 ;
errrsv("3") = +r0e("5")*1.5 ;

remin = 1 ;
remax = 1.05 ;

display remax, remin;


parameters
gammaeboot(boot)
deltaeboot(boot)
alphakeboot(boot)
betaeboot(boot)
rhozeboot(boot)
sigmazeboot(boot)
sigmareboot(boot)
alphar1eboot(boot)
alphar2eboot(boot)
modelboot(boot)
aeboot(boot)
;

parameter
vfderivkss
vfderivbss
ref;



loop(boot$(titi(boot) EQ 0) ,
cs(t) = css(t,boot);
ys(t) = yss(t,boot);
is(t) = iss(t,boot);

* we reinitialize endogenous variables
gammae.l  = gammasv("2") ;
deltae.l  = deltasv("2") ;
alphake.l = alphaksv("2") ;
betae.l   = betasv("2") ;
rhoze.l    = rhozsv("2") ;
sigmaze.l  = sigmazsv("2") ;
sigmare.l  = sigmarsv("2") ;
alphar1e.l  = alphar1sv("2") ;
alphar2e.l  = alphar2sv("2") ;
ae.l        = asv("2") ;

gammae.lo  = gammasv("1") ;
deltae.lo  = deltasv("1") ;
alphake.lo = alphaksv("1") ;
betae.lo   = betasv("1") ;
rhoze.lo    = rhozsv("1") ;
alphar1e.lo  = alphar1sv("1") ;
alphar2e.lo  = alphar2sv("1") ;
ae.lo        = asv("1") ;

gammae.up  = gammasv("3") ;
deltae.up  = deltasv("3") ;
alphake.up = alphaksv("3") ;
betae.up   = betasv("3") ;
rhoze.up    = rhozsv("3") ;
sigmaze.up  = sigmazsv("3") ;
sigmare.up  = sigmarsv("3") ;
alphar1e.up  = alphar1sv("3") ;
alphar2e.up  = alphar2sv("3") ;
ae.up        = asv("3") ;

rhoze.fx = 0;
*betae.fx=beta_true;
*gammae.fx=gamma_true;
*alphake.fx=alphak_true;
*deltae.fx=delta_true;
*sigmaze.fx=sigmaz_true;
*alphar1e.fx  = 1.012 ;
*alphar2e.fx = alphar2 ;
*ae.fx       = a ;

kbar = 1;
bbar = 1;
*cbar   = (1/betae.l+deltae.l-1)/alphake.l - deltae.l ;
kemax     = kbar*2;
kemin     = kbar*0.5;
bemax     = bbar*2;
bemin     = bbar*0.5;
ke.l(t)           = kbar ;
kfin.l            = kbar ;
be.l(t)           = bbar ;
bfin.l            = bbar ;
ke.lo(t)          = kemin ;
ke.up(t)          = kemax ;
kfin.lo           = kemin ;
kfin.up           = kemax ;
be.lo(t)          = bemin ;
be.up(t)          = bemax ;
bfin.lo          = bemin ;
bfin.up          = bemax ;



ze.l(t)           = 0 ;
ze.up(t)          = zemax ;
ze.lo(t)          = zemin ;
zne.l(t,nzb)       = 0 ;
zne.up(t,nzb)       = zemax ;
zne.lo(t,nzb)       = zemin ;
re.l(t)           = alphar1sv("2") ;
re.up(t)          = remax ;
re.lo(t)          = remin ;
rne.l(t,nrb)       = alphar1sv("2") ;
rne.up(t,nrb)       = remax ;
rne.lo(t,nrb)       = remin ;



errz.l(t)         = 0 ;
errr.l(t)         = 0 ;
erri.l(t)         = 0 ;
eulererk.l(t)      = 0 ;
eulererb.l(t)      = 0 ;
eulererkb.l(t)      = 0 ;
eulererbb.l(t)      = 0 ;

pbetasv.l(k)            = 1/card(k) ;
palphaksv.l(k)          = 1/card(k) ;
pgammasv.l(k)           = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
psigmazsv.l(k)          = 1/card(k) ;
psigmarsv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perrrsv.l(t,k)          = 1/card(k) ;
perreulererk.l(t,k)      = 1/card(k) ;
perreulererb.l(t,k)      = 1/card(k) ;
perreulererkb.l(t,k)      = 1/card(k) ;
perreulererbb.l(t,k)      = 1/card(k) ;

palphar1sv.l(k)          =1/card(k) ;
palphar2sv.l(k)          =1/card(k) ;
pasv.l(k)                = 1/card(k) ;


vf_derivk.l(t) =  cs(t)**(-gammae.l)*
            (1-deltae.l +  (1/betae.l+deltae.l-1)/alphake.l*alphake.l*kbar**(alphake.l-1));
vfne_derivk.l(t,nzb,nrb) =  cs(t)**(-gammae.l)*
            (1-deltae.l + (1/betae.l+deltae.l-1)/alphake.l*alphake.l*kbar**(alphake.l-1));

vf_derivb.l(t) =  -cs(t)**(-gammae.l) * (alphar1e.l+3*alphar2e.l*bbar*bbar) ;
vfne_derivb.l(t,nzb,nrb) =  -cs(t)**(-gammae.l) * (alphar1e.l+3*alphar2e.l*bbar*bbar)  ;

vfderivkss = ( sum(t,cs(t))/card(t))**(-gammae.l)*
                       (1-deltae.l +  (1/betae.l+deltae.l-1)/alphake.l*alphake.l*kbar**(alphake.l-1));
vfderivbss = -(sum(t,cs(t))/card(t))**(-gammae.l) * (alphar1e.l+3*alphar2e.l*bbar*bbar)  ;

coeffe.l(ncoeff) = 0 ;
coeffe.l("2")= vfderivkss ;
coeffe.l("3")=  vfderivbss  ;
predict = 0.1;
vf_derivk.lo(t) = 0 ;
vf_derivb.up(t) = 0 ;


errisv("1")         = -0.0*(sum(t,is(t))/card(t));
errisv("3")         = 0.0*(sum(t,is(t))/card(t));
erreulerksv("1")         = -1;
erreulerksv("3")         = +1;
erreulerbsv("1")         = -1;
erreulerbsv("3")         = +1;



*estimation.optfile = 1 ;

solve estimation using nlp maximising entropie;

erreulerksv("1")         = -0.1;
erreulerksv("3")         = +0.1;
erreulerbsv("1")         = -0.1;
erreulerbsv("3")         = +0.1;



*estimation.optfile = 1 ;

solve estimation using nlp maximising entropie;



gammaeboot(boot) = gammae.l ;
deltaeboot(boot) = deltae.l ;
alphakeboot(boot) = alphake.l ;
betaeboot(boot) = betae.l   ;
rhozeboot(boot) = rhoze.l   ;
sigmazeboot(boot) = sigmaze.l  ;
sigmareboot(boot) = sigmare.l  ;
alphar1eboot(boot) = alphar1e.l  ;
alphar2eboot(boot) = alphar2e.l  ;
aeboot(boot)       = ae.l ;
modelboot(boot) = estimation.modelstat ;


) ;





parameters
alphakem
gammaem
deltaem
betaem
rhozem
sigmazem
sigmarem
alphakestd
gammaestd
deltaestd
betaestd
rhozestd
sigmazestd
sigmarestd
alphar1m
alphar1std
alphar2m
alphar2std
;

alphakem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), alphakeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), 1) ;
alphakestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (alphakeboot(boot)-alphakem)*(alphakeboot(boot)-alphakem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
gammaem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), gammaeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
gammaestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (gammaeboot(boot)-gammaem)*(gammaeboot(boot)-gammaem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
deltaem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), deltaeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
deltaestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (deltaeboot(boot)-deltaem)*(deltaeboot(boot)-deltaem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
betaem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), betaeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
betaestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (betaeboot(boot)-betaem)*(betaeboot(boot)-betaem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
rhozem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), rhozeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
rhozestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (rhozeboot(boot)-rhozem)*(rhozeboot(boot)-rhozem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
sigmazem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), sigmazeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
sigmazestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (sigmazeboot(boot)-sigmazem)*(sigmazeboot(boot)-sigmazem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
sigmarem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), sigmareboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
sigmarestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (sigmareboot(boot)-sigmarem)*(sigmareboot(boot)-sigmarem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
alphar1m = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), alphar1eboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
alphar1std = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (alphar1eboot(boot)-alphar1m)*(alphar1eboot(boot)-alphar1m) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;
alphar2m = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), alphar2eboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
alphar2std = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (alphar2eboot(boot)-alphar2m)*(alphar2eboot(boot)-alphar2m) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;

parameters
alphakemse
gammaemse
betaemse
deltaemse
rhozemse
sigmazemse
sigmaremse
alphar1mse
alphar2mse
;

alphakemse = (power((alphakem-alphak_true),2)+alphakestd**2)**0.5;
gammaemse = (power((gammaem-gamma_true),2)+gammaestd**2)**0.5;
deltaemse = (power((deltaem-delta_true),2)+deltaestd**2)**0.5;
betaemse = (power((betaem-beta_true),2)+betaestd**2)**0.5;
rhozemse = (power((rhozem-rhoz_true),2)+rhozestd**2)**0.5;
sigmazemse = (power((sigmazem-sigmaz_true),2)+sigmazestd**2)**0.5;
sigmaremse = (power((sigmarem-sigmar_true),2)+sigmarestd**2)**0.5;
alphar1mse = (power((alphar1m-alphar1),2)+alphar1std**2)**0.5;
alphar2mse = (power((alphar2m-alphar2),2)+alphar2std**2)**0.5;

display alphakem, alphakestd, alphakemse ;
display gammaem, gammaestd, gammaemse ;
display deltaem, deltaestd, deltaemse ;
display betaem, betaestd, betaemse ;
display rhozem, rhozestd, rhozemse ;
display sigmazem, sigmazestd, sigmazemse ;
display sigmarem, sigmarestd, sigmaremse ;
display alphar1m, alphar1std, alphar1mse ;
display alphar2m, alphar2std, alphar2mse ;

display gammaeboot, deltaeboot, alphakeboot, betaeboot, rhozeboot, sigmazeboot, sigmareboot,alphar1eboot, alphar2eboot,aeboot, modelboot ;

*execute_unload 'res_smolyak_100_true.gdx';









display erreulerksv ;
