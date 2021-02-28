scalar starttime; starttime = jnow;

$ontext
*----------------
* data simulation
*----------------
Set ttot    / 1*150/
t(ttot)     / 51*150 / ;
Set boot  /1*20 /;

display coeff.l;

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
dss(ttot,boot)
bc_verif(ttot,boot)
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
zv(boot)
rv(boot)
kv(boot)
bv(boot)
;
variables
yv(boot)
uv(boot)
cv(boot)
k1v(boot)
b1v(boot)
vfs_derivk(boot)
vfs_derivb(boot)
;

equations
eqyv(boot)
equv(boot)
eqcv(boot)
eqk1v(boot)
eqb1v(boot)
eqvfs_derivk(boot)
eqvfs_derivb(boot)
;

eqyv(boot)..   yv(boot) =E= ((1/beta+delta-1)/alphak)*exp(zv(boot))*kv(boot)**alphak  ;
equv(boot)..   uv(boot) =E= vfs_derivk(boot)/ (1-delta + ((1/beta+delta-1)/alphak)*exp(zv(boot))*alphak*kv(boot)**(alphak-1)
                                                                +a*(k1v(boot)-kv(boot))/kv(boot) + a/2*(k1v(boot)-kv(boot))*(k1v(boot)-kv(boot))/kv(boot)/kv(boot)) ;
eqcv(boot)..   cv(boot) =E= uv(boot)**(-1/gamma) ;
eqb1v(boot)..  b1v(boot) =E= -(vfs_derivb(boot)/uv(boot) +1)*(rv(boot)+alphar2*bv(boot)*bv(boot))*(rv(boot)+alphar2*bv(boot)*bv(boot))
                 /(2*alphar2*bv(boot));
eqk1v(boot)..  yv(boot)+b1v(boot)/(rv(boot)+alphar2*bv(boot)*bv(boot)) - bv(boot) - cv(boot) - k1v(boot) + (1-delta)*kv(boot)
                  -a/2*(k1v(boot)-kv(boot))*(k1v(boot)-kv(boot))/kv(boot) =E= 0;

eqvfs_derivk(boot)..    vfs_derivk(boot) =E= coeffini("2")
                          + coeffini("6")*2*kv(boot) + coeffini("7")*bv(boot)+ coeffini("8")*exp(zv(boot))+ coeffini("9")*rv(boot)
*                          + coeffini("16")*3*kv(boot)*kv(boot)+ coeffini("17")*2*kv(boot)*bv(boot)+ coeffini("18")*2*kv(boot)*exp(zv(boot))+ coeffini("19")*2*kv(boot)*rv(boot)
*                          + coeffini("20")*bv(boot)*bv(boot)+ coeffini("21")*bv(boot)*exp(zv(boot))+ coeffini("22")*bv(boot)*rv(boot)
*                          + coeffini("23")*exp(zv(boot))*exp(zv(boot))+ coeffini("24")*exp(zv(boot))*rv(boot)
*                          + coeffini("25")*rv(boot)*rv(boot)
                          ;
eqvfs_derivb(boot)..    vfs_derivb(boot) =E= coeffini("3")
                          + coeffini("7")*kv(boot)+ coeffini("10")*2*bv(boot) + coeffini("11")*exp(zv(boot))+ coeffini("12")*rv(boot)
*                          + coeffini("17")*kv(boot)*kv(boot)
*                          + coeffini("20")*kv(boot)*2*bv(boot)+ coeffini("21")*kv(boot)*exp(zv(boot))+ coeffini("22")*kv(boot)*rv(boot)
*                          + coeffini("26")*3*bv(boot)*bv(boot)+ coeffini("27")*2*bv(boot)*exp(zv(boot))+ coeffini("28")*2*bv(boot)*rv(boot)
*                          + coeffini("29")*exp(zv(boot))*exp(zv(boot))+ coeffini("30")*exp(zv(boot))*rv(boot)
*                          + coeffini("31")*rv(boot)*rv(boot)
                          ;

model simulation /
eqyv
equv
eqcv
eqk1v
eqb1v
eqvfs_derivk
eqvfs_derivb
/;

kv(boot) = kss("1",boot) ;
bv(boot) = bss("1",boot);
k1v.l(boot) = kbar ;
b1v.l(boot) = bbar ;
yv.l(boot)= ybar;
cv.l(boot)= cbar;
uv.l(boot)=cbar**(-gamma);

loop(ttot$(ord(ttot) LE card(ttot)),
zv(boot) = zss(ttot,boot);
rv(boot) = rss(ttot,boot);
solve simulation using mcp;
kv(boot) = k1v.l(boot);
bv(boot) = b1v.l(boot);

css(ttot,boot)=cv.l(boot);
yss(ttot,boot)=yv.l(boot);
bss(ttot+1,boot)=b1v.l(boot);
kss(ttot+1,boot)=k1v.l(boot);
iss(ttot,boot) = kss(ttot+1,boot)-(1-delta)*kss(ttot,boot);
dss(ttot,boot) = bss(ttot+1,boot)/(rss(ttot,boot)+alphar2*bss(ttot,boot)*bss(ttot,boot))-bss(ttot,boot);
bc_verif(ttot,boot)= yss(ttot,boot)+dss(ttot,boot)-css(ttot,boot)-iss(ttot,boot)-a/2*(kss(ttot+1,boot)-kss(ttot,boot))*(kss(ttot+1,boot)-kss(ttot,boot))/kss(ttot,boot);
);

display css,kss, bss, rss, yss, zss,iss, vfs_derivb.l,vfs_derivk.l,bc_verif;


Parameters
boundmax
boundmin
titi(boot) ;
boundmax = 1.25  ;
boundmin  = 0.75 ;
titi(boot) = 1$ ( (smax(t, kss(t,boot)) GE kmax*boundmax)  or (smax(t, bss(t,boot)) GE bmax*boundmax)  or (smin(t, kss(t,boot)) LE kmin*boundmin)  or (smin(t, bss(t,boot)) LE bmin*boundmin)) ;
$offtext

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
sigmarsv(k)
alphar1sv(k)
alphar2sv(k)
asv(k)
errzsv(k)
errrsv(k)
errisv(k)
errysv(k)
erreulerksv(k)
erreulerneksv(k)
erreulerbsv(k)
erreulernebsv(k)
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
erry(t)
eulererk(t)
eulerernek(t)
eulererb(t)
eulererneb(t)

coeffe(ncoeff)

vf_derivk(t)
vf_derivb(t)
vfne(t,nzb)
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
perrysv(t,k)
perrrsv(t,k)
perreulererk(t,k)
perreulerneerk(t,k)
perreulererb(t,k)
perreulerneerb(t,k)
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
eqae
eqpalphaksv
eqpgammasv
eqpdeltasv
eqpbetasv
eqprhozsv
eqpsigmazsv
eqpsigmarsv
eqpalphar1sv
eqpalphar2sv
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
eqerry(t)
eqperrysv(t)

eqvfe(t)
eqvf_derivk(t)
eqvf_derivb(t)
eqvfne(t,nzb,nrb)
eqvfne_derivk(t,nzb,nrb)
eqvfne_derivb(t,nzb,nrb)
eqvf_FOCk(t)
eqvf_FOCb(t)

eqeulererk(t)
eqpeulererksv(t)
eqeulerneerk(t)
eqpeulerneerksv(t)

eqeulererb(t)
eqpeulererbsv(t)
eqeulerneerb(t)
eqpeulerneerbsv(t)
eqeulererkm
eqeulererbm
eqeulerneerkm
eqeulerneerbm
;


eqentropie..     entropie =E=     - predict*sum(k, palphaksv(k)*LOG(1.e-5+palphaksv(k)))
                                  - predict*sum(k, pgammasv(k)*LOG(1.e-5+pgammasv(k)))
                                  - predict*sum(k, pdeltasv(k)*LOG(1.e-5+pdeltasv(k)))
                                  - predict*sum(k, prhozsv(k)*LOG(1.e-5+prhozsv(k)))
                                  - predict*sum(k, psigmazsv(k)*LOG(1.e-5+psigmazsv(k)))
                                  - predict*sum(k, psigmarsv(k)*LOG(1.e-5+psigmarsv(k)))
                                  - predict*sum(k, palphar1sv(k)*LOG(1.e-5+palphar1sv(k)))
                                  - predict*sum(k, palphar2sv(k)*LOG(1.e-5+palphar2sv(k)))
                                  - predict*sum(k, pasv(k)*LOG(1.e-5+pasv(k)))
                                  - predict*sum(k, pbetasv(k)*LOG(1.e-5+pbetasv(k)))
                                  - predict*sum((k,t), perrzsv(t,k)*log(1.e-5+perrzsv(t,k)))
                                  - predict*sum((k,t), perrrsv(t,k)*log(1.e-5+perrrsv(t,k)))
                                  - 0*predict*sum((k,t), perrisv(t,k)*log(1.e-5+perrisv(t,k)))
                                  - 0*predict*sum((k,t), perrysv(t,k)*log(1.e-5+perrysv(t,k)))
*                                  - 100*sum((k,t), perreulererk(t,k)*log(1.e-5+perreulererk(t,k)))
                                  - 100*sum((k,t), perreulerneerk(t,k)*log(1.e-5+perreulerneerk(t,k)))
*                                  - 100*sum((k,t), perreulererb(t,k)*log(1.e-5+perreulererb(t,k)))
                                  - 100*sum((k,t), perreulerneerb(t,k)*log(1.e-5+perreulerneerb(t,k)))
                                  ;



eqvf_derivk(t).. vf_derivk(t)+ 0*eulererk(t)  =E= coeffe("2")
                  + coeffe("6")*2*ke(t)+ coeffe("7")*be(t) + coeffe("8")*exp(ze(t)) +coeffe("9")*re(t)
*                  + coeffe("16")*3*ke(t)*ke(t)+ coeffe("17")*2*ke(t)*be(t)+ coeffe("18")*2*ke(t)*exp(ze(t))+ coeffe("19")*2*ke(t)*re(t)
*                  + coeffe("20")*be(t)*be(t)+ coeffe("21")*be(t)*exp(ze(t))+ coeffe("22")*be(t)*re(t)
*                  + coeffe("23")*exp(ze(t))*exp(ze(t))+ coeffe("24")*exp(ze(t))*re(t)
*                  + coeffe("25")*re(t)*re(t)

;

eqvf_derivb(t)..  vf_derivb(t)+ 0*eulererb(t) =E= coeffe("3")
                  + coeffe("7")*ke(t) + coeffe("10")*2*be(t) + coeffe("11")*exp(ze(t)) + coeffe("12")*re(t)
*                  + coeffe("17")*ke(t)*ke(t)
*                  + coeffe("20")*ke(t)*2*be(t)+ coeffe("21")*ke(t)*exp(ze(t))+ coeffe("22")*ke(t)*re(t)
*                  + coeffe("26")*3*be(t)*be(t)+ coeffe("27")*2*be(t)*exp(ze(t))+ coeffe("28")*2*be(t)*re(t)
*                  + coeffe("29")*exp(ze(t))*exp(ze(t))+ coeffe("30")*exp(ze(t))*re(t)
*                  + coeffe("31")*re(t)*re(t)
;

eqcs(t)..          (cs(t)-0*eulererk(t))**(-gammae) =E= vf_derivk(t)/(1-deltae + ((1/betae+deltae-1)/alphake)*exp(ze(t))*alphake*ke(t)**(alphake-1)
                                                     +ae*(ke(t+1)-ke(t))/ke(t)+ae/2*(ke(t+1)-ke(t))*(ke(t+1)-ke(t))/ke(t)/ke(t)) ;

eqbe(t)..          (cs(t)-0*eulererb(t))**(-gammae)  =E=  ( vf_derivb(t) / (-2*alphar2e*be(t)*be(t+1)/((re(t)+alphar2e*be(t)*be(t))*(re(t)+alphar2e*be(t)*be(t)))-1))$(ord(t) LT card(t))
                                  + ( vf_derivb(t) / (-2*alphar2e*be(t)*bfin/((re(t)+alphar2e*be(t)*be(t))*(re(t)+alphar2e*be(t)*be(t)))-1))$(ord(t) EQ card(t)) ;



eqvfne_derivk(t,nzb,nrb)..  vfne_derivk(t,nzb,nrb) + 0*eulerernek(t)  =E= (coeffe("2")
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

eqvfne_derivb(t,nzb,nrb)..  vfne_derivb(t,nzb,nrb) + 0*eulererneb(t)  =E= (coeffe("3")
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

eqvf_FOCk(t)..       (cs(t)- eulerernek(t))**(-gammae)  =E= betae*sum((nzb,nrb),probz(nzb)*probr(nrb)*(vfne_derivk(t,nzb,nrb)))  ;

eqvf_FOCb(t)..       (cs(t)- eulererneb(t))**(-gammae)  =E= -(re(t)+alphar2e*be(t)*be(t))* betae*sum((nzb,nrb),probz(nzb)*probr(nrb)*(vfne_derivb(t,nzb,nrb)))  ;

* Modeling next period TFP expectation
eqzne(t,nzb)..        zne(t,nzb) =E=  rhoze * ze(t) + sigmaze*err0(nzb);

* Modeling TFP evolution process
eqze(t)..            ze(t+1)   =E=  rhoze * ze(t) + errz(t);

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

eqys(t)..             ys(t) =E= ((1/betae+deltae-1)/alphake)*exp(ze(t))*ke(t)**alphake + erry(t) ;

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
eqae..             ae   =E= sum(k, pasv(k)*asv(k)) ;
eqpasv..           1         =E= sum(k, pasv(k)) ;

eqerrz(t)..           errz(t)    =E= sum(k, perrzsv(t,k)*errzsv(k)) ;
eqperrzsv(t)..        1          =E= sum(k, perrzsv(t,k)) ;

eqerrr(t)..           errr(t)    =E= sum(k, perrrsv(t,k)*errrsv(k)) ;
eqperrrsv(t)..        1          =E= sum(k, perrrsv(t,k)) ;

eqerri(t)..           erri(t)    =E= sum(k, perrisv(t,k)*errisv(k)) ;
eqperrisv(t)..        1          =E= sum(k, perrisv(t,k)) ;

eqerry(t)..           erry(t)    =E= sum(k, perrysv(t,k)*errysv(k)) ;
eqperrysv(t)..        1          =E= sum(k, perrysv(t,k)) ;

eqeulererk(t)..          eulererk(t) =E= sum(k, perreulererk(t,k)*erreulerksv(k)) ;
eqpeulererksv(t)..       1          =E= sum(k, perreulererk(t,k)) ;
eqeulererkm..            sum(t,eulererk(t)) =E= 0 ;

eqeulerneerk(t)..          eulerernek(t) =E= sum(k, perreulerneerk(t,k)*erreulerneksv(k)) ;
eqpeulerneerksv(t)..       1          =E= sum(k, perreulerneerk(t,k)) ;
eqeulerneerkm..            sum(t,eulerernek(t)) =E= 0 ;

eqeulererb(t)..          eulererb(t) =E= sum(k, perreulererb(t,k)*erreulerbsv(k)) ;
eqpeulererbsv(t)..       1          =E= sum(k, perreulererb(t,k)) ;
eqeulererbm..            sum(t,eulererb(t)) =E= 0 ;

eqeulerneerb(t)..          eulererneb(t) =E= sum(k, perreulerneerb(t,k)*erreulernebsv(k)) ;
eqpeulerneerbsv(t)..       1          =E= sum(k, perreulerneerb(t,k)) ;
eqeulerneerbm..            sum(t,eulererneb(t)) =E= 0 ;


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
eqsigmare
eqpsigmarsv
eqbetae
eqpbetasv
eqalphare1
eqpalphar1sv
eqalphare2
eqpalphar2sv
eqae
eqpasv
eqke
eqre
eqze
eqcs
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
eqerry
eqperrysv
*eqvf
eqvf_derivk
eqvf_derivb
*eqvfne
eqvfne_derivk
eqvfne_derivb
eqvf_FOCk
eqvf_FOCb
*eqeulererk
*eqpeulererksv
eqeulerneerk
eqpeulerneerksv
*eqeulererb
*eqpeulererbsv
eqeulerneerb
eqpeulerneerbsv
*eqeulererkm
*eqeulererbm
eqeulerneerkm
eqeulerneerbm
/
;

*initiate chebyshev coefficient for estimation
*$ontext
betasv("1")   = 0.97 ;
betasv("2")   = 0.98 ;
betasv("3")   = 0.99 ;
alphaksv("1")   = 0.2 ;
alphaksv("2")   = 0.36 ;
alphaksv("3")   = 0.5 ;
gammasv("1")    = 0.1 ;
gammasv("2")    = 0.75 ;
gammasv("3")    = 1.5 ;
deltasv("1")    = 0.001 ;
deltasv("2")    = 0.025 ;
deltasv("3")    = 0.05 ;
rhozsv("1")    = 0 ;
rhozsv("2")    = 0 ;
rhozsv("3")    = 0 ;
sigmazsv("1")   = 0.001 ;
sigmazsv("2")   = 0.04 ;
sigmazsv("3")   = 0.08 ;
sigmarsv("1")   = 0.0001 ;
sigmarsv("2")   = 0.003 ;
sigmarsv("3")   = 0.006;
alphar1sv("1")   = 1.001 ;
alphar1sv("2")   = 1.012 ;
alphar1sv("3")   = 1.025;
alphar2sv("1")   = 0.001 ;
alphar2sv("2")   = 0.003;
alphar2sv("3")   = 0.005;
asv("1")   = 0 ;
asv("2")   = 0.1;
asv("3")   = 0.2;
*$offtext

$ontext
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
aeboot(boot)
modelboot(boot)
;

parameter
vfderivkss
vfderivbss
ref;

loop(boot$( (titi(boot) EQ 0) and (ord(boot)   )),
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
ae.l  = asv("2") ;

gammae.lo  = gammasv("1") ;
deltae.lo  = deltasv("1") ;
alphake.lo = alphaksv("1") ;
betae.lo   = betasv("1") ;
rhoze.lo    = rhozsv("1") ;
alphar1e.lo  = alphar1sv("1") ;
alphar2e.lo  = alphar2sv("1") ;
ae.lo  = asv("1") ;

gammae.up  = gammasv("3") ;
deltae.up  = deltasv("3") ;
alphake.up = alphaksv("3") ;
betae.up   = betasv("3") ;
rhoze.up    = rhozsv("3") ;
sigmaze.up  = sigmazsv("3") ;
sigmare.up  = sigmarsv("3") ;
alphar1e.up  = alphar1sv("3") ;
alphar2e.up  = alphar2sv("3") ;
ae.up  = asv("3") ;

rhoze.fx = 0;
*betae.fx=beta_true;
*gammae.fx=gamma_true;
*alphake.fx=alphak_true;
*deltae.fx=delta_true;
*sigmaze.fx=sigmaz_true;
*alphar1e.fx  = 1.012 ;
*alphar2e.fx = 0.003;

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
eulerernek.l(t) = 0 ;
eulererb.l(t)      = 0 ;
eulererneb.l(t) = 0 ;
pbetasv.l(k)            = 1/card(k) ;
palphaksv.l(k)          = 1/card(k) ;
pgammasv.l(k)           = 1/card(k) ;
pdeltasv.l(k)           = 1/card(k) ;
prhozsv.l(k)            = 1/card(k) ;
palphar2sv.l(k)          = 1/card(k) ;
psigmarsv.l(k)          = 1/card(k) ;
palphar1sv.l(k)          = 1/card(k) ;
palphar2sv.l(k)          = 1/card(k) ;
pasv.l(k)          = 1/card(k) ;
perrzsv.l(t,k)          = 1/card(k) ;
perrrsv.l(t,k)          = 1/card(k) ;
perrisv.l(t,k)          = 1/card(k) ;
perrysv.l(t,k)          = 1/card(k) ;
perreulererk.l(t,k)      = 1/card(k) ;
perreulerneerk.l(t,k) = 1/card(k) ;
perreulererb.l(t,k)      = 1/card(k) ;
perreulerneerb.l(t,k) = 1/card(k) ;


vf_derivk.l(t) =  cs(t)**(-gammae.l)*
            (1-deltae.l + ((1/betae.l+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));
vfne_derivk.l(t,nzb,nrb) =  cs(t)**(-gammae.l)*
            (1-deltae.l + ((1/betae.l+deltae.l-1)/alphake.l)*alphake.l*kbar**(alphake.l-1));

vf_derivb.l(t) =  (cs(t)**(-gammae.l) * (-2*alphar2e.l*bbar*bbar/((alphar1e.l+alphar2e.l*bbar*bbar)*(alphar1e.l+alphar2e.l*bbar*bbar))-1));
vfne_derivb.l(t,nzb,nrb) =  (cs(t)**(-gammae.l) * (-2*alphar2e.l*bbar*bbar/((alphar1e.l+alphar2e.l*bbar*bbar)*(alphar1e.l+alphar2e.l*bbar*bbar))-1));

display   vf_derivk.l,vf_derivb.l;

vfderivkss = (sum(t,is(t))/card(t))**(-gammasv("2"))*
            (1-deltasv("2") + ((1/betasv("2")+deltasv("2")-1)/alphaksv("2"))*alphaksv("2")*kbar**(alphaksv("2")-1));
vfderivbss = ((sum(t,is(t))/card(t))**(-gammasv("2")) * (-2*alphar2sv("2")*bbar*bbar/((alphar1sv("2")+alphar2sv("2")*bbar*bbar)*(alphar1sv("2")+alphar2sv("2")*bbar*bbar))-1));

*coeffe.l(ncoeff) = coeff.l(ncoeff);
coeffe.l("2")= vfderivkss ;
coeffe.l("3")=  vfderivbss  ;

predict = 0.1;

*erreulersv("1")         = -0.1*vfderivss ;
*erreulersv("3")         = 0.1*vfderivss ;
*erreulernesv("1")       = -0.1*vfderivss ;
*erreulernesv("3")       = 0.1*vfderivss ;

ref= (sum(t,cs(t))/card(t))**(-gammasv("2"));

display ref;
errisv("1")         = -0.0*(sum(t,is(t))/card(t));
errisv("3")         = 0.0*(sum(t,is(t))/card(t));
errysv("1")         = -0*(sum(t,is(t))/card(t));
errysv("3")         = 0*(sum(t,is(t))/card(t));
*erreulerksv("1")         = -0.2*sum(t,cs(t))/card(t) ;
*erreulerksv("3")         = 0.2*sum(t,cs(t))/card(t);
*erreulerbsv("1")         = -0.2*sum(t,cs(t))/card(t);
*erreulerbsv("3")         = 0.2*sum(t,cs(t))/card(t);
erreulerneksv("1")       = -0.2*sum(t,cs(t))/card(t);
erreulerneksv("3")       = 0.2*sum(t,cs(t))/card(t);
erreulernebsv("1")       = -0.2*sum(t,cs(t))/card(t);
erreulernebsv("3")       = 0.2*sum(t,cs(t))/card(t);

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
aeboot(boot) = ae.l  ;
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
aem
aestd
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
aem = sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), aeboot(boot))/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) ;
aestd = (sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)), (aeboot(boot)-aem)*(aeboot(boot)-aem) )/sum(boot$((modelboot(boot) LE 2) and (titi(boot) EQ 0)),1) )**0.5 ;

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
aemse
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
aemse = (power((aem-a_true),2)+aestd**2)**0.5;

parameters
alphakebias
gammaebias
betaebias
deltaebias
rhozebias
sigmazebias
sigmarebias
alphar1bias
alphar2bias
abias
;

alphakebias = (alphakem-alphak_true)/alphak_true;
gammaebias = (gammaem-gamma_true)/gamma_true;
betaebias = (betaem-beta_true)/beta_true;
deltaebias = (deltaem-delta_true)/delta_true;
rhozebias = (rhozem-rhoz_true)/rhoz_true;
sigmazebias = (sigmazem-sigmaz_true)/sigmaz_true;
sigmarebias = (sigmarem-sigmar_true)/sigmar_true;
alphar1bias = (alphar1m-alphar1)/alphar1;
alphar2bias = (alphar2m-alphar2)/alphar2;
abias = (aem-a_true)/a_true;


display alphakem, alphakestd, alphakemse ;
display gammaem, gammaestd, gammaemse ;
display deltaem, deltaestd, deltaemse ;
display betaem, betaestd, betaemse ;
display rhozem, rhozestd, rhozemse ;
display sigmazem, sigmazestd, sigmazemse ;
display sigmarem, sigmarestd, sigmaremse ;
display alphar1m, alphar1std, alphar1mse ;
display alphar2m, alphar2std, alphar2mse ;
display aem,aestd,aemse;

display gammaeboot, deltaeboot, alphakeboot, betaeboot, rhozeboot, sigmazeboot, sigmareboot,alphar1eboot, alphar2eboot,aeboot, modelboot ;

*execute_unload 'res_smolyak_100_true.gdx';

Parameters
res(boot,*);

res(boot,"betae") = betaeboot(boot) ;
res(boot,"gammae") = gammaeboot(boot) ;
res(boot,"alphae") = alphakeboot(boot) ;
res(boot,"deltae") = deltaeboot(boot) ;
res(boot,"rhoe") = rhozeboot(boot);
res(boot,"sigmaze") = sigmazeboot(boot);
res(boot,"sigmare") = sigmareboot(boot);
res(boot,"alphar1") = alphar1eboot(boot);
res(boot,"alphar2") = alphar2eboot(boot);
res(boot,"a") = aeboot(boot);
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

res_table("rho","True")=rhoz_true;
res_table("rho","Mean")=rhozem;
res_table("rho","S.D.")=rhozestd;
res_table("rho","Bias")=rhozebias;
res_table("rho","MSE")=rhozemse;

res_table("sigmaz","True")=sigmaz_true;
res_table("sigmaz","Mean")=sigmazem;
res_table("sigmaz","S.D.")=sigmazestd;
res_table("sigmaz","Bias")=sigmazebias;
res_table("sigmaz","MSE")=sigmazemse;

res_table("sigmar","True")=sigmar_true;
res_table("sigmar","Mean")=sigmarem;
res_table("sigmar","S.D.")=sigmarestd;
res_table("sigmar","Bias")=sigmarebias;
res_table("sigmar","MSE")=sigmaremse;

res_table("alphar1","True")=alphar1;
res_table("alphar1","Mean")=alphar1m;
res_table("alphar1","S.D.")=alphar1std;
res_table("alphar1","Bias")=alphar1bias;
res_table("alphar1","MSE")=alphar1mse;

res_table("alphar2","True")=alphar2;
res_table("alphar2","Mean")=alphar2m;
res_table("alphar2","S.D.")=alphar2std;
res_table("alphar2","Bias")=alphar2bias;
res_table("alphar2","MSE")=alphar2mse;

res_table("a","True")=a_true;
res_table("a","Mean")=aem;
res_table("a","S.D.")=aestd;
res_table("a","Bias")=abias;
res_table("a","MSE")=aemse;


execute_unload 'rbc-borrowing-kcost-true-yu.gdx',res,res_table,elapsed;








