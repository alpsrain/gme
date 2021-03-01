sets nz  number of nodes for productivity shock / 1*5/ ;

Parameters
err0(nz)
prob(nz)
prob_student(nz)
weight(nz)
f_student(nz)
f_normal(nz)
;

** Gaussian quadrature
* std=1
*gaussion quadrature std=1
err0("1")= -2.85697001387281; err0("2")= -1.35562617997427; err0("3")= 0; err0("4")= 1.35562617997427; err0("5")= 2.85697001387281;
prob("1")=0.0112574113277207;prob("2")=0.222075922005613; prob("3")=0.533333333333333;prob("4")=0.222075922005613 ;prob("5")=0.0112574113277207 ;

* student distribution with DF=5
* pdf student distribution with DF=5 and normal distribution N(0,1)
f_student(nz)=4*2/(pi*5**0.5*3)*(1+err0(nz)*err0(nz)/5)**(-(5+1)/2);
f_normal(nz)=1/((2*pi)**0.5) * exp(-err0(nz)*err0(nz)/2);
weight(nz)=f_student(nz)/f_normal(nz);
prob_student(nz)=weight(nz)*prob(nz);

display weight, prob, prob_student;

execute_unload 'student_quadrature.gdx',prob_student ;  