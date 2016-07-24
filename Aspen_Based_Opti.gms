$ TITLE  BATCHREACT
$ OFFSYMXREF
$ OFFSYMLIST

OPTION LIMROW=0;
OPTION LIMCOL=0;

Parameters
a Price of Salicylic Acid /83.488/
b Price of Acetic Anhydride /30.13/
c Price of Aspirin /133.17/
d Price of Acetic Acid /0.64/
e Price of Toluene /4.7/
f Price of Heating Water /0.02/
C1 Heat Capacity of Salicylic Acid  /161/
C2 Heat Capacity of Acetic Anhydride /168.2/
C3 Heat Capacity of Aspirin /227/
C4 Heat Capacity of Acetic Acid /139.7/
C5 Heat Capacity of Toluene /155.96/
C6 Heat Capacity of Water /75.6/
Tin Inlet temperature /302/
Tout Outlet temperature /333/
delH Heat of Reaction /23500/
K Rate Constant /0.002/
M1 Mole mass of Salicylic Acid /0.1381/
M2 Mole mass of Acetic Anhydride /0.1021/
M3 Mole mass of Aspirin /0.1802/
M4 Mole mass of Acetic Acid /0.0601/
M5 Mole mass of toluene /0.0921/
Mw Mole mass of water /0.018/

;

Equations E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11;

Positive Variables x1, x2, x3, x4, conv, xw, t
x1 Salicylic Acid Mass Flow
x2 Acetic Anhydride Mass Flow
x3 Aspirin Mass Flow
x4 Acetic Acid Mass Flow
x5 Toluene Mass Flow
xw Heating Water Mass Flow
conv Conversion
t Time Per Batch
;

Variable Obj ;

* MAXIMUM REACTOR CAPACITY (BASIS OF 10000KG)
E1.. x1 + x2 + x5 =l= 1e4;

* Acetic Anhydride Excess 50%
E2.. x2 =e= x1/M1 * 1.5 * M2;

* Stoichemistry
E3.. x3 =e= x1/M1 * conv * M3;

E4.. x4 =e= x1/M1 * conv * M4;

* ENERGY BALANCE
*E2.. XW*C6*(333-302) =E= (X1*C1+X2*C2+X5*C5)*Tin-(X3*C3+X4*C4+X5*C5+0.1*X1*C1+0.408*X2*C2)*Tout+(-delH*K*Ca*((0.1*X1+0.48*X2+X3+X4+X5)/(1.1821)));
E5.. xW * C6 / Mw * (Tout - Tin) =e= -(x1 * C1 / M1+ x2 * C2 / M2+ x5 * C5 / M5) * Tin +
(x3 *C3 / M3 + x4 * C4 / M4+ x5 * C5 / M5+ (1-conv) * x1 * C1 * M1 + (x2 - x1*conv*M2/M1) * C2 / M2)* Tout
+delH*conv*x1/M1;

* Conversion vs Time
E6.. Conv =e= 1 - exp(-k*t) ;

E7.. x3/t*60*24 =g= 3526.4;

E8.. x5 =e= 0.868*x1;
E9.. Conv =l= 0.948;
E11.. Conv =g= 0.8;

E10.. Obj =e= (x3*c + x4 * d - x1 * a - x2 * b - f*xW)*24*60/t;

t.l = 20*60;
MODEL BATCHREACT /ALL/;

SOLVE BATCHREACT USING NLP MAXIMISING OBJ;
* Use SNOPT solver to optimize
*OPTION NLP = NILP;
Option NLP = ipopt
