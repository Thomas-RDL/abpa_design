* battery.gms

* BATTERY BARON GAMS PYMOO
* Created on: jan 2022
* Author: Thomas Richard de Latour

* constants
$onEcho >cells_data.csv
n,Cc,Uc,Tc,Wc,Hc,Mc,Type,Form,Cprice,Rint,Ucmax,Ucmin,Crate,Crate_peak
0,3.5,3.65,1,0.0184,0.065,0.049,0,1,2.02,6,4.2,2.5,2.8,5
1,3.35,3.6,1,0.018,0.065,0.0474,0,1,1.8,40,4.2,2.5,2,3
2,3.45,3.6,1,0.0183,0.065,0.0495,0,1,2,5,4.2,2.5,2.5,3
3,65.6,3.6,0.0115,0.125,0.325,0.964,0,0,30,1.1,4.2,2.5,0.5,3
4,60,3.6,0.015,0.114,0.31,0.887,0,0,31,1.1,4.2,2.5,1.18,3
5,4.2,3.6,1,0.02035,0.0703,0.063,0,1,2,20,4.2,2.6,3.5,4.7
6,62.5,3.7,0.029,0.103,0.148,0.975,0,0,29,0.7,4.2,2.7,0.2,3
7,113,3.7,0.033,0.22,0.105,1.825,0,0,79.1,0.5,4.35,2.8,1,4
8,50,3.7,0.027,0.148,0.092,0.85,0,0,35,0.8,4.2,2.8,2,3
9,230,3.2,0.054,0.174,0.207,4.1,1,0,101,0.3,3.65,2.5,1,2
10,280,3.2,0.072,0.1725,0.2,5.22,1,0,116.2,0.25,3.65,2.5,1,2
11,33,3.2,0.0143,0.1445,0.1545,0.625,1,0,12,1.5,3.65,2.5,0.5,2
12,100,3.2,0.0367,0.1305,0.2005,1.92,1,0,51,0.5,3.65,2.5,1,3
13,94,3.7,0.045,0.173,0.125,2.1,0,0,38,0.75,4.15,2.7,1.6,2
14,163,3.2,0.0366,0.174,0.2305,3.19,1,0,84,0.8,3.65,2.5,1,3
15,25,3.2,0.0082,0.141,0.241,0.535,1,0,9,1.2,3.65,2.5,0.5,3
16,4,3.2,1,0.026,0.0655,0.088,1,1,1.3,20,3.65,2.5,3,3


$offEcho
Set
    i
    j;

$call csv2gdx cells_data.csv id=dataPar useHeader=y index=1 values=1,2,3,4,5,6,7,8,9,10,11,13,15 trace =0
$gdxIn cells_data.gdx
$load i = dim1
$load j = dim2
Parameter
    dataPar(i,j) 'cell data sheet'
    Cc(i)
    Uc(i)
    Mc(i)
    n(i)
    Type(i)
    Cprice(i)
    Rint(i)
    Ucmin(i)
    Ucmax(i)
    Crate_peak(i);
    
$onUNDF
$load dataPar
$offUndf

display dataPar;

n(i) = dataPar(i,"n");
Cc(i) = dataPar(i,"Cc");
Uc(i) = dataPar(i,"Uc");
Mc(i) = dataPar(i,"Mc");
Type(i) = dataPar(i,"Type");
Cprice(i) = dataPar(i,"Cprice");
Rint(i) = dataPar(i,"Rint");
Ucmin(i) = dataPar(i,"Ucmin");
Crate_peak(i) = dataPar(i,"Crate_peak");

Scalar

Eavg 'Wh.km-1 Avearage energy consumption' /160/
Tcy '45°C working temperature of the cells K' /338.15/
DOD	'Depth of discharge on average' /50/
SOC ' State of charge on average' /80/
Crate 'C-rate' /1/
PtW 'w/kg Power to weight ratio of the vehicule' /66.7/
M_base ' Mass car without battery' /1200/

*Ageing coefficients
* Qca NMC
a0 /21500/
R /8.314/

*Qca LFP
b0 /165400/
b1 /-4148/
b2 /0.01/

*Qcy NMC
c0 '(Ah-1.K-2)' /8.6772e-6/ 
c1 ' Ah-1.K-1' /-5.1613e-3/ 
c2 ' Ah-1' /7.6788e-1/ 
c3 ' crate-1' /-6.7e-3/ 

*Qcy LFP
d0 /27789/


*GWP
PtT_ef 'Plug-to-Tank efficienc q_PtT' /90/
Int_ef 'Internal efficiency n_ri' /90/
GWP_elec 'GWP value for the production of 1 kW.h of electricity' /0.0855/
LCD ' km Life cycle distance' /200000/
LCD_year 'number of km done in year' /14700/

MANUF_nmc 'GWP value for the manufacturing of 1 kg of nmc cell' /17.12/
MANUF_lfp 'GWP value for the manufacturing of 1 kg of lfp cell' /8.33/

Vef_cyl  ' Efficiency rate of cell level to system level volumetric cylindrical' /0.295/ 
Vef_p  'Efficiency rate of cell level to system level volumetric prismatic and pouch' /0.353/

Mef_cyl 'Efficiency rate of cell level to system level volumetric cylindrical' /0.552/
Mef_p 'Efficiency rate of cell level to system level volumetric prismatic and pouch' /0.575/

A ' coef GWP ' /1/
B ' coef REV ' /0/
C ' coef Cost' /1/
D ' coef Mbatt'/0/
E ' coef Vbatt'/1/
* Variables:
*    Cc  'Capacity of the cell'
*    Mc  'Mass of the cell'
*    Type    '0 = NMC, 1 = LFP*'

Variable
    X(j)
    Y(i)
    Ns 'Number of cells in series'
    Np 'Number of cells in parallels'
	Ubp	'Battery voltage'
	Ebp	'Useful energy'
	REV	'Autonomy'
	E_batt	'Energy to move the battery'
	E_ri	
	GWP	'Global Warming Potential'
    Nsy 'Cyclic Life duration'
    M_batt ' Battery mass'
    V_batt ' Battery volume'
    Cost    'Estimated cost of the battery'
    Ppeak '/ W Power produce by the motor (66.67 Power to weight ratio for peugeot e208)'
    obj 'objective function variable'

	
Binary Variable		Y;
Integer Variable    Ns, Np 'Number of cell serie/parallel';



* Objective to minimize:
Equation 
	fobj 'Objective fuNstion';
	
*fobj.. obj =e= M_batt;
*fobj.. obj =e= (A*((GWP-GWP.lo)/(GWP.up-GWP.lo))+(B*(REV.up-REV)/(REV.up-REV.lo))+(C*(Cost-Cost.lo)/(Cost.up-Cost.lo)));
fobj.. obj =e= (A*((GWP-GWP.lo)/(GWP.up-GWP.lo))+(B*(REV.up-REV)/(REV.up-REV.lo))+(C*(Cost-Cost.lo)/(Cost.up-Cost.lo))+(D*(M_batt-M_batt.lo)/(M_batt.up-M_batt.lo))+(E*(V_batt-V_batt.lo)/(V_batt.up-V_batt.lo)));


* Constraints:
Equation
    g0
	g1(j)
    g2
	g3
	g4
	g6
	g7
	g8
	g9
	g10
	g11
	g12
	g13
	g14
;
	
g0..    sum(i,Y(i))=e=1;    
g1(j).. X(j) - sum(i,Y(i)*dataPar(i,j)) =e= 0;
g2..	(Ebp*0.9 - (Eavg*REV)) =e= 0;
g3..	(Ubp-(X("Uc")*Ns)) =e= 0;
g4..	(Np-(Ebp/(Ubp*X("Cc")))) =e= 0;
g6.. 	Nsy*REV*DOD/100 - LCD =e= 0;
g7..	(E_batt-0.3*M_batt*(Ebp/1000)*LCD/(REV*(PtT_ef/100)*(M_base + M_batt))) =e= 0;
g8..	(E_ri-(1- Int_ef/100)*(Ebp/1000)*LCD/REV) =e= 0;
g9..	GWP-(((GWP_elec*E_batt)+(GWP_elec*E_ri))+Ns*Np*X("Mc")*(X("Type")*MANUF_lfp+(1-X("Type"))*MANUF_nmc)) =e= 0;
g10..   Cost =e= X("Cprice")*Ns*Np;
g11..   X("Crate_peak")*Ebp - (X("Rint")*(Ns/Np)/1000)*((Ppeak/Ubp)**2)=g= Ppeak;
g12..   M_batt =e= Ns*Np*X("Mc")*(X("Form")/Mef_cyl+(1-X("Form"))/Mef_p);
g13..   V_batt =e= (Ns*Np*X("Wc")*X("Hc")*X("Tc"))*(1-X("Form"))/Vef_p + Ns*Np*(pi*X("Hc")*(X("Wc")**2)/4)*X("Form")/Vef_cyl;
g14..   Ppeak/(M_base+M_batt) =e= PtW;

* Bounds for Variables:
X.lo("Cc") = 1;
X.up("Cc") = 280;
X.lo("Mc") = 0.01;
X.up("Mc") = 10;
X.lo("Uc")=3.2;
X.up("Uc")=3.7;
X.lo("Ucmin")=2.5;
X.up("Ucmin")=3;
X.lo("Type")=0;
X.up("Type")=1;
Ns.lo = 1.0;
Ns.up = 1000.0;
Np.lo = 1.0;
Np.up = 100.0;
Ubp.lo = 220;
Ubp.up = 880.0;
Ebp.lo = 40000.0;
Ebp.up = 100000.0;
REV.lo = 300.0;
REV.up = 552.9435;
E_batt.lo = 0.0;
E_batt.up = 10000.0;
E_ri.lo = 0.0;
E_ri.up = 10000.0;
GWP.lo = 3296.4257;
GWP.up = 7500.0;
M_batt.lo = 370.6;
M_batt.up = 700;
Nsy.lo = 0.0;
Nsy.up = 1500.0;
Cost.lo = 5421;
Cost.up = 15000;
Ppeak.lo = 50000;
Ppeak.up = 150000;
V_batt.lo = 0.245;
V_batt.up = 0.525;

* Starting point (mid of initial box):
X.l("Cc") = 10;
X.l("Mc") = 1.3;
X.l("Uc") = 3.2;
X.l("Type") = 1;
X.l("Ucmin") = 2.5;
Ns.l = 100;
Np.l = 10;
Ubp.l = 400.0;
Ebp.l = 50025.0;
REV.l = 550.0;
E_batt.l = 1000.0;
E_ri.l = 1000.0;
GWP.l = 3000.0;
Nsy.l = 1174;
Cost.l = 5000;

Model batt / all / ;

option minlp = BARON;

solve batt using minlp  minimizing obj;
Set
    k
    l;
Parameter
    alldata(k,l) 'distaNse in thousands of miles'
    M_string ' Mass of a module'
    Nsell;
Nsell = X.l("n");
  
display Ns.l, Np.l, Nsell, X.l;



$call csv2gdx cells_data.csv id=alldata useHeader=y index=2 values=1..lastcol trace =0
$gdxIn cells_data.gdx
$load k = dim1
$load l = dim2

M_string = X.l("Mc")*Np.l;

display M_string;

