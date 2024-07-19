$PLUGIN autodec

$PARAM  @ annotated
// alcohol pars - supplied in external file
ke      : log(15)  : stomach emptying (/hr)
ks      : log(0.54)  : absorption from stomach (/hr)
ka      : log(4.5)  : absorption from gut (/hr)
Vmax    : log(20)   : MM Vmax for EtOH (mmol/L/hr)
Km      :  log(10)  : MM Km for EtOH (mM)
CLd     : log(0.9)   : compartmental clearance (L/hr)
Vcent   : log(3.77)  : volume of central compartment (L)
Vtissue : log(7.4)  : volume of peripheral compartment (L)
Pstart  : 0.2 : starting PEth concentration

// PEth pars - fit with this model
kform   : log(0.68)   : formation of PEth from alcohol (umol/hr/L / mM EtOH)
CLp     : log(0.009)   : inter-reservoir exchange  (/hr)
CLrapid  : log(10.7)   : rate of rapid disappearance from reservoir 1 (/hr)
CLterm   : log(0.004)   : rate of terminal elimination from reservoir 2 (/hr)
VR1     : log(1.1)     : volume of reservoir 1
VR2     : log(0.23)     : volume of reservoir 2
VBlood  : log(3.9)    : blood volume (L)

$MAIN
D_Stomach = 0.25;
F_Stomach = 1.0;
F_Gut = 1.0;

// alcohol thetas
tke = exp(ke);
tks = exp(ks);
tka = exp(ks);
tVmax = exp(Vmax);
tKm = exp(Km);
tCLd = exp(CLd);
tVcent = exp(Vcent);
tVtissue = exp(Vtissue);

// PEth thetas
tkform = exp(kform);
tCLp = exp(CLp);
tCLrapid = exp(CLrapid);
tCLterm = exp(CLterm);
tVR1 = exp(VR1);
tVR2 = exp(VR2);
tVBlood = exp(VBlood);

PC1_0 = 0.8*Pstart*tVBlood/tVR1;
PC2_0 = 0.2*Pstart*tVBlood/tVR2;
Pconc_0 = Pstart;

$CMT
Stomach Gut Central Tissue PC1 PC2 Pconc

$ODE
// alcohol portion
double kemp = tke/(1+tks*0.615*0.005*Stomach);
double dxdt_Stomach =  -tks*Stomach - kemp*Stomach;

double dxdt_Gut = kemp*Stomach - tka*Gut;
double met = tVmax*Central/(tKm + Central);
double dxdt_Central = (tks*Stomach + tka*Gut)/tVcent + tCLd*Tissue/tVcent - tCLd*Central - met;
double dxdt_Tissue = tCLd*Central/tVtissue - tCLd*Tissue;

// PEth portion
newP = tkform*Central*0.8826;
double dxdt_PC1 = (newP - tCLp*PC1 + tCLp*PC2 - tCLterm*PC1) / tVR1;
double dxdt_PC2 = (tCLp*PC1 - tCLp*PC2 - tCLrapid*PC2) / tVR2;
double dxdt_Pconc = (dxdt_PC1*tVR1 + dxdt_PC2*tVR2)/tVBlood;

$TABLE
if(self.time==0) { Pconc2 = Pstart; }
//double Pconc2 = (PC1*tVR1 + PC2*tVBlood)/tVR2;

$CAPTURE
newP 
met
//Pconc2
