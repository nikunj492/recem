* DICE
* William D. Nordhaus, August 25, 1993

* Annotated by Rob Dellink
* Note, some changes to the code have been made by Rob Dellink;
* this does NOT affect the parameter values or equations.
* Note: 'abatement' is the same as 'emission control' so 'abatement rate' is 'control rate'.

SETS
T               time periods                                    /1*60/
* TF(T) and TL(T) are subsets of T with obvious interpretation.
TF(T)           first period
TL(T)           last period

* Identify four scenarios to be solved (not all reported by Nordhaus)
VERSIONS        /BASE no climate change no abatement,
                 MARKET with climate change no abatement,
                 OPT_CONT with climate change and abatement,
                 CONCENT concentrations below doubling rate/
* A subset 'scenario' is declared, which will contain the active scenario.
* The set element(s) will be identified later.
SCENARIO(VERSIONS)


SCALARS
BASE            dummy for base scenario                         /1.0/
R               rate of social time preference per year         /.03/
GL0             growth rate of population per year              /.223/
DLAB            decline rate of population growth per decade    /.195/
DELTAM          removal rate carbon per decade                  /.0833/
GA0             initial growth rate for technology per decade   /.15/
DELA            decline rate of technology growth per decade    /.11/
SIG0            CO2EQ-GWP ratio                                 /.519/
GSIGMA          growth of SIGMA per decade                      /-.1168/
DK              depreciation rate on capital per year           /.10/
GAMA            capital elasticity in output                    /.25/
M0              CO2EQ concentration 1965 bln tons of carbon     /677/
TL0             lower stratum temperature (C) 1965              /.10/
T0              atmospheric temperature (C) 1965                /.20/
ATRET           marginal atmospheric retention rate             /.64/
Q0              1965 gross world output tln 1989 US$            /8.519/
LL0             1965 world population millions                  /3369/
K0              1965 value capital billions 1989 US$            /16.03/
C1              coefficient for upper level                     /.226/
LAM             climate feedback factor                         /1.41/
C3              coefficient trans upper to lower stratum        /.440/
C4              coefficient of transfer for lower level         /.02/
A0              initial level of total factor productivity      /.00963/
A1              damage coeff for 2xCO2 (fraction GWP)           /0.0133/
B1              intercept control cost function                 /.0686/
B2              exponent of control cost function               /2.887/
PHIK            transversality coefficient capital              /140/
PHIM            transversality coefficient carbon ($ per ton)   /-9/
PHITE           transversality coeff temp (billion $ per dC)    /-7000/

PARAMETERS
Results         Results from the various scenarios
Results2        More results from the various scenarios
Q(T)            level of gross production
L(T)            level of population and labor
AL(T)           level of total factor productivity (TFP)
SIGMA(T)        emissions-output ration
RR(T)           discount factor
GA(T)           growth rate of TFP from 0 to T
FORCOTH(T)      exogenous forcings from other greenhouse gases
GL(T)           growth rate labor 0 to T
GSIG(T)         cumulative improvement of energy efficiency;

* Only the first element in T goes into TF.
TF(T)      = YES$(ORD(T) EQ 1);
* Only the last element in T goes into TL.
TL(T)      = YES$(ORD(T) EQ CARD(T));

DISPLAY TF, TL;

* Labour supply in period T grows with rate GL(T); this is calculated with a specific function:
* In the first periods, the growth rate increases rapidly but the growth of the growth rate
* declines over time so that eventually the growth rate stabilises (at around 1.14).
GL(T)      = (GL0/DLAB)*(1-EXP(-DLAB*(ORD(T)-1)));
L(T)       = LL0*EXP(GL(T));
* The growth rate of technology (total factor productivity) follows a similar function.
GA(T)      = (GA0/DELA)*(1-EXP(-DELA*(ORD(T)-1)));
AL(T)      = A0*EXP(GA(T));
* The growth rate of emission intensity (emissions per unit production) is reversed:
* GSIGMA is negative, so GSIG is initally high but declines over time.
GSIG(T)    = (GSIGMA/DELA)*(1-EXP(-DELA*(ORD(T)-1)));
SIGMA(T)   = SIG0*EXP(GSIG(T));
* As the model uses decades as periods,
* the discount factor is based on ten times the time preference
RR(T)      = (1+R)**(10*(1-ORD(t)));
* The influence of other greenhouse gasses is exogenous;
* for the first 15 periods a function depending on time is used; after that the value is fixed.
FORCOTH(T) = 1.42;
FORCOTH(T)$(ORD(T) LT 15) = .2604+.125*ORD(T)-.0034*ORD(T)**2;

VARIABLES
MIU(T)          emission control rate GHGs
FORC(T)         radiative forcing W per m2
TE(T)           temperature atmosphere C
TLO(T)          temperature lower ocean C
M(T)            CO2EQ concentration billion ton
E(T)            CO2EQ emissions billion ton
C(T)            consumption trillion US dollar
K(T)            capital stock trillion US dollar
CPC(T)          per capita consumption thousand US dollar
PCY(T)          per capita income thousand US dollar
I(T)            investment trillion US dollar
S(T)            savings rate fraction GDP
TRANS           transversality variable last period
ABCOSTS(T)      tangible relative costs related to abatement
TECOSTS(T)      tangible relative costs related to temp rise
Y(T)            output
UTILPC(T)       utility per capita
UTILITY;

POSITIVE VARIABLES
MIU, E, TE, M, Y, C, K, I;

EQUATIONS
UTIL            objective function
PCUTIL(T)       utility per capita
ABCO(T)         tangible relative costs related to abetement
TECO(T)         tangible relative costs related to temperature rise
YY(T)           output
CC(T)           consumption
KK(T)           capital balance
KK0(T)          initial condition of K
KC(T)           terminal condition of K
CPCE(T)         per capita consumption
PCYE(T)         per capita income
EE(T)           emissions process
SEQ(T)          savings rate
FORCE(T)        radiative forcing
MM(T)           CO2 distribution
MM0(T)          initial condition for M
TTE(T)          atmospheric temperature
TTE0(T)         initial condition for atmospheric temperature
TLE(T)          lower oceanic temperature
TLE0(T)         initial condition for lower oceanic temperature
TRANSE(T)       transversality condition;

* Capital stock equals the old capital stock net of ten years of depreciation
* plus 10 years of (annual) investments.
KK(T+1)..    K(T+1)    =L= (1-DK)**10*K(T)+10*I(T);

* In the first period, the capital stock is fixed.
KK0(TF)..    K(TF)     =E= K0;

* To avoid a zero investments in the last period they are required to be
* at least equal to the return on existing capital.
* This so-called transversality condition ensures that in the last period
* capital grows at the ('normal') steady-state rate.
KC(TL)..     R*K(TL)   =L= I(TL);

* Emissions are 10 (years) times emissions per unit of output
* times (1 minus the abatement rate) times output.
* For output, the CES production function is used.
EE(T)..      E(T)      =G= 10*SIGMA(T)*(1-MIU(T))*AL(T)*L(T)**(1-GAMA)
                           *K(T)**GAMA;

* Radiative forcings depend on CO2 concentrations and forcings of other greenhouse gasses.
FORCE(T)..   FORC(T)   =E= 4.1*(LOG(M(T)/590)/LOG(2))+FORCOTH(T);

* In the first period, CO2 concentrations are given.
MM0(TF)..    M(TF)     =E= M0;

* Similar to the build-up of capital, CO2 concentrations equal
* previous period concentration net of 'depreciation'
* plus the portion of emissions that remains in the atmosphere.
MM(T+1)..    M(T+1)    =E= 590+ATRET*E(T)+(1-DELTAM)*(M(T)-590);

* In the first period, atmospheric temperature is given.
TTE0(TF)..   TE(TF)    =E= T0;

* Each period, atmospheric temperature increases due to radiative forcing,
* decreases due to climate feedback and
* increases (decreases) if atmospheric temperature is lower (higher) than ocean temperature
TTE(T+1)..   TE(T+1)   =E= TE(T)+C1*(FORC(T)-LAM*TE(T)-C3*(TE(T)-TLO(T)));

* In the first period, the temperature of the ocean is given.
TLE0(TF)..   TLO(TF)   =E= TL0;

* The temperature difference between the atmosphere and the ocean is the only factor
* influencing ocean temperature changes.
TLE(T+1)..   TLO(T+1)  =E= TLO(T)+C4*(TE(T)-TLO(T));

* Abatement costs are an exponential function of the abatement rate MIU
* The dummy BASE is used to let abatement costs be zero in the base simulation
ABCO(T)..    ABCOSTS(T)=G= B1*(MIU(T)**B2);

* Damage (temperature) costs are a function of atmospheric temperature
* The dummy BASE is used to let temperature costs be zero in the base simulation
TECO(T)..    TECOSTS(T)=G= BASE*(1-1/(1+(A1/9)*SQR(TE(T))));

* Available production equals output (the CES production function is used again)
* minus abatement and temperature costs.
YY(T)..      Y(T)      =E= AL(T)*L(T)**(1-GAMA)*K(T)**GAMA
                                *(1-ABCOSTS(T))*(1-TECOSTS(T));

* Total savings equal investments, so the savings rate equals investments divided by production.
* Note that 0.001 is used to prevent the "division by zero" error if Y is zero.
* This may be the case if no starting values are provided for Y.
SEQ(T)..     S(T)      =E= I(T)/(.001+Y(T));

* The standard material balance: Consumption plus investments equal production.
CC(T)..      C(T)      =E= Y(T)-I(T);

* Consumption per capita is total consumption divided by population;
* the factor 1000 is used as scaling.
CPCE(T)..    CPC(T)    =E= C(T)*1000/L(T);

* Similarly, per capita production can be calculated.
PCYE(T)..    PCY(T)    =E= Y(T)*1000/L(T);

* There will be utility derived from consumption after the model horizon.
* The parameters determining 'after-horizon-utility are calculated outside the model.
* Last period values are used to calculate utulity in the periods after the horizon.
TRANSE(TL).. TRANS     =E= RR(TL)*(PHIK*K(TL)+PHIM*M(TL)+PHITE*TE(TL));

* The utility in period T is the log of per capita consumption in that period.
* This is a quite common utility function.
PCUTIL(T)..  UTILPC(T) =L= LOG(C(T)/L(T))/.55;

* Total present value utility is the sum of all future per capita utility levels
* times the population, times the discount factor times ten years per period,
* plus the after-horizon-utility.
UTIL..       UTILITY   =E= SUM(T, 10*RR(T)*L(T)*UTILPC(T))+TRANS;

* Note that there is no equation to determine the optimal value of MIU(T).
* This variable is determined endogenously by GAMS as a 'free variable'.

* In the base, there is no abatement, so MIU has to be equal to zero.
MIU.FX(T)   = 0;

* Some lower and upper bounds are necessary for GAMS.
K.LO(T)     = 1;
TE.UP(T)    = 20;
M.LO(T)     = 600;
C.LO(T)     = 2;

* Some options are provided for the maximum number of iterations.
option iterlim  = 99900;
option reslim   = 99999;

* Options reducing the size of the listing file.
*option solprint = off;
*option limrow   = 0;
*option limcol   = 0;
option nlp=ipopt;
MODEL DICE /ALL/;


* //////////===== BASE SCENARIO =====\\\\\\\\\\
* In the base scenario, there is no abatement and hence no abatement costs
* Moreover, the temperature costs are also assumed to be zero (implemented through BASE=0).

SCENARIO(VERSIONS) = NO;
SCENARIO('BASE') = YES;
BASE = 0;
MIU.FX(T) = 0;

SOLVE DICE MAXIMIZING UTILITY USING NLP;

* Calculate an additional parameter.
Q(T) = Y.L(T)/((1-abcosts.L(T))*(1-tecosts.L(T)));
Results('BASE',T,'control rate') = MIU.L(T);
Results('BASE',T,'abat costs') = Y.L(T)*Abcosts.L(T);
Results('BASE',T,'damage costs') = Y.L(T)*Tecosts.L(T);
Results2('BASE','disc cons') = SUM(T, 10*RR(T)*C.L(T));

*Display the results:
DISPLAY SCENARIO;
DISPLAY RR, L, K.L, MIU.L, abcosts.L, tecosts.L, Y.L, Q, I.L, E.L, M.L, C.L, Forc.l;

*DISPLAY GL,L,GA,AL,GSIG,SIGMA,RR,FORCOTH;

*$ontext
* //////////===== MARKET SCENARIO =====\\\\\\\\\\
* In this scenario, there are temperature costs (BASE = 1)
* but as abatement is still fixed at zero, there are no abatement costs.
* Nordhaus has labelled this scenario 'the market solution'.

SCENARIO(VERSIONS) = NO;
SCENARIO('MARKET') = YES;
BASE        = 1;

SOLVE DICE MAXIMIZING UTILITY USING NLP;

Q(T) = Y.L(T)/((1-abcosts.L(T))*(1-tecosts.L(T)));
Results('MARKET',T,'control rate') = MIU.L(T);
Results('MARKET',T,'abat costs') = Y.L(T)*Abcosts.L(T);
Results('MARKET',T,'damage costs') = Y.L(T)*Tecosts.L(T);
Results2('MARKET','disc cons') = SUM(T, 10*RR(T)*C.L(T));

DISPLAY SCENARIO;
DISPLAY RR, L, K.L, MIU.L, abcosts.L, tecosts.L, Y.L, Q, I.L, E.L, M.L, C.L, te.l, Forc.l;

*//////////===== OPT_CONT SCENARIO =====\\\\\\\\\\
* In this third scenario, the optimal abatement rate is calculated.
* This optimal point is where marginal abatement costs equal marginal temperature costs.


SCENARIO(VERSIONS) = NO;
SCENARIO('OPT_CONT') = YES;

* The variable MIU is 'freed' by the .UP and .LO statements.
* Abatement in the (historical) first three periods is fixed at zero.
MIU.UP(T)   = 0.99;     MIU.LO(T)   = 0.01;
MIU.FX('1') = 0;        MIU.FX('2') = 0;        MIU.FX('3') = 0;

SOLVE DICE MAXIMIZING UTILITY USING NLP;

Q(T) = Y.L(T)/((1-abcosts.L(T))*(1-tecosts.L(T)));
Results('OPT_CONT',T,'control rate') = MIU.L(T);
Results('OPT_CONT',T,'abat costs') = Y.L(T)*Abcosts.L(T);
Results('OPT_CONT',T,'damage costs') = Y.L(T)*Tecosts.L(T);
Results2('OPT_CONT','disc cons') = SUM(T, 10*RR(T)*C.L(T));

DISPLAY SCENARIO;
DISPLAY RR, L, K.L, MIU.L, abcosts.L, tecosts.L, Y.L, Q, I.L, E.L, M.L, C.L, Forc.l;



*//////////===== CONCENT SCENARIO =====\\\\\\\\\\
* In the last scenario, an upper bound is placed on CO2 concentrations.
* Note that the abatement rate is still free.

SCENARIO(VERSIONS) = NO;
SCENARIO('CONCENT') = YES;
M.UP(T) = 1180;

SOLVE DICE MAXIMIZING UTILITY USING NLP;

Q(T) = Y.L(T)/((1-abcosts.L(T))*(1-tecosts.L(T)));
Results('CONCENT',T,'control rate') = MIU.L(T);
Results('CONCENT',T,'abat costs') = Y.L(T)*Abcosts.L(T);
Results('CONCENT',T,'damage costs') = Y.L(T)*Tecosts.L(T);
Results2('CONCENT','disc cons') = SUM(T, 10*RR(T)*C.L(T));

DISPLAY SCENARIO;
DISPLAY RR, L, K.L, MIU.L, abcosts.L, tecosts.L, Y.L, Q, I.L, E.L, M.L, C.L, Forc.l;

DISPLAY RESULTS, RESULTS2;
*$offtext
