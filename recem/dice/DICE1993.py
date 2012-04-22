#==============================================================================
# RECEM: Repository of Energy, Climate, and Economics Models
#
#   NOTE: DICE and all related programs in this package [recem.dice] are in
#   the public domain. DICE and various versions of it have been in the public
#   domain since the first version of the models were published by
#   William Nordhaus in 1979.
#
#   You may use this program at your own risk!
#==============================================================================
#
# This version of DICE1993 is from the The Environmental Economics
# and Natural Resources Group of Wageningen University
# See the webpage for GAMS for environmental modeling:
#    http://www.enr.wur.nl/UK/gams/
#
# NOTE THE DIFFERENCE:
#   COMMENTED PROGRAM FROM  dice.gms
# vs.
#COMMENTS FOR THIS PROGRAM
# or
# """
# COMMENTS OF THIS TYPE
# """
#==============================================================================
#
#    * DICE
#    * William D. Nordhaus, August 25, 1993
#
#    * Annotated by Rob Dellink
#    * Note, some changes to the code have been made by Rob Dellink;
#    * this does NOT affect the parameter values or equations.
#    * Note: 'abatement' is the same as 'emission control' so 'abatement rate' is 'control rate'.
#
#    * Identify four scenarios to be solved (not all reported by Nordhaus)
#    VERSIONS        /BASE no climate change no abatement,
#                     MARKET with climate change no abatement,
#                     OPT_CONT with climate change and abatement,
#                     CONCENT concentrations below doubling rate/
#    * A subset 'scenario' is declared, which will contain the active scenario.
#    * The set element(s) will be identified later.
#    SCENARIO(VERSIONS)

from coopr.pyomo import *

def createDICE1993(name='DICE1993', BASE=1.0):
    """
     This version of DICE1993 is from 'The Environmental Economics
     and Natural Resources Group' of Wageningen University

     See the webpage for 'GAMS for environmental modeling' for more:
        http://www.enr.wur.nl/UK/gams/

    See the DICE website for a complete description of DICE:
    http://nordhaus.econ.yale.edu/DICE2007.htm

    To create an abstract optimization model of DICE1993  call:

    createDICE1993(name=NAME,BASE=value)

    where NAME can be any string and is used to give the model a short
    description. The default value is name='DICE1993'.

    BASE=value controls whether TECOSTS [Damage from temperature rise] is set to
    zero or not.
    value=0.0 [no damages], value=1.0 [damages], default is 1.0

    """

    dice = AbstractModel()
    dice.name = name
    if BASE==0.0 or BASE==1.0:
        pass
    else:
        raise ValueError("BASE should be 0.0 [no temp. damages] or 1.0 [include damages].")

#    SETS
#    T               time periods                                    /1*60/
#    * TF(T) and TL(T) are subsets of T with obvious interpretation.
#    TF(T)           first period
#    TL(T)           last period
# We don't need TF and TL as we handle these using other features of Pyomo and Python.

    dice.L0 = Param(initialize=60, doc='L is the number of 10 yr time periods')
    dice.L1 = Param(initialize=15, doc='L1: time period when some parameters\
                    saturate. eg., FORCOTH')
    dice.T = RangeSet(1, dice.L0, doc='Time Periods')

#    SCALARS
#    BASE            dummy for base scenario                         /1.0/
#    R               rate of social time preference per year         /.03/
#    GL0             growth rate of population per year              /.223/
#    DLAB            decline rate of population growth per decade    /.195/
#    DELTAM          removal rate carbon per decade                  /.0833/
#    GA0             initial growth rate for technology per decade   /.15/
#    DELA            decline rate of technology growth per decade    /.11/
#    SIG0            CO2EQ-GWP ratio                                 /.519/
#    GSIGMA          growth of SIGMA per decade                      /-.1168/
#    DK              depreciation rate on capital per year           /.10/
#    GAMA            capital elasticity in output                    /.25/
#    M0              CO2EQ concentration 1965 bln tons of carbon     /677/
#    TL0             lower stratum temperature (C) 1965              /.10/
#    T0              atmospheric temperature (C) 1965                /.20/
#    ATRET           marginal atmospheric retention rate             /.64/
#    Q0              1965 gross world output tln 1989 US$            /8.519/
#    LL0             1965 world population millions                  /3369/
#    K0              1965 value capital billions 1989 US$            /16.03/
#    C1              coefficient for upper level                     /.226/
#    LAM             climate feedback factor                         /1.41/
#    C3              coefficient trans upper to lower stratum        /.440/
#    C4              coefficient of transfer for lower level         /.02/
#    A0              initial level of total factor productivity      /.00963/
#    A1              damage coeff for 2xCO2 (fraction GWP)           /0.0133/
#    B1              intercept control cost function                 /.0686/
#    B2              exponent of control cost function               /2.887/
#    PHIK            transversality coefficient capital              /140/
#    PHIM            transversality coefficient carbon ($ per ton)   /-9/
#    PHITE           transversality coeff temp (billion $ per dC)    /-7000/

    dice.BASE = Param(initialize=BASE, doc='dummy for base scenario')
    dice.R = Param(initialize=0.03, doc='rate of social time preference per year')
    dice.GL0 = Param(initialize=0.223, doc='growth rate of population per year')
    dice.DLAB = Param(initialize=0.195, doc=' decline rate of population growth per decade')
    dice.DELTAM = Param(initialize=0.0833, doc='removal rate carbon per decade')
    dice.GA0 = Param(initialize=0.15, doc='initial growth rate for technology per decade')
    dice.DELA = Param(initialize=0.11, doc='decline rate of technology growth per decade')
    dice.SIG0 = Param(initialize=0.519, doc='CO2EQ-GWP ratio')
    dice.GSIGMA = Param(initialize=-0.1168, doc='growth of SIGMA per decade')
    dice.DK = Param(initialize=0.10, doc='depreciation rate on capital per year')
    dice.GAMA = Param(initialize=0.25, doc='capital elasticity in output')
    dice.M0 = 677 # CO2EQ concentration 1965 bln tons of carbon'
    dice.TL0 = Param(initialize=0.10, doc='lower stratum temperature (C) 1965')
    dice.T0 = Param(initialize=0.20, doc='atmospheric temperature (C) 1965')
    dice.ATRET = Param(initialize=0.64, doc='marginal atmospheric retention rate')
    dice.Q0 = Param(initialize=8.519, doc='1965 gross world output trillion 1989 US$')
    dice.LL0 = Param(initialize=3369, doc='1965 world population millions')
    dice.K0 = Param(initialize=16.03, doc='1965 value capital trillions? 1989 US$')
    dice.C1 = Param(initialize=0.226, doc='coefficient for upper level')
    dice.LAM = Param(initialize=1.41, doc='climate feedback factor')
    dice.C3 = Param(initialize=0.440, doc='coefficient trans upper to lower stratum')
    dice.C4 = Param(initialize=0.02, doc='coefficient of transfer for lower level')
    dice.A0 = Param(initialize=0.00963, doc='initial level of total factor productivity')
    dice.A1 = Param(initialize=0.0133, doc='damage coeff for 2xCO2 (fraction GWP)')
    dice.B1 = Param(initialize=0.0686, doc='intercept control cost function')
    dice.B2 = Param(initialize=2.887, doc='exponent of control cost function')
    dice.PHIK = Param(initialize=140.0, doc='transversality coefficient capital')
    dice.PHIM = Param(initialize=-9.0, doc='transversality coefficient carbon ($ per ton)')
    dice.PHITE = Param(initialize=-7000.0, doc='transversality coeff temp (billion $ per dC)')

#    PARAMETERS
#    Results         Results from the various scenarios
#    Results2        More results from the various scenarios
#    Q(T)            level of gross production
#    L(T)            level of population and labor
#    AL(T)           level of total factor productivity (TFP)
#    SIGMA(T)        emissions-output ration
#    RR(T)           discount factor
#    GA(T)           growth rate of TFP from 0 to T
#    FORCOTH(T)      exogenous forcings from other greenhouse gases
#    GL(T)           growth rate labor 0 to T
#    GSIG(T)         cumulative improvement of energy efficiency;
#
#    * Only the first element in T goes into TF.
#    TF(T)      = YES$(ORD(T) EQ 1);
#    * Only the last element in T goes into TL.
#    TL(T)      = YES$(ORD(T) EQ CARD(T));
#
#    DISPLAY TF, TL;
#
# We don't need TF and TL as we handle these using other features of Pyomo and Python.

#    * Labour supply in period T grows with rate GL(T); this is calculated with a specific function:
#    * In the first periods, the growth rate increases rapidly but the growth of the growth rate
#    * declines over time so that eventually the growth rate stabilises (at around 1.14).
#    GL(T)      = (GL0/DLAB)*(1-EXP(-DLAB*(ORD(T)-1)));
    def GL_define(dice, t):
        expression = (dice.GL0 / dice.DLAB) * (1 - exp(-dice.DLAB * (t - 1)))
        return value(expression)
    dice.GL = Param(dice.T, rule = GL_define, doc="growth rate of labor")

#    L(T)       = LL0*EXP(GL(T));
    def L_define(dice, t):
        expression = dice.LL0 * exp(dice.GL[t])
        return value(expression)
    dice.L = Param(dice.T, rule = L_define, doc="level of population and labor")


#    * The growth rate of technology (total factor productivity) follows a similar function.
#    GA(T)      = (GA0/DELA)*(1-EXP(-DELA*(ORD(T)-1)));
    def GA_define(dice, t):
        expression = (dice.GA0 / dice.DELA) * (1 - exp(-dice.DELA * (t - 1)))
        return value(expression)
    dice.GA = Param(dice.T, rule = GA_define, doc="growth rate of TFP from 0 to T")

#    AL(T)      = A0*EXP(GA(T));
    def AL_define(dice, t):
        expression = dice.A0 * exp(dice.GA[t])
        return value(expression)
    dice.AL = Param(dice.T, rule = AL_define, doc="level of total factor productivity (TFP)")

#    * The growth rate of emission intensity (emissions per unit production) is reversed:
#    * GSIGMA is negative, so GSIG is initally high but declines over time.
#    GSIG(T)    = (GSIGMA/DELA)*(1-EXP(-DELA*(ORD(T)-1)));
    def GSIG_define(dice, t):
        expression =  (dice.GSIGMA / dice.DELA) * (1 - exp(-dice.DELA * (t - 1)))
        return value(expression)
    dice.GSIG = Param(dice.T, rule = GSIG_define, doc="cumulative improvement of energy efficiency")

#    SIGMA(T)   = SIG0*EXP(GSIG(T));
    def SIGMA_define(dice, t):
        expression = dice.SIG0 * exp(dice.GSIG[t])
        return value(expression)
    dice.SIGMA = Param(dice.T, rule = SIGMA_define, doc="emissions-output ratio")

#    * As the model uses decades as periods,
#    * the discount factor is based on ten times the time preference
#    RR(T)      = (1+R)**(10*(1-ORD(t)));
    dice.RR = Param(dice.T, \
    rule = lambda dice, t: value((1 + dice.R) ** (10 * (1 - t))),\
    doc="discount factor")

#    * The influence of other greenhouse gasses is exogenous;
#    * for the first 15 periods a function depending on time is used; after that the value is fixed.
#    FORCOTH(T) = 1.42;
#    FORCOTH(T)$(ORD(T) LT 15) = .2604+.125*ORD(T)-.0034*ORD(T)**2;
    def FORCOTH_define(dice, t):
        return t < value(dice.L1) and (0.2604 + 0.125 * t - 0.0034 * t * t) or 1.42
    dice.FORCOTH = Param(dice.T, rule = FORCOTH_define)

#    VARIABLES
#
#    MIU(T)          emission control rate GHGs
#    FORC(T)         radiative forcing W per m2
#    TE(T)           temperature atmosphere C
#    TLO(T)          temperature lower ocean C
#    M(T)            CO2EQ concentration billion ton
#    E(T)            CO2EQ emissions billion ton
#    C(T)            consumption trillion US dollar
#    K(T)            capital stock trillion US dollar
#    CPC(T)          per capita consumption thousand US dollar
#    PCY(T)          per capita income thousand US dollar
#    I(T)            investment trillion US dollar
#    S(T)            savings rate fraction GDP
#    TRANS           transversality variable last period
#    ABCOSTS(T)      tangible relative costs related to abatement
#    TECOSTS(T)      tangible relative costs related to temp rise
#    Y(T)            output
#    UTILPC(T)       utility per capita
#    UTILITY;
#
#    POSITIVE VARIABLES
#    MIU, E, TE, M, Y, C, K, I;
#
#    * Some lower and upper bounds are necessary for GAMS.
#    K.LO(T)     = 1;
#    TE.UP(T)    = 20;
#    M.LO(T)     = 600;
#    C.LO(T)     = 2;


    dice.MIU = Var(dice.T, within=NonNegativeReals, bounds = (0.0, 1.0), doc="emission control rate GHGs")
    dice.FORC = Var(dice.T, bounds = (0.0, None), doc="radiative forcing W per m2")
    dice.TE = Var(dice.T, within=NonNegativeReals, initialize = 1.0,  bounds = (0.0, 20.0), doc="temperature atmosphere C")
    dice.TLO = Var(dice.T, doc="temperature lower ocean C")
    dice.M = Var(dice.T, within=NonNegativeReals, bounds = (0.0, None), initialize = dice.M0, doc="CO2EQ concentration billion ton")
    dice.E = Var(dice.T, bounds = (0.0, 1000), initialize=30.0, doc="CO2EQ emissions billion ton")
    dice.C = Var(dice.T, within=NonNegativeReals, initialize=100.0, bounds = (2.0, None), doc="consumption trillion US dollar")
    dice.K = Var(dice.T, within=NonNegativeReals, bounds = (1.0, None), initialize=10.0, doc="capital stock trillion US dollar")
    dice.CPC = Var(dice.T,  doc="per capita consumption thousand US dollar")
    dice.PCY = Var(dice.T, doc="per capita income thousand US dollar")
    dice.I = Var(dice.T, within=NonNegativeReals, bounds = (0.0, None), initialize = 1.0, doc="investment trillion US dollar")
    dice.S = Var(dice.T, initialize = 0.1, doc="savings rate fraction GDP")
    dice.TRANS = Var(doc="transversality variable last period")
    dice.ABCOSTS = Var(dice.T, initialize=0.0, doc="tangible relative costs related to abatement")
    dice.TECOSTS = Var(dice.T, initialize=0.0, doc="tangible relative costs related to temp rise")
    dice.Y = Var(dice.T, initialize = 20.0, within=NonNegativeReals, bounds = (0.0, None), doc="output")
    dice.UTILPC = Var(dice.T, doc="utility per capita")

#    EQUATIONS
#    UTIL            objective function
#    PCUTIL(T)       utility per capita
#    ABCO(T)         tangible relative costs related to abetement
#    TECO(T)         tangible relative costs related to temperature rise
#    YY(T)           output
#    CC(T)           consumption
#    KK(T)           capital balance
#    KK0(T)          initial condition of K
#    KC(T)           terminal condition of K
#    CPCE(T)         per capita consumption
#    PCYE(T)         per capita income
#    EE(T)           emissions process
#    SEQ(T)          savings rate
#    FORCE(T)        radiative forcing
#    MM(T)           CO2 distribution
#    MM0(T)          initial condition for M
#    TTE(T)          atmospheric temperature
#    TTE0(T)         initial condition for atmospheric temperature
#    TLE(T)          lower oceanic temperature
#    TLE0(T)         initial condition for lower oceanic temperature
#    TRANSE(T)       transversality condition;


#    * Capital stock equals the old capital stock net of ten years of depreciation
#    * plus 10 years of (annual) investments.
#    KK(T+1)..    K(T+1)    =L= (1-DK)**10*K(T)+10*I(T);
#
#    * In the first period, the capital stock is fixed.
#    KK0(TF)..    K(TF)     =E= K0;
#
#    * To avoid a zero investments in the last period they are required to be
#    * at least equal to the return on existing capital.
#    * This so-called transversality condition ensures that in the last period
#    * capital grows at the ('normal') steady-state rate.
#    KC(TL)..     R*K(TL)   =L= I(TL);

    def KK_define(dice, t):
        if t == 1:
            return (dice.K[t] == dice.K0)
        else:
            return (dice.K[t] <= (1 - dice.DK) ** 10.0 * dice.K[t - 1] + 10 * dice.I[t - 1])
    dice.KK = Constraint(dice.T, rule = KK_define, doc="capital balance equation with special initial conditions")

    def KC_define(dice, t):
        if t == value(dice.L0):
            return (dice.R * dice.K[value(dice.L0)] <= dice.I[value(dice.L0)])
        else:
            return Constraint.Skip
    dice.KC = Constraint(dice.T, rule = KC_define, doc="transversality condition on capital")

#    * Emissions are 10 (years) times emissions per unit of output
#    * times (1 minus the abatement rate) times output.
#    * For output, the CES production function is used.
#    EE(T)..      E(T)      =G= 10*SIGMA(T)*(1-MIU(T))*AL(T)*L(T)**(1-GAMA)
#                               *K(T)**GAMA;
    def EE_define(dice, t):
        return (dice.E[t] == 10.0 * dice.SIGMA[t] * (1 - dice.MIU[t]) \
                * dice.AL[t] * dice.L[t] ** (1 - dice.GAMA) * dice.K[t] ** dice.GAMA)
    dice.EE = Constraint(dice.T, rule = EE_define, doc="emissions  process")

#    * Radiative forcings depend on CO2 concentrations and forcings of other greenhouse gasses.
#    FORCE(T)..   FORC(T)   =E= 4.1*(LOG(M(T)/590)/LOG(2))+FORCOTH(T);

    def FORCE_define(dice, t):
        return (dice.FORC[t] == 4.1 * (log(dice.M[t] / 590.0) / log(2.0)) + dice.FORCOTH[t])
    dice.FORCE = Constraint(dice.T, rule = FORCE_define, doc="radiative forcing")


#    * In the first period, CO2 concentrations are given.
#    MM0(TF)..    M(TF)     =E= M0;
#
#    * Similar to the build-up of capital, CO2 concentrations equal
#    * previous period concentration net of 'depreciation'
#    * plus the portion of emissions that remains in the atmosphere.
#    MM(T+1)..    M(T+1)    =E= 590+ATRET*E(T)+(1-DELTAM)*(M(T)-590);

    def MM_define(dice, t):
        if t == 1:
            return (dice.M[t] == dice.M0)
        else:
            return (dice.M[t] == 590.0 + dice.ATRET * dice.E[t - 1] \
                    + (1 - dice.DELTAM) * (dice.M[t - 1] - 590.0))
    dice.MM = Constraint(dice.T, rule = MM_define, doc="CO2 concentration in the atmosphere")

#    * In the first period, atmospheric temperature is given.
#    TTE0(TF)..   TE(TF)    =E= T0;
#
#    * Each period, atmospheric temperature increases due to radiative forcing,
#    * decreases due to climate feedback and
#    * increases (decreases) if atmospheric temperature is lower (higher) than ocean temperature
#    TTE(T+1)..   TE(T+1)   =E= TE(T)+C1*(FORC(T)-LAM*TE(T)-C3*(TE(T)-TLO(T)));

    def TTE_define(dice, t):
        if t == 1:
            return (dice.TE[t] == dice.T0)
        else:
            return(dice.TE[t] == dice.TE[t - 1] + dice.C1 * (dice.FORC[t - 1] \
                - dice.LAM * dice.TE[t - 1] - dice.C3 * (dice.TE[t - 1] - dice.TLO[t - 1])))
    dice.TTE = Constraint(dice.T, rule = TTE_define, doc="temperature in the atmosphere")

#    * In the first period, the temperature of the ocean is given.
#    TLE0(TF)..   TLO(TF)   =E= TL0;
#
#    * The temperature difference between the atmosphere and the ocean is the only factor
#    * influencing ocean temperature changes.
#    TLE(T+1)..   TLO(T+1)  =E= TLO(T)+C4*(TE(T)-TLO(T));

    def TLE_define(dice, t):
        if t == 1:
            return (dice.TLO[t] == dice.TL0)
        else:
            return(dice.TLO[t] == dice.TLO[t - 1] + dice.C4 * (dice.TE[t - 1] - dice.TLO[t - 1]))
    dice.TLE = Constraint(dice.T, rule = TLE_define, doc="temperature in the ocean")

#    * Abatement costs are an exponential function of the abatement rate MIU
#    * The dummy BASE is used to let abatement costs be zero in the base simulation
#    ABCO(T)..    ABCOSTS(T)=G= B1*(MIU(T)**B2);

    def ABCO_define(dice, t):
        return (dice.ABCOSTS[t] >= dice.B1 * (dice.MIU[t] ** dice.B2))
    dice.ABCO = Constraint(dice.T, rule = ABCO_define, doc="tangible relative costs related to abatement")

#    * Damage (temperature) costs are a function of atmospheric temperature
#    * The dummy BASE is used to let temperature costs be zero in the base simulation
#    TECO(T)..    TECOSTS(T)=G= BASE*(1-1/(1+(A1/9)*SQR(TE(T))));

    def TECO_define(dice, t):
        return (dice.TECOSTS[t] >= dice.BASE * (1 - 1 / (1 + (dice.A1 / 9) * sqrt(dice.TE[t]))))
    dice.TECO = Constraint(dice.T, rule = TECO_define, doc="tangible relative costs related to temperature rise")

#    * Available production equals output (the CES production function is used again)
#    * minus abatement and temperature costs.
#    YY(T)..      Y(T)      =E= AL(T)*L(T)**(1-GAMA)*K(T)**GAMA
#                                    *(1-ABCOSTS(T))*(1-TECOSTS(T));

    def YY_define(dice, t):
        return (dice.Y[t] == dice.AL[t] * dice.L[t] ** (1 - dice.GAMA) \
                * dice.K[t] ** dice.GAMA * (1 - dice.ABCOSTS[t]) * (1 - dice.TECOSTS[t]))
    dice.YY = Constraint(dice.T, rule = YY_define, doc="output")

#    * Total savings equal investments, so the savings rate equals investments divided by production.
#    * Note that 0.001 is used to prevent the "division by zero" error if Y is zero.
#    * This may be the case if no starting values are provided for Y.
#    SEQ(T)..     S(T)      =E= I(T)/(.001+Y(T));

    def SEQ_define(dice, t):
        return (dice.S[t] == dice.I[t]/(0.001 + dice.Y[t]))
    dice.SEQ = Constraint(dice.T, rule = SEQ_define, doc="savings rate")

#    * The standard material balance: Consumption plus investments equal production.
#    CC(T)..      C(T)      =E= Y(T)-I(T);

    def CC_define(dice, t):
        return (dice.C[t] == dice.Y[t] - dice.I[t])
    dice.CC = Constraint(dice.T, rule = CC_define, doc="consumption")

#    * Consumption per capita is total consumption divided by population;
#    * the factor 1000 is used as scaling.
#    CPCE(T)..    CPC(T)    =E= C(T)*1000/L(T);

    def CPCE_define(dice, t):
        return (dice.CPC[t] == dice.C[t] * 1000 / dice.L[t])
    dice.CPCE = Constraint(dice.T, rule = CPCE_define, doc="per capita consumption")

#    * Similarly, per capita production can be calculated.
#    PCYE(T)..    PCY(T)    =E= Y(T)*1000/L(T);

    def PCYE_define(dice, t):
        return (dice.PCY[t] == dice.Y[t] * 1000 / dice.L[t])
    dice.PCYE = Constraint(dice.T, rule = PCYE_define, doc="per capita production")

#    * There will be utility derived from consumption after the model horizon.
#    * The parameters determining 'after-horizon-utility are calculated outside the model.
#    * Last period values are used to calculate utility in the periods after the horizon.
#    TRANSE(TL).. TRANS     =E= RR(TL)*(PHIK*K(TL)+PHIM*M(TL)+PHITE*TE(TL));

    def TRANSE_define(dice):
        return (dice.TRANS == dice.RR[value(dice.L0)] * (dice.PHIK * dice.K[value(dice.L0)] \
                + dice.PHIM * dice.M[value(dice.L0)] + dice.PHITE * dice.TE[value(dice.L0)]))
    dice.TRANSE = Constraint(rule = TRANSE_define, doc="transversality condition")

#    * The utility in period T is the log of per capita consumption in that period.
#    * This is a quite common utility function.
#    PCUTIL(T)..  UTILPC(T) =L= LOG(C(T)/L(T))/.55;

    def PCUTIL_define(dice, t):
        return (dice.UTILPC[t] <= log(dice.C[t] / dice.L[t]) / 0.55)
    dice.PCUTIL = Constraint(dice.T, rule = PCUTIL_define, doc="per capita utility")
#
#    * Total present value utility is the sum of all future per capita utility levels
#    * times the population, times the discount factor times ten years per period,
#    * plus the after-horizon-utility.
#    UTIL..       UTILITY   =E= SUM(T, 10*RR(T)*L(T)*UTILPC(T))+TRANS;

    def UTIL_define(dice):
        return (dice.TRANS + sum(10.0 * dice.RR[t] * dice.L[t] * dice.UTILPC[t] for t in dice.T))
    dice.UTIL = Objective(rule = UTIL_define, sense = maximize, doc="Objective function: total present value utility")

# Note that there is no equation to determine the optimal value of MIU(T).
# This variable is determined endogenously by GAMS(AMPL) as a 'free variable'.
    return dice