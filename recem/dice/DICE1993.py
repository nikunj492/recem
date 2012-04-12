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
############################################################################
# DICE
# William D. Nordhaus, August 25, 1993

# Annotated by Rob Dellink
# Note, some changes to the code have been made by Rob Dellink;
# this does NOT affect the parameter values or equations.
# Note: 'abatement' is the same as 'emission control' so 'abatement rate' is 'control rate'.

import sys
import os
from coopr.pyomo import *


#
# Setup
#
def createDICE1993(name='DICE1993'):
    """
     This version of DICE1993 is from 'The Environmental Economics
     and Natural Resources Group' of Wageningen University

     See the webpage for 'GAMS for environmental modeling' for more:
        http://www.enr.wur.nl/UK/gams/

    See the DICE website for a complete description of DICE:
    http://nordhaus.econ.yale.edu/DICE2007.htm

    To create an abstract optimization model of DICE1993  call:

    createDICE1993(name)

    where name can be any string and is used to give the model a short
    description. The default value is:

    name='DICE1993'
    """

    dice = AbstractModel()
    dice.name = 'DICE 93'
# SETS

    dice.L0 = Param(initialize=40)                          # L is the number of 10 yr timeperiods
    dice.L1 = Param(initialize=15)                         # L1: some parameters saturate
    dice.T = RangeSet(1, dice.L0)

#SCALAR PARAMETERS

    dice.BASE = Param(initialize=0)             # dummy for base scenario
    dice.R = Param(initialize=0.03)                # rate of social time preference per year
    dice.GL0 = Param(initialize=0.223)              # growth rate of population per year
    dice.DLAB = Param(initialize=0.195)             # decline rate of population growth per decade
    dice.DELTAM = Param(initialize=0.0833)           # removal rate carbon per decade
    dice.GA0 = Param(initialize=0.15)              # initial growth rate for technology per decade
    dice.DELA = Param(initialize=0.11)             # decline rate of technology growth per decade
    dice.SIG0 = Param(initialize=0.519)             # CO2EQ-GWP ratio
    dice.GSIGMA = Param(initialize=-0.1168)           # growth of SIGMA per decade
    dice.DK = Param(initialize=0.10)               # depreciation rate on capital per year
    dice.GAMA = Param(initialize=0.25)             # capital elasticity in output
    #dice.M0 = Param(initialize=)               # CO2EQ concentration 1965 bln tons of carbon
    dice.M0 = 677                                   # CO2EQ concentration 1965 bln tons of carbon
    dice.TL0 = Param(initialize=0.10)              # lower stratum temperature (C) 1965
    dice.T0 = Param(initialize=0.10)               # atmospheric temperature (C) 1965
    dice.ATRET = Param(initialize=0.64)            # marginal atmospheric retention rate
    dice.Q0 = Param(initialize=8.519)               # 1965 gross world output tln 1989 US$
    dice.LL0 = Param(initialize=3369)              # 1965 world population millions
    dice.K0 = Param(initialize=16.03)               # 1965 value capital billions 1989 US$
    dice.C1 = Param(initialize=0.226)               # coefficient for upper level
    dice.LAM = Param(initialize=1.41)              # climate feedback factor
    dice.C3 = Param(initialize=0.440)               # coefficient trans upper to lower stratum
    dice.C4 = Param(initialize=0.02)               # coefficient of transfer for lower level
    dice.A0 = Param(initialize=0.00963)               # initial level of total factor productivity
    dice.A1 = Param(initialize=0.0133)               # damage coeff for 2xCO2 (fraction GWP)
    dice.B1 = Param(initialize=0.0686)               # intercept control cost function
    dice.B2 = Param(initialize=2.887)                  # exponent of control cost function
    dice.PHIK = Param(initialize=140)                   # transversality coefficient capital
    dice.PHIM = Param(initialize=-9)                    # transversality coefficient carbon ($ per ton)
    dice.PHITE = Param(initialize=-7000)                # transversality coeff temp (billion $ per dC)




# COMPUTED PARAMETERS

# Labour supply in period T grows with rate GL[t]; this is calculated with a specific function:
# In the first periods, the growth rate increases rapidly but the growth of the growth rate
# declines over time so that eventually the growth rate stabilises (at around 1.14).

#param GL{t in T}      := (GL0/DLAB)*(1-exp(-DLAB*(t-1)));	#	growth rate labor 0 to T
#def GL_define(dice, t):
#    return (dice.GL0.value / dice.DLAB.value) * (1 - exp(-dice.DLAB.value * (dice.T[t] - 1)))
    def GL_define(dice, t):
        express = (dice.GL0 / dice.DLAB) * (1 - exp(-dice.DLAB * (dice.T[t] - 1)))
        return value(express)
    dice.GL = Param(dice.T, rule = GL_define)

#param L{t in T}       := LL0*exp(GL[t]);       #	level of population and labor
#def L_define(dice, t):
#    return dice.LL0.value * exp(dice.GL[t].value)
    def L_define(dice, t):
        express = dice.LL0 * exp(dice.GL[t])
        return value(express)
    dice.L = Param(dice.T, rule = L_define)


## The growth rate of technology (total factor productivity) follows a similar function.

#param GA{t in T}      = (GA0/DELA)*(1-exp(-DELA*(t-1)));   #	growth rate of TFP from 0 to T
    def GA_define(dice, t):
    #    return (dice.GA0.value / dice.DELA.value) * (1 - exp(-dice.DELA.value * (dice.T[t] - 1)))
        express = (dice.GA0 / dice.DELA) * (1 - exp(-dice.DELA * (dice.T[t] - 1)))
        return value(express)
    dice.GA = Param(dice.T, rule = GA_define)

#param AL{t in T}      = A0*exp(GA[t]);     # level of total factor productivity (TFP)
    def AL_define(dice, t):
        express = dice.A0 * exp(dice.GA[t])
        return value(express)
    dice.AL = Param(dice.T, rule = AL_define)

## The growth rate of emission intensity (emissions per unit production) is reversed:
## GSIGMA is negative, so GSIG is initally high but declines over time.

#param GSIG{t in T}    = (GSIGMA/DELA)*(1-exp(-DELA*(t-1)));	#	cumulative improvement of energy efficiency
    def GSIG_define(dice, t):
        express =  (dice.GSIGMA / dice.DELA) * (1 - exp(-dice.DELA * (dice.T[t] - 1)))
        return value(express)
    dice.GSIG = Param(dice.T, rule = GSIG_define)

#param SIGMA{t in T}   = SIG0*exp(GSIG[t]);         #	emissions-output ration
    def SIGMA_define(dice, t):
        express = dice.SIG0 * exp(dice.GSIG[t])
        return value(express)
    dice.SIGMA = Param(dice.T, rule = SIGMA_define)

## As the model uses decades as periods,
## the discount factor is based on ten times the time preference

#param RR{t in T}      = (1+R)**(10*(1-t));         #	discount factor
    dice.RR = Param(dice.T, rule = lambda dice, t: value((1 + dice.R) ** (10 * (1 - dice.T[t]))))

## The influence of other greenhouse gases is exogenous;
## for the first 15 periods a function depending on time is used; after that the value is fixed.

#param FORCOTH{t in T} =
#                if t<n_early
#                    then .2604+.125*t-.0034*t**2dice.
#                    else 1.42;
    def FORCOTH_define(dice, t):
#    express =  t < dice.L1.value and (0.2604 + 0.125 * dice.T[t] - 0.0034 * dice.T[t] * dice.T[t]) or 1.42
        return t < value(dice.L1) and (0.2604 + 0.125 * dice.T[t] - 0.0034 * dice.T[t] * dice.T[t]) or 1.42
    dice.FORCOTH = Param(dice.T, rule = FORCOTH_define)


#dice.MIU = Param(dice.T, \
#		initialize = 0.0)                     # emission control rate GHGs



#VARIABLES

    dice.MIU = Var(dice.T, bounds = (0.0, 0.999))                     # emission control rate GHGs
    dice.FORC = Var(dice.T, bounds = (0.0, None))                   # radiative forcing W per m2
    dice.TE = Var(dice.T, bounds = (0.0, 20.0))                     # temperature atmosphere C
    dice.TLO = Var(dice.T)                                            # temperature lower ocean C
    dice.M = Var(dice.T, bounds = (0.0, None), \
                  initialize = dice.M0)                      # CO2EQ concentration billion ton
    dice.E = Var(dice.T, bounds = (0.0, 1000), initialize=30)                      # CO2EQ emissions billion ton

    dice.C = Var(dice.T, initialize = 2.0, bounds = (2.0, None))                    # consumption trillion US dollar
    dice.K = Var(dice.T, bounds = (1.0, None), initialize=10)                    # capital stock trillion US dollar
    #dice.CPC = Var(dice.T)                                      #per capita consumption thousand US dollar
    #dice.PCY = Var(dice.T)                                       #per capita income thousand USdollar
    dice.I = Var(dice.T, bounds = (0.0, None), \
                  initialize = 1.0)                               #investment trillion US dollar
    #dice.S = Var(dice.T)                                         #savings rate fraction GDP
    dice.TRANS = Var()                                           #transversality variable last period

    dice.ABCOSTS = Var(dice.T)                                   #tangible relative costs related to abatement
    dice.TECOSTS = Var(dice.T)                                   #tangible relative costs related to temp rise

    dice.Y = Var(dice.T, bounds = (0.0, None), \
                  initialize = 5.0)                                #output
    dice.UTILPC = Var(dice.T)                                    #utility per capita

# CONSTRAINTS


#KK{T}           capital balance
# Capital stock equals the old capital stock net of ten years of depreciation plus 10 years of (annual) investments.
#subject to KK{t in T: 1<t<n_periods}:
#            K[t] <= (1-DK)**10*K[t-1]+10*I[t-1];

#KK0              initial condition of K
# In the first period, the capital stock is fixed.
#subject to KK0:
#            K[1] = K0;

#KC            terminal condition of K
# To avoid a zero investments in the last period they are required to be
# at least equal to the return on existing capital.
# This so-called transversality condition ensures that in the last period
# capital grows at the ('normal') steady-state rate.
#subject to KC:
#            R*K[n_periods] <= I[n_periods];

    def KK_define(dice, t):
        if t == 1:
            return (dice.K[t] == dice.K0)
        else:
            return (dice.K[t] <= (1 - dice.DK) ** 10.0 * dice.K[t - 1] + 10 * dice.I[t - 1])
    dice.KK = Constraint(dice.T, rule = KK_define)

    def KC_define(dice, t):
        if t == dice.L0.value:
            return (dice.R * dice.K[dice.L0.value] <= dice.I[dice.L0.value])
        else:
            return Constraint.Skip
    dice.KC = Constraint(dice.T, rule = KC_define)

#EE{T}           emissions process
# Emissions are 10 (years) times emissions per unit of output
# times (1 minus the abatement rate) times output.
# For output, the CES production function is used.
#subject to EE{t in T}:
#            E[t] >= 10*SIGMA[t]*(1-MIU[t])*AL[t]*L[t]**(1-GAMA)*K[t]**GAMA;

    def EE_define(dice, t):
        return (dice.E[t] >= 10.0 * dice.SIGMA[t] * (1 - dice.MIU[t]) \
                * dice.AL[t] * dice.L[t] ** (1 - dice.GAMA) * dice.K[t] ** dice.GAMA)
    dice.EE = Constraint(dice.T, rule = EE_define)

#FORCE{T}        radiative forcing
# Radiative forcings depend on CO2 concentrations and forcings of other greenhouse gases.
#subject to FORCE{t in T}:
#            FORC[t] = 4.1*(log(M[t]/590)/log(2))+FORCOTH[t];

    def FORCE_define(dice, t):
        return (dice.FORC[t] == 4.1 * (log(dice.M[t] / 590) / log(2.0)) + dice.FORCOTH[t])
    dice.FORCE = Constraint(dice.T, rule = FORCE_define)


#MM0              initial condition for M
# In the first period, CO2 concentrations are given.
#subject to MM0:
#            M[1] = M0;
#MM{T}           CO2 distribution
# Similar to the build-up of capital, CO2 concentrations equal
# previous period concentration net of 'depreciation'
# plus the portion of emissions that remains in the atmosphere.
#subject to MM{t in T: t>1}:
#            M[t] = 590+ATRET*E[t-1]+(1-DELTAM)*(M[t-1]-590);

    def MM_define(dice, t):
        if t == 1:
            return (dice.M[t] == dice.M0)
        else:
            return (dice.M[t] == 590.0 + dice.ATRET * dice.E[t - 1] \
                    + (1 - dice.DELTAM) * (dice.M[t - 1] - 590.0))
    dice.MM = Constraint(dice.T, rule = MM_define)

#TTE0            initial condition for atmospheric temperature
# In the first period, atmospheric temperature is given.
#subject to TTE0:
#          TE[1] = T0;
#TTE{T}          atmospheric temperature
# Each period, atmospheric temperature increases due to radiative forcing,
# decreases due to climate feedback and
# increases (decreases) if atmospheric temperature is lower (higher) than ocean temperature
#subject to TTE{t in T: t>1}:
#            TE[t] = TE[t-1]+C1*(FORC[t-1]-LAM*TE[t-1]-C3*(TE[t-1]-TLO[t-1]));

    def TTE_define(dice, t):
        if t == 1:
            return (dice.TE[t] == dice.T0)
        else:
            return(dice.TE[t] == dice.TE[t - 1] + dice.C1 * (dice.FORC[t - 1] \
                - dice.LAM * dice.TE[t - 1] - dice.C3 * (dice.TE[t - 1] - dice.TLO[t - 1])))
    dice.TTE = Constraint(dice.T, rule = TTE_define)

#TLE0            initial condition for lower oceanic temperature
# In the first period, the temperature of the ocean is given.
#subject to TLE0:
#           TLO[1] = TL0;
#TLE{T}          lower oceanic temperature
# The temperature difference between the atmosphere and the ocean is the only factor
# influencing ocean temperature changes.
#subject to TLE{t in T: t>1}:
#            TLO[t] = TLO[t-1]+C4*(TE[t-1]-TLO[t-1]);

    def TLE_define(dice, t):
        if t == 1:
            return (dice.TLO[t] == dice.TL0)
        else:
            return(dice.TLO[t] == dice.TLO[t - 1] + dice.C4 * (dice.TE[t - 1] - dice.TLO[t - 1]))
    dice.TLE = Constraint(dice.T, rule = TLE_define)

#ABCO{T}         tangible relative costs related to abetement
# Abatement costs are an exponential function of the abatement rate MIU
# The dummy BASE is used to let abatement costs be zero in the base simulation
#subject to ABCO{t in T}:
#            ABCOSTS[t] >=  B1*(MIU[t]**B2);

    def ABCO_define(dice, t):
        return (dice.ABCOSTS[t] >= dice.B1 * (dice.MIU[t] ** dice.B2))
    dice.ABCO = Constraint(dice.T, rule = ABCO_define)

#TECO{T}         tangible relative costs related to temperature rise
# Damage (temperature) costs are a function of atmospheric temperature
# The dummy BASE is used to let temperature costs be zero in the base simulation
#subject to TECO{t in T}:
#            TECOSTS[t] >= BASE*(1-1/(1+(A1/9)*sqrt(TE[t])));

    def TECO_define(dice, t):
        return (dice.TECOSTS[t] >= dice.BASE * (1 - 1 / (1 + (dice.A1 / 9) * sqrt(dice.TE[t]))))
    dice.TECO = Constraint(dice.T, rule = TECO_define)

#YY{T}           output
# Available production equals output (the CES production function is used again)
# minus abatement and temperature costs.
#subject to YY{t in T}:
#            Y[t] = AL[t]*L[t]**(1-GAMA)*K[t]**GAMA*(1-ABCOSTS[t])*(1-TECOSTS[t]);
#            Y[t] = AL[t]*L[t]**(1-GAMA)*K[t]**GAMA;

    def YY_define(dice, t):
        return (dice.Y[t] == dice.AL[t] * dice.L[t] ** (1 - dice.GAMA) \
                * dice.K[t] ** dice.GAMA * (1 - dice.ABCOSTS[t]) * (1 - dice.TECOSTS[t]))
    dice.YY = Constraint(dice.T, rule = YY_define)

#SEQ{T}          savings rate
# Total savings equal investments, so the savings rate equals investments divided by production.
# Note that 0.001 is used to prevent the "division by zero" error if Y is zero.
# This may be the case if no starting values are provided for Y.
#subject to SEQ{t in T}:
#            S[t] = I[t]/(.001+Y[t]);


#CC{T}           consumption
# The standard material balance: Consumption plus investments equal production.
#subject to CC{t in T}:
#            C[t] = Y[t]-I[t];

    def CC_define(dice, t):
        return (dice.C[t] == dice.Y[t] - dice.I[t])
    dice.CC = Constraint(dice.T, rule = CC_define)

#CPCE{T}         per capita consumption
# Consumption per capita is total consumption divided by population;
# the factor 1000 is used as scaling.
#subject to CPCE{t in T}:
#            CPC[t] = C[t]*1000/L[t];


#PCYE{T}         per capita income
# Similarly, per capita production can be calculated.
#subject to PCYE{t in T}:
#            PCY[t] = Y[t]*1000/L[t];


#TRANSE{T}       transversality condition;
# There will be utility derived from consumption after the model horizon.
# The parameters determining 'after-horizon-utility are calculated outside the model.
# Last period values are used to calculate utulity in the periods after the horizon.
#subject to TRANSE:
#            TRANS = RR[n_periods]*(PHIK*K[n_periods]+PHIM*M[n_periods]
#                +PHITE*TE[n_periods]);

    def TRANSE_define(dice):
        return (dice.TRANS == dice.RR[dice.L0.value] * (dice.PHIK * dice.K[dice.L0.value] \
                        + dice.PHIM * dice.M[dice.L0.value]))
    dice.TRANSE = Constraint(rule = TRANSE_define)

#PCUTIL{T}       utility per capita
# The utility in period T is the log of per capita consumption in that period.
# This is a quite common utility function.
#subject to PCUTIL{t in T}:
#            UTILPC[t] <= log(C[t]/L[t])/.55;

    def PCUTIL_define(dice, t):
        return (dice.UTILPC[t] <= log(dice.C[t] / dice.L[t]) / 0.55)
    dice.PCUTIL = Constraint(dice.T, rule = PCUTIL_define)

#OBJECTIVE            objective function
# Total present value utility is the sum of all future per capita utility levels
# times the population, times the discount factor times ten years per period,
# plus the after-horizon-utility.
#maximize OBJECTIVE:
#            TRANS + sum{ t in T} 10*RR[t]*L[t]*UTILPC[t];

    def OBJ_define(dice):
        return (dice.TRANS + sum(10.0 * dice.RR[t] * dice.L[t].value * dice.UTILPC[t] for t in dice.T))
    dice.OBJ = Objective(rule = OBJ_define, sense = maximize)

# Note that there is no equation to determine the optimal value of MIU(T).
# This variable is determined endogenously by GAMS(AMPL) as a 'free variable'.
    return dice