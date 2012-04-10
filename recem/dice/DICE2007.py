#==============================================================================
# This is a translation of the GAMS version of  William Nordhaus's DICE:
# http://nordhaus.econ.yale.edu/DICE2007_short.gms
#
# See the DICE website
# http://nordhaus.econ.yale.edu/DICE2007.htm
#==============================================================================


#  DICE delta version 8
#  July 17, 2008.
#  This version is used for the DICE book, A Question of Balance (YUP, 2008).
#  We have included only the base, Hotelling, and optimal runs.
#  Exclude statements are removed so that it can run as a self-contained program.
#  Created September 5, 2008.
#
#  Note that this can be loaded into a data reading program,

import sys
import os
from coopr.pyomo import *


##
## Setup
##

M = AbstractModel()
M.name = 'DICE2007'

# SCALARS FOR SET DEFINITIONS
M.L0 = Param(initialize = 60, doc = 'Number of Time Periods')               # Time periods
M.L1 = Param(initialize = 12, doc = 'Time Period at which some parameters saturate')
M.L2 = Param(initialize = 25, doc = 'Time Period at which some parameters saturate')

#SETS
M.T = RangeSet(1, M.L0, doc = 'Time Periods')                # Time periods

#SCALARS

#** Preferences
M.B_ELASMU = Param(initialize = 2.0, doc = 'Elasticity of marginal utility of consumption')
M.B_PRSTP = Param(initialize = 0.015, doc = 'Initial rate of social time preference per year')

#  DON'T DECLARE AS PARAMS BUT AS NUMBERS IN CLASS M (PYTHON HACK)
#M.B_ELASMU = 2.0 # Elasticity of marginal utility of consumption
#M.B_PRSTP = 0.015 # Initial rate of social time preference per year

#** Population and technology
M.POP0 = Param(initialize = 6514, doc = '2005 world population millions')
M.GPOP0 = Param(initialize = 0.35, doc = 'Growth rate of population per decade')
M.POPASYM = Param(initialize = 8600, doc = 'Asymptotic population')
M.A0 = Param(initialize = 0.02722, doc = 'Initial level of total factor productivity')
M.GA0 = Param(initialize = 0.092, doc = 'Initial growth rate for technology per decade')
M.DELA = Param(initialize = 0.001, doc = 'Decline rate of technol change per decade')
M.DK = Param(initialize = 0.1, doc = 'Depreciation rate on capital per year')
M.GAMA = Param(initialize = 0.3, doc = 'Capital elasticity in production function')
M.Q0 = Param(initialize = 61.1, doc = '2005 world gross output trill 2005 US dollars')
M.K0 = Param(initialize = 137, doc = '2005 value capital trill 2005 US dollars')

#** Emissions
M.SIG0 = Param(initialize = 0.13418, doc = 'CO2-equivalent emissions-GNP ratio 2005')
M.GSIGMA = Param(initialize = -0.073, doc = 'Initial growth of sigma per decade')
M.DSIG = Param(initialize = 0.003, doc = 'Decline rate of decarbonization per decade')
M.DSIG2 = Param(initialize = 0.000, doc = 'Quadratic term in decarbonization')
M.ELAND0 = Param(initialize = 11.000, doc = 'Carbon emissions from land 2005(GtC per decade')

#** Carbon cycle
M.MAT2000 = Param(initialize = 808.9, doc = 'Concentration in atmosphere 2005 (GtC)')
M.MU2000 = Param(initialize = 1255, doc = 'Concentration in upper strata 2005 (GtC)')
M.ML2000 = Param(initialize = 18365, doc = 'Concentration in lower strata 2005 (GtC)')
M.B11 = Param(initialize = 0.810712, doc = 'Carbon cycle transition matrix')
M.B12 = Param(initialize = 0.189288, doc = 'Carbon cycle transition matrix')
M.B21 = Param(initialize = 0.097213, doc = 'Carbon cycle transition matrix')
M.B22 = Param(initialize = 0.852787, doc = 'Carbon cycle transition matrix')
M.B23 = Param(initialize = 0.05, doc = 'Carbon cycle transition matrix')
M.B32 = Param(initialize = 0.003119, doc = 'Carbon cycle transition matrix')
M.B33 = Param(initialize = 0.996881, doc = 'Carbon cycle transition matrix')

#** Climate model
M.T2XCO2 = Param(initialize = 3, doc = 'Equilibrium temp impact of CO2 doubling oC')
#M.T2XCO2 = 3 # Equilibrium temp impact of CO2 doubling oC

M.FEX0 = Param(initialize = -0.06, doc = 'Estimate of 2000 forcings of non-CO2 GHG')
M.FEX1 = Param(initialize = 0.3, doc = 'Estimate of 2100 forcings of non-CO2 GHG')
M.TOCEAN0 = Param(initialize = 0.0068, doc = '2000 lower strat. temp change (C) from 1900')
M.TATM0 = Param(initialize = 0.7307, doc = '2000 atmospheric temp change (C)from 1900')
M.C1 = Param(initialize = 0.22, doc = 'Climate-equation coefficient for upper level')
M.C3 = Param(initialize = 0.3, doc = 'Transfer coeffic upper to lower stratum')
M.C4 = Param(initialize = 0.05, doc = 'Transfer coeffic for lower level')
M.FCO22X = Param(initialize = 3.8, doc = 'Estimated forcings of equilibrium CO2 doubling')
#M.FCO22X = 3.8 # Estimated forcings of equilibrium co2 doubling

#** Climate damage parameters calibrated for quadratic at 2.5 C for 2105
M.A1 = Param(initialize = 0, doc = 'Damage intercept')
M.A2 = Param(initialize = 0.0028388, doc = 'Damage quadratic term')
M.A3 = Param(initialize = 2, doc = 'Damage exponent')

#  DON'T DECLARE AS PARAMS BUT AS NUMBERS IN CLASS M (PYTHON HACK)
#M.A1 = 0 # Damage intercept
#M.A2 = 0.0028388 # Damage quadratic term
#M.A3 = 2 # Damage exponent

#** Abatement cost
M.EXPCOST2 = Param(initialize = 2.8, doc = 'Exponent of control cost function')
M.PBACK = Param(initialize = 1.17, doc = 'Cost of backstop 2005 000$ per tC 2005')
M.BACKRAT = Param(initialize = 2, doc = 'Ratio initial to final backstop cost')
M.GBACK = Param(initialize = 0.05, doc = 'Initial cost decline backstop pc per decade')
#M.LIMMIU = Param(initialize = 1, doc = 'Upper limit on control rate')
M.LIMMIU = 1 # ANOTHER PYTHON PYOMO HACK

#** Participation
#M.PARTFRACT1 = Param(initialize = 1, doc = 'Fraction of emissions under control regime 2005')
M.PARTFRACT1 = Param(initialize = 0.25372, doc = 'Fraction of emissions under control regime 2005')
M.PARTFRACT2 = Param(initialize = 1, doc = 'Fraction of emissions under control regime 2015')
M.PARTFRACT21 = Param(initialize = 1, doc = 'Fraction of emissions under control regime 2205')
M.DPARTFRACT = Param(initialize = 0, doc = 'Decline rate of participation')

#** Availability of fossil fuels
#M.FOSSLIM = Param(initialize = 6000, doc = 'Maximum cumulative extraction fossil fuels')
M.FOSSLIM = 6000 # ANOTHER PYTHON PYOMO HACK

#** Scaling and inessential parameters')
M.SCALE1 = Param(initialize = 194, doc = 'Scaling coefficient in the objective function')
M.SCALE2 = Param(initialize = 381800, doc = 'Scaling coefficient in the objective function')


#* Definitions for outputs of no economic interest
#SETS
#TODO
#M.TFIRST = Set(M.T)
#M.TLAST = Set(M.T)
#M.TEARLY = Set(M.T)
#M.TLATE = Set(M.T)

# COMPUTED PARAMETERS

#M.AA1 = M.A1;
#M.AA2 = M.A2;
#M.AA3 = M.A3;
#M.ELASMU = M.B_ELASMU;
#M.PRSTP  = M.B_PRSTP;

# two ways of doing the same thing: see above

#M.AA1 = Param(initialize = M.A1, doc = 'Variable A1')
#M.AA2 = Param(initialize = M.A2, doc = 'Variable A2')
#M.AA3 = Param(initialize = M.A3, doc = 'Variable A3')
#M.ELASMU = Param(initialize = M.B_ELASMU, doc = 'Variable elasticity of marginal utility of consumption')
#M.PRSTP = Param(initialize = M.B_PRSTP, doc = 'Variable initial rate of social time preference per year'
#M.LAM = Param(initialize = M.FCO22X/M.T2XCO2, doc = 'Climate model parameter')

M.AA1 = Param(rule = lambda M: M.A1.value, doc = 'Variable A1')
M.AA2 = Param(rule = lambda M: M.A2.value, doc = 'Variable A2')
M.AA3 = Param(rule = lambda M: M.A3.value, doc = 'Variable A3')
M.ELASMU = Param(rule = lambda M: M.B_ELASMU.value, doc = 'Variable elasticity of marginal utility of consumption')
M.PRSTP = Param(rule = lambda M: M.B_PRSTP.value, doc = 'Variable initial rate of social time preference per year')
M.LAM = Param(rule = lambda M: M.FCO22X.value/M.T2XCO2.value, doc = 'Climate model parameter')

#
#* Unimportant definitions to reset runs
#M.TFIRST(T) = YES$(ORD(T) EQ 1);
#M.TLAST(T)  = YES$(ORD(T) EQ CARD(T));
#M.TEARLY(T) = YES$(ORD(T) LE 20);
#M.TLATE(T)  = YES$(ORD(T) GE 21);

#TODO
#M.B11 = 1 - M.B12;
#M.B21 = 587.473*M.B12/1143.894;
#M.B22 = 1 - M.B21 - M.B23;
#M.B32 = 1143.894*M.B23/18340;
#M.B33 = 1 - M.B32 ;


#* Important parameters for the model
#
#M.R = Param(M.T, doc = 'Instantaeous rate of social time preference')
#M.GL = Param(M.T, doc = 'Growth rate of labor 0 to T')
#M.GCOST1 = Param(M.T, doc = 'Growth of cost factor')


#GFACPOP(T) =   (exp(GPOP0*(ORD(T)-1))-1)/exp(GPOP0*(ORD(T)-1));
def GFACPOP_rule(M, t):
    GPOP0 = M.GPOP0.value
    return ((exp(GPOP0*(t-1)) -1)/exp(GPOP0*(t-1)))
M.GFACPOP = Param(M.T, rule = GFACPOP_rule, doc = 'Growth factor population')

#L(T)=POP0* (1- GFACPOP(T))+GFACPOP(T)*POPASYM;
def L_rule(M, t):
    return(M.POP0.value*(1 - M.GFACPOP[t].value) + M.GFACPOP[t].value*M.POPASYM.value)
M.L = Param(M.T, rule = L_rule, doc = 'Level of population and labor')

#GA(T)=GA0*EXP(-DELA*10*(ORD(T)-1));
def GA_rule(M, t):
    return (M.GA0.value*exp(-M.DELA.value*10*(t - 1)))
M.GA = Param(M.T, rule = GA_rule, doc = 'Growth rate of productivity from 0 to T')

#AL("1") = A0;
#LOOP(T, AL(T+1)=AL(T)/((1-GA(T))););
def AL_rule(M, t):
    if t == 1:
        return M.A0.value
    else:
        return M.AL[t-1].value/(1 - M.GA[t-1].value)
M.AL = Param(M.T, rule = AL_rule, doc = 'Level of total factor productivity')


#GSIG(T)=GSIGMA*EXP(-DSIG*10*(ORD(T)-1)-DSIG2*10*((ORD(T)-1)**2));
def GSIG_rule(M, t):
    return (M.GSIGMA.value*exp(-M.DSIG.value*10*(t-1) - M.DSIG2.value*10*(t-1)**2))
M.GSIG = Param(M.T, rule = GSIG_rule, doc = 'Cumulative improvement of energy efficiency')

#SIGMA("1")=SIG0;
#LOOP(T,SIGMA(T+1)=(SIGMA(T)/((1-GSIG(T+1)))););
def SIGMA_rule(M,t):
    if t == 1:
        return M.SIG0.value
    else:
        return M.SIGMA[t-1].value/(1 - M.GSIG[t].value)
M.SIGMA = Param(M.T, rule = SIGMA_rule, doc = 'CO2-equivalent-emissions output ratio')


#COST1(T) = (PBACK*SIGMA(T)/EXPCOST2)* ( (BACKRAT-1+ EXP (-GBACK* (ORD(T)-1) ) )/BACKRAT);
def COST1_rule(M, t):
    return ((M.PBACK.value*M.SIGMA[t].value/M.EXPCOST2.value)*((M.BACKRAT.value - 1 + exp(-M.GBACK.value*(t-1)))/M.BACKRAT.value))
M.COST1 = Param(M.T, rule = COST1_rule, doc = 'Adjusted cost for backstop')

#ETREE(T) = ELAND0*(1-0.1)**(ORD(T)-1);
def ETREE_rule(M, t):
    return M.ELAND0.value*(1 - 0.1)**(t-1)
M.ETREE = Param(M.T, rule = ETREE_rule, doc = 'Emissions from deforestation')

#RR(T)=1/((1+PRSTP)**(10*(ORD(T)-1)));
def RR_rule(M, t):
    return (1 + M.PRSTP.value)**(10*(1 - t))
M.RR = Param(M.T, rule = RR_rule, doc = 'Average utility social discount rate')

#FORCOTH(T)= FEX0+ .1*(FEX1-FEX0)*(ORD(T)-1)$(ORD(T) LT 12)+ 0.36$(ORD(T) GE 12);
def FORCOTH_rule(M, t):
    return (t < M.L1.value  and (M.FEX0.value + 0.1*(M.FEX1.value - M.FEX0.value)*(t-1)) or (M.FEX0.value + 0.36))
M.FORCOTH = Param(M.T, rule = FORCOTH_rule, doc = 'Exogenous forcing for other greenhouse gases')

#PARTFRACT(T) = PARTFRACT21;
#PARTFRACT(T)$(ORD(T)<25) = PARTFRACT21 + (PARTFRACT2-PARTFRACT21)*EXP(-DPARTFRACT*(ORD(T)-2));
#PARTFRACT("1")= PARTFRACT1;
def PARTFRACT_rule(M, t):
    if t ==1:
        return M.PARTFRACT1.value
    elif 1 < t < M.L2.value:
        return M.PARTFRACT21.value + (M.PARTFRACT2.value - M.PARTFRACT21.value)*exp(-M.DPARTFRACT.value*(t-2))
    else:
        return M.PARTFRACT21.value
M.PARTFRACT = Param(M.T, rule = PARTFRACT_rule, doc = 'Fraction of emissions in control regime')


#VARIABLES
M.MIU = Var(M.T, within = NonNegativeReals, bounds = (0.0, M.LIMMIU), doc = 'Emission control rate GHGs')
M.FORC = Var(M.T, doc = 'Radiative forcing in watts per m2')
M.TATM = Var(M.T, within = NonNegativeReals, doc = 'Temperature of atmosphere in degrees C')
M.TOCEAN = Var(M.T, within = NonNegativeReals, doc = 'Temperature of lower oceans degrees C')
M.MAT = Var(M.T, within = NonNegativeReals, doc = 'Carbon concentration in atmosphere GtC')
M.MATAV = Var(M.T, within = NonNegativeReals, doc = 'Average concentrations')
M.MU = Var(M.T, within = NonNegativeReals, doc = 'Carbon concentration in shallow oceans Gtc')
M.ML = Var(M.T, within = NonNegativeReals, doc = 'Carbon concentration in lower oceans GtC')
M.E = Var(M.T, within = NonNegativeReals, doc = 'CO2-equivalent emissions GtC')
M.C = Var(M.T, initialize = 5, within = NonNegativeReals, doc = 'Consumption trillions US dollars')
M.K = Var(M.T, initialize = 10, within = NonNegativeReals, doc = 'Capital stock trillions US dollars')
M.CPC = Var(M.T, doc = 'Per capita consumption thousands US dollars')
M.PCY = Var(M.T, doc = 'Per capita income thousands US dollars')
M.I = Var(M.T, within = NonNegativeReals, doc = 'Investment trillions US dollars')
M.S = Var(M.T, doc = 'Gross savings rate as fraction of gross world product')
M.RI = Var(M.T, doc = 'Real interest rate per annum')
M.Y = Var(M.T, within = NonNegativeReals, doc = 'Gross world product net of abatement and damages')
M.YGROSS = Var(M.T, within = NonNegativeReals, doc = 'Gross world product GROSS of abatement and damages')
M.YNET = Var(M.T, doc = 'Output net of damages equation')
M.DAMAGES = Var(M.T, doc = 'Damages')
M.ABATECOST = Var(M.T, doc = 'Cost of emissions reductions')
M.CCA = Var(M.T, within = NonNegativeReals, doc = 'Cumulative industrial carbon emissions GTC')
M.PERIODU = Var(M.T, doc = 'One period utility function')
M.UTILITY = Var(M.T, doc = 'Total utility')



##**  UPPER AND LOWER BOUNDS: GENERAL CONDITIONS FOR STABILITY

#K.LO(T)         = 100;
M.K.bounds = (100, None)
#MAT.LO(T)       = 10;
M.MAT.bounds = (10, None)
#MU.LO(T)        = 100;
M.MU.bounds = (100, None)
#ML.LO(T)        = 1000;
M.ML.bounds = (1000, None)
#C.LO(T)         = 20;
M.C.bounds = (20, None)
#TOCEAN.UP(T)    = 20;
#TOCEAN.LO(T)    = -1;
M.TOCEAN.bounds = (-1, 20)
#TATM.UP(T)      = 20;
M.TATM.bounds = (0, 20)
#MIU.UP(T)       = LIMMIU;
M.MIU.bounds = (0, M.LIMMIU)
#PARTFRACT("1")= 0.25372;
#M.PARTFRACT.initialize = 0.25372

##* First period predetermined by Kyoto Protocol
#MIU.FX("1")     = 0.005;
def FIXMIU_rule(M, t):
    if t==1:
        return (0.005, M.MIU[t], 0.005)
    else:
        return Constraint.Skip
M.FIXMIU = Constraint(M.T, rule = FIXMIU_rule)

##** Fix savings assumption for standardization if needed
#S.FX(T)=.22;
M.S.bounds = (0.22,0.22)

##** Cumulative limits on carbon use at 6000 GtC
#CCA.UP(T) = FOSSLIM;
M.CCA.bounds = (None, M.FOSSLIM)



#EQUATIONS


#** Equations of the model

#CCTFIRST(TFIRST).. CCA(TFIRST)=E=0;
#CCACCA(T+1)..      CCA(T+1)=E=CCA(T)+ E(T);
def CCACCA_rule(M, t):
    if t == 1:
        return (M.CCA[t] == 0)
    else:
        return (M.CCA[t] == M.CCA[t-1] + M.E[t-1])
#M.CCTFIRST = Constraint(M.T, rule = CCTFIRST_rule, doc = 'First period cumulative carbon')
M.CCACCA = Constraint(M.T, rule = CCACCA_rule, doc = 'Cumulative carbon emissions')

#KK0(TFIRST)..      K(TFIRST) =E= K0;
#KK(T)..            K(T+1) =L= (1-DK)**10 *K(T)+10*I(T);
#KC(TLAST)..        .02*K(TLAST) =L= I(TLAST);
#M.KK0 = Constraint(M.T, rule = KK0_rule, doc = 'Initial condition for capital')
#M.KC = Constraint(M.T, rule = KC_rule, doc = 'Terminal condition for capital')
def KK_rule(M,t):
    if t == 1:
        return (M.K[t] == M.K0)
    else:
        return (M.K[t] <= (1 - M.DK)**10 *M.K[t-1] + 10*M.I[t-1])
M.KK = Constraint(M.T, rule = KK_rule, doc = 'Capital balance equation with special initial conditions')

def KC_rule(M,t):
    if t == M.L0.value:
        return (0.02*M.K[t] <= M.I[t])
    else:
        return Constraint.Skip
M.KC = Constraint(M.T, rule = KC_rule, doc = 'Capital balance terminal conditions')

#EE(T)..            E(T)=E=10*SIGMA(T)*(1-MIU(T))*AL(T)*L(T)**(1-GAMA)*K(T)**GAMA + ETREE(T);
def EE_rule(M, t):
    return (M.E[t] == 10*M.SIGMA[t]*(1 - M.MIU[t])*M.AL[t]*M.L[t]**(1-M.GAMA)*M.K[t]**M.GAMA + M.ETREE[t])
M.EE = Constraint(M.T, rule = EE_rule, doc = 'Emissions equation')

#FORCE(T)..         FORC(T) =E=  FCO22X*((log((Matav(T)+.000001)/596.4)/log(2)))+FORCOTH(T);
def FORCE_rule(M, t):
    return (M.FORC[t] == M.FCO22X*((log((M.MATAV[t] + 0.000001)/596.4)/log(2.0))) + M.FORCOTH[t])
M.FORCE = Constraint(M.T, rule = FORCE_rule, doc = 'Radiative forcing equation')

#MMAT0(TFIRST)..    MAT(TFIRST) =E= MAT2000;
#MMAT(T+1)..        MAT(T+1)    =E= MAT(T)*b11+MU(T)*b21 + E(T);
#M.MMAT0 = Constraint(M.T, rule = MMAT0_rule, doc = 'Starting atmospheric concentration')
def MMAT_rule(M, t):
    if t==1:
        return (M.MAT[t] == M.MAT2000)
    else:
        return (M.MAT[t] == M.MAT[t-1]*M.B11 + M.MU[t-1]*M.B21 + M.E[t-1])
M.MMAT = Constraint(M.T, rule = MMAT_rule, doc = 'Atmospheric concentration equation')

#MMU0(TFIRST)..     MU(TFIRST)  =E= MU2000;
#MMU(T+1)..         MU(T+1)     =E= MAT(T)*b12+MU(T)*b22+ML(T)*b32;
#M.MMU0 = Constraint(M.T, rule = MMU0_rule, doc = 'Initial shallow ocean concentration')
def MMU_rule(M, t):
    if t==1:
        return (M.MU[t] == M.MU2000)
    else:
        return (M.MU[t] == M.MAT[t-1]*M.B12 + M.MU[t-1]*M.B22 + M.ML[t-1]*M.B32)
M.MMU = Constraint(M.T, rule = MMU_rule, doc = 'Shallow ocean concentration')

#MML0(TFIRST)..     ML(TFIRST)  =E= ML2000;
#MML(T+1)..         ML(T+1)     =E= ML(T)*b33+b23*MU(T);
#M.MML0 = Constraint(M.T, rule = MML0_rule, doc = 'Initial lower ocean concentration')
def MML_rule(M, t):
    if t==1:
        return (M.ML[t] == M.ML2000)
    else:
        return (M.ML[t] == M.ML[t-1]*M.B33 + M.MU[t-1]*M.B23)
M.MML = Constraint(M.T, rule = MML_rule, doc = 'Lower ocean concentration')

#MMATAVEQ(t)..      MATAV(T)    =E= (MAT(T)+MAT(T+1))/2
def MMATAVEQ_rule(M, t):
    if t < M.L0:
        return (M.MATAV[t] == 0.5*(M.MAT[t] + M.MAT[t+1]))
    else:
        return (M.MATAV[t] == M.MAT[t])
M.MMATAVEQ = Constraint(M.T, rule = MMATAVEQ_rule, doc = 'Average concentrations equation')


#TATM0EQ(TFIRST)..  TATM(TFIRST) =E= TATM0;
#TATMEQ(T+1)..      TATM(T+1) =E= TATM(t)+C1*(FORC(t+1)-LAM*TATM(t)-C3*(TATM(t)-TOCEAN(t)));
#M.TATM0EQ = Constraint(M.T, rule = TATM0EQ_rule, doc = 'Initial condition for atmospheric temperature')
def TATMEQ_rule(M, t):
    if t == 1:
        return (M.TATM[t] == M.TATM0)
    else:
        return (M.TATM[t] == M.TATM[t-1] + M.C1*(M.FORC[t] - M.LAM*M.TATM[t-1] - M.C3*(M.TATM[t-1] - M.TOCEAN[t-1])))
M.TATMEQ = Constraint(M.T, rule = TATMEQ_rule, doc = 'Temperature-climate equation for atmosphere')

#TOCEAN0EQ(TFIRST)..  TOCEAN(TFIRST) =E= TOCEAN0;
#TOCEANEQ(T+1)..    TOCEAN(T+1) =E= TOCEAN(T)+C4*(TATM(T)-TOCEAN(T));
#M.TOCEAN0EQ = Constraint(M.T, rule = TOCEAN0EQ_rule, doc = 'Initial condition for lower ocean temperature')
def TOCEANEQ_rule(M, t):
    if t == 1:
        return (M.TOCEAN[t] == M.TOCEAN0)
    else:
        return (M.TOCEAN[t] == M.TOCEAN[t] + M.C4*(M.TATM[t]-M.TOCEAN[t]))
M.TOCEANEQ = Constraint(M.T, rule = TOCEANEQ_rule, doc = 'Temperature-climate equation for lower oceans')

#YGROSSEQ(T)..   YGROSS(T) =e= AL(T)*L(T)**(1-GAMA)*K(T)**GAMA;
def YGROSSEQ_rule(M, t):
    return (M.YGROSS[t] == M.AL[t]*M.L[t]**(1-M.GAMA)*M.K[t]**M.GAMA)
M.YGROSSEQ = Constraint(M.T, rule = YGROSSEQ_rule, doc = 'Output gross equation')

#DAMEQ(T)..      DAMAGES(t) =E= YGROSS(T)- YGROSS(T)/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
def DAMEQ_rule(M, t):
    return (M.DAMAGES[t] == M.YGROSS[t]- M.YGROSS[t]/(1 + M.AA1*M.TATM[t]+ M.AA2*M.TATM[t]**M.AA3))
M.DAMEQ = Constraint(M.T, rule = DAMEQ_rule, doc = 'Damage equation')

#YNETEQ(T)..     YNET(T) =E=  YGROSS(T)/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
def YNETEQ_rule(M, t):
    return (M.YNET[t] ==  M.YGROSS[t]/(1 + M.AA1*M.TATM[t] + M.AA2*M.TATM[t]**M.AA3))
M.YNETEQ = Constraint(M.T, rule = YNETEQ_rule, doc = 'Output net of damages equation')

#ABATEEQ(T)..    ABATECOST(T) =E= (PARTFRACT(T)**(1-expcost2))*YGROSS(T)*(cost1(t)*(MIU(T)**EXPcost2));
def ABATEEQ_rule(M, t):
    return (M.ABATECOST[t] == (M.PARTFRACT[t]**(1 - M.EXPCOST2))*M.YGROSS[t]*(M.COST1[t]*(M.MIU[t]**M.EXPCOST2)))
M.ABATEEQ = Constraint(M.T, rule = ABATEEQ_rule, doc = 'Cost of emissions reductions equation')

#YY(T)..         Y(T) =E= YGROSS(T)*((1-(PARTFRACT(T)**(1-expcost2))*cost1(t)*(MIU(T)**EXPcost2)))/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
def YY_rule(M, t):
    return (M.Y[t] == M.YGROSS[t]*((1-(M.PARTFRACT[t]**(1-M.EXPCOST2))*M.COST1[t]*(M.MIU[t]**M.EXPCOST2)))/(1+M.AA1*M.TATM[t]+ M.AA2*M.TATM[t]**M.AA3))
M.YY = Constraint(M.T, rule = YY_rule, doc = 'Output net equation')

#SEQ(T)..        S(T)    =E= I(T)/(.001+Y(T));
def SEQ_rule(M, t):
    return (M.S[t] == M.I[t]/(0.001+M.Y[t]))
M.SEQ = Constraint(M.T, rule = SEQ_rule, doc = 'Savings rate equation')

#RIEQ(T)..       RI(T)   =E= GAMA*Y(T)/K(T)- (1-(1-DK)**10)/10  ;
def RIEQ_rule(M, t):
    return (M.RI[t] == M.GAMA*M.Y[t]/M.K[t]- (1-(1-M.DK)**10)/10)
M.RIEQ = Constraint(M.T, rule = RIEQ_rule, doc = 'Interest rate equation')

#CC(T)..         C(T)    =E= Y(T)-I(T);
def CC_rule(M, t):
    return (M.C[t] == M.Y[t] - M.I[t])
M.CC = Constraint(M.T, rule = CC_rule, doc = 'Consumption equation')

#CPCE(T)..       CPC(T)  =E= C(T)*1000/L(T);
def CPCE_rule(M, t):
    return (M.CPC[t] == M.C[t]*1000/M.L[t])
M.CPCE = Constraint(M.T, rule = CPCE_rule, doc = 'Per capita consumption definition')

#PCYE(T)..       PCY(T)  =E= Y(T)*1000/L(T);
def PCYE_rule(M, t):
    return (M.PCY[t] == M.Y[t]*1000/M.L[t])
M.PCYE = Constraint(M.T, rule = PCYE_rule, doc = 'Per capita income definition')

#PERIODUEQ(T)..  PERIODU(T)  =E=   ((C(T)/L(T))**(1-ELASMU)-1)/(1-ELASMU);
def PERIODUEQ_rule(M, t):
    return (M.PERIODU[t] == ((M.C[t]/M.L[t])**(1 - M.ELASMU) - 1)/(1 - M.ELASMU))
M.PERIODUEQ = Constraint(M.T, rule = PERIODUEQ_rule, doc = 'Instantaneous utility function equation')

#UTIL..          UTILITY =E= SUM(T, 10 *RR(T)*L(T)*(PERIODU(T))/scale1)+ scale2 ;
def UTIL_rule(M):
    return (sum(10*M.RR[t]*M.L[t]*M.PERIODU[t]/M.SCALE1 for t in M.T) + M.SCALE2)
M.UTIL = Objective(rule = UTIL_rule, sense = maximize, doc = 'Objective function')


#* Optimal run
#* Solution for optimal run

#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
##########################################################################################
##########################################################################################

##* Definition of opt results

##Parameters
#Year(t)         Date
#opt_y(t)
#opt_cpc(t)
#opt_s(t)
#opt_indem(t)
#opt_sigma(t)
#opt_tatm(t)
#opt_mat(t)
#opt_tax(t)
#opt_ri(t)
#opt_rr(t)
#opt_al(t)
#opt_forcoth(t)
#opt_l(t)
#opt_etree(t)
#opt_yy(t)
#opt_cc(t)
#opt_miu(t)
#opt_wem(t)
#opt_ri(t)
#opt_dam(t)
#opt_abate(t)
#opt_mcemis(t)
#opt_utility   ;

#Year(t)         = 2005 +10*(ord(t)-1);
#opt_y(t)=y.l(t);
#opt_cpc(t)=cpc.l(t);
#opt_s(t)=s.l(t)     ;
#opt_indem(t)= e.l(t)-etree(t);;
#opt_sigma(t)=sigma(t) ;
#opt_tatm(t)=tatm.l(t)  ;
#opt_mat(t)=mat.l(t)     ;
#opt_tax(t)=-1*ee.m(t)*1000/(kk.m(t)+.00000000001)       ;
#opt_ri(t)=ri.l(t);
#opt_rr(t)=rr(t)   ;
#opt_al(t)=al(t)    ;
#opt_forcoth(t)=forcoth(t);
#opt_l(t)=l(t);
#opt_etree(t)=etree(t);
#opt_yy(t)=yy.m(t)     ;
#opt_cc(t)=cc.m(t)      ;
#opt_miu(t)=miu.l(t)     ;
#opt_wem(t)= e.l(t);
#opt_ri(t)=ri.l(t)         ;
#opt_dam(t)= damages.l(t);
#opt_abate(t) = abatecost.l(t);
#opt_mcemis(t)= expcost2*cost1(t)*miu.l(t)**(expcost2-1)/sigma(t)*1000;
#opt_utility=utility.l        ;

##* Reset for initial conditions

#aa1 = a1;
#aa2 = a2;
#aa3 = a3;

#PBACK =  1.17 ;
#PARTFRACT1 = 1;
#PARTFRACT2 = 1;
#PARTFRACT21 = 1;
#partfract(t) = partfract21;
#cost1(T) = (PBACK*SIGMA(T)/EXPCOST2)* ( (BACKRAT-1+ EXP (-gback* (ORD(T)-1) ) )/BACKRAT);
#PARTFRACT(T)$(ord(T)<25) = Partfract21 + (PARTFRACT2-Partfract21)*exp(-DPARTFRACT*(ORD(T)-2));
#partfract("1")= PARTFRACT1;

#TATM.up(t)  = 10 ;
#mat.up(T)= 4000;
#miu.up(t)= 1;
#miu.lo(t)= .001;
#k.lo(t) = 1;
#k.up(t) = 1000000;
#K0 = 137;
#miu.fx("1")=.005;
#partfract("1")=        0.25372  ;

##* Estimate Hoteling rents
##* parameter estimates

#aa1 = 0;
#aa2 = 0;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;

##parameters
#miuhotel(t)    estimate of Hoteling rents;
#miuhotel(t)=miu.l(t);

##* Definition of hotelling results

##Parameters
#Year(t)         Date
#hotel_y(t)
#hotel_cpc(t)
#hotel_s(t)
#hotel_indem(t)
#hotel_sigma(t)
#hotel_tatm(t)
#hotel_mat(t)
#hotel_tax(t)
#hotel_ri(t)
#hotel_rr(t)
#hotel_al(t)
#hotel_forcoth(t)
#hotel_l(t)
#hotel_etree(t)
#hotel_yy(t)
#hotel_cc(t)
#hotel_miu(t)
#hotel_wem(t)
#hotel_ri(t)
#hotel_dam(t)
#hotel_abate(t)
#hotel_mcemis(t)
#hotel_utility   ;

#Year(t)         = 2005 +10*(ord(t)-1);
#hotel_y(t)=y.l(t);
#hotel_cpc(t)=cpc.l(t);
#hotel_s(t)=s.l(t)     ;
#hotel_indem(t)= e.l(t)-etree(t);;
#hotel_sigma(t)=sigma(t) ;
#hotel_tatm(t)=tatm.l(t)  ;
#hotel_mat(t)=mat.l(t)     ;
#hotel_tax(t)=-1*ee.m(t)*1000/(kk.m(t)+.0000001)       ;
#hotel_ri(t)=ri.l(t);
#hotel_rr(t)=rr(t)   ;
#hotel_al(t)=al(t)    ;
#hotel_forcoth(t)=forcoth(t);
#hotel_l(t)=l(t);
#hotel_etree(t)=etree(t);
#hotel_yy(t)=yy.m(t)     ;
#hotel_cc(t)=cc.m(t)      ;
#hotel_miu(t)=miu.l(t)     ;
#hotel_wem(t)= e.l(t);
#hotel_ri(t)=ri.l(t)         ;
#hotel_dam(t)= damages.l(t);
#hotel_abate(t) = abatecost.l(t);
#hotel_mcemis(t)= expcost2*cost1(t)*miu.l(t)**(expcost2-1)/sigma(t)*1000;
#hotel_utility=utility.l        ;
#* Reset for initial conditions

#aa1 = a1;
#aa2 = a2;
#aa3 = a3;

#PBACK =  1.17 ;
#PARTFRACT1 = 1;
#PARTFRACT2 = 1;
#PARTFRACT21 = 1;
#partfract(t) = partfract21;
#cost1(T) = (PBACK*SIGMA(T)/EXPCOST2)* ( (BACKRAT-1+ EXP (-gback* (ORD(T)-1) ) )/BACKRAT);
#PARTFRACT(T)$(ord(T)<25) = Partfract21 + (PARTFRACT2-Partfract21)*exp(-DPARTFRACT*(ORD(T)-2));
#partfract("1")= PARTFRACT1;

#TATM.up(t)  = 10 ;
#mat.up(T)= 4000;
#miu.up(t)= 1;
#miu.lo(t)= .001;
#k.lo(t) = 1;
#k.up(t) = 1000000;
#K0 = 137;
#miu.fx("1")=.005;
#partfract("1")=        0.25372  ;


##* Base-25per defined as 250 years of no action with miu at Hotelling control rates
##* Definition of base_250yr results

##* Control statements
#MIU.lo("1")=miuhotel("1");
#MIU.lo("2")=miuhotel("2");
#MIU.lo("3")=miuhotel("3");
#MIU.lo("4")=miuhotel("4");
#MIU.lo("5")=miuhotel("5");
#MIU.lo("6")=miuhotel("6");
#MIU.lo("7")=miuhotel("7");
#MIU.lo("8")=miuhotel("8");
#MIU.lo("9")=miuhotel("9");
#MIU.lo("10")=miuhotel("10");
#MIU.lo("11")=miuhotel("11");
#MIU.lo("12")=miuhotel("12");
#MIU.lo("13")=miuhotel("13");
#MIU.lo("14")=miuhotel("14");
#MIU.lo("15")=miuhotel("15");
#MIU.lo("16")=miuhotel("16");
#MIU.lo("17")=miuhotel("17");
#MIU.lo("18")=miuhotel("18");
#MIU.lo("19")=miuhotel("19");
#MIU.lo("20")=miuhotel("20");
#MIU.lo("21")=miuhotel("21");
#MIU.lo("22")=miuhotel("22");
#MIU.lo("23")=miuhotel("23");
#MIU.lo("24")=miuhotel("24");
#MIU.lo("25")=miuhotel("25");

#MIU.up("1")=miuhotel("1");
#MIU.up("2")=miuhotel("2");
#MIU.up("3")=miuhotel("3");
#MIU.up("4")=miuhotel("4");
#MIU.up("5")=miuhotel("5");
#MIU.up("6")=miuhotel("6");
#MIU.up("7")=miuhotel("7");
#MIU.up("8")=miuhotel("8");
#MIU.up("9")=miuhotel("9");
#MIU.up("10")=miuhotel("10");
#MIU.up("11")=miuhotel("11");
#MIU.up("12")=miuhotel("12");
#MIU.up("13")=miuhotel("13");
#MIU.up("14")=miuhotel("14");
#MIU.up("15")=miuhotel("15");
#MIU.up("16")=miuhotel("16");
#MIU.up("17")=miuhotel("17");
#MIU.up("18")=miuhotel("18");
#MIU.up("19")=miuhotel("19");
#MIU.up("20")=miuhotel("20");
#MIU.up("21")=miuhotel("21");
#MIU.up("22")=miuhotel("22");
#MIU.up("23")=miuhotel("23");
#MIU.up("24")=miuhotel("24");
#MIU.up("25")=miuhotel("25");

#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;
#solve CO2 maximizing UTILITY using nlp ;

##*Output
##Parameters
#Year(t)         Date
#base25_y(t)
#base25_cpc(t)
#base25_s(t)
#base25_indem(t)
#base25_sigma(t)
#base25_tatm(t)
#base25_mat(t)
#base25_tax(t)
#base25_ri(t)
#base25_rr(t)
#base25_al(t)
#base25_forcoth(t)
#base25_l(t)
#base25_etree(t)
#base25_yy(t)
#base25_cc(t)
#base25_miu(t)
#base25_wem(t)
#base25_ri(t)
#base25_dam(t)
#base25_abate(t)
#base25_mcemis(t)
#base25_mcemis(t)
#base25_utility
#base25_k(t) ;

#Year(t)         = 2005 +10*(ord(t)-1);
#base25_y(t)=y.l(t);
#base25_cpc(t)=cpc.l(t);
#base25_s(t)=s.l(t)     ;
#base25_indem(t)= e.l(t)-etree(t);;
#base25_sigma(t)=sigma(t) ;
#base25_tatm(t)=tatm.l(t)  ;
#base25_mat(t)=mat.l(t)     ;
#base25_tax(t)=-1*ee.m(t)*1000/(kk.m(t)+.00000000001)       ;
#base25_ri(t)=ri.l(t);
#base25_rr(t)=rr(t)   ;
#base25_al(t)=al(t)    ;
#base25_forcoth(t)=forcoth(t);
#base25_l(t)=l(t);
#base25_etree(t)=etree(t);
#base25_yy(t)=yy.m(t)     ;
#base25_cc(t)=cc.m(t)      ;
#base25_miu(t)=miu.l(t)     ;
#base25_wem(t)= e.l(t);
#base25_ri(t)=ri.l(t)         ;
#base25_dam(t)= damages.l(t);
#base25_abate(t) = abatecost.l(t);
#base25_mcemis(t)= expcost2*cost1(t)*miu.l(t)**(expcost2-1)/sigma(t)*1000;
#base25_utility=utility.l        ;
#base25_mcemis(t) = expcost2*cost1(t)*miu.l(t)**(expcost2-1)/sigma(t)*1000;
#base25_k(t) = k.l(t);



