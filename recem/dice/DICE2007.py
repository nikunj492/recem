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
#==============================================================================
# This is a translation of the GAMS version of  William Nordhaus's DICE:
# http://nordhaus.econ.yale.edu/DICE2007_short.gms
#
# See the DICE website
# http://nordhaus.econ.yale.edu/DICE2007.htm
#
# NOTE THE DIFFERENCE:
#   COMMENTED PROGRAM FROM  DICE2007_short.gms
# vs.
#COMMENTS FOR THIS PROGRAM
# or
# """
# COMMENTS OF THIS TYPE
# """
#==============================================================================


#    $ontext
#    DICE delta version 8
#    July 17, 2008.
#    This version is used for the DICE book, A Question of Balance (YUP, 2008).
#    We have included only the base, Hotelling, and optimal runs.
#    Exclude statements are removed so that it can run as a self-contained program.
#    Created September 5, 2008.
#
#    Note that this can be loaded into a data reading program,
#    $offtext

__all__ = ['createDICE2007']

from coopr.pyomo import *

"""
we define a an abstract optimization model using the AbstractModel
class of the Pyomo modeling language of the Coopr optimization tools
packages. See

https://software.sandia.gov/trac/coopr

The model is an abstract optimization model which is filled with data from
an appropriate data file (a *.dat file in the AMPL data format). This is then
translated by Pyomo to an instance file (another special data format file) and
sent to specialized solver programs like GLPK, CBC, IPOPT. The solution is
read by Pyomo and reported in human readable format.
"""
def createDICE2007(name='DICE2007', LIMMIU=1.0, FOSSLIM=6000.0):
    """
    This is a translation of the GAMS version of  William Nordhaus's DICE:
    http://nordhaus.econ.yale.edu/DICE2007_short.gms

    See the DICE website
    http://nordhaus.econ.yale.edu/DICE2007.htm

    To create an abstract optimization model of DICE2007  call:

    createDICE2007(name=NAME, LIMMIU=value, FOSSLIM=value)

    where NAME can be any string (used a short description of the model), the
    LIMMIU value is the upper limit on the fraction of CO2 emission under abatement,
    and the FOSSLIM value is the maximum cumulative extraction of fossil fuels (in GtC).
    The default values are:

    name='DICE2007 Optimal'
    LIMMIU=1
    FOSSLIM=6000

    This version of DICE2007 calculates the Optimum trajectory. Compare with
    the book: A Question of Balance: Economic Modeling of Global Warming on the
    DICE website.
    """

    if 0 <= LIMMIU <= 1:
        pass
    else:
        raise ValueError('LIMMIU should be between 0 and 1! Default is 1.')
    if FOSSLIM < 0:
        raise ValueError('FOSSLIM should be positive! Default value is 6000.')

    M = AbstractModel()
    M.name = name


    #   SETS  T                 Time periods                     /1*60/ ;

    M.L0 = Param(initialize = 60, doc = 'Number of Time Periods')
    M.L1 = Param(initialize = 12, doc = 'Time Period at which some parameters saturate')
    M.L2 = Param(initialize = 25, doc = 'Time Period at which some parameters saturate')

    #Equivalent to SET T of DICE
    M.T = RangeSet(1, M.L0, doc = 'Time Periods')


#   SCALARS

#    ** Preferences
#     B_ELASMU   Elasticity of marginal utility of consumption     /  2.0    /
#     B_PRSTP    Initial rate of social time preference per year   / .015    /

    M.B_ELASMU = Param(initialize = 2.0, doc = 'Elasticity of marginal utility of consumption')
    M.B_PRSTP = Param(initialize = 0.015, doc = 'Initial rate of social time preference per year')


#    ** Population and technology
#     POP0     2005 world population millions                  /6514     /
#     GPOP0    Growth rate of population per decade            /.35      /
#     POPASYM  Asymptotic population                           / 8600    /
#     A0       Initial level of total factor productivity      /.02722   /
#     GA0      Initial growth rate for technology per decade   /.092      /
#     DELA     Decline rate of technol change per decade       /.001     /
#     DK       Depreciation rate on capital per year           /.100     /
#     GAMA     Capital elasticity in production function       /.300     /
#     Q0       2005 world gross output trill 2005 US dollars   /61.1     /
#     K0       2005 value capital trill 2005 US dollars        /137.     /

    M.POP0 = Param(initialize = 6514.0, doc = '2005 world population millions')
    M.GPOP0 = Param(initialize = 0.35, doc = 'Growth rate of population per decade')
    M.POPASYM = Param(initialize = 8600.0, doc = 'Asymptotic population')
    M.A0 = Param(initialize = 0.02722, doc = 'Initial level of total factor productivity')
    M.GA0 = Param(initialize = 0.092, doc = 'Initial growth rate for technology per decade')
    M.DELA = Param(initialize = 0.001, doc = 'Decline rate of technol change per decade')
    M.DK = Param(initialize = 0.1, doc = 'Depreciation rate on capital per year')
    M.GAMA = Param(initialize = 0.3, doc = 'Capital elasticity in production function')
    M.Q0 = Param(initialize = 61.1, doc = '2005 world gross output trill 2005 US dollars')
    M.K0 = Param(initialize = 137.0, doc = '2005 value capital trill 2005 US dollars')


#    ** Emissions
#     SIG0     CO2-equivalent emissions-GNP ratio 2005         /.13418    /
#     GSIGMA   Initial growth of sigma per decade              /-.0730    /
#     DSIG     Decline rate of decarbonization per decade      /.003   /
#     DSIG2    Quadratic term in decarbonization               / .000   /
#     ELAND0   Carbon emissions from land 2005(GtC per decade) / 11.000  /

    M.SIG0 = Param(initialize = 0.13418, doc = 'CO2-equivalent emissions-GNP ratio 2005')
    M.GSIGMA = Param(initialize = -0.073, doc = 'Initial growth of sigma per decade')
    M.DSIG = Param(initialize = 0.003, doc = 'Decline rate of decarbonization per decade')
    M.DSIG2 = Param(initialize = 0.000, doc = 'Quadratic term in decarbonization')
    M.ELAND0 = Param(initialize = 11.000, doc = 'Carbon emissions from land 2005(GtC per decade')


#    ** Carbon cycle
#     MAT2000  Concentration in atmosphere 2005 (GtC)          /808.9   /
#     MU2000   Concentration in upper strata 2005 (GtC)        /1255     /
#     ML2000   Concentration in lower strata 2005 (GtC)        /18365    /
#     b11      Carbon cycle transition matrix                  /0.810712 /
#     b12      Carbon cycle transition matrix                  /0.189288 /
#     b21      Carbon cycle transition matrix                  /0.097213 /
#     b22      Carbon cycle transition matrix                  /0.852787 /
#     b23      Carbon cycle transition matrix                  /0.05     /
#     b32      Carbon cycle transition matrix                  /0.003119 /
#     b33      Carbon cycle transition matrix                  /0.996881 /

    M.MAT2000 = Param(initialize = 808.9, doc = 'Concentration in atmosphere 2005 (GtC)')
    M.MU2000 = Param(initialize = 1255.0, doc = 'Concentration in upper strata 2005 (GtC)')
    M.ML2000 = Param(initialize = 18365.0, doc = 'Concentration in lower strata 2005 (GtC)')
    M.B11 = Param(initialize = 0.810712, doc = 'Carbon cycle transition matrix')
    M.B12 = Param(initialize = 0.189288, doc = 'Carbon cycle transition matrix')
    M.B21 = Param(initialize = 0.097213, doc = 'Carbon cycle transition matrix')
    M.B22 = Param(initialize = 0.852787, doc = 'Carbon cycle transition matrix')
    M.B23 = Param(initialize = 0.05, doc = 'Carbon cycle transition matrix')
    M.B32 = Param(initialize = 0.003119, doc = 'Carbon cycle transition matrix')
    M.B33 = Param(initialize = 0.996881, doc = 'Carbon cycle transition matrix')


#    ** Climate model
#     T2XCO2   Equilibrium temp impact of CO2 doubling oC      / 3 /
#     FEX0     Estimate of 2000 forcings of non-CO2 GHG        / -.06   /
#     FEX1     Estimate of 2100 forcings of non-CO2 GHG        / 0.30   /
#     TOCEAN0  2000 lower strat. temp change (C) from 1900     /.0068   /
#     TATM0    2000 atmospheric temp change (C)from 1900       /.7307   /
#     C1       Climate-equation coefficient for upper level    /.220    /
#     C3       Transfer coeffic upper to lower stratum         /.300    /
#     C4       Transfer coeffic for lower level                /.050    /
#     FCO22X   Estimated forcings of equilibrium co2 doubling  /3.8     /

    M.T2XCO2 = Param(initialize = 3.0, doc = 'Equilibrium temp impact of CO2 doubling oC')
    M.FEX0 = Param(initialize = -0.06, doc = 'Estimate of 2000 forcings of non-CO2 GHG')
    M.FEX1 = Param(initialize = 0.3, doc = 'Estimate of 2100 forcings of non-CO2 GHG')
    M.TOCEAN0 = Param(initialize = 0.0068, doc = '2000 lower strat. temp change (C) from 1900')
    M.TATM0 = Param(initialize = 0.7307, doc = '2000 atmospheric temp change (C)from 1900')
    M.C1 = Param(initialize = 0.22, doc = 'Climate-equation coefficient for upper level')
    M.C3 = Param(initialize = 0.3, doc = 'Transfer coeffic upper to lower stratum')
    M.C4 = Param(initialize = 0.05, doc = 'Transfer coeffic for lower level')
    M.FCO22X = Param(initialize = 3.8, doc = 'Estimated forcings of equilibrium CO2 doubling')


#    ** Climate damage parameters calibrated for quadratic at 2.5 C for 2105
#     A1       Damage intercept                                / 0.00000    /
#    A2       Damage quadratic term                           /  0.0028388 /
#     A3       Damage exponent                                 / 2.00       /

    M.A1 = Param(initialize = 0.0, doc = 'Damage intercept')
    M.A2 = Param(initialize = 0.0028388, doc = 'Damage quadratic term')
    M.A3 = Param(initialize = 2.0, doc = 'Damage exponent')


#    ** Abatement cost
#     EXPCOST2   Exponent of control cost function               /2.8   /
#     PBACK      Cost of backstop 2005 000$ per tC 2005          /1.17  /
#     BACKRAT    Ratio initial to final backstop cost            / 2    /
#     GBACK      Initial cost decline backstop pc per decade     /.05   /
#     LIMMIU     Upper limit on control rate                     / 1    /

    M.EXPCOST2 = Param(initialize = 2.8, doc = 'Exponent of control cost function')
    M.PBACK = Param(initialize = 1.17, doc = 'Cost of backstop 2005 000$ per tC 2005')
    M.BACKRAT = Param(initialize = 2.0, doc = 'Ratio initial to final backstop cost')
    M.GBACK = Param(initialize = 0.05, doc = 'Initial cost decline backstop pc per decade')
    M.LIMMIU = LIMMIU # ANOTHER PYTHON PYOMO HACK (WE READ/INITIALIZE THIS PARAM WHEN
#                       WE CALL THIS FUNCTION TO CREATE A DICE2007 MODEL)


#    ** Participation
#     PARTFRACT1  Fraction of emissions under control regime 2005 /1      /
#     PARTFRACT2  Fraction of emissions under control regime 2015 /1      /
#     PARTFRACT21 Fraction of emissions under control regime 2205 /1      /
#     DPARTFRACT  Decline rate of participation                   /0      /
#
#    ** Availability of fossil fuels
#     FOSSLIM  Maximum cumulative extraction fossil fuels         / 6000  /
#
#    ** Scaling and inessential parameters
#      scale1 Scaling coefficient in the objective function       /194    /
#      scale2 Scaling coefficient in the objective function       /381800 / ;

    M.PARTFRACT1 = Param(initialize = 1.0, doc = 'Fraction of emissions under control regime 2005')
    M.PARTFRACT2 = Param(initialize = 1.0, doc = 'Fraction of emissions under control regime 2015')
    M.PARTFRACT21 = Param(initialize = 1.0, doc = 'Fraction of emissions under control regime 2205')
    M.DPARTFRACT = Param(initialize = 0.0, doc = 'Decline rate of participation')

    M.FOSSLIM = FOSSLIM # ANOTHER PYTHON PYOMO HACK

    M.SCALE1 = Param(initialize = 194, doc = 'Scaling coefficient in the objective function')
    M.SCALE2 = Param(initialize = 381800, doc = 'Scaling coefficient in the objective function')

#    * Definitions for outputs of no economic interest
#    SETS
#          TFIRST(T)
#          TLAST(T)
#          TEARLY(T)
#          TLATE(T);

#NOTE: We don't use the above SETS in the OPTIMUM CASE.


#      AA1           Variable A1
#      AA2           Variable A2
#      AA3           Variable A3
#      ELASMU        Variable elasticity of marginal utility of consumption
#      PRSTP         Variable nitial rate of social time preference per year
#      LAM           Climate model parameter
#    * Unimportant definitions to reset runs
#    TFIRST(T) = YES$(ORD(T) EQ 1);
#    TLAST(T)  = YES$(ORD(T) EQ CARD(T));
#    TEARLY(T) = YES$(ORD(T) LE 20);
#    TLATE(T)  = YES$(ORD(T) GE 21);
#    AA1 = A1;
#    AA2 = A2;
#    AA3 = A3;
#    ELASMU = B_ELASMU;
#    PRSTP  = B_PRSTP;

    M.AA1 = Param(rule = lambda M: value(M.A1), doc = 'Variable A1')
    M.AA2 = Param(rule = lambda M: value(M.A2), doc = 'Variable A2')
    M.AA3 = Param(rule = lambda M: value(M.A3), doc = 'Variable A3')
    M.ELASMU = Param(rule = lambda M: value(M.B_ELASMU), doc = 'Variable elasticity of marginal utility of consumption')
    M.PRSTP = Param(rule = lambda M: value(M.B_PRSTP), doc = 'Variable initial rate of social time preference per year')




#    M.B11 = 1 - M.B12;
#    M.B21 = 587.473*M.B12/1143.894;
#    M.B22 = 1 - M.B21 - M.B23;
#    M.B32 = 1143.894*M.B23/18340;
#    M.B33 = 1 - M.B32 ;

    M.B11 = Param(rule = lambda M: 1 - value(M.B12), doc = 'Carbon cycle transition matrix')
    M.B21 = Param(rule = lambda M: 587.473*value(M.B12)/1143.894, doc = 'Carbon cycle transition matrix')
    M.B22 = Param(rule = lambda M: 1 - value(M.B21) - value(M.B23), doc = 'Carbon cycle transition matrix')
    M.B32 = Param(rule = lambda M: 1143.894*value(M.B23)/18340, doc = 'Carbon cycle transition matrix')
    M.B33 = Param(rule = lambda M: 1 - value(M.B32), doc = 'Carbon cycle transition matrix')



#    PARAMETERS
#      L(T)          Level of population and labor
#      AL(T)         Level of total factor productivity
#      SIGMA(T)      CO2-equivalent-emissions output ratio
#      R(T)          Instantaeous rate of social time preference
#      RR(T)         Average utility social discount rate
#      GA(T)         Growth rate of productivity from 0 to T
#      FORCOTH(T)    Exogenous forcing for other greenhouse gases
#      GL(T)         Growth rate of labor 0 to T
#      GCOST1        Growth of cost factor
#      GSIG(T)       Cumulative improvement of energy efficiency
#      ETREE(T)      Emissions from deforestation
#      COST1(t)      Adjusted cost for backstop
#      PARTFRACT(T)  Fraction of emissions in control regime
#      AA1           Variable A1
#      AA2           Variable A2
#      AA3           Variable A3
#      ELASMU        Variable elasticity of marginal utility of consumption
#      PRSTP         Variable nitial rate of social time preference per year
#      LAM           Climate model parameter
#      Gfacpop(T)    Growth factor population ;

#   * Important parameters for the model

#   LAM     = FCO22X/ T2XCO2;
    M.LAM = Param(rule = lambda M: value(M.FCO22X)/value(M.T2XCO2), doc = 'Climate model parameter')

#   GFACPOP(T) =   (exp(GPOP0*(ORD(T)-1))-1)/exp(GPOP0*(ORD(T)-1));
    def GFACPOP_rule(M, t):
        paramvalue = (exp(M.GPOP0*(t-1)) -1)/exp(M.GPOP0*(t-1))
        return value(paramvalue)
    M.GFACPOP = Param(M.T, rule = GFACPOP_rule, doc = 'Growth factor population')

#   L(T)=POP0* (1- GFACPOP(T))+GFACPOP(T)*POPASYM;
    def L_rule(M, t):
        paramvalue = M.POP0 * (1 - M.GFACPOP[t]) + M.GFACPOP[t] * M.POPASYM
        return value(paramvalue)
    M.L = Param(M.T, rule = L_rule, doc = 'Level of population and labor')

#   GA(T)=GA0*EXP(-DELA*10*(ORD(T)-1));
    def GA_rule(M, t):
        paramvalue = M.GA0 * exp(- M.DELA * 10 *(t - 1))
        return value(paramvalue)
    M.GA = Param(M.T, rule = GA_rule, doc = 'Growth rate of productivity from 0 to T')

#   AL("1") = A0;
#   LOOP(T, AL(T+1)=AL(T)/((1-GA(T))););
    def AL_rule(M, t):
        if t == 1:
            return value(M.A0)
        else:
            return value(M.AL[t-1])/(1 - value(M.GA[t-1]))
    M.AL = Param(M.T, rule = AL_rule, doc = 'Level of total factor productivity')


#   GSIG(T)=GSIGMA*EXP(-DSIG*10*(ORD(T)-1)-DSIG2*10*((ORD(T)-1)**2));
    def GSIG_rule(M, t):
        paramvalue = M.GSIGMA*exp(-M.DSIG*10*(t-1) - M.DSIG2*10*(t-1)**2)
        return value(paramvalue)
    M.GSIG = Param(M.T, rule = GSIG_rule, doc = 'Cumulative improvement of energy efficiency')

#   SIGMA("1")=SIG0;
#   LOOP(T,SIGMA(T+1)=(SIGMA(T)/((1-GSIG(T+1)))););
    def SIGMA_rule(M,t):
        if t == 1:
            return value(M.SIG0)
        else:
            return value(M.SIGMA[t-1])/(1 - value(M.GSIG[t]))
    M.SIGMA = Param(M.T, rule = SIGMA_rule, doc = 'CO2-equivalent-emissions output ratio')


#   COST1(T) = (PBACK*SIGMA(T)/EXPCOST2)* ( (BACKRAT-1+ EXP (-GBACK* (ORD(T)-1) ) )/BACKRAT);
    def COST1_rule(M, t):
        paramvalue = (M.PBACK*M.SIGMA[t]/M.EXPCOST2)*((M.BACKRAT - 1 + exp(-M.GBACK*(t-1)))/M.BACKRAT)
        return value(paramvalue)
    M.COST1 = Param(M.T, rule = COST1_rule, doc = 'Adjusted cost for backstop')

#   ETREE(T) = ELAND0*(1-0.1)**(ORD(T)-1);
    def ETREE_rule(M, t):
        return value(M.ELAND0)*(1 - 0.1)**(t-1)
    M.ETREE = Param(M.T, rule = ETREE_rule, doc = 'Emissions from deforestation')

#   RR(T)=1/((1+PRSTP)**(10*(ORD(T)-1)));
    def RR_rule(M, t):
        return (1 + value(M.PRSTP))**(10*(1 - t))
    M.RR = Param(M.T, rule = RR_rule, doc = 'Average utility social discount rate')

#   FORCOTH(T)= FEX0+ .1*(FEX1-FEX0)*(ORD(T)-1)$(ORD(T) LT 12)+ 0.36$(ORD(T) GE 12);
    def FORCOTH_rule(M, t):
        return (t < value(M.L1)  and (value(M.FEX0) + 0.1*(value(M.FEX1) - value(M.FEX0))*(t-1)) or (value(M.FEX0) + 0.36))
    M.FORCOTH = Param(M.T, rule = FORCOTH_rule, doc = 'Exogenous forcing for other greenhouse gases')

#   PARTFRACT(T) = PARTFRACT21;
#   PARTFRACT(T)$(ORD(T)<25) = PARTFRACT21 + (PARTFRACT2-PARTFRACT21)*EXP(-DPARTFRACT*(ORD(T)-2));
#   PARTFRACT("1")= PARTFRACT1;
#   PARTFRACT("1")= 0.25372;
    def PARTFRACT_rule(M, t):
        if t==1:
#            return value(M.PARTFRACT1)
            return 0.25372
        elif 1 < t < value(M.L2):
#            return value(M.PARTFRACT21) + (value(M.PARTFRACT2) - value(M.PARTFRACT21))*exp(-value(M.DPARTFRACT)*(t-2))
            paramvalue = M.PARTFRACT21 + (M.PARTFRACT2 - M.PARTFRACT21)*exp(-M.DPARTFRACT*(t-2))
            return value(paramvalue)
        else:
            return value(M.PARTFRACT21)
    M.PARTFRACT = Param(M.T, rule = PARTFRACT_rule, doc = 'Fraction of emissions in control regime')


#    VARIABLES
#     MIU(T)          Emission control rate GHGs
#     FORC(T)         Radiative forcing in watts per m2
#     TATM(T)         Temperature of atmosphere in degrees C
#     TOCEAN(T)       Temperatureof lower oceans degrees C
#     MAT(T)          Carbon concentration in atmosphere GtC
#     MATAV(T)        Average concentrations
#     MU(T)           Carbon concentration in shallow oceans Gtc
#     ML(T)           Carbon concentration in lower oceans GtC
#     E(T)            CO2-equivalent emissions GtC
#     C(T)            Consumption trillions US dollars
#     K(T)            Capital stock trillions US dollars
#     CPC(T)          Per capita consumption thousands US dollars
#     PCY(t)          Per capita income thousands US dollars
#     I(T)            Investment trillions US dollars
#     S(T)            Gross savings rate as fraction of gross world product
#     RI(T)           Real interest rate per annum
#     Y(T)            Gross world product net of abatement and damages
#     YGROSS(T)       Gross world product GROSS of abatement and damages
#     YNET(T)         Output net of damages equation
#     DAMAGES(T)      Damages
#     ABATECOST(T)    Cost of emissions reductions
#     CCA(T)          Cumulative industrial carbon emissions GTC
#     PERIODU(t)      One period utility function
#     UTILITY;
#
#    POSITIVE VARIABLES MIU, TATM, TOCE, E, MAT, MATAV, MU, ML, Y, YGROSS, C, K, I, CCA ;

    M.MIU = Var(M.T, within = NonNegativeReals, bounds = (0.0, M.LIMMIU), doc = 'Emission control rate GHGs')
    M.FORC = Var(M.T, doc = 'Radiative forcing in watts per m2')
    M.TATM = Var(M.T, initialize = 0, within = NonNegativeReals, doc = 'Temperature of atmosphere in degrees C')
    M.TOCEAN = Var(M.T, initialize = 0, within = NonNegativeReals, doc = 'Temperature of lower oceans degrees C')
    M.MAT = Var(M.T, within = NonNegativeReals, doc = 'Carbon concentration in atmosphere GtC')
    M.MATAV = Var(M.T, within = NonNegativeReals, doc = 'Average concentrations')
    M.MU = Var(M.T, within = NonNegativeReals, doc = 'Carbon concentration in shallow oceans Gtc')
    M.ML = Var(M.T, within = NonNegativeReals, doc = 'Carbon concentration in lower oceans GtC')
    M.E = Var(M.T, within = NonNegativeReals, doc = 'CO2-equivalent emissions GtC')
    M.C = Var(M.T, initialize = 20, within = NonNegativeReals, doc = 'Consumption trillions US dollars')
    M.K = Var(M.T, initialize = 100, within = NonNegativeReals, doc = 'Capital stock trillions US dollars')
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



#   **  UPPER AND LOWER BOUNDS: GENERAL CONDITIONS FOR STABILITY

#   K.LO(T)         = 100;
    M.K.bounds = (100, None)
#   MAT.LO(T)       = 10;
    M.MAT.bounds = (10, None)
#   MU.LO(T)        = 100;
    M.MU.bounds = (100, None)
#   ML.LO(T)        = 1000;
    M.ML.bounds = (1000, None)
#   C.LO(T)         = 20;
    M.C.bounds = (20, None)
#   TOCEAN.UP(T)    = 20;
#   TOCEAN.LO(T)    = -1;
    M.TOCEAN.bounds = (-1, 20)
#   TATM.UP(T)      = 20;
    M.TATM.bounds = (0, 20)
#   MIU.UP(T)       = LIMMIU;
    M.MIU.bounds = (0, M.LIMMIU)
#   PARTFRACT("1")= 0.25372;
#See above in PARTFRACT Param definition

#   * First period predetermined by Kyoto Protocol
#   MIU.FX("1")     = 0.005;
    def FIXMIU_rule(M, t):
        if t==1:
            return (0.005, M.MIU[t], 0.005)
        else:
            return Constraint.Skip
    M.FIXMIU = Constraint(M.T, rule = FIXMIU_rule)

#   ** Fix savings assumption for standardization if needed
#   S.FX(T)=.22;
    M.S.bounds = (0.22,0.22)

#   ** Cumulative limits on carbon use at 6000 GtC
#   CCA.UP(T) = FOSSLIM;
    M.CCA.bounds = (None, M.FOSSLIM)



#   EQUATIONS


#   ** Equations of the model

#   CCTFIRST(TFIRST).. CCA(TFIRST)=E=0;
#   CCACCA(T+1)..      CCA(T+1)=E=CCA(T)+ E(T);
    def CCACCA_rule(M, t):
        if t == 1:
            return (M.CCA[t] == 0)
        else:
            return (M.CCA[t] == M.CCA[t-1] + M.E[t-1])
#   M.CCTFIRST = Constraint(M.T, rule = CCTFIRST_rule, doc = 'First period cumulative carbon')
    M.CCACCA = Constraint(M.T, rule = CCACCA_rule, doc = 'Cumulative carbon emissions')

#   KK0(TFIRST)..      K(TFIRST) =E= K0;
#   KK(T)..            K(T+1) =L= (1-DK)**10 *K(T)+10*I(T);
#   KC(TLAST)..        .02*K(TLAST) =L= I(TLAST);
#   M.KK0 = Constraint(M.T, rule = KK0_rule, doc = 'Initial condition for capital')
#   M.KC = Constraint(M.T, rule = KC_rule, doc = 'Terminal condition for capital')
    def KK_rule(M,t):
        if t == 1:
            return (M.K[t] == M.K0)
        else:
            return (M.K[t] <= (1 - M.DK)**10 *M.K[t-1] + 10*M.I[t-1])
    M.KK = Constraint(M.T, rule = KK_rule, doc = 'Capital balance equation with special initial conditions')

    def KC_rule(M,t):
        if t == value(M.L0):
            return (0.02*M.K[t] <= M.I[t])
        else:
            return Constraint.Skip
    M.KC = Constraint(M.T, rule = KC_rule, doc = 'Capital balance terminal conditions')

#   EE(T)..            E(T)=E=10*SIGMA(T)*(1-MIU(T))*AL(T)*L(T)**(1-GAMA)*K(T)**GAMA + ETREE(T);
    def EE_rule(M, t):
        return (M.E[t] == 10*M.SIGMA[t]*(1 - M.MIU[t])*M.AL[t]*M.L[t]**(1-M.GAMA)*M.K[t]**M.GAMA + M.ETREE[t])
    M.EE = Constraint(M.T, rule = EE_rule, doc = 'Emissions equation')

#   FORCE(T)..         FORC(T) =E=  FCO22X*((log((Matav(T)+.000001)/596.4)/log(2)))+FORCOTH(T);
    def FORCE_rule(M, t):
        return (M.FORC[t] == M.FCO22X*((log((M.MATAV[t] + 0.000001)/596.4)/log(2.0))) + M.FORCOTH[t])
    M.FORCE = Constraint(M.T, rule = FORCE_rule, doc = 'Radiative forcing equation')

#   MMAT0(TFIRST)..    MAT(TFIRST) =E= MAT2000;
#   MMAT(T+1)..        MAT(T+1)    =E= MAT(T)*b11+MU(T)*b21 + E(T);
#   M.MMAT0 = Constraint(M.T, rule = MMAT0_rule, doc = 'Starting atmospheric concentration')
    def MMAT_rule(M, t):
        if t==1:
            return (M.MAT[t] == M.MAT2000)
        else:
            return (M.MAT[t] == M.MAT[t-1]*M.B11 + M.MU[t-1]*M.B21 + M.E[t-1])
    M.MMAT = Constraint(M.T, rule = MMAT_rule, doc = 'Atmospheric concentration equation')

#   MMU0(TFIRST)..     MU(TFIRST)  =E= MU2000;
#   MMU(T+1)..         MU(T+1)     =E= MAT(T)*b12+MU(T)*b22+ML(T)*b32;
#   M.MMU0 = Constraint(M.T, rule = MMU0_rule, doc = 'Initial shallow ocean concentration')
    def MMU_rule(M, t):
        if t==1:
            return (M.MU[t] == M.MU2000)
        else:
            return (M.MU[t] == M.MAT[t-1]*M.B12 + M.MU[t-1]*M.B22 + M.ML[t-1]*M.B32)
    M.MMU = Constraint(M.T, rule = MMU_rule, doc = 'Shallow ocean concentration')

#   MML0(TFIRST)..     ML(TFIRST)  =E= ML2000;
#   MML(T+1)..         ML(T+1)     =E= ML(T)*b33+b23*MU(T);
#   M.MML0 = Constraint(M.T, rule = MML0_rule, doc = 'Initial lower ocean concentration')
    def MML_rule(M, t):
        if t==1:
            return (M.ML[t] == M.ML2000)
        else:
            return (M.ML[t] == M.ML[t-1]*M.B33 + M.MU[t-1]*M.B23)
    M.MML = Constraint(M.T, rule = MML_rule, doc = 'Lower ocean concentration')

#   MMATAVEQ(t)..      MATAV(T)    =E= (MAT(T)+MAT(T+1))/2
    def MMATAVEQ_rule(M, t):
        if t < M.L0:
            return (M.MATAV[t] == 0.5*(M.MAT[t] + M.MAT[t+1]))
        else:
            return (M.MATAV[t] == M.MAT[t])
    M.MMATAVEQ = Constraint(M.T, rule = MMATAVEQ_rule, doc = 'Average concentrations equation')


#   TATM0EQ(TFIRST)..  TATM(TFIRST) =E= TATM0;
#   TATMEQ(T+1)..      TATM(T+1) =E= TATM(t)+C1*(FORC(t+1)-LAM*TATM(t)-C3*(TATM(t)-TOCEAN(t)));
#   M.TATM0EQ = Constraint(M.T, rule = TATM0EQ_rule, doc = 'Initial condition for atmospheric temperature')
    def TATMEQ_rule(M, t):
        if t == 1:
            return (M.TATM[t] == M.TATM0)
        else:
            return (M.TATM[t] == M.TATM[t-1] + M.C1*(M.FORC[t] - M.LAM*M.TATM[t-1] - M.C3*(M.TATM[t-1] - M.TOCEAN[t-1])))
    M.TATMEQ = Constraint(M.T, rule = TATMEQ_rule, doc = 'Temperature-climate equation for atmosphere')

#   TOCEAN0EQ(TFIRST)..  TOCEAN(TFIRST) =E= TOCEAN0;
#   TOCEANEQ(T+1)..    TOCEAN(T+1) =E= TOCEAN(T)+C4*(TATM(T)-TOCEAN(T));
#   M.TOCEAN0EQ = Constraint(M.T, rule = TOCEAN0EQ_rule, doc = 'Initial condition for lower ocean temperature')
    def TOCEANEQ_rule(M, t):
        if t == 1:
            return (M.TOCEAN[t] == M.TOCEAN0)
        else:
            return (M.TOCEAN[t] == M.TOCEAN[t-1] + M.C4*(M.TATM[t-1]-M.TOCEAN[t-1]))
    M.TOCEANEQ = Constraint(M.T, rule = TOCEANEQ_rule, doc = 'Temperature-climate equation for lower oceans')

#   YGROSSEQ(T)..   YGROSS(T) =e= AL(T)*L(T)**(1-GAMA)*K(T)**GAMA;
    def YGROSSEQ_rule(M, t):
        return (M.YGROSS[t] == M.AL[t]*M.L[t]**(1-M.GAMA)*M.K[t]**M.GAMA)
    M.YGROSSEQ = Constraint(M.T, rule = YGROSSEQ_rule, doc = 'Output gross equation')

#   DAMEQ(T)..      DAMAGES(t) =E= YGROSS(T)- YGROSS(T)/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
    def DAMEQ_rule(M, t):
        return (M.DAMAGES[t] == M.YGROSS[t]- M.YGROSS[t]/(1 + M.AA1*M.TATM[t]+ M.AA2*M.TATM[t]**M.AA3))
    M.DAMEQ = Constraint(M.T, rule = DAMEQ_rule, doc = 'Damage equation')

#   YNETEQ(T)..     YNET(T) =E=  YGROSS(T)/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
    def YNETEQ_rule(M, t):
        return (M.YNET[t] ==  M.YGROSS[t]/(1 + M.AA1*M.TATM[t] + M.AA2*M.TATM[t]**M.AA3))
    M.YNETEQ = Constraint(M.T, rule = YNETEQ_rule, doc = 'Output net of damages equation')

#   ABATEEQ(T)..    ABATECOST(T) =E= (PARTFRACT(T)**(1-expcost2))*YGROSS(T)*(cost1(t)*(MIU(T)**EXPcost2));
    def ABATEEQ_rule(M, t):
        return (M.ABATECOST[t] == (M.PARTFRACT[t]**(1 - M.EXPCOST2))*M.YGROSS[t]*(M.COST1[t]*(M.MIU[t]**M.EXPCOST2)))
    M.ABATEEQ = Constraint(M.T, rule = ABATEEQ_rule, doc = 'Cost of emissions reductions equation')

#   YY(T)..         Y(T) =E= YGROSS(T)*((1-(PARTFRACT(T)**(1-expcost2))*cost1(t)*(MIU(T)**EXPcost2)))/(1+aa1*TATM(T)+ aa2*TATM(T)**aa3);
    def YY_rule(M, t):
        return (M.Y[t] == M.YGROSS[t]*((1-(M.PARTFRACT[t]**(1-M.EXPCOST2))*M.COST1[t]*(M.MIU[t]**M.EXPCOST2)))/(1+M.AA1*M.TATM[t]+ M.AA2*M.TATM[t]**M.AA3))
    M.YY = Constraint(M.T, rule = YY_rule, doc = 'Output net equation')

#   SEQ(T)..        S(T)    =E= I(T)/(.001+Y(T));
    def SEQ_rule(M, t):
        return (M.S[t] == M.I[t]/(0.001+M.Y[t]))
    M.SEQ = Constraint(M.T, rule = SEQ_rule, doc = 'Savings rate equation')

#   RIEQ(T)..       RI(T)   =E= GAMA*Y(T)/K(T)- (1-(1-DK)**10)/10  ;
    def RIEQ_rule(M, t):
        return (M.RI[t] == M.GAMA*M.Y[t]/M.K[t]- (1-(1-M.DK)**10)/10)
    M.RIEQ = Constraint(M.T, rule = RIEQ_rule, doc = 'Interest rate equation')

#   CC(T)..         C(T)    =E= Y(T)-I(T);
    def CC_rule(M, t):
        return (M.C[t] == M.Y[t] - M.I[t])
    M.CC = Constraint(M.T, rule = CC_rule, doc = 'Consumption equation')

#   CPCE(T)..       CPC(T)  =E= C(T)*1000/L(T);
    def CPCE_rule(M, t):
        return (M.CPC[t] == M.C[t]*1000/M.L[t])
    M.CPCE = Constraint(M.T, rule = CPCE_rule, doc = 'Per capita consumption definition')

#   PCYE(T)..       PCY(T)  =E= Y(T)*1000/L(T);
    def PCYE_rule(M, t):
        return (M.PCY[t] == M.Y[t]*1000/M.L[t])
    M.PCYE = Constraint(M.T, rule = PCYE_rule, doc = 'Per capita income definition')

#   PERIODUEQ(T)..  PERIODU(T)  =E=   ((C(T)/L(T))**(1-ELASMU)-1)/(1-ELASMU);
    def PERIODUEQ_rule(M, t):
        return (M.PERIODU[t] == ((M.C[t]/M.L[t])**(1 - M.ELASMU) - 1)/(1 - M.ELASMU))
    M.PERIODUEQ = Constraint(M.T, rule = PERIODUEQ_rule, doc = 'Instantaneous utility function equation')

#   UTIL..          UTILITY =E= SUM(T, 10 *RR(T)*L(T)*(PERIODU(T))/scale1)+ scale2 ;
    def UTIL_rule(M):
        return (sum(10*M.RR[t]*M.L[t]*M.PERIODU[t]/M.SCALE1 for t in M.T) + M.SCALE2)
    M.UTIL = Objective(rule = UTIL_rule, sense = maximize, doc = 'Objective function: total discounted utility')

    return M
