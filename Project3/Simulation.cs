using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Oblig2
{
    class Simulation
    {
        // file name for opacity matrix
        public const string OPM_FILE = @"opacity.txt";

        // constants used in calculations
        //      1) physical constants
        public const double G = 0.00000000006673; // Gravitational constant: m^3 / (kg*s^2)
        public const double A_RAD = 0.00000000000000075657; // Radiation constant: J*m^(-3)*K^(-4)
        public const double K = 0.000000000000000000000013806488; // Boltzmann constant: J/K
        public const double SIGMA = 0.00000005670373; // Stefan-Boltzmann constant: W*m^(-2)*K
        public const double ATOMIC_MASS_UNIT = 0.000000000000000000000000001660538922; // kg
        public const double N_A = 602214129000000000000000.0; // Avogadro constant: mol^(-1)
        public const double L_SUN = 384600000000000000000000000.0; // W
        public const double M_SUN = 1989000000000000000000000000000.0; // kg
        public const double R_SUN = 696000000.0; // m

        //      2) mass fractions (X: hydrogen, Y3/Y: helium-3/helium-4, Z: metals)
        public const double X = 0.7;
        public const double Y3 = 2E-5; // changed from the given value: 1E-10
        public const double Y = 0.29;
        public const double Z = 0.01;
        public const double Z7Li = 1E-12;
        public const double Z7Be = 1E-12;

        //      3) element masses (u)
        public const double M_e = 0.0005486;
        public const double M_p = 1.007; // since the elements are ionized
        public const double M_3He = 3.016;
        public const double M_4He = 4.0026;
        public const double M_7Be = 7.017;
        public const double M_7Li = 7.016;

        //      4) constants for dynamic step size scaling
        public const double LOWER_DM_LIMIT = -1E-3; // if the step size goes below this value (1 g), stop shrinking the step size and abort

        //      5) flags for which processes in the PP chain to use
        public const bool USE_PP = true;
        public const bool USE_33 = true;
        public const bool USE_34 = true;
        public const bool USE_E7 = true;
        public const bool USE_17X = true;
        public const bool USE_17 = true;

        //      6) flag to switch convection on and off
        public const bool CONVECTION = true;

        // variable parameters (updated every step)
        public double L; // luminosity (W)
        public double R; // radius (m)
        public double M; // mass (kg)
        public double rho; // density (kg/m^(3))
        public double T; // temperature (K)
        public double P; // pressure (Pa)
        public double dm; // resolution

        // use static or dynamic dm?
        public bool dynamicDm;
        public double upperLimit; // max change (%/100 of current parameter value)
        public double lowerLimit; // min change (%/100 of current parameter value)
        public double scaleDownFactor; // scale down by this factor
        public double scaleUpFactor; // scale up by this factor
        protected bool forceStop; // true = break out of loop (end simulation prematurely)
        protected bool dmChanged; // set to true if we need to change the step size and do the current step over again

        // other variables
        protected double mu; // average molecular weight (kg)
        protected double C_EOS; // constant to be used in equation of state
        public double ALPHA; // constant used in calculating the mixing length (in reality a free parameter)
        protected int currentStep; // tracks the current step (for debugging purposes)
        public int plotStep; // keep data from 1 in ... steps (to save on memory)

        // store initial values of some parameters here
        public double initialR; // unit: solar radii
        public double initialP;
        public double initialL;
        public double initialT;
        public double initialRho;
        public double initialM;
        public double initialDm; // initial step size

        // store final values of some parameters here
        public double finalR; // unit: solar radii
        public double finalP;
        public double finalL;
        public double finalT;
        public double finalRho;
        public double finalM; // point where execution halts (unit: solar masses)
        public double finalEpsilon;

        // arrays to store results (collected at regular intervals)
        public List<double> LValues;
        public List<double> RValues;
        public List<double> MValues;
        public List<double> rhoValues;
        public List<double> TValues;
        public List<double> PValues;
        public List<double> epsilonValues;

        // variables for epsilons (see arrays below for explanation)
        public double epsilon;
        public double e_pp;
        public double e_33;
        public double e_34;
        public double e_e7;
        public double e_17x;
        public double e_17;

        // arrays to store epsilon values for individual processes
        public List<double> e_pp_values; // PPI - PPIII
        public List<double> e_33_values; // PPI
        public List<double> e_34_values; // PPII & PPIII
        public List<double> e_e7_values; // PPII
        public List<double> e_17x_values; // 7Li --> alpha (PP II)
        public List<double> e_17_values; // 7Be --> 8B (PP III)

        // variables for PP chain contributions
        public double PPI; // fraction of total energy
        public double PPIIandIII; // fraction of total energy

        // arrays to store PP chain contributions
        public List<double> PPI_values;
        public List<double> PPIIandIII_values;

        // variables for fluxes
        public double FC; // convective flux
        public double FR; // radiative flux
        public double F; // total flux

        // arrays to store fluxes as fractions of total flux
        public List<double> FC_values;
        public List<double> FR_values;

        // variables for gradients
        public double nabla; // actual gradient (= nabla_rad when there is no convection)
        public double nabla_rad; // radiation-only gradient

        // arrays to store gradients
        public List<double> nabla_values;
        public List<double> nabla_rad_values;

        // opacity matrix (from file)
        public MatrixFromFile om;

        public Simulation()
        {
            this.om = new MatrixFromFile(OPM_FILE); // get opacity matrix from file

            // calculate average molecular weight using X, Y, Z, assuming full ionization
            // (3.20) where we assume all metals are fully ionized
            // and use that almost all the He is 4He:
            double mu0 = (1 / (X + 0.25 * Y + 0.5 * Z));
            // since H, He also fully ionized, we apply the correction (3.24)
            // with nH+/nH = 1, nHe+/nHe = 0, nHe++/nHe = 1
            // (3.27) cannot be used directly, since it is derived from (3.15) and not (3.20)
            this.mu = (mu0 / (1.0 + mu0 * (X + 0.5 * Y)));
            // this.mu = 0.6; // DEBUG

            // C_EOS is here a constant to be used in equation of state: C = k / (mu * mass_unit)
            // mu is defined relative to the hydrogen mass (see 3.15), so mass_unit is this quantity
            this.C_EOS = (K / (this.mu * M_p * ATOMIC_MASS_UNIT));
        }

        // calculate opacity (using interpolation on table read from .txt file)
        // returns opacity in SI units (m^2 / kg)
        public double GetKappa(double T, double rho)
        {
            double logT = Math.Log10(T);
            double rhoGC = (rho / 1000.0); // unit conversion: kg*m^(-3) --> g*cm^(-3)
            double T6 = (T * Math.Pow(10, -6)); // T6 has unit 10*(6) K
            double logR = Math.Log10(rhoGC / T6); // r = P/T6 where T6 is T in units of 10^6 K

            // find closest matches to logT in the opacity table
            // such that logTPrev <= logT <= logTNext and the three are as close as possible
            // assume the column keys are sorted in ascending order
            double logTPrev = 0;
            double logTNext = 0;
            int n = om.RowKeys.Length;
            for (int i = 0; i < (n - 1); i++) // loop through each PAIR of indices (hence n-1)
            {
                if ((om.RowKeys[i] == logT) || (om.RowKeys[i + 1] == logT)) // exact match
                {
                    logTPrev = logT;
                    logTNext = logT;
                    break;
                }
                else if ((om.RowKeys[i] < logT) && (om.RowKeys[i + 1] > logT)) // between two values
                {
                    logTPrev = om.RowKeys[i];
                    logTNext = om.RowKeys[i + 1];
                    break;
                }
                else if ((i == 0) && (om.RowKeys[i] > logT)) // only happens if we need to extrapolate at the low end
                {
                    logTPrev = om.RowKeys[0];
                    logTNext = om.RowKeys[1];
                    break;
                }
                else if ((i == (n - 2)) && (om.RowKeys[i + 1] < logT)) // only happens if we need to extrapolate at the high end
                {
                    logTPrev = om.RowKeys[n - 2];
                    logTNext = om.RowKeys[n - 1];
                }
                // else keep looping
            }

            // find closest matches to logR in the opacity table
            // such that logRPrev <= logR <= logRNext and the three are as close as possible
            // assume the row keys are sorted in ascending order
            double logRPrev = 0;
            double logRNext = 0;
            n = om.ColKeys.Length;
            for (int i = 0; i < (n - 1); i++) // loop through each PAIR of indices (hence n-1)
            {
                if ((om.ColKeys[i] == logR) || (om.ColKeys[i + 1] == logR)) // exact match
                {
                    logRPrev = logR;
                    logRNext = logR;
                    break;
                }
                else if ((om.ColKeys[i] < logR) && (om.ColKeys[i + 1] > logR)) // between two values
                {
                    logRPrev = om.ColKeys[i];
                    logRNext = om.ColKeys[i + 1];
                    break;
                }
                else if ((i == 0) && (om.ColKeys[i] > logR)) // only happens if we need to extrapolate at the low end
                {
                    logRPrev = om.ColKeys[0];
                    logRNext = om.ColKeys[1];
                    break;
                }
                else if ((i == (n - 2)) && (om.ColKeys[i + 1] < logR)) // only happens if we need to extrapolate at the high end
                {
                    logRPrev = om.ColKeys[n - 2];
                    logRNext = om.ColKeys[n - 1];
                }
                // else keep looping
            }

            // intrapolate (2D) - 4 points in xy-space, each with a z value.
            // using a plane in xyz-space (linear interpolation) with x = logR, y = logT, z = logK
            // since we know little about what the function logK(log r, log T) looks like
            //
            // *2              |              *1
            //                 |
            //                 |
            // --------------------------------
            //                 |
            //        o        |
            // *3              |              *4
            //
            // Use the closest point to calculate the partial derivatives.
            // (If we are at point o, closest to 3, we use points 3 and 4 to calculate dz/dx and 3 and 2 to calculate dz/dy.)
            double dzdx;
            double dzdy;
            double dx;
            double dy;

            try
            {
                if (logRPrev == logRNext) // exact match
                {
                    dx = 0;
                    dzdx = 0;
                }
                else
                {
                    // delta x (log r)
                    dx = logR - logRPrev;

                    // we are closer to logTPrev than we are to logTNext
                    if (Math.Abs(logT - logTPrev) < Math.Abs(logTNext - logT))
                    {
                        // use definition of partial derivative at y = logTPrev
                        dzdx = ((om.Get(logRNext, logTPrev) - om.Get(logRPrev, logTPrev)) / (logRNext - logRPrev));
                    }
                    else
                    {
                        // use definition of partial derivative at y = logTNext
                        dzdx = ((om.Get(logRNext, logTNext) - om.Get(logRPrev, logTNext)) / (logRNext - logRPrev));
                    }
                }

                if (logTPrev == logTNext) // exact match
                {
                    dy = 0;
                    dzdy = 0;
                }
                else
                {
                    // delta y (log T)
                    dy = logT - logTPrev;

                    // we are closer to logRPrev than we are to logRNext
                    if (Math.Abs(logR - logRPrev) < Math.Abs(logRNext - logR))
                    {
                        // use definition of partial derivative at x = logRPrev
                        dzdy = ((om.Get(logRPrev, logTNext) - om.Get(logRPrev, logTPrev)) / (logTNext - logTPrev));
                    }
                    else
                    {
                        // use definition of partial derivative at x = logRNext
                        dzdy = ((om.Get(logRNext, logTNext) - om.Get(logRNext, logTPrev)) / (logTNext - logTPrev));
                    }
                }

                double dz = (dzdx * dx) + (dzdy * dy); // total derivative
                double z = om.Get(logRPrev, logTPrev);
                double logK = (z + dz);

                return ((Math.Pow(10, logK) / 10.0)); // the intrapolated opacity value (unit conversion: cm^2/g ---> m^2/kg)
            }
            catch (KeyNotFoundException) // happens when logR = logT = 0, an indicator of some numerical problem
            {
                this.forceStop = true;
                return 1;
            }
        }

        // calculate radiation pressure
        public double GetPrad(double T)
        {
            return (A_RAD * Math.Pow(T, 4) / 3.0);
        }

        // calculate P from T, rho
        public double GetP(double T, double rho)
        {
            double Pg; // gas pressure

            // calculate P using the ideal gas law on the from of (3.6)
            // NOTE: We must calculate and add/subtract the radiation pressure here, since the ideal gas law uses gas pressure only

            Pg = (this.C_EOS * T * rho); // gas pressure
            return (Pg + GetPrad(T)); // total pressure = gas pressure + radiation pressure
        }

        // calculate rho from T, P
        public double GetRho(double T, double P)
        {
            double Pg; // gas pressure

            // calculate rho using the ideal gas law on the from of (3.6)
            // NOTE: We must calculate and add/subtract the radiation pressure here, since the ideal gas law uses gas pressure only

            Pg = this.P - GetPrad(T); // gas pressure = total pressure - radiation pressure

            // C is here a constant to be used in equation of state: C = k / (mu * mass_unit)
            return (Pg / (this.C_EOS * T));
        }

        // calculate n_i or n_k for an element (used in calculating r_ik)
        // massFraction is X/Y/Z
        // m is the mass of the molecule in atomic masses
        // returns n in SI units
        public double n(double massFraction, double rho, double m)
        {
            return (massFraction * rho / (m * ATOMIC_MASS_UNIT));
        }

        // calculate reaction rate (# of reactions / second) for a certain step in one of the pp chains
        // n_i and n_k are number densities in m^(-3)
        // lambda_ik is the reaction rate in m^3 / second (based on data from table 2.1)
        // d_ik is the Kronecker delta function. Its value is 1 if i and k are the same type of element, 0 otherwise
        public double r_ik(double n_i, double n_k, double rho, double lambda_ik, double d_ik)
        {
            return (n_i * n_k * lambda_ik / (rho * (1 + d_ik)));
        }

        // convert NA_lambda-values from table 2.1 (where NA is Avogadro's constant and
        // lambda is given in CGS units) to pure lambda values in SI units (m^3/s)
        public double ConvertLambda(double N_A_lambda)
        {
            // 10^-6 is the conversion factor: CGS --> SI
            return ((N_A_lambda / N_A) * Math.Pow(10, -6));
        }

        // computes the reaction rate (m^3/s) of a given step in the pp chain
        // (using T9 as temperature in 10^9 K)
        // n_e is the electron density (m^-3) used for the e7 step (default value: 0)
        public double lambda_ik(string stepName, double T9, double n_e = 0)
        {
            double lambda = 0;

            // T9-based values
            double T916 = Math.Pow(T9, (1.0 / 6.0)); // T9^(1/6)
            double T913 = Math.Pow(T916, 2); // T9^(1/3)
            double T923 = Math.Pow(T913, 2); // T9^(2/3)
            double T943 = T9 * T913; // T9^(4/3)
            double T9m13 = Math.Pow(T913, -1); // T9^-(1/3)
            double T9m12 = Math.Pow((T913 * T916), -1); // T9^-(1/2) 
            double T9m23 = Math.Pow(T9m13, 2); // T9^-(2/3)
            double T9m1 = Math.Pow(T9m13, 3); // T9^-1
            double T9m32 = T9m1 * T9m12; // T9^-(3/2)
            double T9m53 = T9m1 * T9m23; // T9^-(5/3)
            double T9_34_56 = Math.Pow((T9 / (1 + 0.0495 * T9)), (5.0 / 6.0)); // used in lambda_34
            double T9_34_m13 = Math.Pow((T9 / (1 + 0.0495 * T9)), -(1.0 / 3.0)); // used in lambda_34
            double T9_17_56 = Math.Pow((T9 / (1 + 0.0759 * T9)), (5.0 / 6.0)); // used in lambda_17'
            double T9_17_m13 = Math.Pow((T9 / (1 + 0.0759 * T9)), -(1.0 / 3.0)); // used in lambda_17'

            switch (stepName)
            {
                case "pp":
                    lambda = ConvertLambda(4.01 * Math.Pow(10, -15) * T9m23 * Math.Exp(-3.380 * T9m13) * (1 + 0.123 * T913 + 1.09 * T923 + 0.938 * T9));
                    break;

                case "33":
                    lambda = ConvertLambda(6.04 * Math.Pow(10, 10) * T9m23 * Math.Exp(-12.276 * T9m13) * (1 + 0.034 * T913 - 0.522 * T923 - 0.124 * T9 + 0.353 * T943 + 0.213 * T9m53));
                    break;

                case "34":
                    lambda = ConvertLambda(5.61 * Math.Pow(10, 6) * T9_34_56 * T9m32 * Math.Exp(-12.826 * T9_34_m13));
                    break;

                case "e7":
                    double e7lambda = 1.34 * Math.Pow(10, -10) * T9m12 * (1 - 0.537 * T913 + 3.86 * T923 + 0.0027 * T9m1 * Math.Exp(2.515 * Math.Pow(10, -3) * T9m1));
                    // RATE MUST NOT EXCEED l.51e-07/n_e F0R T9 LESS THAN 0.001.
                    double treshold = 1.51 * Math.Pow(10, -7) / n_e;
                    if ((T9 < 0.001) && (e7lambda > treshold))
                    {
                        e7lambda = treshold;
                    }
                    lambda = ConvertLambda(e7lambda);
                    break;

                case "17'": // 7Li
                    lambda = ConvertLambda(1.096 * Math.Pow(10, 9) * T9m23 * Math.Exp(-8.472 * T9m13) - 4.830 * Math.Pow(10, 8) * T9_17_56 * T9m23 * Math.Exp(-8.472 * T9_17_m13) + 1.06 * Math.Pow(10, 10) * T9m32 * Math.Exp(-30.442 * T9m1));
                    break;

                case "17": // 7Be
                    lambda = ConvertLambda(3.11 * Math.Pow(10, 5) * T9m23 * Math.Exp(-10.262 * T9m13) + 2.53 * Math.Pow(10, 3) * T9m32 * Math.Exp(-7.306 * T9m1));
                    break;

                default: // in case of a typo in the code
                    throw new Exception("Unknown fusion process (lambda_ik)");
            }

            return lambda;
        }

        // converts an amount of energy in MeV to J
        public double ConvertToJ(double energyMeV)
        {
            // 1 eV = 1.6 * 10^-19 J
            // 1 MeV = 10^6 * 1.6 * 10^-19 J = 1.6 * 10^-13 J
            return (energyMeV * 1.6 * Math.Pow(10, -13));
        }

        // calculare energy/second from pp chain reactions
        public double GetEpsilon(double rho, double T)
        {
            double T9 = (T / (Math.Pow(10, 9))); // T in units of 10^9 kelvin

            // we add to this quantity for each reaction:
            // epsilon = sum( r_ik, Q_ik )
            // where r_ik is the reaction rate and Q_ik is the energy produced per reaction
            // (reactions / second) * (energy / reaction) = (energy / second)
            double epsilon = 0;

            // Note: We do not need to know n_1H and n_2H, just n_H.
            // We assume that we get all the 2H we need from the first step in the pp chain
            // (rate of production = rate of destruction), and treat step 1) and 2) as one
            // in terms of reaction rate. As such we assume that ALL the hydrogen is 1H
            // to begin with.

            // number densities (calculated from mass fractions)
            double n_H = this.n(X, rho, M_p); // assumes no deuterium initially
            double n_3He = this.n(Y3, rho, M_3He);
            double n_4He = this.n((Y - Y3), rho, M_4He);
            double n_7Be = this.n(Z7Be, rho, M_7Be);
            double n_7Li = this.n(Z7Li, rho, M_7Li);

            // calculate number density of electrons as a sum of:
            double n_e_H = n_H; // fully ionized hydrogen
            double n_e_He = 2 * (n_3He + n_4He); // fully ionized helium
            double n_e_Z = this.n(Z, rho, (2 * M_p)); // using (3.19), accurate since we only look at electrons and not total

            double n_e = n_e_H + n_e_He + n_e_Z;

            // variable to keep track of the energy produced in each step
            // (excepting the energy carried by neutrinos as these do not
            // contribute to the luminosity/temperature/photon pressure)
            // we will measure this energy in J
            double Q_ik;

            // other variables for easier debugging
            double r;

            // lambdas for individual steps (used in calculating branch contributions)
            double lambda_pp = 0;
            double lambda_33 = 0;
            double lambda_34 = 0;
            double lambda_e7 = 0;
            double lambda_17x = 0;
            double lambda_17 = 0;

            // reset the energy contributions so they are set to 0 if the flag for the step is set to false
            this.e_pp = 0;
            this.e_33 = 0;
            this.e_34 = 0;
            this.e_e7 = 0;
            this.e_17x = 0;
            this.e_17 = 0;

            // Reaction                                         Delta function      Lambda function
            // ====================================================================================
            // 1)   1H + 1H --> 2H (all variants)               d_ik = 1            lambda_pp
            // 2)   2H + 1H --> 3He (all variants)              (treat these steps as one step)
            if (USE_PP)
            {
                Q_ik = ConvertToJ(0.15 + 1.02 + 5.49); // (2.1) + (2.2): 1.02 MeV from e+
                lambda_pp = lambda_ik("pp", T9);
                r = r_ik(n_H, n_H, rho, lambda_pp, 1);
                this.e_pp = r * Q_ik;
                epsilon += this.e_pp;
            }

            // 3a)  3He + 3He --> 4He + 1H + 1H (PP I)          d_ik = 1            lambda_33
            if (USE_33)
            {
                Q_ik = ConvertToJ(12.86); // (2.3)
                lambda_33 = lambda_ik("33", T9); // using this to calculate the PPI contibution
                r = r_ik(n_3He, n_3He, rho, lambda_33, 1);
                this.e_33 = r * Q_ik;
                epsilon += this.e_33;
            }

            // 3bc) 3He + 4He --> 7Be (PP II && PP III)         d_ik = 0            lambda_34
            if (USE_34)
            {
                Q_ik = ConvertToJ(1.59); // (2.4) & (2.5)
                lambda_34 = lambda_ik("34", T9); // using this to calculate the PPII (and III) contibution
                r = r_ik(n_3He, n_4He, rho, lambda_34, 0);
                this.e_34 = r * Q_ik;
                epsilon += this.e_34;
            }

            // 4b)  7Be + e- --> 7Li (PP II)                    d_ik = 0            lambda_e7
            if (USE_E7)
            {
                Q_ik = ConvertToJ(0.81); // (2.4): weighted average of dE
                lambda_e7 = lambda_ik("e7", T9, n_e); // use electron density here!
                r = r_ik(n_7Be, n_e, rho, lambda_e7, 0);
                this.e_e7 = r * Q_ik;
                epsilon += this.e_e7;
            }

            // 5b)  7Li + 1H --> 4He + 4He (PP II)              d_ik = 0            lambda_17'
            if (USE_17X)
            {
                Q_ik = ConvertToJ(17.35); // (2.4)
                lambda_17x = lambda_ik("17'", T9);
                r = r_ik(n_7Li, n_H, rho, lambda_17x, 0);
                this.e_17x = r * Q_ik;
                epsilon += this.e_17x;
            }

            // 4c)  7Be + 1H --> 8B (PP III)                    d_ik = 0            lambda_17
            // 5c)  8B --> 8Be + e+ (PP III)                    (instant decay)
            // 6c)  8Be --> 4He + 4 He (PP III)                 (instant decay)
            if (USE_17)
            {
                Q_ik = ConvertToJ(0.14 + 1.02 + 6.88 + 3.00); // (2.5): 1.02 MeV from e+
                lambda_17 = lambda_ik("17", T9);
                r = r_ik(n_7Be, n_H, rho, lambda_17, 0);
                this.e_17 = r * Q_ik;
                epsilon += this.e_17;
            }

            // variables used in calculating the fraction a branch of the chain contributes to the total energy
            double lambda_PPI = 0; // branching into PP I
            double lambda_PPIIandIII = 0; // branching into PPII and III (assuming PPIII contribution is negligible)

            // check flags so we choose the right steps (and not disabled steps)
            // PP I
            if (USE_33)
            {
                lambda_PPI += lambda_33;
            }
            // else 0 (no subsequent step)

            // PP II & III
            if (USE_34)
            {
                lambda_PPIIandIII += lambda_34;
            }
            else
            {
                // PP II contributions
                if (USE_E7)
                {
                    lambda_PPIIandIII += lambda_e7;
                }
                else if (USE_17X)
                {
                    lambda_PPIIandIII += lambda_17x;
                }

                // PPIII contributions
                if (USE_17)
                {
                    lambda_PPIIandIII += lambda_17;
                }

                // if none of these steps were enabled then 0
            }

            // determine how much each chain contributed to the total energy using fractions
            // we assume that these lambdas cannot both be zero
            double X_PPI = (lambda_PPI / (lambda_PPI + lambda_PPIIandIII));
            double X_PPIIandIII = (lambda_PPIIandIII / (lambda_PPI + lambda_PPIIandIII));
            
            // take differing energies from different steps into account (using fractions for common steps)
            this.PPI = (X_PPI * this.e_pp) + this.e_33;
            this.PPIIandIII = (X_PPIIandIII * this.e_pp) + this.e_34 + this.e_e7 + this.e_17x + this.e_17;

            return epsilon;
        }

        // calculate partial derivative for r with respect to m
        public double Get_dr_dm(double r, double rho)
        {
            double dr_dm = (1.0 / (4 * Math.PI * Math.Pow(r, 2) * rho));

            this.CheckParameter(r, dr_dm); // check for numerical errors and adjust step size if possible

            return dr_dm;
        }

        // calculate partial derivative for P with respect to m
        public double Get_dP_dm(double m, double r)
        {
            double dP_dm = -(G * m / (4 * Math.PI * Math.Pow(r, 4)));

            this.CheckParameter(this.P, dP_dm); // check for numerical errors and adjust step size if possible

            return dP_dm;
        }

        // calculate partial derivative for L with respect to m
        public double Get_dL_dm(double rho, double T)
        {
            double dL_dm = this.GetEpsilon(rho, T);

            this.CheckParameter(this.L, dL_dm); // check for numerical errors and adjust step size if possible

            return dL_dm;
        }

        // calculate partial derivative for T with respect to m (taking the possibility of convection into account)
        public double Get_dT_dm(double T, double rho, double P, double m, double r, double L)
        {
            double dT_dm;
            double kappa = this.GetKappa(T, rho);            
            double nabla_rad = this.GetNabla_rad(kappa, L, T, P, m);

            if (this.UseConvection(nabla_rad)) // convection zone
            {
                double g = this.Get_g(m, r);
                double a = this.GetA(T, kappa, rho, g, P);
                double b = this.GetB(T, rho, P);
                double c = this.GetC(kappa, P, g, T);
                double d = this.GetD(nabla_rad);
                double xi = this.GetXi(a, b, c, d);
                double nabla = this.GetNabla(xi, a);

                dT_dm = -((g * T * nabla) / (4 * Math.PI * Math.Pow(r, 2) * P)); // dT/dm in the convective case
            }
            else // no convection: radiation zone
            {
                dT_dm = this.Get_dT_dm_rad(kappa, L, r, T); // use the radiation-only partial derivative
            }

            this.CheckParameter(T, dT_dm); // check for numerical errors and adjust step size if possible

            return dT_dm;
        }

        // calculate the pressure scale height (used for debugging)
        public double GetHP(double P, double g, double rho)
        {
            return (P / (g * rho));
        }

        // calculate the gradient when there is only radiation, no convection (unitless)
        public double GetNabla_rad(double kappa, double L, double T, double P, double m)
        {
            return ((3 * kappa * P * L) / (64 * SIGMA * G * m * Math.PI * Math.Pow(T, 4)));
        }

        // calculate the adiabatic gradient (unitless)
        public double nabla_ad
        {
            get
            {
                return (0.4); // see report for why this is 2/5
            }
        }

        // determine whether the criteria for convection to occur are met: the instability criterion, eq. (4.56)
        public bool UseConvection(double nabla_rad)
        {
            // Note that the nabla in (4.56) is nabla_rad - we assume radiation only, and check for stability
            // If this part of the star is not stable with regard to convection, we know that nabla > nabla_ad
            // and have to proceed using the convective equations instead of the radiative ones.

            if (CONVECTION) // flag to enable/disable convection entirely
            {
                return (nabla_rad > this.nabla_ad); // true: convection. false: no convection.
            }
            else
            {
                return false; // no convection
            }
        }

        // calculate partial derivative for T with respect to m (assuming only radiation)
        public double Get_dT_dm_rad(double kappa, double L, double r, double T)
        {
            return -((3 * kappa * L) / (256 * Math.Pow(Math.PI, 2) * SIGMA * Math.Pow(r, 4) * Math.Pow(T, 3)));
        }

        // calculate the heat capacity at constant pressure (per unit mass), unit: J/(K*kg)
        public double cP
        {
            get
            {
                return ((5 * K) / (2 * this.mu * ATOMIC_MASS_UNIT));
            }
        }

        // calculate the local gravitational acceleration (absolute value), using Newton's law of gravitation
        // (only mass inside the spherical shell contributes)
        public double Get_g(double m, double r)
        {
            return ((G * m) / (Math.Pow(r, 2)));
        }
        
        // calculate the "constant" of eq. (4.77)
        public double GetA(double T, double kappa, double rho, double g, double P)
        {
            return (256 * SIGMA * Math.Pow(T, 3) * g / (3 * kappa * this.cP * Math.Pow(ALPHA, 2) * Math.Sqrt(rho) * Math.Pow(P, 1.5)));
        }

        // calculate the "constant" of eq. (4.78)
        public double GetB(double T, double rho, double P)
        {
            return (0.25 * this.cP * T * Math.Pow(ALPHA, 2) * Math.Sqrt(rho * P));
        }

        // calculate the "constant" of eq. (4.79) and (4.80)
        public double GetC(double kappa, double P, double g, double T)
        {
            return ((16 * SIGMA * g * Math.Pow(T, 4)) / (3 * kappa * P));
        }

        // the final constant in the cubic equation for xi
        public double GetD(double nabla_rad)
        {
            return (this.nabla_ad - nabla_rad);
        }
        
        // calculate the convective gradient (unitless)
        public double GetNabla(double xi, double a)
        {
            return ((xi + a) * xi + this.nabla_ad); // nested form
        }

        // calculate xi when nabla is known (currently not used)
        public double GetXi(double nabla, double a)
        {
            return ((-a + Math.Sqrt(Math.Pow(a, 2) + 4 * (nabla - this.nabla_ad))) / 2); // = 0 when nabla = nabla_ad
        }

        // solve cubic equation to obtain xi = (nabla - nabla*)^(0.5)
        // note that a, b, c and d are constants from the original equations, not the coefficients
        public double GetXi(double a, double b, double c, double d)
        {
            // Solve the cubic equation for xi: Note that Solve.Cubic((b / c), 1, a, d) should always give the same result.
            double[] roots = Solve.Cubic((4 / a), 1, a, d).ToArray<double>();

            if (roots.Length != 1) // we can expect that there will always be exactly one solution (see report)
            {
                throw new Exception("Unexpected number of solutions for xi: " + roots.Length);
            }
            else if (roots[0] < 0) // additionally, we know that the solution will be positive (see report)
            {
                throw new Exception("Negative xi from cubic equation: " + roots[0]);
            }

            return roots[0];
        }

        // calculate the radiative flux
        public double GetFR(double F, double FC)
        {
            return (F - FC);
        }

        // calculate the total flux (from luminosity)
        public double GetF(double L, double r)
        {
            return (L / (4 * Math.PI * Math.Pow(r, 2)));
        }

        // calculate the convective flux
        public double GetFC(double b, double xi, double nabla_rad)
        {
            if (this.UseConvection(nabla_rad)) // convection zone
            {
                return (b * Math.Pow(xi, 3)); // by eq. (4.78)
            }
            else // radiative zone
            {
                return 0; // no convective flux
            }
        }

        // calculate other things need for plots and for the next step (call this before updating the parameters)
        public void StepCalculations(double m, double r, double P, double rho, double T, double L)
        {
            // calculate properties from the current parameters (that we got in the previous step)
            double kappa = this.GetKappa(T, rho);
            this.nabla_rad = this.GetNabla_rad(kappa, L, T, P, m); // we will plot this

            // total flux
            this.F = this.GetF(L, r);

            // if there is convection
            if (this.UseConvection(nabla_rad))
            {
                // calculate the following additional properties
                double g = this.Get_g(m, r);
                double a = this.GetA(T, kappa, rho, g, P);
                double b = this.GetB(T, rho, P);
                double c = this.GetC(kappa, P, g, T);
                double d = this.GetD(this.nabla_rad);
                double xi = this.GetXi(a, b, c, d);
                this.nabla = this.GetNabla(xi, a); // we will plot this

                this.FC = this.GetFC(b, xi, this.nabla_rad); // convective flux
            }
            else // no convection
            {
                this.FC = 0; // no convective flux
                this.nabla = this.nabla_rad;
            }

            this.FR = this.GetFR(this.F, this.FC);
        }

        // check for numerical problems before and after the step
        // especially useful if the step size is dynamic, so we can avoid problems by scaling the step size down
        // NOTE: We do not scale UP the step size here, to do that we check all relevant parameters at once
        public void CheckParameter(double x, double dx_dm)
        {
            // calculate change in parameter from partial derivative
            double delta_x = (dx_dm * this.dm);

            // do this check only if the step size is not static
            if (this.dynamicDm)
            {
                // check for numerical difficulties that can arise (values that are out of bounds)
                if (Double.IsNaN(x) || Double.IsNaN(delta_x))
                {
                    if (this.dm < LOWER_DM_LIMIT) // if the step size is not already too small
                    {
                        // scale down the step size
                        this.dm = this.scaleDownFactor * dm;
                        this.dmChanged = true;
                    }
                    else
                    {
                        this.forceStop = true; // halt execution
                    }
                }
                else if (Double.IsInfinity(x) || Double.IsInfinity(delta_x))
                {
                    if (this.dm < LOWER_DM_LIMIT) // if the step size is not already too small
                    {
                        // scale down the step size
                        this.dm = this.scaleDownFactor * dm;
                        this.dmChanged = true;
                    }
                    else
                    {
                        this.forceStop = true; // halt execution
                    }
                }
                // check whether we get too large a step (more than the given % of reference value for ANY parameter)
                else if (Math.Abs(delta_x / x) > upperLimit)
                {
                    if (this.dm < LOWER_DM_LIMIT) // if the step size is not already too small
                    {
                        // scale down the step size
                        this.dm = this.scaleDownFactor * dm;
                        this.dmChanged = true;
                    }
                    else
                    {
                        this.forceStop = true; // halt execution
                    }
                }
            }
            else
            {
                // check for numerical difficulties that can arise (values that are out of bounds)
                if (Double.IsNaN(x) || Double.IsNaN(delta_x))
                {
                    this.forceStop = true; // halt execution
                }
                else if (Double.IsInfinity(x) || Double.IsInfinity(delta_x))
                {
                    this.forceStop = true; // halt execution
                }
            }
        }

        // call if the step is not too large: check if it was too small
        // (less than the given % of the current value for ALL parameters in the arrays)
        public void ScaleUp(double[] x, double[] delta_x)
        {
            if (this.dynamicDm)
            {
                bool scaleUp = true; // set this to false if ANY parameter changes too much

                // check all parameters
                for (int i = 0; i < Math.Max(x.Length, delta_x.Length); i++)
                {
                    if (!(Math.Abs(delta_x[i] / x[i]) < this.lowerLimit))
                    {
                        scaleUp = false;
                    }
                }

                if (scaleUp)
                {
                    // scale up the step size (at a different rate than the scaling down so we don't create an infinite loop)
                    this.dm = this.scaleUpFactor * dm;
                    this.dmChanged = true;
                }
            }
        }

        // calculate values for initial step
        public void FirstStep(double dm)
        {
            /* DEBUGGING ONLY
            double T = 0.9E6;
            double rho = 55.9;
            double M = 0.99 * M_SUN;
            double R = 0.84 * R_SUN;
            double L = L_SUN;
            double g = this.Get_g(M, R);
            double kappaDebug = this.GetKappa(T, rho);
            double kappa = 3.9;
            double nabla_rad = this.GetNabla_rad(kappa, L, T, P, M);
            if (this.UseConvection(nabla_rad))
            {
                double HP = this.GetHP(P, g, rho);
                double a = this.GetA(T, kappa, rho, g, P);
                double b = this.GetB(T, rho, P);
                double c = this.GetC(kappa, P, g, T);
                double d = this.GetD(nabla_rad);
                double xi = this.GetXi(a, b, c, d);
                double nabla = this.GetNabla(xi, a);
            } */

            this.Step(dm);
        }

        // calculate new values for the next step
        public void Step(double dm)
        {
            if (dm > 0.0) // check that dm is actually negative or 0
            {
                throw new Exception("dm must be a negative number or 0.");
            }

            // variables to contain each step in the Runge-Kutta method
            // See: https://en.wikipedia.org/wiki/Runge-kutta
            // with t --> m, h --> dm, y --> partial derivative with respect to m (for r, P, L, T respectively)
            double k1;
            double k2;
            double k3;
            double k4;

            // calculate mass steps for Runge-Kutta method
            double m1 = this.M;
            double m2 = this.M + 0.5 * dm; // half-step
            double m3 = m2; // half-step
            double m4 = this.M + dm; // full step

            // check the mass parameter (the others will be checked when calculated)
            CheckParameter(this.M, 1.0); // dm/dm = 1

            // calculate change in radius (delta_R): dr/dm depends on r, but not m
            k1 = Get_dr_dm(this.R, this.rho);
            k2 = Get_dr_dm((this.R + 0.5 * dm * k1), this.rho);
            k3 = Get_dr_dm((this.R + 0.5 * dm * k2), this.rho);
            k4 = Get_dr_dm((this.R + dm * k3), this.rho);
            double delta_R = (dm / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

            // calculate change in pressure (delta_P): dP/dm depends on m, but not P
            k1 = Get_dP_dm(this.M, this.R);
            k2 = Get_dP_dm((this.M + 0.5 * dm), this.R);
            k3 = Get_dP_dm((this.M + 0.5 * dm), this.R);
            k4 = Get_dP_dm((this.M + dm), this.R);
            double delta_P = (dm / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

            // calculate change in luminosity (delta_L): dL/dm depends on NEITHER m or L
            double epsilon = Get_dL_dm(this.rho, this.T);
            k1 = epsilon;
            k2 = epsilon;
            k3 = epsilon;
            k4 = epsilon;
            double delta_L = (dm / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

            // calculate change in temperature (delta_T): dT/dm depends on both T and m (in the convective zone; only T in the radiative zone)
            k1 = Get_dT_dm(this.T, this.rho, this.P, this.M, this.R, this.L);
            k2 = Get_dT_dm((this.T + 0.5 * dm * k1), this.rho, this.P, (this.M + 0.5 * dm), this.R, this.L);
            k3 = Get_dT_dm((this.T + 0.5 * dm * k2), this.rho, this.P, (this.M + 0.5 * dm), this.R, this.L);
            k4 = Get_dT_dm((this.T + dm * k3), this.rho, this.P, (this.M + dm), this.R, this.L);
            double delta_T = (dm / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

            // finally, check whether the step size was too small if there were no other problems
            if (!this.dmChanged)
            {
                double[] parameters = { this.M, this.R, this.P, this.L, this.T };
                double[] deltas = { this.dm, delta_R, delta_P, delta_L, delta_T };
                ScaleUp(parameters, deltas);
            }

            if (!this.dmChanged) // either we use static dm, or we did not need to change it this step (so we can continue)
            {
                // calculate variables for plotting (NOTE: Must be done before we update the parameters!)
                this.StepCalculations(M, R, P, rho, T, L);

                // update parameters
                this.R += delta_R;
                this.P += delta_P;
                this.L += delta_L;
                this.T += delta_T;
                this.epsilon = epsilon;

                // finally calculate rho (using the updated T, P values) since drho/dm is not given
                this.rho = this.GetRho(this.T, this.P);

                // proceed to next step: dm < 0
                this.M += dm;

                // only now do we increment the step counter
                this.currentStep++;
            }
            else // we DID change the step size and will do this step over again with the new step size
            {
                this.dmChanged = false; // reset the flag
            }
        }

        // save current data values to arrays (to be used in plots)
        public void SaveCurrentValues()
        {
            this.LValues.Add(this.L / this.initialL); // plot luminosity (relative to initial value)
            this.RValues.Add(this.R / R_SUN); // plot radius in solar radii (x axis in our plots)
            this.MValues.Add(this.M / (this.initialM * M_SUN)); // plot mass (relative to initial value, given in solar masses)
            this.rhoValues.Add(this.rho);
            this.TValues.Add(this.T);
            this.PValues.Add(this.P);
            this.epsilonValues.Add(this.epsilon);
            this.e_pp_values.Add(this.e_pp / this.epsilon); // relative value
            this.e_33_values.Add(this.e_33 / this.epsilon); // relative value
            this.e_34_values.Add(this.e_34 / this.epsilon); // relative value
            this.e_e7_values.Add(this.e_e7 / this.epsilon); // relative value
            this.e_17x_values.Add(this.e_17x / this.epsilon); // relative value
            this.e_17_values.Add(this.e_17 / this.epsilon); // relative value
            this.PPI_values.Add(this.PPI / this.epsilon); // relative value
            this.PPIIandIII_values.Add(this.PPIIandIII / this.epsilon); // relative value
            this.FC_values.Add(this.FC / this.F); // relative value
            this.FR_values.Add(this.FR / this.F); // relative value
            this.nabla_values.Add(this.nabla);
            this.nabla_rad_values.Add(this.nabla_rad);
        }

        // prepares the simulation object for a new simulation
        public void ClearData()
        {
            // reset the step counter
            this.currentStep = 0;

            // reset flags
            this.forceStop = false;

            // clear data values in arrays (for a new round of calculations)
            this.LValues = new List<double>();
            this.RValues = new List<double>();
            this.MValues = new List<double>();
            this.rhoValues = new List<double>();
            this.TValues = new List<double>();
            this.PValues = new List<double>();
            this.epsilonValues = new List<double>();
            this.e_pp_values = new List<double>();
            this.e_33_values = new List<double>();
            this.e_34_values = new List<double>();
            this.e_e7_values = new List<double>();
            this.e_17x_values = new List<double>();
            this.e_17_values = new List<double>();
            this.PPI_values = new List<double>();
            this.PPIIandIII_values = new List<double>();
            this.FC_values = new List<double>();
            this.FR_values = new List<double>();
            this.nabla_values = new List<double>();
            this.nabla_rad_values = new List<double>();
        }

        // starts the loop that runs the simulation with the given parameters
        public void Run()
        {
            this.ClearData(); // clear data from previous simulations (if any) -- also resets this.currentStep
            this.SaveInitialValues(); // save initial values
            this.FirstStep(this.dm); // perform first step

            // loop until we reach the core or we're forced to stop short
            // this.currentStep will increment only if the step size did not change (see Step() for details)
            while ((this.M > 0) && (!forceStop))
            {
                this.Step(this.dm); // perform next step

                if (this.currentStep % this.plotStep == 0) // for every ...th step
                {
                    this.SaveCurrentValues(); // save values at this point
                }
            }

            this.SaveFinalValues(); // save final values
        }

        // record the initial values of the parameters
        public void SaveInitialValues()
        {
            // store initial parameter values
            this.initialM = (this.M / M_SUN); // unit --> solar masses
            this.initialDm = this.dm;
            this.initialL = this.L;
            this.initialP = this.P;
            this.initialR = (this.R / R_SUN); // unit --> solar radii
            this.initialRho = this.rho;
            this.initialT = this.T;

            this.SaveCurrentValues(); // make sure initial values are plotted
        }

        // record the final values of the parameters when execution halts
        public void SaveFinalValues()
        {
            this.SaveCurrentValues(); // make sure final values are plotted

            // store final parameter values: it's very important that this happens BEFORE the values are scaled down!
            this.finalM = (this.M / M_SUN); // unit --> solar masses
            this.finalL = this.L;
            this.finalP = this.P;
            this.finalR = (this.R / R_SUN); // unit --> solar radii
            this.finalRho = this.rho;
            this.finalT = this.T;
            this.finalEpsilon = this.epsilon;

            // finally scale all non-relative values that are not strictly decreasing to between 0 (min) and 1 (max)
            double maxT = this.TValues.Last<double>();
            for (int i = 0; i < this.TValues.Count; i++)
            {
                this.TValues[i] = (this.TValues[i] / maxT);
            }

            double maxRho = this.rhoValues.Last<double>();
            for (int i = 0; i < this.rhoValues.Count; i++)
            {
                this.rhoValues[i] = (this.rhoValues[i] / maxRho);
            }

            double maxP = this.PValues.Last<double>();
            for (int i = 0; i < this.PValues.Count; i++)
            {
                this.PValues[i] = (this.PValues[i] / maxP);
            }

            double maxEpsilon = this.epsilonValues.Last<double>();
            for (int i = 0; i < this.epsilonValues.Count; i++)
            {
                this.epsilonValues[i] = (this.epsilonValues[i] / maxEpsilon);
            }

            // reset flags
            this.forceStop = false;
        }
    }
}