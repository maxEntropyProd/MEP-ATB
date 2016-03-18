module globalVars
   implicit none
   ! Specify version
   character(len=80),parameter:: version = 'V4.1 (17 Mar 2012)' ! This same as V3, but replace rxn 9 with dN remineralization. This version differs from V 4 only in 
                                                                ! the calculation of delGR for reactions.  Use of exponents was minimized to improve speed.  Results should be the same
                                                                ! and they were checked.  Slight differeces, but nothing important.  So, results from runs with 4.0 and 4.1 can be intermixed.
   ! Defined parameters
   integer, parameter:: nc = 12 ! number of concentration state variables
   integer, parameter:: nbs = 4 ! number of bioloigcal structure state variables
   integer, parameter:: nSys = 2 ! number of other state variables, probably related to sigma calculation
   integer, parameter:: nfeed = 2 ! number of different feed conditions (i.e., stepping ch4 on and off)
   integer, parameter:: nrxn = 9 ! total number of reactions
   integer, parameter:: nocv = 6 ! number of optimal control variables at a single time point.
   integer, parameter:: nMaxKnot = 50 ! maximum number of spline knots NOT including the   
                                      ! one at the start of the interval (at ti).  Note, just using linear interpolation (not splines from TSPACK).
                                      
   ! derived parameter values from above                                      
   integer, parameter:: nsv = nSys + nc + nbs ! total number of state variables (first state vars are entropy related, sigma and sigmaW in this version)

   ! Global variables
   logical allocateVarsCalled ! used to see if system varialbes allocated.
   real(8), allocatable:: tknots(:), ocvMat(:,:) ! Matrix of OCV for spline fitting.  First row is from previous interval.
   
   real(8) ti, tii ! start and end times of current interval.
   real(8), allocatable:: solTpts(:) ! time points to save ODE and other variable solutions at.
   logical saveSoln ! whether ODE and other var are saved
   integer nODEpts
   real(8) state0(nsv) ! initial conditions of state for current interval stored here.
   integer fEval ! number of evaluations of objective function of an interval
   character basefilename*80   ! all files opened or created use this basefilename
   character*100 cnames(nc), bsnames(nbs), ocvnames(nocv), rxnnames(nrxn) ! names for Tecplot output files
   
   ! Unit filename
   logical unitsOpen ! whether units have been open yet.
   integer TCunit ! holds the concentration state variables, c
   integer TBSunit ! holds the Biological structure variables, bs
   integer TOCVunit ! holds the optimal control variables, ocv
   integer TRunit ! unit for reaction output
   integer TGunit ! output for reaction delta G's
   integer TSunit ! output of internal entropy and production rate sigma, sigmaDot, aveSigmaDot
   integer TSINVunit ! entropy output data only at end of interval.
   integer BOBunit ! output for BOBYQA info   
   integer DEBUGunit ! when debuging is needed.
   
   real(8) sigma, sigmaDot, aveSigmaDot ! entropy from internal production, internal entropy production, and averagted entropy prod.
   real(8) sigmaW, aveSigmaWdot ! weighted sigma and sigmaDot-average value over the tInfinity interval
   real(8), pointer:: ch4, ch3oh, h2co3, dC, hno3, nh3, dN, o2, pch4, pch3oh, pco2, po2, bs1, bs2, bs3, bs4
   real(8) ch4_F, ch3oh_F, h2co3_F, dC_F, hno3_F, nh3_F, dN_F, o2_F, pch4_F, pch3oh_F, pco2_F, po2_F, bs1_F, bs2_F, bs3_F, bs4_F !for feed
   real(8), pointer:: epsilon1, omega1, epsilon2, omega2, epsilon3, epsilon4        
   
   ! BOBYQA allocatable variables
   REAL(8), DIMENSION(:), ALLOCATABLE :: bl, bu ! Lower and upper bounds
   REAL(8), DIMENSION(:), ALLOCATABLE :: u ! Minimizing vector.
   REAL(8), DIMENSION(:), ALLOCATABLE :: w ! work vector.
      
   ! declarations for variables in namelist
   real(8) t0 ! start of simulation time (d)
   namelist /params/  t0
   integer nKnots ! number of spline points for time integration NOT including the first point, which comes from the previous interval.
   namelist /params/ nKnots
   integer nSolTpts ! number of time points per interval for saving solution in files, not including initial point of interval.
   namelist /params/ nSolTpts
   integer nIntervals ! number of intervals to run simulation for
   namelist /params/ nIntervals
   real(8) tIntv ! time (in days) of an interval
   real(8) tInfinity ! time (in days) of infinate horizon
   real(8) wInfinity ! value of weighting (discounting) function at tInfinity; range: (0,1].  Value of 1 means no weighting
   namelist /params/ tIntv, tInfinity, wInfinity
   real(8), target:: c(nc,nfeed+1), bs(nbs,nfeed+1) ! concentration and biological state variables at a specified time point and their feed conc.
   namelist /params/ c, bs
   real(8), target:: ocv(nocv,3) ! optimal control variables at time t, and their lower and upper bounds.
   namelist /params/ ocv
   real(8), target:: chon(nbs,3) ! CHON stoichiometry for each BS, unit carbon (i.e., #H, #O, #N, or alpha, beta, gamma)
   namelist /params/ chon  
   real(8) nu(nbs), kappa(nbs) ! reaction kinetic terms
   namelist /params/ nu, kappa
   real(8) tFeedCyclingOn, tFeedCyclingOff ! times (d) when feed cycling should be on and off
   real(8) tFeed2OnPeriod ! Duration (d) for feed 2 to be on (feed 1 is set to same duration for now).
   namelist /params/ tFeedCyclingOn, tFeedCyclingOff, tFeed2OnPeriod
   real(8) fLoss ! fraction of POM stat is lost via dilution from chemostat.  fLoss=0 means batch-like operation
   namelist /params/ fLoss
   
   ! variables associated with experiment
   real(8) Temp ! temperature (K)
   real(8) is ! ionic strength (M)
   real(8) pH
   real(8) volL, volG ! volume of the liquid and gas spaces in the reactor (m3).
   real(8) fLiq, fGas ! liquid and gas flow rates into reactor (m3/d)
   real(8) kLa ! Liquid-side mass transfer velocity, kL (m/d), times the air-water interface area, a (m2), including bubbles. 
               ! Overall kLa units m3/d
   namelist /params/ Temp, is, pH, volL, volG, fLiq, fGas, kLa

   ! BOBYQA parameters    
   integer npt  ! is the number of interpolation conditions. Its value must be in
                ! the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not recommended.
   real(8) rhobeg ! around one tenth of the greatest expected change to a variable
   real(8) rhoend ! accuracy that is required in the final values of the variables
   integer iprint ! print output (0,1,2,3).  0 for no output
   integer maxfun ! maximum function evaluations.                
   namelist /params/ npt, rhobeg, rhoend, iprint, maxfun
   
   ! BiM options
   real(8) rtol ! relative tolerance
   real(8) atol ! absolute tolerance
   integer BiMmaxIter ! max iterations given by iwk(1), default: 100000
   integer maxAttempts ! maximum times retry integration of interval
   real(8) h0 ! initial set size.  making this very small 1e-10 cans solve integration problems
   real(8) hmax !  HMAX. MAXIMUM INTEGRATION STEP. (DEFAULT = (TEND-T0)/8)   
   namelist /params/ rtol, atol, BiMmaxIter, maxAttempts, h0, hmax
   integer BiMretry, BiMfail ! counts on ODE integration problems   
   
   contains
   subroutine point2State ()
      implicit none
      ! This routine sets the pointers to the state vectors.  Note, when the BiM (ODE solver) passes the
      ! total state vector, x, to the state evaluating routine, the following can be used to set values
      ! of these pointers:
      ! sigma = x(1); sigmaW = x(2); c(1:nc,1) = x(nSys+1:nc+nSys); bs(1:nbs,1) = x(nc+nSys+1:nSys+nc+nbs) 
      ! Also setup the names for tecplot output
      ch4    => c(1,1);  cnames(1)  = 'CH<sub>4</sub>'
      ch3oh  => c(2,1);  cnames(2)  = 'CH<sub>3</sub>OH'
      h2co3  => c(3,1);  cnames(3)  = 'H<sub>2</sub>CO<sub>3</sub>'
      dC     => c(4,1);  cnames(4)  = 'dC'
      hno3   => c(5,1);  cnames(5)  = 'HNO<sub>3</sub>'
      nh3    => c(6,1);  cnames(6)  = 'NH<sub>3</sub>'
      dN     => c(7,1);  cnames(7)  = 'dN'
      o2     => c(8,1);  cnames(8)  = 'O<sub>2</sub>'
      pch4   => c(9,1);  cnames(9)  = 'pCH<sub>4</sub>'
      pch3oh => c(10,1); cnames(10) = 'pCH<sub>3</sub>OH'
      pco2   => c(11,1); cnames(11) = 'pCO<sub>2</sub>'
      po2    => c(12,1); cnames(12) = 'pO<sub>2</sub>'
      
      bs1 => bs(1,1); bsnames(1) = '<userdef>s</userdef><sub>1</sub>'
      bs2 => bs(2,1); bsnames(2) = '<userdef>s</userdef><sub>2</sub>'
      bs3 => bs(3,1); bsnames(3) = '<userdef>s</userdef><sub>3</sub>'
      bs4 => bs(4,1); bsnames(4) = '<userdef>s</userdef><sub>4</sub>'
      
      ! point to control variables
      epsilon1 => ocv(1,1); ocvnames(1) = '<greek>e</greek><sub>1</sub>' 
      omega1   => ocv(2,1); ocvnames(2) = '<greek>w</greek><sub>1</sub>'
      epsilon2 => ocv(3,1); ocvnames(3) = '<greek>e</greek><sub>2</sub>'
      omega2   => ocv(4,1); ocvnames(4) = '<greek>w</greek><sub>2</sub>'
      epsilon3 => ocv(5,1); ocvnames(5) = '<greek>e</greek><sub>3</sub>'
      epsilon4 => ocv(6,1); ocvnames(6) = '<greek>e</greek><sub>4</sub>'
      
      ! Set the reaction names
      rxnnames(1)  = 'rCH<sub>4</sub>+HNO<sub>3</sub>->CH<sub>3</sub>OH'
      rxnnames(2)  = 'rCH<sub>4</sub>+NH<sub>3</sub>->CH<sub>3</sub>OH'
      rxnnames(3)  = 'rCH<sub>3</sub>OH+HNO<sub>3</sub>->CO<sub>2</sub>'
      rxnnames(4)  = 'rCH<sub>3</sub>OH+NH<sub>3</sub>->CO<sub>2</sub>'
      rxnnames(5)  = 'r<userdef>s</userdef><sub>1</sub>'
      rxnnames(6)  = 'r<userdef>s</userdef><sub>2</sub>'
      rxnnames(7)  = 'r<userdef>s</userdef><sub>3</sub>'
      rxnnames(8)  = 'r<userdef>s</userdef><sub>4</sub>'
      rxnnames(9)  = 'r<userdef>s</userdef>-decomp'
      return
   end subroutine point2State
   
   subroutine initialize ()
      implicit none
      integer i
      ti = t0 ! start time of first interval
      tknots(1) = t0
      
      ! copy the state variables from PRM file to state0 (IC for ODE soln)
      state0(1)            = 0. ! value of sigma
      state0(2)            = 0. ! value of sigmaW
      state0(nSys+1:nc+nSys)       = c(1:nc,1) ! concentration variables
      state0(nc+nSys+1:nSys+nc+nbs) = bs(1:nbs,1) ! biological strucures.
      
      ! copy ocv to first column of ocvMat, which is needed for spline interpolation.
      ! Since reactions could be calculated, just fill whole ocvMat with IC
      do i=1,nknots+1
         ocvMat(i,1:nocv) = ocv(1:nocv,1)
      end do
                  
      ! copy the bounds on the OCV to the vectors for pVTDIRECT, this only
      ! needs to be done once, since pVTDIRECt does not change them (i think)
      do i=1,nknots
         bl( (i-1)*nocv+1:i*nocv ) = ocv(1:nocv,2)
         bu( (i-1)*nocv+1:i*nocv ) = ocv(1:nocv,3)
         ! initialize vector for BOBYQA (this is not needed for pVTDIRECT)
         u( (i-1)*nocv+1:i*nocv )  = ocv(1:nocv,1)
      end do
            
      return
   end subroutine initialize
   
   subroutine allocateVars ()
      implicit none
      ! allocate arrays
      allocate( tknots(nknots+1), ocvMat(nknots+1, nocv) )
      allocate( solTpts(nSolTpts) )
      ! BOBYQA space
      allocate( bl(nocv*nknots) )
      allocate( bu(nocv*nknots) )
      allocate( u(nocv*nknots) )
      allocate( w((npt+5)*(npt+nocv*nknots)+3*nocv*nknots*(nocv*nknots+5)/2) )
      
      allocateVarsCalled = .true.
      return      
   end subroutine allocateVars
   
   subroutine deallocateVars ()
      implicit none
      ! allocate arrays
      deallocate( tknots, ocvMat )
      !deallocate( TSP_yp, TSP_sigma )
      deallocate( solTpts )
      ! VTDIRECT space
      deallocate( bl )
      deallocate( bu )
      deallocate( u )
      deallocate( w )
      allocateVarsCalled = .false.
      return      
   end subroutine deallocateVars      
end module globalVars

module reactionSystem
   ! User based module to calculate reaction related variables.
   use globalVars   
   implicit none
   ! pointers to CHON elements
   real(8), pointer:: alpha1, beta1, gamma1
   real(8), pointer:: alpha2, beta2, gamma2
   real(8), pointer:: alpha3, beta3, gamma3
   real(8), pointer:: alpha4, beta4, gamma4
   ! Standard free energies of formation
   real(8) delGF0_bs(nbs)
   real(8) delGF0_ch4, delGF0_ch3oh, delGF0_h2co3, delGF0_dC, delGF0_hno3, delGF0_nh3, &
           delGF0_dN, delGF0_o2, delGF0_h2o, delGF0_ch4G, delGF0_ch3ohG, delGF0_o2G, delGF0_co2G 
   real(8) kH_ch4, kH_ch3oh, kH_h2co3, kH_o2
   real(8) a1, b11, d1, e11
   real(8)     b12, e12
   real(8) a2, b21, d2, e21
   real(8)     b22,     e22
   real(8) a3(nbs), b3(nbs), fp3(nbs)
   real(8) b4, d4, e4, f4, gammadCN
   
   real(8) delGR(nrxn), rxn(nrxn)
   
   contains
   
   subroutine initializeRS ()
      use thermoData
      ! For this version, I am not calculating  pH or temp changes, so delG of formation need only
      ! be calculated once, which saves some overhead.
      
      call getDelGF0 ()
      ! These are for vapor phase calculations (i.e., Henery's constants) in uM/atm. These are calculated
      ! once since this version does not allow temp, is, or pH to vary during the simulation.
      kH_ch4   = solCH4(temp,is,pH)
      kH_ch3oh = solCH3OH(temp,is,pH)
      kH_h2co3 = solH2CO3(temp,is,pH)
      kH_o2    = solO2(temp,is,pH)

      return
   end subroutine initializeRS

   subroutine stateEqns (m,t,x,dxdt,ierr,rpar,ipar)
      ! State equations governing sigma and sigmaDot.
      use globalVars
      use thermoData, only: Rpm3 ! this is gas constant in (atm m^3/(g-mmol K)) (note, mmol not mol)
      real(8) t,x(m),dxdt(m),rpar(*)
      integer m,ierr,ipar(*)
      
      integer i, j
      real(8) tempV, weightS

      ! Copy state to names via pointer association
      sigma = x(1); sigmaW = x(2); c(1:nc,1) = x(nSys+1:nc+nSys); bs(1:nbs,1) = x(nc+nSys+1:nSys+nc+nbs) 
      
      ! Don't allow any state variables to be less than zero
      forall(i=1:nc , c(i,1) < 0) c(i,1) = 0.
      forall(i=1:nbs, bs(i,1) < 0) bs(i,1) = 0.
     
      ! Set the feed concentrations
      call feedControl (t)
      
      ! get OCV at time t from previous spine fit to ocvMat
      call getOCV (t)  
      call transformOCV () ! do any problem specific transformations to OCV before proceeding.

      call evalReactionSystem (t) ! evaluate all components of the reaction system
      
      ! Construct and evaluate the state equations    
      ! dxdt(1) integrated value of sigmadot  (i.e., sigma)
      dxdt(1) = sigmaDot ! (J/d/deg-K) and aveSigmaDot = sigma/(t-ti) averaged value of sigmaDot at time t.
      
      ! weighted sigma: sigmaW      
      dxdt(2) = sigmaDot*weightS(t)
      
      ! ch4 balance
      dxdt(3) = -rxn(1) - rxn(2) + fLiq*(ch4_F - ch4)/volL + kLa*(pch4*kH_ch4 - ch4)/volL
 
      ! ch3oh balance
      dxdt(4) = d1*(rxn(1) + rxn(2)) - (rxn(3) + rxn(4)) + fLiq*(ch3oh_F - ch3oh)/volL + kLa*(pch3oh*kH_ch3oh - ch3oh)/volL

      ! h2co3 balance
      dxdt(5) = d2*(rxn(3) + rxn(4)) + epsilon3*(1. - epsilon3)*sum(rxn(5:4+nbs)) + d4*rxn(9) &
              + fLiq*(h2co3_F - h2co3)/volL + kLa*(pco2*kH_h2co3 - h2co3)/volL
              
      ! dC
      dxdt(6) = (1. - epsilon3)**2*sum(rxn(5:4+nbs)) - rxn(9) + fLiq*(dC_F - fLoss*dC)/volL
      
      ! hno3
      dxdt(7) = -a1*rxn(1) - a2*rxn(3) + fLiq*(hno3_F - hno3)/volL
      
      ! nh4
      tempV = 0.
      do i=1,nbs
         tempV = tempV + ( (1. - epsilon3)*chon(i,3) + (chon(i,3) - chon(3,3)) )*epsilon3*rxn(i+4)
      end do
      dxdt(8) = -a1*rxn(2) - a2*rxn(4) + tempV  + f4*rxn(9) + fLiq*(nh3_F - nh3)/volL
      
      ! dN
      dxdt(9) = (1. - epsilon3)**2*dot_product(chon(1:nbs,3),rxn(5:4+nbs)) - gammadCN*rxn(9) + fLiq*(dN_F - fLoss*dN)/volL
      
      ! O2 
      dxdt(10) = -b11*rxn(1) - b12*rxn(2) - b21*rxn(3) - b22*rxn(4) - dot_product(a3(1:nbs),rxn(5:4+nbs)) - b4*rxn(9) &
               + fLiq*(o2_F - o2)/volL + kLa*(po2*kH_o2 - o2)/volL
              
      ! pch4
      dxdt(11) = ( fGas*(pch4_F - pch4) + kLa*Rpm3*temp*(ch4 - pch4*kH_ch4) )/volG
 
      ! pch3oh
      dxdt(12) = ( fGas*(pch3oh_F - pch3oh) + kLa*Rpm3*temp*(ch3oh - pch3oh*kH_ch3oh) )/volG
    
      ! co2
      dxdt(13) = ( fGas*(pco2_F - pco2) + kLa*Rpm3*temp*(h2co3 - pco2*kH_h2co3) )/volG

      ! o2
      dxdt(14) = ( fGas*(po2_F - po2) + kLa*Rpm3*temp*(o2 - po2*kH_o2) )/volG
      
      ! bs1
      dxdt(15) = epsilon1*(rxn(1) + rxn(2)) - rxn(5) + fLiq*(bs1_F - fLoss*bs1)/volL

      ! bs2
      dxdt(16) = epsilon2*(rxn(3) + rxn(4)) - rxn(6) + fLiq*(bs2_F - fLoss*bs2)/volL
      
      ! bs3
      dxdt(17) = epsilon3*sum( rxn(5:4+nbs) ) - rxn(7) + fLiq*(bs3_F - fLoss*bs3)/volL
      
      ! bs4
      dxdt(18) = epsilon4*rxn(9) - rxn(8) + fLiq*(bs4_F - fLoss*bs4)/volL
      
      ierr = 0 ! no problems in evaluating state eqns.
      return      
   end subroutine stateEqns

   subroutine dummyJac(m,t,y,jac,ldjac,ierr,rpar,ipar)
      ! Dummy Jacobian routine needed for BiM
      implicit none
      integer m, ldjac, ierr, ipar(*)
      real(8) t,y(m),jac,rpar(*)
      return
   end subroutine dummyJac

   subroutine evalReactionSystem (t)
      real(8) t
      ! This routine evaluates all aspects of the reaction system for other routines to use
      ! it assumes that the ocv vector has already been set by a call to the getOCV routine and
      ! that the state variable concentrations are specified.
      
      ! call getDelGF) () since Temp, pH and IS are contant in this model, this is called once in initialize routine
      call getStoichiometry () ! sets the reaction stiochiometric coefficients
      call rxnDelG () ! calculates reaction delta G.
      call reactions () ! calculates reaction rates.
      call calcSigmaDot () ! calculates internal entropy production
      return
   end subroutine evalReactionSystem
   
   subroutine getDelGF0 ()
      use thermoData
      ! For this run, I will assume the pH, Temp and IS do not change, so that
      ! Gibbs free energies of formation only need to be calculated once.  If this
      ! is not so, then these need to be calculated at each iteration, which is
      ! computationally significant.
      ! Note, these data are in kJ/mol (not J)
      delGF0_ch4 = dGf0(methaneaqsp,Temp,is,pH)
      delGF0_ch3oh = dGf0(methanolsp,Temp,is,pH)
      delGF0_h2co3 = dGf0(co2totsp,Temp,is,pH)
      delGF0_dC = dGf0(ch2osp,Temp,is,pH) ! same as CH2O
      delGF0_hno3 = dGf0(nitratesp,Temp,is,pH)
      delGF0_nh3 = dGf0(ammoniasp,Temp,is,pH)
      delGF0_dN = dGf0(ammoniasp,Temp,is,pH) ! same as NH3
      delGF0_o2 = dGf0(o2aqsp,Temp,is,pH)

      delGF0_bs(1:nbs) = dGf0(yeastsp,Temp,is,pH)
      
      delGF0_h2o = dGf0(h2osp,Temp,is,pH)
      
      return
   end subroutine getDelGF0 
   
   subroutine getStoichiometry ()
      integer i
      real(8) bsT
      real(8) sml
      ! this routine solves for the reaction stiochiometric coefficients based on ocv values.
      ! point to Biological structure compositions
      alpha1 => chon(1,1); beta1 => chon(1,2); gamma1 => chon(1,3); 
      alpha2 => chon(2,1); beta2 => chon(2,2); gamma2 => chon(2,3); 
      alpha3 => chon(3,1); beta3 => chon(3,2); gamma3 => chon(3,3); 
      alpha4 => chon(4,1); beta4 => chon(4,2); gamma4 => chon(4,3);  
      
      ! Reactions 1 and 2 catalyzed by bs1 for ch4 uptake and CH3OH production (partial oxidation)
      ! omaga1       [ ch4 + a1 hno3 + b11 o2 -> epsilon1 bs1 + d1 ch3oh + e11 h2o ]
      ! (1 - omaga1) [ ch4 + a1 nh3  + b12 o2 -> epsilon1 bs1 + d1 ch3oh + e12 h2o ]
      ! 
      ! Omega1 determines how bs1 is allocated to r1 versus r2.
      ! elemental balances round each equation gives:
      a1  = epsilon1*gamma1
      d1  = 1. - epsilon1
      b11 = (4. - 5.*a1 - 2.*d1 - alpha1*epsilon1 + 2.*beta1*epsilon1)/4.
      e11 = (4. + a1 - 4.*d1 - alpha1*epsilon1)/2.

      b12 = (4. + 3.*a1 - 2.*d1 - alpha1*epsilon1 + 2.*beta1*epsilon1)/4.
      e12 = (4. + 3.*a1 - 4.*d1 - alpha1*epsilon1)/2.
      
      ! Reactions 3 and 4 catalyzed by bs2 for ch3oh updake
      ! omega2      [ ch3oh + a2 hno3 + b21 o2 -> epsilon2 bs2 + d2 h2co3 + e21 h2o ]
      ! (1 - omega2)[ ch3oh + a2 nh3  + b22 o2 -> epsilon2 bs2 + d2 h2co3 + e22 h2o ]
      a2  = epsilon2*gamma2
      d2  = 1. - epsilon2
      b21 = (6. - epsilon2*(4. + alpha2 - 2.*beta2 + 5.*gamma2))/4.
      e21 = 1. - (epsilon2*(-2. + alpha2 - gamma2))/2.

      b22 = (6. - epsilon2*(4. + alpha2 - 2.*beta2 - 3.*gamma2))/4.
      e22 = 1. - (epsilon2*(-2. + alpha2 - 3.*gamma2))/2.
      
      ! Reactions 5 - 8, descruction of biological structures
      ! bsi + a3i o2 -> espsilon3 bs3 + (1-epsilon3)[gammai(epsilon3 nh3 + (1-epsilon3)dN) + (epsilon3 h2co3 + (1-epsilon3) dC)] + epsilon3(gammai-gamma3)nh3 + b3i h2o
      ! where i =1,4.  For this version of the model, CHON composition of all four biological structures
      ! is the same; so that, a31 = a32, ... a34; b31 = b32 ...b34, but keep code general anyway
      bsT = max(sum(bs(1:nbs,1)), epsilon(bsT))
      do i=1,nbs
         a3(i) = (chon(i,1) - 2.*chon(i,2) - 3.*chon(i,3) + 4.*epsilon3 - alpha3*epsilon3 + 2.*beta3*epsilon3 &
            + 3.*gamma3*epsilon3 - 4.*epsilon3**2)/4.0
         b3(i) = (-2. + chon(i,1) - 3.*chon(i,3) + 2.*epsilon3 - alpha3*epsilon3 + 3.*gamma3*epsilon3)/2.0
         fp3(i) = bs(i,1)/bsT ! used indescriminant bs decomposition
      end do
      
      ! Reaction 9, consumption of dC and dN
      ! dCN + b4 o2 -> espsilon4 bs4 + d4 h2co3 + e4 h2o + f4 NH3
      ! where dCN = dN+dC and [dCN] = [dC], gammadCN = [dN]/[dC]
      sml = epsilon(sml) ! small compared to 1. to avoid divide by zero.
      gammadCN = dN/max(dC,sml)   
      b4 = (4. - epsilon4*(4. + alpha4 - 2.*beta4 - 3.*gamma4))/4.
      d4 = 1. - epsilon4
      e4 = -(epsilon4*(-2. + alpha4 - 3.*gamma4))/2.
      f4 = -(epsilon4*gamma4) + gammadCN
      
      return
   end subroutine getStoichiometry
   
   subroutine rxnDelG ()
      ! Calculates reaction gibbs free energies of reaction at given state conditions.
      ! Note, RkJ is the gas constant in kJ (not J)
      use thermoData, only: RkJ
      integer i
      real(8) delGRconc, xpon1, xpon2, xpon3, xpon4
      real(8), parameter:: ln106 = log(1.0d-6) ! this is used to convert micromolar to molar      
      real(8), parameter:: sml = epsilon(sml) ! small compared to 1. to avoid log of zero.
      
      ! Rxn 1
      delGR(1) = (epsilon1*delGF0_bs(1) + d1*delGF0_ch3oh + e11*delGF0_h2o) &
               - (delGF0_ch4 + a1*delGF0_hno3 + b11*delGF0_o2) ! This is the standard rxn free energy
      ! account for concentrations not at 1 M (state variables here are in microMolar, or mmol/m3)
      !delGRconc = RkJ*temp*log( max( (bs1*1.d-6)**epsilon1 * (ch3oh*1.d-6)**d1, sml ) &
      !                          /max( ch4*1.d-6 * (hno3*1.d-6)**a1 * (o2*1.d-6)**b11, sml ) )
      delGRconc = RkJ*temp*( epsilon1*(log(bs1+sml)+ln106) + d1*(log(ch3oh+sml)+ln106) - (log(ch4+sml)+ln106) - a1*(log(hno3+sml)+ln106) - b11*(log(o2+sml)+ln106) )
      delGR(1) = delGR(1) + delGRconc ! correct for concentration of species
      
      ! Rxn 2
      delGR(2) = (epsilon1*delGF0_bs(1) + d1*delGF0_ch3oh + e12*delGF0_h2o) &
               - (delGF0_ch4 + a1*delGF0_nh3 + b12*delGF0_o2) ! This is the standard rxn free energy
      ! account for concentrations not at 1 M (state variables here are in microMolar, or mmol/m3)
!      delGRconc = RkJ*temp*log( max( (bs1*1.d-6)**epsilon1 * (ch3oh*1.d-6)**d1, sml ) &
!                                /max( ch4*1.d-6 * (nh3*1.d-6)**a1 * (o2*1.d-6)**b12, sml ) )
      delGRconc = RkJ*temp*(epsilon1*(log(bs1+sml)+ln106) + d1*(log(ch3oh+sml)+ln106) - (log(ch4+sml)+ln106) - a1*(log(nh3+sml)+ln106) - b12*(log(o2+sml)+ln106) )
      delGR(2) = delGR(2) + delGRconc ! correct for concentration of species

      ! Rxn 3
      delGR(3) = (epsilon2*delGF0_bs(2) + d2*delGF0_h2co3 + e21*delGF0_h2o) &
               - (delGF0_ch3oh + a2*delGF0_hno3 + b21*delGF0_o2)
      !delGRconc =  RkJ*temp*log( max( (bs2*1.d-6)**epsilon2 * (h2co3*1.d-6)**d2, sml )  &
      !                           /max( ch3oh*1.d-6 * (hno3*1.d-6)**a2 * (o2*1.d-6)**b21, sml ) )
      delGRconc = RkJ*temp*(epsilon2*(log(bs2+sml)+ln106) + d2*(log(h2co3+sml)+ln106) - (log(ch3oh+sml)+ln106) - a2*(log(hno3+sml)+ln106) - b21*(log(o2+sml)+ln106) )
      delGR(3) = delGR(3) + delGRconc ! correct for concentration of species                                 
                                 
      ! Rxn 4
      delGR(4) = (epsilon2*delGF0_bs(2) + d2*delGF0_h2co3 + e22*delGF0_h2o) &
               - (delGF0_ch3oh + a2*delGF0_nh3 + b22*delGF0_o2)
      !delGRconc =  RkJ*temp*log( max( (bs2*1.d-6)**epsilon2 * (h2co3*1.d-6)**d2, sml )  &
      !                           /max( ch3oh*1.d-6 * (nh3*1.d-6)**a2 * (o2*1.d-6)**b22, sml ) )
      delGRconc = RkJ*temp*( epsilon2*(log(bs2+sml)+ln106) + d2*(log(h2co3+sml)+ln106) - (log(ch3oh+sml)+ln106) - a2*(log(nh3+sml)+ln106) - b22*(log(o2+sml)+ln106) )
      delGR(4) = delGR(4) + delGRconc ! correct for concentration of species      
      
      ! Rxns 5 - 8 
      do i=1, nbs
         delGR(4+i) = epsilon3*delGF0_bs(3) + (1.-epsilon3)*( chon(i,3)*(epsilon3*delGF0_nh3 + (1.-epsilon3)*delGF0_dN) &
                  + (epsilon3*delGF0_h2co3 + (1.-epsilon3)*delGF0_dC) ) &
                  + epsilon2*(chon(i,3)-gamma3)*delGF0_nh3 + b3(i)*delGF0_h2o - (delGF0_bs(i) + a3(i)*delGF0_o2)
         xpon1 = (1.-epsilon3)*chon(i,3)*epsilon3 + epsilon3*(chon(i,3) - gamma3) !nh3 coef
         xpon2 = (1.-epsilon3)**2*chon(i,3) ! dN coef
         xpon3 = (1.-epsilon3)*epsilon3 ! h2co3 coef
         xpon4 = (1.-epsilon3)**2 ! dC coef
         !delGRconc = RkJ*temp*log( max( (bs3*1.d-6)**epsilon3 * (nh3*1.d-6)**xpon1 * (dN*1.d-6)**xpon2 * (h2co3*1.d-6)**xpon3 * (dC*1.d-6)**xpon4, sml ) )
         delGRconc = RkJ*temp*( epsilon3*(log(bs3+sml)+ln106) + xpon1*(log(nh3+sml)+ln106) + xpon2*(log(dN+sml)+ln106) + xpon3*(log(h2co3+sml)+ln106) + xpon4*(log(dC+sml)+ln106) &
                   - (log(bs(i,1)+sml)+ln106) - a3(i)*(log(o2+sml)+ln106) )
         delGR(4+i) = delGR(4+i) + delGRconc ! - RkJ*temp*log( max( bs(i,1)*1.d-6 * (o2*1.d-6)**a3(i), sml ) ) ! bs1 conc.
      end do
      
      ! Rxn 9
      delGR(9) = epsilon4*delGF0_bs(4) + d4*delGF0_h2co3 + e4*delGF0_h2o + f4*delGF0_nh3 &
               - (delGF0_dC + gammadCN*delGF0_dN + b4*delGF0_o2 )
      !delGRconc = RkJ*temp*log( max( (bs4*1.d-6)**epsilon4 * (h2co3*1.d-6)**d4 * (nh3*1.d-6)**f4 , sml ) &
      !                         /max( (dC*1.d-6) * (o2*1.d-6)**b4 , sml ) )
      delGRconc = RkJ*temp*( epsilon4*(log(bs4+sml)+ln106) + d4*(log(h2co3+sml)+ln106) + f4*(log(nh3+sml)+ln106) - (log(dC+sml)+ln106) - b4*(log(o2+sml)+ln106) )
      delGR(9) = delGR(9) + delGRconc         
                            
      return
   end subroutine rxnDelG
   
   real(8) function fG(delG) 
      real(8) delG
      ! Local def
      real(8) xpon
      ! define the following function to shut reactions down as delGR approaches 0.
      xpon = 10.
      if (delG >= 0.) then
         fG = 0.0
         return
      else
         fG = 1.0 - exp(xpon*delG)
         return
      end if
   end function fG

   subroutine reactions ()
      ! calculates reaction rates at current state. It is assumed that the state 
      ! and ocv vectors have been set before calling this routine, as well as delGR for each reaction
      
      ! Rxns 1 and 2
      rxn(1) = nu(1)*epsilon1**2*(1.-epsilon1**2)*( ch4/(ch4 + kappa(1)*epsilon1**4) )*( hno3/(hno3 + kappa(1)*epsilon1**4) ) &
             *( o2/(o2 + kappa(1)*epsilon1**4) )*fG(delGR(1))*omega1*bs1
      rxn(2) = nu(1)*epsilon1**2*(1.-epsilon1**2)*( ch4/(ch4 + kappa(1)*epsilon1**4) )*( nh3/(nh3 + kappa(1)*epsilon1**4) ) &
             *( o2/(o2 + kappa(1)*epsilon1**4) )*fG(delGR(2))*(1.-omega1)*bs1
             
      ! Rxns 3 and 4
      rxn(3) = nu(2)*epsilon2**2*(1.-epsilon2**2)*( ch3oh/(ch3oh + kappa(2)*epsilon2**4) )*( hno3/(hno3 + kappa(2)*epsilon2**4) ) &
             *( o2/(o2 + kappa(2)*epsilon2**4) )*fG(delGR(3))*omega1*bs2
      rxn(4) = nu(2)*epsilon2**2*(1.-epsilon2**2)*( ch3oh/(ch3oh + kappa(2)*epsilon2**4) )*( nh3/(nh3 + kappa(2)*epsilon2**4) ) &
             *( o2/(o2 + kappa(2)*epsilon2**4) )*fG(delGR(4))*(1.-omega1)*bs2
             
      ! Rxns 5-8
      rxn(5) = nu(3)*epsilon3**2*(1.-epsilon3**2)*( bs1/(bs1 + kappa(3)*epsilon3**4) )*( o2/(o2 + kappa(3)*epsilon3**4) ) &
             *fG(delGR(5))*fp3(1)*bs3
      rxn(6) = nu(3)*epsilon3**2*(1.-epsilon3**2)*( bs2/(bs2 + kappa(3)*epsilon3**4) )*( o2/(o2 + kappa(3)*epsilon3**4) ) &
             *fG(delGR(6))*fp3(2)*bs3
      rxn(7) = nu(3)*epsilon3**2*(1.-epsilon3**2)*( bs3/(bs3 + kappa(3)*epsilon3**4) )*( o2/(o2 + kappa(3)*epsilon3**4) ) &
             *fG(delGR(7))*fp3(3)*bs3
      rxn(8) = nu(3)*epsilon3**2*(1.-epsilon3**2)*( bs4/(bs4 + kappa(3)*epsilon3**4) )*( o2/(o2 + kappa(3)*epsilon3**4) ) &
             *fG(delGR(8))*fp3(4)*bs3
             
      ! Rxn 9
      rxn(9) = nu(4)*epsilon4**2*(1.-epsilon4**2)*( dC/(dC + kappa(4)*epsilon4**4) ) &
             *( o2/(o2 + kappa(4)*epsilon4**4) )*fG(delGR(9))*bs4
     
      return
   end subroutine reactions  
       
   subroutine calcSigmaDot ()
      ! this calculates the internal entropy production rate, sigmaDot.  Since destruction of chemical potential
      ! dominates, do not bother calculating entropy of mixing.
      sigmaDot = - volL*dot_product(rxn(1:nrxn),delGR(1:nrxn))/temp ! units: J/(d K)
      return
   end subroutine calcSigmaDot     
end module reactionSystem

program methanotroph_BOBYQA_2_V2
   ! This is the first version of the methanotrophic distributed metabolic
   ! network based on whole reactions (not half reactions) and using BOBYQA.  
   ! This is the serial version of the code
   !  
   ! Versions:
   !  1.0   10 Mar 2012 -  This version is based on methanotroph_BOBYQA_V7, but the equation for complete CH4
   !                       oxidation has been removed.  Rxns 11 and 12, that were the partial oxidation of CH4 to CH3OH
   !                       have been moved up to Rxns 1 and 2.  Total number of rxns in now 9, and only 4 BS are needed.
   !  2.0   11 Mar 2012 -  This version differs from V1 in that CH3OH is not transported out of the reactor, but
   !                       can act as a storage compound that is not subject to dilution or gas loss.
   !  3.0   13 Mar 2012 -  Adding parameter that acts to reduce POM from leaving the chemostat.  Simple test to simulate
   !                       some aspects of the biofilm and flocs that develop w/o adding more state var.
   !  4.0   14 Mar 2012 -  This is same as V3, but rxn 9 has been modified so that dN is taken up with dC, and excess N is
   !                       remineralized as NH3.  This was done because dN was accumulating over time, which is not consistent
   !                       with observations.
   !  4.1   17 Mar 2012 -  Change calculation of delGR in rxnDelG to reduce the number of exponents, as VTune showed high computation load.
   !                       Simulation with V4_Int10-20_2_cntl reduce CPU time from 45 min to 30 min, so significant speed increase.
   
   use globalVars ! Global variables.
   use reactionSystem
   implicit none
   integer mpierr

   ! BOBYQA declarations
   integer nocvIntv ! number of optimization variables
   real(8) objFcnValue
   
   ! Other declarations
   integer narg, intv, i, iflag
   character char80*80, fmtstr*80
   integer ierr, ipar(1) ! these are just dummy variables for call to stateEqns
   real(8) rpar(1), xdot(nsv) ! these are just dummy variables for call to stateEqns
   integer(8) startT, endT, startTi, endTi, clockCnt ! used for timing runs
   real(8) aveSigmaDotTinfinity, aveSigmaWdotTinfinity
   integer nScreenCNT
   character sDate*9, sTime*8
   
   ! Nothing allocated or open yet
   allocateVarsCalled = .false.
   unitsOpen = .false. 
    
   ! read in parameter file.  First try getting filename from commandline
   ! if that fails, then have processes zero get name and send it along 
   ! to the other processes.
   narg = command_argument_count () ! see if a filename has been included in command line
   if (narg == 1) then
      call get_command_argument(narg,value=basefilename)
   else 
      write(6,'(a,$)') 'Enter parameter filename base, no extension: '
      read(5,*) basefilename
   end if
   
   call point2State () ! this associates state variable names to vectors, must be called before PRM read.   

   ! readin parameters file basefilename.prm
   call readPRM (mpierr)
   ! check for error in opening the PRM file
   if (mpierr /= 0) then
      ! an error occured reading the prm file by one or more processes
      write(6,'(/a/)') 'Error opening or reading '//trim(basefilename)//'.prm file. Aborting.'
      call finishIt
   end if
   
   ! allocate space for problem
   nocvIntv = nocv*nknots ! total number of optimal control variables over interval
   call allocateVars () 
   
   ! initialize state at starting time, t0 and other things
   call initialize()    
   call initializeRS () ! initializes the reactive system module.
   call feedControl (ti) ! set the feed composition
         
   ! Open output file (only process 0 actually does this, all others are kicked back).
   call openOutputs ()
   
   ! setup some timing info
   call system_clock(startT, clockCnt)
   
   ! Save initial conditions
   ! for splineFit tknots must be set, but it can be any monotonically increasing
   ! function, as ocvMat has all the same values at this point
   do i=1,nKnots
      tknots(i+1) = tknots(i) + 1.
   end do
   !call splineFit () ! fit splines to data in ocvMat (this is not needed in Ver 7)
   call stateEqns (nsv,ti,state0,xdot,ierr,rpar,ipar) 
   call saveSolnPoint (ti)
   call date(sDate); call time(sTime)
   write(6,'(a)') 'Version: '//trim(version)
   write(6,'(a)') 'Start Time: '//sDate//' '//sTime
   write(6,'(a)') ' File: '//trim(basefilename)
   ! *** Main loop to find OCV values over interval.  Loop to cover desired simulation time
   nScreenCNT = 0
   do intv = 1, nIntervals ! number of intervals to run
      BiMretry = 0; BiMfail = 0 ! track BiM integration problems for each process over the inverval only.      
      call system_clock(startTi)
      tii = ti + tInfinity ! end time (d) of current interval.
      ! set times where spline knots will be set.  Use exponential spacing if wInfinity is <1  
      call spaceKnots ()

      ! Specify ODE integration outpout time points. Since points are not being save to file at
      ! this point, just specify only the end point. In other models that are sensitive to small
      ! integration errors, changing ODE output points could produce a different result.  For now,
      ! assume the models is robust enought that the ODE solution (and the optimization) do
      ! not depend on the ODE integration output points, as this should speed up integration.
      saveSoln = .false.
      nODEpts = 1
      solTpts(1) = tii 
         
      ! Call BOBYQA
      fEval = 0
      call BOBYQA (nocvIntv,npt,u,bl,bu,rhobeg,rhoend,iprint,maxfun,w)
      ! BOBYQA does not provide any variables regarding output, so nothing to store.
      call calfun (nocvIntv,u,objFcnValue) ! insures system evaluated at ocv solution
      aveSigmaDotTinfinity = aveSigmaDot  ! values at end of "infinite" interval
      aveSigmaWdotTinfinity = aveSigmaWdot
      
      ! Reset the interval to tIntv
      tii = ti + tIntv
      ! Specify the ODE and other variable output times
      nODEpts = nSolTpts
      do i=1,nODEpts
         solTpts(i) = ti + real(i)*(tii - ti)/real(nODEpts)
      end do
      ! with saveSoln set to true, this stores all solution outputs. The solution at tii is not saved here
      ! because it could differ slightly from that during the optimum search.  Very small changes in the
      ! initial conditions for the next interval can cause large differences.
      saveSoln = .true.
      call calfun (nocvIntv,u,objFcnValue)      
      call system_clock(endTi) 
      ! Save sigma data at interval end only
       write(TSINVunit,'(f10.4,1x,f6.2,3(1x,i7),4(1x,g13.5))') tii, real(endTi-startTi)/real(clockCnt)/60.0, fEval, BiMretry, BiMfail, aveSigmaDot, aveSigmaDotTinfinity, aveSigmaWdot, aveSigmaWdotTinfinity
      
      ! copy final state for IC for next interval, zero out sigma.
      state0(1) = sigma; state0(2) = sigmaW; state0(nSys+1:nc+nSys) = c(1:nc,1); state0(nc+nSys+1:nSys+nc+nbs) = bs(1:nbs,1) 
      ! copy final control state to IC for next interval, which require interpolation
      ! because tIntv may not be fall on a knot within tInfinity
      call getOCV (tii)
      ocvMat(1,1:nocv) = ocv(1:nocv,1)
      if (mod(nScreenCNT,80) == 0) write(6,'(/a)') 'Time (d) Intv  <SigDot>Inv  <SigDot>Inf    fEval  CPU (min) BiM-Retry BiM Fails'
      write(6,'(1x,f6.2,2x,i5,2x,g11.5,2x,g11.5,1x,i7,3x,f6.2,5x,i6,4x,i6)') tii, intv, aveSigmaDot, aveSigmaDotTinfinity, fEval, real(endTi-startTi)/real(clockCnt)/60.0, BiMretry, BiMfail
      nScreenCNT = nScreenCNT + 1
      
      ! update initial time for next interval.
      tknots(1) = tii
      ti = tii
   end do
   
   ! *** end main loop
   call system_clock(endT)
   
   write(6,'(/a,f8.2,a)') 'Total time: ', real(endT-startT)/60.d0/real(clockCnt), ' (min).'
   write(6,'(a)') 'Finished simulation'
   
   ! deallocate and finish
   call finishIt ()
end program methanotroph_BOBYQA_2_V2

subroutine readPRM (ioE)
   ! This routine reads the parameters in the namelist stored in basefilename.prm
   use globalVars
   implicit none
   integer ioE
   ! local declarations
   integer PRMunit, ioerr, mpierr
   
   ! open a readonly unit to the prm file. 
   open(newunit=PRMunit,file=trim(basefilename)//'.prm',action='read',status='old',iostat=ioerr)   
   if (ioerr /= 0) then 
      write(6,'(a,i5,a)') 'Could not open PRM file' 
   else
      read(PRMunit,nml=params)
      close(unit=PRMunit)
   end if
   
   ioE = abs(ioerr) ! remove negative sign, if present
   return
end subroutine readPRM

subroutine openOutputs()
   ! This file opens the output files and write a Tecplot header to each.
   use globalVars
   implicit none

   ! local declarations
   integer i
   character strg*1000, tmp*80

   ! Begin body
   ! Output for state variables, C
   open(newunit=TCunit,file=trim(basefilename)//'.tc')
   ! write Tecplot header
   strg = 'Variables = "Time (d)" '
   do i=1,nc
      strg = trim(strg)//' "'//trim(cnames(i))//'"'
   end do
   write(TCunit,'(a)') strg
   write(TCunit,'(a)') 'Zone T="State Var. C" '

   ! Output for state variables, BS
   open(newunit=TBSunit,file=trim(basefilename)//'.tbs')
   ! write Tecplot header
   strg = 'Variables = "Time (d)"'
   do i=1,nbs
     strg = trim(strg)//' "'//trim(bsnames(i))//'"'
   end do
   write(TBSunit,'(a)') strg
   write(TBSunit,'(a)') 'Zone T="State Var. BS" '

   ! Optimal control variables   
   open(newunit=TOCVunit,file=trim(basefilename)//'.tocv')
   ! write Tecplot header
   strg = 'Variables = "Time (d)"'
   do i=1,nocv
     strg = trim(strg)//' "'//trim(ocvnames(i))//'"'
   end do
   write(TOCVunit,'(a)') strg
   write(TOCVunit,'(a)') 'Zone T="OCV" '

   ! Output for rxn rate variables
   open(newunit=TRunit,file=trim(basefilename)//'.tr')
   ! write Tecplot header
   strg = 'Variables = "Time (d)"'
   do i=1,nrxn
     strg = trim(strg)//' "'//trim(rxnnames(i))//'"'
   end do
   write(TRunit,'(a)') strg
   write(TRunit,'(a)') 'Zone T="Rxn Rates" '

   ! Output for rxn delta G variables, use reaction names
   open(newunit=TGunit,file=trim(basefilename)//'.tG')
   ! write Tecplot header
   strg = 'Variables = "Time (d)"'
   do i=1,nrxn
     strg = trim(strg)//' "'//trim(rxnnames(i))//'"'
   end do
   write(TGunit,'(a)') strg
   write(TGunit,'(a)') 'Zone T="Rxn delG" '
   
   ! Output for Entropy Production rate at each grid point.
   open(newunit=TSunit,file=trim(basefilename)//'.tS')
   ! write Tecplot header
   strg = 'Variables = "Time (d)" "sigma" "sigmaDot" "sigmaW" "sigmaWdot"'
   write(TSunit,'(a)') strg
   write(TSunit,'(a)') 'Zone T="<<greek>sigma</greek>> (t)"'

   ! Output for Entropy Production rate for output at interval End.
   open(newunit=TSINVunit,file=trim(basefilename)//'.tSinv')
   ! write Tecplot header
   strg = 'Variables = "Time (d)" "CPU (min)" "fEval" "BiM Retries" "BiM Fails" "aveSigmaDot" "aveSigmaDotTinfinity" "aveSigmaWdot" "aveSigmaWdotTinfinity"'
   write(TSINVunit,'(a)') strg
   write(TSINVunit,'(a)') 'Zone T="<<greek>sigma</greek>>Intv (t)"'
   
   ! open unit for debuging
   open(newunit=DEBUGunit,file=trim(basefilename)//'.debug')
   ! write Tecplot header
   strg = 'Variables = "Time (d)" "N total (mmol)" "C total (mmol)" '
   write(DEBUGunit,'(a)') strg
   write(DEBUGunit,'(a)') 'Zone T="Total N and C"'

   unitsOpen = .true. 
   return
end subroutine openOutputs

subroutine closeOutputs()
   ! This routine closes all output files, return if not process 0   
   use globalVars
   implicit none

   close(unit=TCunit)
   close(unit=TBSunit)
   close(unit=TOCVunit)
   close(unit=TRunit)
   close(unit=TGunit)
   close(unit=TSunit)
   close(unit=TSINVunit)
   close(unit=DEBUGunit)
   unitsOpen = .false.
   return
end subroutine closeOutputs

subroutine calfun (nBOB,uu,J)
   ! This is an interface routine to SSsoln for BOBYQA.
   use globalVars
   implicit none
   ! Dummy variables.
   integer nBOB
   real(8) uu(nBOB), J
   
   ! Local variables.
   integer i
   integer infeasibleSoln

   ! First copy the optimal control variables to a matrix.
   ! Note, the first ROW of ocv contains the value of the ocv from the previous interval (i.e., last row).
   ! the time points for the interval should have already been set in tknots prior to the VTDIRECT call.
   do i=1,nknots
      ocvMat(i+1,1:nocv) = uu(1+(i-1)*nocv:i*nocv)
   end do
   
   ! fit the SPLINE functions to the control variable being manipulated by DIRECT
   ! which covers the optimization interval (ti, tii)
   !call splineFit () ! fit splines to data in ocvMat (Not needed in Ver 7)
   call integrateState (infeasibleSoln) ! note, infeasibleSoln not used for BOBYQA 
   
   !J = -aveSigmaDot ! average sigmaDot (take negative, because BOBYQA finds minimum)
   ! base optimum on weighted (discounted) sigmaDot
   J = -aveSigmaWdot ! average sigmaWdot (take negative, because BOBYQA finds minimum)
   fEval = fEval + 1
   if (fEval > maxfun+2) write(6,*) 'fEval > maxfun ?'
   return
end subroutine calfun

subroutine getOCV(t)
   use globalVars, only: ocv, ocvMat, nocv, nknots, tknots
   implicit none
   real(8) t
   ! This routine extract u from umat at time t.  Linear interpolation is used. It is a replacement for TSPACK.
   ! Note dimensions: ocvMat(1:nknots+1,1:nocv), ocv(1:nocv,1:3)
   ! Local declarations
   integer nlow
   real(8) delt
       
   ! get location of t within tknots
   call interp1D (t, nknots+1, tknots, nlow, delt)
   ocv(1:nocv,1) = ocvMat(nlow,1:nocv) + delt*(ocvMat(nlow+1,1:nocv) - ocvMat(nlow,1:nocv) )
   ! Assume optimization routines never pass ocv to ocvMat that are beyond ocv(1:nocv,2) and ocv(1:nocv,3)
   !forall(i=1:nocv, ocv(i,1) < ocv(i,2)) ocv(i,1) = ocv(i,2) ! lower bound constraint 
   !forall(i=1:nocv, ocv(i,1) > ocv(i,3)) ocv(i,1) = ocv(i,3) ! upper bound constraint    
   return
end subroutine getOCV

Subroutine interp1D (t, ndim_t, tvec, nlow, delt)
    Implicit None              
    Integer ndim_t, nlow
    real(8) t, tvec(ndim_t), delt
    !     This routine determines nlow such that:
    !        tvec(nlow) <= t <= tvec(nlow+1)
    !     Note, if t = tvec(ndim_t), then nlow is set to ndim_t-1 
    !     delt is given by:
    !        delt = [t - tvec(nlow)]/[tvec(nlow+1)-tvec(nlow)]
    !     that can be used for linear interpolation

    !     Input
    !        t        time that is to be found in tvec
    !        ndim_t   Number of elements USED in tvec (not it's declared dimension)
    !        tvec     Vector containing increasing values of t
    !
    !     Output
    !        nlow     Value such that: tvec(nlow) <= t <= tvec(nlow+1)
    !        delt     value such that: [t - tvec(nlow)]/[tvec(nlow+1)-tvec(nlow)]

    !     Local Declarations:
    ! Check to see if t is < tvec(1) or >= tvec(ndim_t)
    if (t < tvec(1)) then
        nlow = 1
        delt = 0.d0
        return
    end if
    if (t >= tvec(ndim_t)) then
        nlow = ndim_t - 1
        delt = 1.0d0
        return
    end if
    nlow = maxloc(tvec, dim = 1, mask = tvec .le. t)
    delt = (t - tvec(nlow))/( tvec(nlow+1)-tvec(nlow) )
    return
end subroutine interp1D

subroutine transformOCV ()
   ! This routine applies any necessary transformations to the OCV
   ! This routine is **PROBLEM SPECIFIC**
   use globalVars
   implicit none
   real(8) denom
   
   ! For this version of the model not transformations are needed.

   return
end subroutine transformOCV

subroutine integrateState (iflag)
   use globalVars
   use reactionSystem
   implicit none
   integer iflag ! set to 1 if infeasible, otherwise 0.
   ! local declarations
   integer i
   logical repeat
   integer attempts
   ! Declarations needed for BiM
   integer, parameter:: lenw = 14+10+8*nsv+4*10*nsv+2*nsv**2, leniw = nsv+37
   integer  ierflg, ijac, mljac, mujac, iout, idid, ierr
   real(8) tcurr, tend, rtolL, atolL, h
   real(8) wk(lenw), x(nsv), xdot(nsv)
   integer iwk(leniw)
   real(8) rpar(1)
   integer ipar(1)
   ! External stateEqns, dummyJac - don't need to specify because of use reactionSystem statement
   
   iflag = 0 ! if set to 1, then ocv is an infeasible point.
   
   ! initialize state values at ti, which needs to be set from end of previous interval.
   x(1:nsv) = state0(1:nsv)
   
   ! Begin integration loop.  solution will only be saved if saveSoln is set to true
   tcurr = ti
   do i=1, nODEpts
      tend = solTpts(i)           
      repeat = .True.
      attempts = 0
      rtolL = rtol
      atolL = atol
      Do While (repeat)
         repeat = .False. 
         h = h0 ! initial BiM stepsize
         ijac = 0 ! BiM: 0=calculate jacobian numerically, otherwise use analytical.
         iout = 0 ! BiM: set to 1 to call solout
         mljac = nsv; mujac = nsv !BiM: banding aspect of jacobian? just set to nsv.
         iwk(1:8) = 0
         wk(1:14) = 0.0d0
         iwk(1) = BiMmaxIter ! default iterations: 100000
         !iwk(2) = 10 ! ORDMIN, 4<=ORDMIN<=12. (DEFAULT = 4).
         wk(1) = epsilon(wk) ! set true precision
         wk(2) = hmax !  HMAX. MAXIMUM INTEGRATION STEP. (DEFAULT = (TEND-T0)/8)
         idid = 0
         call BiM (nsv,stateEqns,tcurr,tend,x,h,rtolL,atolL,dummyJac,ijac,mljac,mujac, &
                   wk,lenw,iwk,leniw,rpar,ipar,iout,idid)
         ! Check for integration errors
         if (idid < 0) Then 
            if (idid == -1) then
               write(6,'(a)') 'Wrong input params for BiM; stopping.'
               call finishIt () ! this is not the way to properly handle this, but this error should not occur.
            end if
            attempts = attempts + 1
            !Write(6,'(a)') 'integrateState:: integration error, retrying.'
            rtolL = rtolL*10.
            atolL = atolL*10.
            repeat = .True.
            if (attempts > maxAttempts) then
               !Write(6,'(a,f9.3,a)') 'integrateState:: integration FAILURE at t = ',tcurr, &
               !   ', returning.'
               iflag = 1 ! mark point as infeasible.
               aveSigmaDot = 0.0
               aveSigmaWdot = 0.0
               BiMfail = BiMfail + 1
               return
            end if
            BiMretry = BiMretry + 1 
         end If
      end do
      ! save solution (if saveSoln is true, and only process 0 will do this).
      ! tcurr will be at tend value, unless error in BiM occured.
      sigma  = x(1)
      sigmaW = x(2)
      aveSigmaDot  = (sigma -state0(1))/(tcurr - ti) ! calculate the average sigmaDot value at current time (J/d/deg-K)
      aveSigmaWdot = (sigmaW-state0(2))/(tcurr - ti) ! calculate the average sigmaWdot value at current time (J/d/deg-K)
      if (saveSoln) then
         ! Since BiM may have interpolated to tend, make sure reaction system is evaluated at current time
         call stateEqns (nsv,tcurr,x,xdot,ierr,rpar,ipar) 
         call saveSolnPoint (tcurr)
      end if
   end do
   return
end subroutine integrateState

subroutine solout(m,t,y,f,k,ord,irtrn)   
   ! needed by BiM, but not used in this case (dummy routine)
   implicit none
   integer m,k,ord,irtrn
   real(8) t(*),y(m,*),f(m,*)
   return
end subroutine solout

subroutine saveSolnPoint (t)
   ! Saves solution at current time point.
   ! This assumes the state has been update at time t and that
   ! optimal control variables, reactions, etc have been defined.     
   use globalVars
   use reactionSystem
   use thermoData, only: Rpm3 !gas constant in m3/mmol
   implicit none
   real(8) t
   integer, parameter:: nTotal = 2
   real(8) total(nTotal)
   real(8) weightS
      
   call writeVec(TCunit,t,nc,c(1:nc,1))     ! state variables, C
   call writeVec(TBSunit,t,nbs,bs(1:nbs,1))  ! biological structure
   call writeVec(TOCVunit,t,nocv,ocv(1:nocv,1))  ! optimal control variables
   call writeVec(TRunit,t,nrxn,rxn(1:nrxn))  ! reaction rates
   call writeVec(TGunit,t,nrxn,delGR(1:nrxn))  ! reaction free energies
   
   ! Output for Entropy Production rate at each grid point.
   write(TSunit,'(f10.4,4(1x,g13.5))') t, sigma, sigmaDot, sigmaW, sigmaDot*weightS(t)
   
   ! debug output
   ! C and N balance
   total(1) = volL*( hno3+nh3+dN+dot_product(chon(1:nbs,3),bs(1:nbs,1)) )
   total(2) = volL*( ch4+ch3oh+h2co3+dC+sum(bs(1:nbs,1)) ) + volG*(pch4+pch3oh+pco2)/(temp*Rpm3)
   call writeVec(DEBUGunit,t,nTotal,total(1:nTotal)) 

   return
end subroutine saveSolnPoint

subroutine writeVec(iunit,t,nx,x)
    ! This routine write a vector to an already open unit
    implicit none
    integer iunit, nx
    real(8) t, x(nx)
    ! Local declarations
    integer i
    character fmt*80
    
    ! Begin Body
    fmt = '(f10.4,##(1x,g13.5))'
    write(fmt(8:9),'(i2)') nx
    write(iunit,fmt) t,(x(i),i=1,nx)
    return
end subroutine writeVec

real(8) function weightS (t)
   ! weighting (or discouting) function on sigma
   use globalVars, only: wInfinity, tInfinity, ti
   implicit none
   real(8) t
   real(8) kW
   kW = -log(wInfinity)/tInfinity
   weightS = exp(-kW*(t-ti))
   return
end function weightS

subroutine spaceKnots ()
   ! this routine spaces the splne knots exponentially; however, if the
   ! sigma weighting function is close to linear, knots are spaced evently.
   ! NOTE tknot(1) must be set prior to calling this routine, and it should equal ti 
   use globalVars
   implicit none
   
   integer i
   real(8) weightS, wS
   
   if (wInfinity >= 0.999) then
      ! just space the knots equally.
      do i=1,nKnots
         tknots(i+1) = tknots(i) + (tii - ti)/real(nknots)
      end do      
      return
   end if
   ! Space knots exponentially (see weightS (t), as this is its inverse)
   wS = 1.0d0
   do i=1,nKnots-1
      wS = wS - (1.d0 - wInfinity)/real(nknots)
      tknots(i+1) = tInfinity*log(wS)/log(wInfinity) + tknots(1)
   end do 
   tknots(nKnots+1) = tknots(1) + tInfinity !just avoids any rounding errors in above expression  
   return
end subroutine spaceKnots

subroutine feedControl (t)
   ! This routine sets which feed to use based on time and parameters set in the input file (*.prm)
   ! This version changes between feed 1 and 2 only. 
   ! The default feed (1) is in c(*,2), the alternate feed (2) is in c(*,3).  
   ! stepFcn = 0 selects for feed 1, while stepFcn = 1 selects for feed 2. 
   use globalVars
   implicit none
   real(8) t
   
   !real(8), parameter:: sig = 10.d0 ! This is not a good value, just for testing.
   real(8), parameter:: sig = 50.d0 ! this will cause a change to start to occur ~ 144 min before step
   !real(8), parameter:: sig = 100.d0 ! this will cause a change to start to occur ~ 72 min before step
   !real(8), parameter:: sig = 300.d0 ! this will cause a change to start to occur ~ 15 min before step
   real(8) tElap, tPhase, tDisc, stepFcn, fOn, fOff, sign 
   
   ! While turning the feed on or off could be done simply, it would introduce
   ! discontinuities in the first derivative of the feed concentration.  This can
   ! cause problems for integration when steps occur.  Consequently, use a 
   ! continuous function to do the changing.
   
   ! First get the step function setup based "centered" on tFeedCyclingOn
   tElap = t - tFeedCyclingOn ! time since tFeedCyclingOn
   tPhase = tElap/tFeed2OnPeriod
   ! determine closest step (up or down)
   tDisc = dble(nint(tphase))
   sign = -1.0d0
   if (modulo(int(tdisc),2) == 1) sign = +1.0d0    
   stepFcn = 1.0d0/( 1.0d0 + exp(sign*sig*(tphase-tdisc)) )
   
   ! now place bounds on the interval over which stepping can occur.
   fOn  = 1.0d0/(1.0d0 + exp(-sig*(t-tFeedCyclingOn )) )
   fOff = 1.0d0/(1.0d0 + exp( sig*(t-tFeedCyclingOff)) )
   stepFcn = stepFcn*fOn*fOff
   
   ! set the feed concentrations
   ch4_F    = c(1,2)  + stepFcn*( c(1,3)  - c(1,2) )
   ch3oh_F  = c(2,2)  + stepFcn*( c(2,3)  - c(2,2) )
   h2co3_F  = c(3,2)  + stepFcn*( c(3,3)  - c(3,2) )
   dC_F     = c(4,2)  + stepFcn*( c(4,3)  - c(4,2) )
   hno3_F   = c(5,2)  + stepFcn*( c(5,3)  - c(5,2) )
   nh3_F    = c(6,2)  + stepFcn*( c(6,3)  - c(6,2) )
   dN_F     = c(7,2)  + stepFcn*( c(7,3)  - c(7,2) )
   o2_F     = c(8,2)  + stepFcn*( c(8,3)  - c(8,2) )
   pch4_F   = c(9,2)  + stepFcn*( c(9,3)  - c(9,2) )
   pch3oh_F = c(10,2) + stepFcn*( c(10,3) - c(10,2) )
   pco2_F   = c(11,2) + stepFcn*( c(11,3) - c(11,2) )
   po2_F    = c(12,2) + stepFcn*( c(12,3) - c(12,2) )

   bs1_F = bs(1,2) + stepFcn*( bs(1,3) - bs(1,2) )
   bs2_F = bs(2,2) + stepFcn*( bs(2,3) - bs(2,2) )
   bs3_F = bs(3,2) + stepFcn*( bs(3,3) - bs(3,2) )
   bs4_F = bs(4,2) + stepFcn*( bs(4,3) - bs(4,2) )      
   
   return
end subroutine feedControl 

subroutine finishIt ()
   ! This routine closes open units, deallocates variables then stops.
   use globalVars
   implicit none
   integer mpierr
   
   if (allocateVarsCalled) call deallocateVars ()
   if (unitsOpen) call closeOutputs () ! close open units
   stop
end subroutine finishIt