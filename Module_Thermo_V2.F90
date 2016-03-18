module thermoData 
    ! This thermodynamic data is from Alberty2003 and his Mathematica 
    ! BiochemThermo`BasicBiochemData3.m data.
    ! For bacterial structure, have used data on yeast from Battley1998.
    ! Energy units are in kJ/mole
    implicit none
    type GHdata
        real(8) dGf
        real(8) dHf
        real(8) zi
        real(8) nH
    end type GHdata
    real(8), parameter:: Rp = 0.0820574587 ! gas constant (atm L/(g-mol K))
    real(8), parameter:: Rpm3 = 0.0820574587d-6 ! gas constant (atm m^3/(g-mmol K)) (note, mmol not mol)
    real(8), parameter:: RkJ = 8.314472d-3 ! gas constant (kJ/(g-mol K))
   
    ! Note, in order for all quantities to be treated as arrays, it is necessary to dimension species
    ! even if they only have one speces, such as h2osp(1).
    type(GHdata) atpsp(3), methaneaqsp(1), adpsp(3), pisp(2), h2osp(1), o2aqsp(1), yeastsp(1), ch2osp(1)
    type(GHdata) nadoxsp(1), nadredsp(1), sulfatesp(2), ammoniasp(2), co2gsp(1), co2totsp(3), nitratesp(2)
    type (GHdata) methanegsp(1), nitritesp(2), n2aqsp(1), h2saqsp(3), h2sgsp(1), methanolsp(1), methanolgsp(1)
    type (GHdata) o2gsp(1)
    
    data atpsp /GHdata(-2768.1,-3619.21,-4,12),GHdata(-2811.48,-3612.91,-3,13),GHdata(-2838.18,-3627.91,-2,14)/
    data methaneaqsp /GHdata(-34.33,-89.04,0,4)/
    data methanegsp /GHdata(-50.72,-74.81,0,4)/
    data methanolsp /GHdata(-175.31,-245.93,0,4)/
    data methanolgsp /GHdata(-162.3,-201.0,0,4)/ ! gas phase data from CRC Handbook of Chemistry and Physics, Vol 91.
    data adpsp /GHdata(-1906.13,-2626.54,-3,12),GHdata(-1947.1,-2620.94,-2,13),GHdata(-1971.98,-2638.54,-1,14)/
    data pisp /GHdata(-1096.1,-1299.,-2,1),GHdata(-1137.3,-1302.6,-1,2)/
    data h2osp /GHdata(-237.19,-285.83,0,2)/  
    data o2gsp /GHdata(0.,0.,0.,0.)/
    data o2aqsp /GHdata(16.4,-11.7,0,0)/
    data yeastsp /GHdata(-87.96,-133.09,0,1.613)/ !Battley1998: CH1.613 O0.557 N0.158 P0.012 S0.003 K0.022 Mg0.003 Ca0.001
    data ch2osp /GHdata(-152.65,-210.365,0,2)/
    data nadoxsp /GHdata(0,0,-1,26)/
    data nadredsp /GHdata(22.65,-31.94,-2,27)/
    data sulfatesp /GHdata(-744.53,-909.27,-2,0),GHdata(-755.91,-887.34,-1,1)/
    data ammoniasp /GHdata(-26.5,-80.29,0,3),GHdata(-79.31,-132.51,1,4)/
    data co2gsp /GHdata(-394.36,-393.5,0,0)/
    data co2totsp /GHdata(-527.81,-677.14,-2,0),GHdata(-586.77,-691.99,-1,1),GHdata(-623.11,-699.63,0,2)/
    data nitratesp /GHdata(-108.74,-205.,-1,0),GHdata(-111.25,-207.36,0,1)/
    data nitritesp /GHdata(-32.2,-104.6,-1,0), GHdata(-50.6,-119.2,0,1)/
    data n2aqsp /GHdata(18.7,-10.54,0,0)/
    data h2saqsp /GHdata(85.8,33.1,-2,0),GHdata(12.08,-17.6,-1,1),GHdata(-27.83,-39.7,0,2)/
    ! These data for H2S in g are from Thauer1977 (for the dGf), and also
    ! from the book "Chemical Principles" by P. Atkins and L. Jones, gotten from google books for dHf.
    ! These data should also be in "The NBS tables of chemical thermodynamic properties : 
    ! selected values for inorganic and C1 and C2 organic substances in SI units" by Donald D Wagman; et al.
    data h2sgsp /GHdata(-33.56,-20.63,0,2)/
    
    contains
    
    real(8) function dGf0(sp,T,is,pH)
        ! This function returns the Gibbs energy of formation for the compound sp at temperature T (K),
        ! ionic strength, is (M), and pH in kJ/mol.  This is based on Alberty2003.
        !   sp:  contains, dGf, dHf, zi and nH at 298.15 and zero is and pH.
        ! 
        ! 30 Mar 2006: First written; Joe Vallino
        implicit none
        type (GHdata) sp(:)
        real(8) T, is, pH
    
        ! Local declarations
        integer i, nsp
        real(8) gibbscoeff, dGzeroT, pHterm, istermG, gpfnsp, R, Tref
        integer, parameter:: maxsp=10
        real(8) eterms(maxsp), maxet
        
        !**** Begin Program *****
        nsp = size(sp)
        if (nsp > maxsp) then
            write(6,*) 'Error::dGf0: a compund has too many species, increase maxsp'
            stop
        end if
        R = 8.31451d0 ! (J/mol/K)
        Tref = 298.15d0
        do i=1,nsp  ! sum for each species of the compound sp
            gibbscoeff=(9.20483d0*T)/1.D3-(1.284668d0*T**2)/1.D5+(4.95199d0*T**3)/1.D8 
            dGzeroT=(sp(i)%dGf*T)/Tref+sp(i)%dHf*(1-T/Tref)
            pHterm=sp(i)%nH*R*T*log(10.d0**(-pH))/1000.d0
            istermG=gibbscoeff*(sp(i)%zi**2-sp(i)%nH)*sqrt(is)/(1.d0+1.6d0*sqrt(is));
            eterms(i)=-1.d3*(dGzeroT-pHterm-istermG)/(R*T)
        end do
        ! Need to bring the exponents down to values that can be calculated
        ! without overflow (underflow ok), so subtract the largest one.
        maxet = maxval(eterms(1:nsp))
        dGf0 = 0.d0
        do i=1,nsp
            dGf0 = dGf0 + exp(eterms(i) - maxet)
        end do
        dGf0 = -R*T*(log(dGf0) + maxet)/1.d3
        return
    end function dGf0
    
   real(8) function solCH4(T, is, pH)
     ! This function returns the solubiility of methane in water (uM/atm) at temperature T (K),
     ! pH and ionic strength is (M).  It is based on Alberty's values for CH4 in gas and aquious phase.
     ! However, these reactions are not affected by pH nor is, so could still have used NIST data,
     ! but this function provides consistency with Alberty's data and approach, so will keep.
     implicit none
     real(8) T, pH, is
     ! Local declarations
     real(8) dGrxn
     real(8), parameter:: R = 8.31451d0 ! (J/mol/K)
     ! Use the following reaction:  CH4(g) -> CH4(aq)
     dGrxn = dGf0(methaneaqsp,T,is,pH) - dGf0(methanegsp,T,is,pH)
     
     solCH4 = exp(-1.d3*dGrxn/(T*R))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
     
     return
   end function solCH4

    real(8) function solCH3OH(T, is, pH)
        ! This function returns the solubiility of methanol in water (uM/atm) at temperature T (K),
        ! pH and ionic strength is (M).  It is based on Alberty's values for CH4 in aquious phase
        ! and CH4 in the gas phase from CRC data.
        ! However, these reactions are not affected by pH nor is, so could still have used NIST data,
        ! but this function provides consistency with Alberty's data and approach, so will keep.
        implicit none
        real(8) T, pH, is
        ! Local declarations
        real(8) dGrxn
        real(8), parameter:: R = 8.31451d0 ! (J/mol/K)
        ! Use the following reaction:  CH3OH(g) -> CH3OH(aq)
        dGrxn = dGf0(methanolsp,T,is,pH) - dGf0(methanolgsp,T,is,pH)
        
        solCH3OH = exp(-1.d3*dGrxn/(T*R))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
        
        return
    end function solCH3OH

    real(8) function solH2CO3(T, is, pH)
        ! This function returns the solubiility of carbon dioxide in water (uM/atm) at temperature T (K)
        ! pH, and ionic strength is (M).  That is, given the CO2 partial pressure in ATM, this returns the 
        ! total concentration of all carbonate species (ie., DIC).  This is from Alberty2003, pg 151.
        ! Rxn is:  H2CO3(aq) -> CO2(g) + H2O
        implicit none
        real(8) T, pH, is
        ! Local declarations
        real(8) dGrxn
        real(8), parameter:: R = 8.31451d0 ! (J/mol/K)
        dGrxn = dGf0(co2gsp,T,is,pH) + dGf0(h2osp,T,is,pH) - dGf0(co2totsp,T,is,pH)

        solH2CO3 = exp(1.d3*dGrxn/(T*R))*1.0d6*1.01325d0 ! Last number converts Bar to ATM

        return
    end function solH2CO3    
    
   real(8) function solO2(T, is, pH)
     ! This function returns the solubiility of oxygen in water (uM/atm) at temperature T (K),
     ! pH and ionic strength is (M).  It is based on Alberty's values for O2 in gas and aquious phase.
     ! However, these reactions are not affected by pH nor is, so could still have used NIST data,
     ! but this function provides consistency with Alberty's data and approach, so will keep.
     implicit none
     real(8) T, pH, is
     ! Local declarations
     real(8) dGrxn
     real(8), parameter:: R = 8.31451d0 ! (J/mol/K)
     ! Use the following reaction:  O2(g) -> O2(aq)
     dGrxn = dGf0(o2aqsp,T,is,pH) - dGf0(o2gsp,T,is,pH)
     
     solO2 = exp(-1.d3*dGrxn/(T*R))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
     
     return
   end function solO2
   
!    real(8) function solCH4_old(T)
!      ! This function returns the solubiility of methane in water (uM/atm) at temperature T (K).
!      ! Data from http://webbook.nist.gov/chemistry/
!        implicit none
!        real(8) T
!        
!        solCH4_old = 1.4d-3*dexp( 1700.*(1.0/T-1.0/298.0) )*1.0d6
!        
!        return
!    end function solCH4_old            
!    real(8) function solO2(T)
!      ! This function returns the solubiility of oxygen in water (uM/atm) at temperature T (K).
!      ! Data from http://webbook.nist.gov/chemistry/
!        implicit none
!        real(8) T
!        
!        solO2 = 1.3d-3*dexp( 1700.*(1.0/T-1.0/298.0) )*1.0d6
!        
!        return
!    end function solO2
!
!    real(8) function solCO2(T)
!        ! This function returns the solubiility of carbon dioxide in water (uM/atm) at temperature T (K).
!        ! Data from http://webbook.nist.gov/chemistry/
!        implicit none
!        real(8) T
!
!        solCO2 = 3.4d-2*dexp( 2400.*(1.0/T-1.0/298.0) )*1.0d6
!
!        return
!    end function solCO2    
    
end module thermoData

