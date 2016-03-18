Module CommParameters
   ! This is version 2, 18 Mar 2010
   ! 03 Aug 2010  Ver 2.1  Adding ability to control WeederTech DOT board.
   ! Added some global variables to allow for correcting for sample loop dilution.
   Implicit none

   !============= Define Which ports go to which equip =================

   integer, parameter:: oxigrafPort = 4
   integer, parameter:: ch4Port     = 5
   integer, parameter:: valcoPort   = 6
   integer, parameter:: mfmPort     = 7
   integer, parameter:: WTDOTport   = 9 ! note port 8 has the MKS 647C controller.

   !===================CHARACTER CONSTANTS===================================

	CHARACTER(LEN=1), PARAMETER :: SOH	= CHAR( 1)
	CHARACTER(LEN=1), PARAMETER :: STX	= CHAR( 2)
	CHARACTER(LEN=1), PARAMETER :: ETX	= CHAR( 3)
	CHARACTER(LEN=1), PARAMETER :: EOT	= CHAR( 4)
	CHARACTER(LEN=1), PARAMETER :: ENQ	= CHAR( 5)
	CHARACTER(LEN=1), PARAMETER :: ACK	= CHAR( 6)
	CHARACTER(LEN=1), PARAMETER :: NAK	= CHAR(15)
	CHARACTER(LEN=1), PARAMETER :: DLE	= CHAR(16)
	CHARACTER(LEN=1), PARAMETER :: ESC	= CHAR(27)
	CHARACTER(LEN=1), PARAMETER :: CR	= CHAR(13)
	CHARACTER(LEN=2), PARAMETER :: crlf	= CHAR(13)//CHAR(10)
	
	! Some WeederTech WTDOT commands.  These are based on which channel the WTDOT board is set to
	! and what channels instruments/valves are connected to.
	! These commands are for turning on or off the gas pump that circulates gases to the detectors.
	character(len=3), parameter:: N2on       = 'ALA' ! Set channel A on Board A to low  (i.e., turn it on )
	character(len=3), parameter:: N2Off      = 'AHA' ! Set channel A on Board A to high (i.e., turn it off)
	character(len=3), parameter:: gasPumpOn  = 'ALB' ! Set channel B on Board A to low  (i.e., turn it on )
	character(len=3), parameter:: gasPumpOff = 'AHB' ! Set channel B on Board A to high (i.e., turn it off)
	character(len=3), parameter:: feed2On    = 'ALC' ! Set channel C on Board A to low  (i.e., turn it on )
	character(len=3), parameter:: feed2Off   = 'AHC' ! Set channel C on Board A to high (i.e., turn it off)
	character(len=3), parameter:: feed2Stat  = 'ARC' ! Read channel C on Board A 
	character(len=3), parameter:: feed2L     = 'ACL' ! On board A, Channel C is Low  (response to ARC)
	character(len=3), parameter:: feed2H     = 'ACH' ! On board A, Channel C is High (response to ARC)
	character(len=3), parameter:: f2SmpOn    = 'ALD' ! Set channel D on Board A to low  (i.e., turn it on )
	character(len=3), parameter:: f2SmpOff   = 'AHD' ! Set channel D on Board A to high (i.e., turn it off)

   !================== Other parameters

   real, parameter:: badRead = -9.99 ! values to use for bad port reads.
   real sampleTime ! Time to spend sampling one reactor (min).  This is needed to purge detectors, etc.
   real loopTime ! How often to intitiate a complete sampling (min).  This is entered as hours, but converted to nearist min.
                    ! loopTime > sampleTime*(nMCports+2)
   integer, parameter:: infoUnit=100, datUnit=101, rcvUnit=102  ! unit numbers for *.info, *.dat and *.rcv files.
   integer, parameter:: onOffUnit = 103 ! file that has times of FeedGas On and Off cycling

   integer, parameter:: maxNoPorts = 10 !Current number of sample ports on Valco Valve.

   character*15, parameter:: okDigits = '+-0123456789. ,'

   real co2Feed, o2Feed, ch4Feed  ! Concentrations of feed gas
   real o2Adj, co2Adj, ch4Adj ! 
   
   integer noteFlag  ! Is set to 1 when new note is entered.
   
   ! Version 3.0 additions
   logical gasCycleRunning ! If true, then the gasCycling of feed 2 is currently running.   
   logical Feed2GasOn ! current state of feed 2 gas valve
   real co2Feed2, o2Feed2, ch4Feed2  ! Concentrations of feed gas 2
   real co2Feed2cur, o2Feed2cur, ch4Feed2cur  ! Current concentrations of feed gas 2 depending on valve position
   integer noteCnt ! moved to here
   real(8) tzero, t0GasCycling
   real(8) tCycleOnOff ! Every tCycleOnOff days, Feed gas is switched between feed gas 1 or 2.
   integer ncycledPorts, cycledPorts(maxNoPorts) ! Which MCs have the feed gas cycled
   character filename*80 
   real(8) MCheadSpaceVol(maxNoPorts), SLvol
   real(8) SL_O2, SL_CO2, SL_CH4, SL_o2_new, SL_co2_new, SL_ch4_new ! concentrations of these gases in the sample loop.
   
   integer mixtankPort ! port that is used to sample mix-tank composition
   integer mixtankVentPort ! Port that is connected to the VENT of the mix-tank  
   integer feed2GasPort ! Port that is used to sample the feed2Gas composition.  Note, must open gas value to do so.
   real(8) SLpurgeTime ! How long sample loop needs to be purged with the mix-tank vent gas (min).
   integer MCports(maxNoPorts), nMCports ! Contains how many MCs to sample (not including mix-tank and
                                            ! feed 2 gas port. 
   character debugStr*1024, tempStr*1024 ! just using for debugging
   
   namelist /recovery/ o2Feed, co2Feed, ch4Feed, sampleTime, loopTime
   namelist /recovery/ tcycleOnOff, gasCycleRunning, ncycledPorts, cycledPorts, o2Feed2, co2Feed2, ch4Feed2
   namelist /recovery/ tzero, noteCnt, t0GasCycling, MCheadSpaceVol, SLvol, mixtankVentPort, SLpurgeTime
   namelist /recovery/ MCports, nMCports, mixtankPort, feed2GasPort
   
   contains
   subroutine writeRecovery
      ! This routine writes recovery info
      open(unit=rcvUnit,file=trim(filename)//'.rcv',status='unknown')
      write(rcvUnit,nml=recovery) 
      close(unit=rcvUnit)
      return
   end subroutine writeRecovery
   
   subroutine readRecovery
      ! This routine reads the recovery info
      open(unit=rcvUnit,file=trim(filename)//'.rcv',status='old')
      read(rcvUnit,nml=recovery) 
      close(unit=rcvUnit)
      return
   end subroutine readRecovery
End Module CommParameters

Program MonitorGases
   ! This Program is used to monitor the gass composition for the bioreactor setup.
   ! Versions
   !  15 May 05 1.0  First writing to read multiple ports for actual equipment
   !  23 Jun 05 1.01 Fixed serial port problem with sbits
   !  18 Aug 05 1.02 Added more purging of serial port data.
   !  23 Sep 05 1.03 Just have purge time now. Also look at menu control more often.
   !  11 Oct 05 1.04 Problem getting valco port to switch following calibration
   !  12 Jan 06 1.05 Modify program so that it also samples feed gas as well as microcosm exit gas
   !                 Also, values are displayed every 10 sec instead of just when storred.
   !  3 Jul 08  1.06 Modify time output format so that > 1000 days can be handelled.
   !                 Note, this is being compilter with Intel Fortran...
   ! 19 Nov 09  1.07 Attempt to sample anayzers repeatively so as to obtain values
   !                 with lower noise.
   ! 20 Mar 10  1.08 Adding ability to record gas analyzer values in a separate file. This
   !                 will allow easier determination of CH4/Air mix concentration that is now
   !                 being controlled by the MKS mass flow controllers (or will be soon).  Also
   !                 will try to improve analyzer readings by averaging reads.
   ! 25 Jun 10  1.09 Modified readCH4 so that an error returned in a read does not cause the
   !                 data to be lost (i.e., reported as -9.99).
   ! 03 Aug 10  2.0  Modifying the code to allow a new gas sampling loop that includes a gas
   !                 pump that requires turning on and off. Also, code will also be added
   !                 to allow the switching on and off a secondary feed for the ATB experiment.
   ! 06 Aug 10  2.01 Added ASCO valve to control N2 on an added code to do so.
   ! 18 Aug 10  3.0  Adding ability to cycle gasses on user specified reactors, also added recovery file.
   ! 01 Sep 10  3.01 Added ability to execute batch file to update web page followind data save.
   ! 10 Sep 10  3.02 Option to turn pump off/on during gas monitoring option 4.
   !                 Set valco value to port that is connected to the exit of the gas-mix tank.
   ! 11 Sep 10  3.03 Accounting for dilution of MC head space gas by sample loop.
   ! 12 Sep 10  3.04 More options regarding sample loop purging
   ! 13 Sep 10  4.0  Moved all routines to one file so that file versions are easier to handle.
   !                 Changed loop and data file to allow for feed gas 2 sampling and saving to dat file
   !                 Microcosm sample order now depends on state of feed gas 2 status.
   !                 Allow feed gas 2 to be sampled with control via ASCO on/off valve.
   ! 30 Sep 10  4.01 Fixed problem with invalid filename causing program to crash in menu option 4.
   ! 15 Mar 11  4.02 Changed gas cycling so that feed gas 2 is on for even days and off for odd.
   !                 This way, feed 2 will be turned on as soon as gas cycling is started if no
   !                 delay is given to the start of cycling.
   ! 15 May 12  4.03 Modified call to "system" for exporting images so that it does not wait for cmd to return
   !
   ! Notes for next version
   !  1. Have gas mixture composition print in notes info (done V. 2.0)
   
   use ifport 
   use ifwin
   use ifWINTY 
   use ifcore
   use CommParameters
   implicit None
   integer iresult, i, j, portNo, byteCnt, ok2read
   integer dt(8), ifail, optNo, cntrl
   real(8) datavec(maxNoPorts,21), beginT, endT, loopStartT, loopCurT, loopMin
   character  formStr*80, key*1, lineRead*80, tecStr*1000, longStr*120, version*20, ptN*2
   logical keyHit, infoExists, datExists, rcvExists, getOpt, testing
   ! Declarations with Version 3.0 and above
   real(8) eltime, t2cycle
   logical saveData, OnOffExists
   character outStr*80, onoff*3
   integer irow
   real(8) purgeTime, purgeTime0
   integer iports(maxNoPorts), nports

   ! Data declarations for console cursor access
   integer fhandle
   logical logstat
	Type(T_CONSOLE_SCREEN_BUFFER_INFO) conbuf
      Type (T_COORD)        dwCursorPosition

   ! Program version
   version = '4.03 (15 May 2012)'

   ! Get the handle of the console
   fhandle = GetStdHandle(STD_OUTPUT_HANDLE)

   ! Get name to label data files
   write(6,'(a)') 'Program Version: '//version 
   write(6,'(a,$)') 'Enter file name for data (no extension): '
   read(5,'(a)') filename

   ! Four files are used:  
   !  filename.info  contains information (notes) about the run
   !  filename.dat   Contains the collected data.
   ! If either files already exists, then data is appended (allows for restarts)
   !  filename.rcv   Contains all the user entered data to restart program.
   !  filename.OnOff Contains the times when feed gas swiching occurs.  It may not be used, but will be created.

   ! See if the files exist
   inquire(file=trim(filename)//'.info' ,exist=infoExists)
   inquire(file=trim(filename)//'.dat'  ,exist=datExists)
   inquire(file=trim(filename)//'.rcv'  ,exist=rcvExists)
   inquire(file=trim(filename)//'.OnOff',exist=OnOffExists)
  
   !open files
   if (infoExists) then
      write(6,*) '  **Warning**: file: ',trim(filename)//'.info',' exits. Will append data'
      open(unit=infoUnit,file=trim(filename)//'.info',status='old',access='append')
   else
      open(unit=infoUnit,file=trim(filename)//'.info',status='new')
   end if

   if (datExists) then
      write(6,*) '  **Warning**: file: ',trim(filename)//'.dat',' exits. Will append data'
      open(unit=datUnit,file=trim(filename)//'.dat',status='old',access='append')
   else
      open(unit=datUnit,file=trim(filename)//'.dat',status='new')
   end if

   if (OnOffExists) then
      write(6,*) '  **Warning**: file: ',trim(filename)//'.OnOff',' exits. Will append data'
      open(unit=onOffUnit,file=trim(filename)//'.OnOff',status='old',access='append')
   else
      open(unit=onOffUnit,file=trim(filename)//'.OnOff',status='new')
   end if

   ! Setup com ports
   call initializePorts()

   ! Initialize the gas analyzer offsets (only needed if not sampling feed gas)
   o2Adj = 0.0; co2Adj = 0.0; ch4Adj = 0.0
   feed2GasPort = 0 ! initialze to zero if not using.

   ! See if a recovery file exists
   if (rcvExists) then
      write(6,'(a,$)') 'Recovery file found! Do you wish to use it (y/n) ?'
      read(5,'(a)') lineRead
      if (scan(lineRead,'yY') /= 0) then
         call readRecovery 
      else
         call getProgramRunParams
      end if
   else
      call getProgramRunParams
   end if
   
   ! Regardless of whether gas cycling is on or not, get status of valve
   call writeWeederT(feed2Stat)
   ! get responce, which should equal feed2L
   call readWeederT(outStr)
   if (outStr == feed2L) then 
      Feed2GasOn = .true.
   else if (outStr == feed2H) then
      Feed2GasOn = .false.
   else
      write(6,'(a)') 'WARNING:: could not obtain status of feed gas 2 valve'
      write(6,'(2a)') 'WTDOT Output response was: ',trim(outStr)
      write(6,'(a)') 'Setting Feed2GasOn state to false'
      Feed2GasOn = .false. ! just set to off then
   end if
   if (feed2GasOn .and. (.not. gasCycleRunning) ) then
      write(6,'(a)') 'WARNING:: Feed 2 gas valve is on, but gasCycleRunning is false'
      write(6,'(a)') ' Will attempt to turn feed gas 2 valve off...'
      call writeWeederT(feed2Off)
      ! get responce, which should just be the return of command
      call readWeederT(outStr)
      if (outStr /= feed2Off) then 
         write(6,'(a)') 'WARNING:: Failed to turn Feed 2 Gas off!'
         write(6,'(2a)') 'WTDOT Output response was: ',trim(outStr)
      end if
   end if

   ! Alow users to enter info about run:
   call DATE_AND_TIME (values=dt)
   write(infoUnit,'(a)') 'Program Version: '//version 
   write(infoUnit,'(2(a,i2.2),a,i4,3(a,i2.2))') 'Program Started: ',dt(2),'/',dt(3),'/',dt(1),' ',dt(5),':',dt(6),':',dt(7)
   write(infoUnit,'(a,i2)')      'Mix-tank sample port: ', mixtankPort
   write(infoUnit,'(a,i2)')      'Mix-tank vent   port: ', mixtankVentPort
   write(infoUnit,'(a,6(x,i1))') 'Microcosm ports:      ',(MCports(i),i=1,nMCports)
   write(infoUnit,'(a,f5.1)') 'Sample loop purge time (min): ', SLpurgeTime
   write(infoUnit,'(a,f5.2)') 'Sample time (min): ', sampleTime
   write(infoUnit,'(a,f5.2)') 'Loop time    (hr): ', loopTime/60.0
   write(infoUnit,'(a,f6.2)') 'O2  feed (%): ', o2Feed
   write(infoUnit,'(a,f6.2)') 'CO2 feed (%): ', co2Feed
   write(infoUnit,'(a,f6.2)') 'CH4 feed (%): ', ch4Feed
   if (gasCycleRunning) then
      write(infoUnit,'(a,6(1x,i1))') 'The following MCs have been configured to have feed gas cycling:', (cycledPorts(i),i=1,ncycledPorts)
      write(infoUnit,'(a,f5.2,a)') 'Feed gas will be cycled every ', tCycleOnOff,' days.'
      write(infoUnit,'(a,i2)') 'Feed gas 2 port: ', feed2GasPort
      write(infoUnit,'(a,f6.2)') 'O2  feed 2 (%): ', o2Feed2
      write(infoUnit,'(a,f6.2)') 'CO2 feed 2 (%): ', co2Feed2
      write(infoUnit,'(a,f6.2)') 'CH4 feed 2 (%): ', ch4Feed2
   end if
   write(infoUnit,'(a)') '--------------------------- Experiment Description ----------------------------'
   write(6,'(a)') 'Enter description of experiment.  Blank line terminates entry.'
   do 
      write(6,'(a,$)') '> '
      read(5,'(a)') lineRead
      if (lineRead == ' ') exit
      write(infoUnit,'(a)') lineRead
   end do
   write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
   write(infoUnit,*) ' '

   ! Write tecplot header if this is a new data file
   if (.not. datExists) then
      tecStr = 'Variables = "M" "D" "Y" "Hr" "Min" "Sec" "O2 Feed" "CO2 Feed" "CH4 Feed" "O2 Feed 2" "CO2 Feed 2" "CH4 Feed 2"'
      longStr = ' "Time ## (d)" "Port ##" "Flow ## (sccm)" "O2r ## (%)" "O2 ## (%)" "CO2r ## (%)" "CO2 ## (%)" "CH4r ## (%)" "CH4 ## (%)"'
      ! Must have space for sampling mix tank and feed gas 2, even if not used.
      do i=1,2
         ptN = 'F'//char(48+i)
         longStr(8:9)     = ptN ! Time
         longStr(22:23)   = ptN ! Port
         longStr(32:33)   = ptN ! Flow
         longStr(48:49)   = ptN ! O2r
         longStr(60:61)   = ptN ! O2
         longStr(74:75)   = ptN ! CO2r
         longStr(87:88)   = ptN ! CO2
         longStr(101:102) = ptN ! CH4r
         longStr(114:115) = ptN ! CH4
         tecStr = trim(tecStr)//longStr
      end do
      do i=1,nMCports
         ptN = 'M'//char(48+i)
         longStr(8:9)     = ptN ! Time
         longStr(22:23)   = ptN ! Port
         longStr(32:33)   = ptN ! Flow
         longStr(48:49)   = ptN ! O2r
         longStr(60:61)   = ptN ! O2
         longStr(74:75)   = ptN ! CO2r
         longStr(87:88)   = ptN ! CO2
         longStr(101:102) = ptN ! CH4r
         longStr(114:115) = ptN ! CH4
         tecStr = trim(tecStr)//longStr
      end do
      tecStr = trim(tecStr)//' "Notes"'
      write(datUnit,'(a)') tecStr
      write(datUnit,'(a)') 'Zone'
   end if  
   
   ! Store the input info in a recovery file
   call writeRecovery
      
   !  ------------------------------------------------------
   !  -------------- Begin Test loop ---------------------
   !  ------------------------------------------------------
   write(6,'(a,$)') 'Do you want to test gas sampling loop hardware/software? (y/n): '
   read(5,'(a)') lineRead
   testing = .false.
   if (scan(lineRead,'yY') /= 0) testing = .true.
   call setIports (iports, nports, Feed2GasOn)
   do while (testing)
      Write(6,'(a)') 'Running test reads...' 
      call pumpOn ()
      write(6,'(a)') 'Gas pump should be on, check flow'
      do i=1,nports
         !  Set the valve to position 1 to nports
         portNo = iports(i) 
         Call valcoControl(portNo)
         call convPort2data (portNo, irow) ! row in datavec that the port correspond to.
         datavec(irow,14) = dble(portNo)

         call readGasAnalyzers (portNo,datavec(irow,:),tzero,.true.)

         longStr =  '(2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,1x,i2,##x,f6.2,6x,f6.2,6x,f7.3,6x,f8.4)'
         write(longStr(43:44),'(i2)') 2*(i-1)+1

         write(6,longStr) (int(datavec(irow,j)),j=1,6), int(datavec(irow,14)), real(datavec(irow,15)), (real(datavec(irow,j)),j=17,21,2)
      end do 
      call pumpOff ()
      write(6,'(a)') 'Gas pump should be off, check flow'
      write(6,'(a)') 'Done with test reads.'
      write(6,'(a$)') 'Run test reads again? (y/n): '
      read(5,'(a)') lineRead
      testing = .false.
      if (scan(lineRead,'yY') /= 0) testing = .true.
   end do

   Write(6,'(a)') 'Sampling begining...'
   write(6,'(a)') 'Hit ESC at anytime to bring up Menu.'

   !  ------------------------------------------------------
   !  -------------- Begin Sample loop ---------------------
   !  ------------------------------------------------------
   noteFlag = 0
   call pumpOff () ! make sure pump is off
   portNo = mixtankVentPort
   Call valcoControl(portNo) ! Move to valvo port connected to mix-tank VENT
   datavec(portNo,14) = dble(portNo)
   purgeTime0 = RTC() ! Start the clock on the purging of the sample loop with vent gas
   
   sampling: do
      call f2SampleGasOff ! make sure feed 2 gas sample valve is closed.
      ! Make sure sample loop has been purged with vent gas for at least SLpurgeTime amount of time.
      purgeTime = (RTC() - purgeTime0)/60.0 ! time that sample loop as been purged in min.
      do while (purgeTime < SLpurgeTime)
         purgeTime = (RTC() - purgeTime0)/60.0
         ! get the cursor position
         logstat = GetConsoleScreenBufferInfo(fhandle, conbuf)
         call readGasAnalyzers (portNo,datavec(portNo,:),tzero,.true.)
         longStr = '(''At: '',2(i2.2,'':''),i2.2,'' Port: '',i3,'' Flow: '',f6.2,'' O2: '',f6.2,'' CO2: '',f7.3,'' CH4: '',f8.4)'
         write(6,longStr)  (int(datavec(portNo,j)),j=4,6), int(datavec(portNo,14)), real(datavec(portNo,15)), (real(datavec(portNo,j)),j=17,21,2)
         write(6,'(a,f5.1,a)') 'Need to wait ',SLpurgeTime-purgeTime,' min. for sample loop purge. ESC to skip.'
         ! put the cursor back to where it was
         logstat = SetConsoleCursorPosition(fhandle, conbuf.dwCursorPosition)
         call SLEEPQQ (10000) ! Sleep for 10 sec.
         ! See if ESC key was hit
         keyHit = PEEKCHARQQ ( )
         if (keyHit) then
            !See if key hit is ESC, if so exit.
            key = GETCHARQQ( )
            if (key == ESC) then
               write(6,'(//a,$)') 'Do you want to skip sample loop purge? (y/n): '
               read(5,'(a)') lineRead
               if (scan(lineRead,'yY') /= 0) exit 
            end if
         end if
      end do       
      
      call cycleFeedGas2 ! set feed 2 gas valve position if using.
      ! Determine which ports to sample depending state of Feed2Gas valve
      call setIports (iports, nports, Feed2GasOn)
      portNo = iports(1)
      Call valcoControl(portNo) ! Move to mix-tank sample port
      ! turn gas pump on
      call pumpOn ()
      saveData = .true.
      loopStartT = RTC() ! marks the time, in sec, that a sampling loop started
      
      readPort: do i=1,nports
         !  Set the valve to position 1 to nports
         portNo = iports(i) 
         Call valcoControl(portNo)
         call convPort2data (portNo, irow) ! row in datavec that the port correspond to.
         if (portNo == feed2GasPort) call f2SampleGasOn ! open valve to sample feed 2 gas.
         datavec(irow,14) = dble(portNo)

         datavec(irow,7) = o2Feed
         datavec(irow,8) = co2Feed
         datavec(irow,9) = ch4Feed

         datavec(irow,10) = o2Feed2cur
         datavec(irow,11) = co2Feed2cur
         datavec(irow,12) = ch4Feed2cur

         ! Now wait sufficient time for system to be flushed (purgeTime)
         ! Get time in seconds since 1970
         beginT = RTC()
         wait: do
            ! Sleepy sleep for 10 sec.
            call SLEEPQQ (10000)
            ! See if ESC key was hit
            keyHit = PEEKCHARQQ ( )
            if (keyHit) then
               !See if key hit is ESC, if is so exit.
               key = GETCHARQQ( )
               if (key == ESC) then
                  call pumpOff () ! turn pump off while in menu.
                  call shortMenu(noteCnt, cntrl)
                  select case (cntrl)
                     case(1) ! just continue sampling port like nothing happened.
                        call pumpOn () ! turn pump on and return to sampling.
                     case(2) ! restart readPort sample loop, do not save anything
                        purgeTime0 = RTC() ! Assume sample loop got contaminated
                        cycle sampling
                     case(3) ! exit readPort sample loop, do not save anything
                        saveData = .false.
                        exit readPort
                     case(9) ! exit program
                        write(6,'(a,$)') 'Exit program? (y/n): '
                        read(5,'(a)') lineRead
                        if (scan(lineRead,'yY') /= 0) exit sampling
                        call pumpOn () ! turn pump on and return to sampling.
                  end select ! if neither case is selected, just continue sampling port
               end if
            end if
            endT = RTC()
            if ( (endT-beginT)/60. >= sampleTime ) exit wait

            ! ok, just print update, but don't store.
            call readGasAnalyzers (portNo,datavec(irow,:),tzero,.false.)
            ! get the cursor position
            logstat = GetConsoleScreenBufferInfo(fhandle, conbuf)
            ! write current reading out at the same location
            longStr = '(''Purge t:'',f4.1,,'' min'','' Port:'',i2,'' Flow: '',f6.2,'' O2:'',f6.2,'' CO2:'',f7.3,'' CH4:'',f8.4)'
            write(6,longStr)  (endT-beginT)/60., int(datavec(irow,14)), real(datavec(irow,15)), (real(datavec(irow,j)),j=17,21,2)
            ! put the cursor back to where it was
            logstat = SetConsoleCursorPosition(fhandle, conbuf.dwCursorPosition)
         end do wait

         call readGasAnalyzers (portNo,datavec(irow,:),tzero,.false.)

         if (portNo == feed2GasPort) call f2SampleGasOff ! close valve to sample feed 2 gas.
      
         longStr =  '(2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,1x,i2,##x,f6.2,6x,f6.2,6x,f7.3,6x,f8.4)'
         write(longStr(43:44),'(i2)') 2*(i-1)+1
         formStr = ' '
         write(formStr,longStr) (int(datavec(irow,j)),j=1,6), int(datavec(irow,14)), real(datavec(irow,15)), (real(datavec(irow,j)),j=17,21,2)
         write(6,'(a)') formStr(1:79)
         
         ! Set the concentration of the gases in the sample loop to the current measured value
         SL_o2 = SL_o2_new; SL_co2 = SL_co2_new; SL_ch4 = SL_ch4_new
      end do readPort
      ! Turn gas pump off
      call pumpOff ()
      call f2SampleGasOff ! make sure feed 2 sample valve is off
      
      ! Store data in *.dat file (designed for Tecplot read) 
      if (saveData) then   
         ! Feed gas 2 was not sampled, then just copy info from main feed sampling
         if (.not. Feed2GasOn) datavec(2,:) = datavec(1,:)
         longStr = '(2(f3.0,1x),f5.0,3(1x,f3.0),6(1x,f7.4),#(1x,f10.5,1x,f4.0,3(1x,f6.2),2(1x,f7.3),2(1x,f8.4),1x),i3)'
         write(longStr(40:40),'(i1)') nMCports+2
         write(datUnit,longStr) (datavec(1,j),j=1,12), ((datavec(i,j),j=13,21),i=1,nMCports+2), noteCnt*noteFlag
         noteFlag = 0
         ! Generate tecplot images for web page by calling dos command batch file.
         ! iresult = system('exportImages.bat > exportImages.txt')
         iresult = system('start "ExportImages" cmd /C exportImages.bat') ! this allows the program to return without waiting.
      end if
      
      ! Now wait loopTime until time to begin next total sampling
      purgeTime0 = RTC() ! Start the clock on the purging of the sample loop with vent gas
      loopWait: do
         ! make sure that the gas-mix tank is purging the sample loop.
         portNo = mixtankVentPort ! port connected to exit of gas mix tank
         Call valcoControl(portNo)
         datavec(portNo,14) = dble(portNo)
         
         call cycleFeedGas2 ! set feed 2 gas valve position if using.
         loopCurT = RTC() ! get current time, in sec, and check loop time.
         loopMin = (loopCurT-loopStartT)/60.0
         if ( loopMin >= loopTime ) exit loopWait
         ! get the cursor position
         logstat = GetConsoleScreenBufferInfo(fhandle, conbuf)
         
         call readGasAnalyzers (portNo,datavec(portNo,:),tzero,.true.)
         longStr = '(''At: '',2(i2.2,'':''),i2.2,'' Port: '',i3,'' Flow: '',f6.2,'' O2: '',f6.2,'' CO2: '',f7.3,'' CH4: '',f8.4)'
         write(6,longStr)  (int(datavec(portNo,j)),j=4,6), int(datavec(portNo,14)), real(datavec(portNo,15)), (real(datavec(portNo,j)),j=17,21,2)
         
         write(6,'(a,f7.2,a)') 'Time to next sampling: ',loopTime-loopMin,  &
                               ' (min). Hit ESC for LONG menu'
         if (gasCycleRunning) then
            if (Feed2GasOn) then 
               onoff = ' on'
            else
               onoff = 'off'
            end if
            call DATE_AND_TIME (values=dt)
            call julday(dt(1),dt(2),dt(3),dt(5),dt(6),real(dt(7)),eltime,ifail)
            eltime = eltime - t0GasCycling ! elapsed time in days since start of gas cycling
            t2cycle = tcycleOnOff - mod(eltime,tcycleOnOff)
            write(6,'(3a,f6.3,a)') 'Gas feed 2 valve is ', onoff,'. Valve cycle in ',t2cycle,' days.'
         end if
         logstat = SetConsoleCursorPosition(fhandle, conbuf.dwCursorPosition)
         
         ! Sleepy sleep for 10 sec.
         call SLEEPQQ (10000)
         ! See if ESC key was hit
         keyHit = PEEKCHARQQ ( )
         if (keyHit) then
            !See if key hit is ESC, if is so exit.
            key = GETCHARQQ( )
            if (key == ESC) then
               call longMenu(noteCnt, cntrl)
               select case (cntrl)
                  case(1) ! just continue sampling port like nothing happened.
                     cycle loopWait
                  case(2) ! restart sample loop BEFORE loopTime has expired
                     exit loopWait
                  case(3) ! Get new feed gas concentrations, but continue to wait
                     cycle loopWait
                  case(4) ! Record gas analyzer data, but continue to wait on return.
                     purgeTime0 = RTC() ! Assume sample loop got contaminated
                     cycle loopWait
                  case(5) ! Turn gas cycling on, but continue to wait.
                     cycle loopWait
                  case(6) ! Turn gas cycling off, but continue to wait.
                     cycle loopWait
                  case(7) ! Test feed gas cycling valve, but continue to wait.
                     cycle loopWait
                  case(8) ! Enter note, but continue to wait.
                     cycle loopWait
                  case(9) ! Get new headspace gas volumes, but continue to wait.
                     cycle loopWait
                  case(10) ! exit program
                     write(6,'(a,$)') 'Exit program? (y/n): '
                     read(5,'(a)') lineRead
                     if (scan(lineRead,'yY') /= 0) exit sampling
                     cycle loopWait
               end select ! if no case is selected, just continue to wait (shouldn't happen)
            end if
         end if
      end do loopWait

   end do sampling
   !  ------------------------------------------------------
   ! -------------- End Main Sampling Loop -----------------
   !  ------------------------------------------------------
   ! make sure pump is off and feed 2 sample valve is closed.
   ! however, leave the feed gas switcing valve unchanged.
   ! place valco port onto purging with vent gas
   call pumpOff ()
   call f2SampleGasOff
   portNo = mixtankVentPort
   Call valcoControl(portNo) ! Move to valvo port connected to mix-tank VENT
   
   ! Let user know status of port before leaving program.
   call writeWeederT(feed2Stat)
   ! get responce, which should just be the return of command
   call readWeederT(outStr)
   if (outStr == feed2L) then 
      write(6,'(/a)') 'NOTE: feed gas 2 valve is ON'
   else if (outStr == feed2H) then
      write(6,'(/a)') 'NOTE: feed gas 2 valve is OFF'
   else
      write(6,'(/a)')  'WARNING:: could not obtain status of feed gas 2 valve'
      write(6,'(2a)') '  WTDOT Output response was: ',trim(outStr)
      write(6,'(a)')  '  You should check valve status manualy!'
      Feed2GasOn = .false. ! just set to off then
   end if

   !  Release all the ports
   iresult = SPORT_RELEASE (oxigrafPort)
   iresult = SPORT_RELEASE (ch4Port    )
   iresult = SPORT_RELEASE (valcoPort  )
   iresult = SPORT_RELEASE (mfmPort    )
   iresult = SPORT_RELEASE (WTDOTport  )
   
   ! Although recovery file should be up to date, write anyway
   call writeRecovery

   close(unit=infoUnit)
   close(unit=datUnit )
   close(unit=onOffUnit )

   stop
end program MonitorGases

subroutine shortMenu(noteCnt, optNo) 
   use DFLIB
   use CommParameters, only: noteFlag, infoUnit, o2Feed, co2Feed, ch4Feed, gasCycleRunning, &
                             o2Feed2, co2Feed2, ch4Feed2, ncycledPorts, cycledPorts, tCycleOnOff, &
                             feed2Off, feed2On
   implicit none 
   integer noteCnt, optNo

   ! Local declarations
   integer dt(8), input
   logical getOpt
   character formStr*80, lineRead*80, outStr*80

   CALL BEEPQQ(1000, 2000)
   getOpt = .true.
   opt: do while (getOpt)
      write(6,'(/)') 
      write(6,'(a)') 'Sampling Paused; Options:'
      write(6,'(a)') '  1. Return to sampling.'
      write(6,'(a)') '  2. Restart sample loop imediately after entering a note.'
      write(6,'(a)') '  3. Exit current sampling w/o saving this data point.'
      write(6,'(a)') '  9. Exit program.'
      read(5,*) optNo
      select case (optNo)
      case(1) ! Just return like nothing happened.
         write(6,*) '  Returning to sampling, press ESC for menu'
         write(6,*) ' '
         getOpt = .false.
         
      case(2) ! Restart sample loop and the entering of a note
         noteCnt = noteCnt + 1
         noteFlag = 1
         call DATE_AND_TIME (values=dt)
         write(formStr,'(2(a,i2.2),a,i4,3(a,i2.2),a)') &
           ' (',int(dt(2)),'/',int(dt(3)),'/',int(dt(1)),' ',int(dt(5)),':',int(dt(6)),':',int(dt(7)),')'
         write(infoUnit,'(a,i3,a)') 'Note Number: ',noteCnt,trim(formStr)
         write(6,'(a,i3)') 'Enter description for Note: ',noteCnt
         write(6,'(a)') 'Enter blank line to return to sampling (calibrate before returning)'
         do 
            write(6,'(a,$)') '> '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            write(infoUnit,'(a)') lineRead
         end do
         write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
         write(infoUnit,*) ' '
         getOpt = .false.
         write(6,*) '  ** Returing to sampling, press ESC for menu'
         
      case(3)
         write(6,*) '  Returning to sampling, press ESC for menu'
         write(6,*) ' '
         getOpt = .false.
       
      case(9)
         write(6,*) '  Exiting program.' 
         getOpt = .false.
         
      case default
         write(6,*) '  ** Option not recognized, try again **'
      end select
   end do opt
   return
end subroutine shortMenu

subroutine longMenu(noteCnt, optNo) 
   use DFLIB
   use CommParameters, only: noteFlag, infoUnit, o2Feed, co2Feed, ch4Feed, gasCycleRunning, &
                             o2Feed2, co2Feed2, ch4Feed2, ncycledPorts, cycledPorts, tCycleOnOff, &
                             feed2Off, feed2On, writeRecovery, o2Feed2Cur, co2Feed2Cur, ch4Feed2Cur, &
                             feed2GasPort, MCheadSpaceVol, SLvol, nMCports
   implicit none 
   integer noteCnt, optNo

   ! Local declarations
   integer dt(8), input, i
   logical getOpt
   character formStr*80, lineRead*80, outStr*80

   CALL BEEPQQ(1000, 2000)
   getOpt = .true.
   opt: do while (getOpt)
      write(6,'(//)') 
      write(6,'(a)') 'Sampling Paused; Options:'
      write(6,'(a)') '  1. Return to sampling.'
      write(6,'(a)') '  2. Restart sample loop IMMEDIATELY after entering a note.'
      write(6,'(a)') '  3. Enter new feed gas conc. and note.'
      write(6,'(a)') '  4. Begin recording gas analyzer readings.'
      write(6,'(a)') '  5. Begin gas cycling with Feed 2.'
      write(6,'(a)') '  6. End gas cycling of Feed 2.'
      write(6,'(a)') '  7. Test feed gas cycling valve.'
      write(6,'(a)') '  8. Enter a note, but continue to wait for next sample event.'
      write(6,'(a)') '  9. Enter new microcosm headspace gas volumes.'
      write(6,'(a)') ' 10. Exit program'
      read(5,*) optNo
      select case (optNo)
      case(1)
         write(6,*) '  Returning to sampling, press ESC for menu'
         write(6,*) ' '
         getOpt = .false.
         
      case(2)
         ! Enter a note and start sampling MCs imediatly on return.
         noteCnt = noteCnt + 1
         noteFlag = 1
         call DATE_AND_TIME (values=dt)
         write(formStr,'(2(a,i2.2),a,i4,3(a,i2.2),a)') &
           ' (',int(dt(2)),'/',int(dt(3)),'/',int(dt(1)),' ',int(dt(5)),':',int(dt(6)),':',int(dt(7)),')'
         write(infoUnit,'(a,i3,a)') 'Note Number: ',noteCnt,trim(formStr)
         write(6,'(a,i3)') 'Enter description for Note: ',noteCnt
         write(6,'(a)') 'Enter blank line to return to sampling (calibrate before returning)'
         do 
            write(6,'(a,$)') '> '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            write(infoUnit,'(a)') lineRead
         end do
         write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
         write(infoUnit,*) ' '
         getOpt = .false.
         write(6,*) '  ** Returing to sampling, press ESC for menu'
         
      case(3)
         noteCnt = noteCnt + 1
         noteFlag = 1
         ! Get concentrations of feed gas
         write(6,'(a)') 'Enter feed gas concentrations:'
         write(6,'(a,$)') 'O2  (%): '
         read(5,*) o2Feed
         write(6,'(a,$)') 'CO2 (%): '
         read(5,*) co2Feed
         write(6,'(a,$)') 'CH4 (%): '
         read(5,*) ch4Feed
         
         if (gasCycleRunning) then
            write(6,'(a)') 'Enter feed 2 gas concentrations:'
            write(6,'(a,$)') 'O2  (%): '
            read(5,*) o2Feed2
            write(6,'(a,$)') 'CO2 (%): '
            read(5,*) co2Feed2
            write(6,'(a,$)') 'CH4 (%): '
            read(5,*) ch4Feed2
         end if
   
         call DATE_AND_TIME (values=dt)
         write(formStr,'(2(a,i2.2),a,i4,3(a,i2.2),a)') &
           ' (',int(dt(2)),'/',int(dt(3)),'/',int(dt(1)),' ',int(dt(5)),':',int(dt(6)),':',int(dt(7)),')'
         write(infoUnit,'(a,i3,a)') 'Note Number: ',noteCnt,trim(formStr)
         write(infoUnit,'(a)'     ) 'New feed gas composition entered:'
         write(infoUnit,'(a,f6.2)') '  O2  feed (%): ', o2Feed
         write(infoUnit,'(a,f6.2)') '  CO2 feed (%): ', co2Feed
         write(infoUnit,'(a,f6.2)') '  CH4 feed (%): ', ch4Feed
         if (gasCycleRunning) then
            write(infoUnit,'(a)')      'New feed 2 gas composition entered:'
            write(infoUnit,'(a,f6.2)') '  O2  feed 2 (%): ', o2Feed2
            write(infoUnit,'(a,f6.2)') '  CO2 feed 2 (%): ', co2Feed2
            write(infoUnit,'(a,f6.2)') '  CH4 feed 2 (%): ', ch4Feed2
         end if
         
         write(6,'(a,i3)') 'Enter description for Note: ',noteCnt
         write(6,'(a)') 'Enter blank line to return to sampling (calibrate before returning)'
         do 
            write(6,'(a,$)') '> '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            write(infoUnit,'(a)') lineRead
         end do
         write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
         write(infoUnit,*) ' '
         getOpt = .false.
         call writeRecovery ! update recovery file
         write(6,*) '  ** Returing to sampling, press ESC for menu'
         
      case(4) ! Save gas anlyzer data to a specified file
         call readAndStoreAnalzers (1)
         write(6,'(/,a)') '  ** Returing to sampling, press ESC for menu'
         getOpt = .false.
         
      case(5) ! Begin cycling of feed gas 2
         noteCnt = noteCnt + 1
         noteFlag = 1
         gasCycleRunning = .true.
         call gasCyclingParams
         call DATE_AND_TIME (values=dt)
         write(formStr,'(2(a,i2.2),a,i4,3(a,i2.2),a)') &
           ' (',int(dt(2)),'/',int(dt(3)),'/',int(dt(1)),' ',int(dt(5)),':',int(dt(6)),':',int(dt(7)),')'
         write(infoUnit,'(a,i3,a)') 'Note Number: ',noteCnt,trim(formStr)
         write(infoUnit,'(a,6(1x,i1))') 'The following MCs have been configured to have feed gas cycling:' &
                                         , (cycledPorts(i),i=1,ncycledPorts)
         write(infoUnit,'(a,f5.2,a)') 'Feed gas will be cycled every ', tCycleOnOff,' days.'
         write(infoUnit,'(a,i2)') 'Feed gas 2 port: ', feed2GasPort
         write(infoUnit,'(a,f6.2)') 'O2  feed 2 (%): ', o2Feed2
         write(infoUnit,'(a,f6.2)') 'CO2 feed 2 (%): ', co2Feed2
         write(infoUnit,'(a,f6.2)') 'CH4 feed 2 (%): ', ch4Feed2
         write(6,'(a,i3)') 'Enter description for Note: ',noteCnt
         write(6,'(a)') 'Enter blank line to return to sampling (calibrate before returning)'
         do 
            write(6,'(a,$)') '> '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            write(infoUnit,'(a)') lineRead
         end do
         write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
         write(infoUnit,*) ' '
         getOpt = .false.
         call writeRecovery ! update recovery file
         write(6,*) '  ** Returing to sampling, press ESC for menu'
         
      case(6) ! Stop gas cycling
         noteCnt = noteCnt + 1
         noteFlag = 1
         gasCycleRunning = .false.
         ! Make sure valve if off
         call forceFeedGas2Valve (0)  
         ! Add a note.       
         call DATE_AND_TIME (values=dt)
         write(formStr,'(2(a,i2.2),a,i4,3(a,i2.2),a)') &
           ' (',int(dt(2)),'/',int(dt(3)),'/',int(dt(1)),' ',int(dt(5)),':',int(dt(6)),':',int(dt(7)),')'
         write(infoUnit,'(a,i3,a)') 'Note Number: ',noteCnt,trim(formStr)
         write(infoUnit,'(a)') ' Gas cycling of feed 2 turned off'        
         write(6,'(a,i3)') 'Enter description for Note: ',noteCnt
         write(6,'(a)') 'Enter blank line to return to sampling (calibrate before returning)'
         do 
            write(6,'(a,$)') '> '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            write(infoUnit,'(a)') lineRead
         end do
         write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
         write(infoUnit,*) ' '
         getOpt = .false. 
         call writeRecovery ! update recovery file
         write(6,*) '  ** Returing to sampling, press ESC for menu'

      case(7) ! Test feed cycling valve
         write(6,'(a)') 'If gas cycling is ON, the valve status will revert to what it was on return.'
         write(6,'(a)') 'But, if gas cycling is OFF, valve retains state on return to sampling!!'
         do                                                              
            write(6,'(a,$)') 'Enter 1 to turn feed gas 2 valve on, 0 for off, or <cr> to exit: '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            read(lineRead,*) input
            call forceFeedGas2Valve (input)
         end do
         getOpt = .false.
         write(6,*) '  ** Returing to sampling, press ESC for menu'
       
     case(8)
         ! Pause to allow entering of a note
         ! Unlike option 2, this just continues to wait until next sampling event.
         noteCnt = noteCnt + 1
         noteFlag = 1
         call DATE_AND_TIME (values=dt)
         write(formStr,'(2(a,i2.2),a,i4,3(a,i2.2),a)') &
           ' (',int(dt(2)),'/',int(dt(3)),'/',int(dt(1)),' ',int(dt(5)),':',int(dt(6)),':',int(dt(7)),')'
         write(infoUnit,'(a,i3,a)') 'Note Number: ',noteCnt,trim(formStr)
         write(6,'(a,i3)') 'Enter description for Note: ',noteCnt
         write(6,'(a)') 'Enter blank line to return to sampling.'
         do 
            write(6,'(a,$)') '> '
            read(5,'(a)') lineRead
            if (lineRead == ' ') exit
            write(infoUnit,'(a)') lineRead
         end do
         write(infoUnit,'(a)') '-------------------------------------------------------------------------------'
         write(infoUnit,*) ' '
         getOpt = .false.
         write(6,*) '  ** Returing to sampling, press ESC for menu'      
         
      case(9) 
         ! Enter new values for MC headspace gas volume (mL)
         write(6,'(a)') 'Enter headspace gas volumes for MCs (mL): '
         read(5,*) (MCheadSpaceVol(i),i=1,nMCports)
         write(6,'(a,$)') 'Enter gas volume of sample loop (mL): '
         read(5,*) SLvol
         call writeRecovery ! update recovery file
         getOpt = .false.
     
      case(10)
         write(6,*) '  Exiting program.' 
         getOpt = .false.
         
      case default
         write(6,*) '  ** Option not recognized, try again **'
      end select
   end do opt
   return
end subroutine longMenu

subroutine convPort2data (iport, idatavec)
   ! This routine maps the port number to the correct datavec row.
   ! The datavec rows correspond as follows, regardless of sampling or use of a second feed gas:
   ! [mixtankPort, feed2gasPort, MCports(1), ..., MCports(nMCports)]
   use CommParameters, only: mixtankPort, feed2GasPort, MCports, nMCports
   implicit none
   integer iport, idatavec
   
   integer i
   
   if (iport == mixtankPort) then
      idatavec = 1
      return
   end if
   
   if (iport == feed2gasPort) then
      idatavec = 2
      return
   end if
   
   do i=1,nMCports
      if (iport == MCports(i)) then
         idatavec = i+2
         return
      end if
   end do
   
   write(6,'(a)')   'ERROR::convPort2data: port value incorrect, Stopping!'
   write(6,'(a,i)') '  iport = ',iport
   stop      
end subroutine convPort2data

Subroutine readGasAnalyzers (iport,datavec,tzero,NoAdjust)
   use CommParameters, only: o2Adj, co2Adj, ch4Adj, o2Feed, co2Feed, ch4Feed, &
                             SL_O2, SL_CO2, SL_CH4, SLvol, MCheadSpaceVol, &
                             SL_o2_new, SL_co2_new, SL_ch4_new, mixtankVentPort, &
                             mixtankPort, feed2GasPort, o2Feed2, co2Feed2, ch4Feed2
   implicit none
   integer iport ! current valco port selected
   logical NoAdjust ! If true, not adjustment is done to the "corrected" gas values
   real(8) datavec(*), tzero

   !local declarations
   integer dt(8), ifail, i
   real flow, o2, co2, ch4

   ! Get time stamp and julian time
   call DATE_AND_TIME (values=dt)
   datavec(1:6) = (/dble(dt(2)),dble(dt(3)),dble(dt(1)),dble(dt(5)),dble(dt(6)),dble(dt(7))/)
   call julday(dt(1),dt(2),dt(3),dt(5),dt(6),real(dt(7)),datavec(13),ifail)
   datavec(13) = datavec(13) - tzero

   ! Sample Mass flow meter
   call readMassFlowMeter (flow)
   datavec(15) = dble(flow)

   ! Sample oxigraph
   call readOxigraf(o2, co2)
   datavec(16) = dble(o2)
   datavec(18) = dble(co2)

   ! Sample CH4
   call readCH4 (ch4) 
   datavec(20) = dble(ch4)
   
   ! Copy this data over to the "adjusted" gas values columns
   datavec(17) = datavec(16)
   datavec(19) = datavec(18)
   datavec(21) = datavec(20)    
   if (NoAdjust) return   

   if (iport == mixtankPort .or. iport == feed2GasPort) then
      ! Either mix tank (feed 1) or feed 2 port is being sampled, so generate offsets
      if (iport == mixtankPort) then
         o2Adj = o2 - o2Feed
         co2Adj = co2 - co2Feed
         ch4Adj = ch4 - ch4Feed
         datavec(17) = dble(o2Adj)
         datavec(19) = dble(co2Adj)
         datavec(21) = dble(ch4Adj)
      else
         ! Note, if the feed2gasPort is being sampled, then the feed 2 gas MUST be on.
         o2Adj = o2 - o2Feed2
         co2Adj = co2 - co2Feed2
         ch4Adj = ch4 - ch4Feed2
         datavec(17) = dble(o2Adj)
         datavec(19) = dble(co2Adj)
         datavec(21) = dble(ch4Adj)
      end if
      ! Store the observed concentration to use for the next SL values
      SL_o2_new  = o2  - o2Adj
      SL_co2_new = co2 - co2Adj
      SL_ch4_new = ch4 - ch4Adj
      
   else if (iport == mixtankVentPort) then
      ! Sample loop is being purged with mix-tank vent gas, just report the values as is
      datavec(17) = datavec(16)
      datavec(19) = datavec(18)
      datavec(21) = datavec(20)    
   else
      ! A microcosm is being sampled, so acount for offset and dilution due to sample loop volume.
      datavec(17) = datavec(16) - dble(o2Adj)
      datavec(19) = datavec(18) - dble(co2Adj)
      datavec(21) = datavec(20) - dble(ch4Adj)
      ! Store the observed concentration to use for the next SL values
      SL_o2_new  = datavec(17)
      SL_co2_new = datavec(19)
      SL_ch4_new = datavec(21)
      ! Correct for dulition of MC headspace by sample loop gas.
      call convPort2data (iport, i); i = i - 2 ! i is the microcosm number, as listed in MCports.
      datavec(17) = datavec(17)*(1.0 + SLvol/MCheadSpaceVol(i)) - SL_O2 *SLvol/MCheadSpaceVol(i)
      datavec(19) = datavec(19)*(1.0 + SLvol/MCheadSpaceVol(i)) - SL_CO2*SLvol/MCheadSpaceVol(i)
      datavec(21) = datavec(21)*(1.0 + SLvol/MCheadSpaceVol(i)) - SL_CH4*SLvol/MCheadSpaceVol(i)
   end if
 
   return
end Subroutine readGasAnalyzers

subroutine getProgramRunParams
   use CommParameters
   implicit none
   ! This routine get the program run parameters from the user
   ! Note, these parameters are latter saved in a recovery file that will
   ! be read in.
   
   ! Local declarations
   character  lineRead*80
   integer dt(8), i, ifail   
   real TinHours   
   integer t0day, t0month, t0year ! date to base tzero on.   
   
   do
      write(6,'(a,$)') 'Enter number of microcosms to sample: '
      read(5,*) nMCports
      if (nMCports+3 > maxNoPorts) then
         write(6,*) 'Error:: max number of ports: ',maxNoPorts
      else
         exit
      end if
   end do
   ! Get which ports to sample
   ! Note, it is now assumed that the gas feed stream will always be available to sample (mix tank)
   write(6,'(a,$)') '   Enter ports connected to MCs: '
   read(5,*) (MCports(i),i=1,nMCports)

   ! Find how ports connected to gas mix-tank.
   write(6,'(a,$)') 'Enter port for sampling mix-tank: '
   read(5,*) mixtankPort      
   write(6,'(a,$)') 'Enter port connected to mix-tank vent: '
   read(5,*) mixtankVentPort   
   write(6,'(a,$)') 'Minimum time to purge sample loop with mix-tank vent gas (min): '
   read(5,*) SLpurgeTime
   
   write(6,'(a)') 'Enter headspace gas volumes for MCs (mL): '
   read(5,*) (MCheadSpaceVol(i),i=1,nMCports)
   write(6,'(a,$)') 'Enter gas volume of sample loop (mL): '
   read(5,*) SLvol

   ! Get concentrations of feed gas
   write(6,'(a)') 'Enter feed gas concentrations:'
   write(6,'(a,$)') 'O2  (%): '
   read(5,*) o2Feed
   write(6,'(a,$)') 'CO2 (%): '
   read(5,*) co2Feed
   write(6,'(a,$)') 'CH4 (%): '
   read(5,*) ch4Feed
   
   ! Get time to collect a sample
   write(6,'(a,$)') 'Enter time to purge and sample a single microcosm (min): '
   read(5,*) sampleTime
   
   ! Get time to wait to purge gas samplers
   write(6,'(a,$)') 'Enter how often to initiate sampling of reactors (hrs): '
   read(5,*) TinHours
   loopTime = TinHours*60.0 ! loopTime is stored as min., but input from user is in hours.
   if (sampleTime*(nMCports+2) > loopTime) write(6,*) 'Warning, total MC sample time exceeds MC sampling loop time!' 
   
   ! Aks user if gas cycling of feed 2 is currently running from previous program run.
   ! This could occur if the program was stopped while running gas cycling and the recovery
   ! file is either not present or was not used by user.
   write(6,'(a,$)') 'Is gas cycling of feed 2 currently running or start it now? (y/n): '
   read(5,'(a)') lineRead
   gasCycleRunning = .false.
   if (scan(lineRead,'yY') /= 0) gasCycleRunning = .true.      
   if (gasCycleRunning) call gasCyclingParams

   ! Get date for zero time
   write(6,'(a,$)') 'Enter month, day, year of Tzero (return to use today): '
   read(5,'(a)') lineRead
   if (lineRead == ' ') then
      call DATE_AND_TIME (values=dt)
      t0day = dt(3)
      t0month = dt(2)
      t0year = dt(1)
   else
      read(lineRead,*) t0month, t0day, t0year
   end if
   ! Set time zero based on above date.  This will be used to ease tecplot graphing
   call julday(t0year,t0month,t0day,0,0,0.0,tzero,ifail)

   ! Get number to start notes
   write(6,'(a$)') 'Number to begin note count (hit return to start at 1): '
   read(5,'(a)') lineRead
   if (lineRead == ' ') then
      noteCnt = 0
   else
      read(lineRead,*) noteCnt
   end if
   
   return
end subroutine getProgramRunParams

subroutine gasCyclingParams
   ! This gets the conditions for gas cycling
   use CommParameters
   implicit none   
   ! local declarations
   character lineRead*80
   integer dt(8), ifail, i
   
   ! Get start time of gas cycling
   write(6,'(a,$)') 'Enter Julian day of start of gas cycling (return to start now): '
   read(5,'(a)') lineRead
   if (lineRead == ' ') then
      call DATE_AND_TIME (values=dt) ! get date and time
      call julday(dt(1),dt(2),dt(3),dt(5),dt(6),real(dt(7)),t0GasCycling,ifail) ! convert to julian day
   else
      read(lineRead,*) t0GasCycling
   end if
   
   write(6,'(a,$)') 'How many MCs will have cycled feed gas? '
   read(5,*) ncycledPorts
   write(6,'(a,$)') 'Which MCs will have cycled feed gas? ' 
   read(5,*) (cycledPorts(i),i=1,ncycledPorts)
   write(6,'(a,$)') 'Enter port for feed gas 2 sampling: ' 
   read(5,*) feed2GasPort
   write(6,'(a,$)') 'Enter number of days to cycle between feed gases: '
   read(5,*) tCycleOnOff
   write(6,'(a)') 'Enter feed 2 gas concentrations:'
   write(6,'(a,$)') 'O2  (%): '
   read(5,*) o2Feed2
   write(6,'(a,$)') 'CO2 (%): '
   read(5,*) co2Feed2
   write(6,'(a,$)') 'CH4 (%): '
   read(5,*) ch4Feed2
   
   return
end subroutine gasCyclingParams        

subroutine setIports (iports, nports, feed2On)
   ! This routine sets up the list of ports to sample depending on whether the feed2 gas
   ! is currently on.
   use CommParameters, only: mixtankPort, feed2GasPort, mixtankVentPort, &
                             MCports,cycledPorts,nMCports,ncycledPorts
   implicit none
   integer nports, iports(*)
   logical feed2On
   
   ! local declarations
   integer i, j, nonCycled
   
   ! executable code
   iports(1) = mixtankPort ! is always sampled first to obtain instrument drift
   
   if (feed2On) then
      ! Set both cycled and non cycled ports, as well as the sampling of feed 2 gas
      nonCycled = 0
      MCscan: do i=1,nMCports
         do j=1,ncycledPorts
            if(MCports(i) == cycledPorts(j)) cycle MCscan
         end do
         nonCycled = nonCycled + 1
         iports(nonCycled+1) = MCports(i)         
      end do MCscan
      iports(nonCycled+2) = feed2GasPort ! this samples the feed gas 2 port
      do j=1,ncycledPorts
         iports(nonCycled+2+j) = cycledPorts(j)
      end do
      nports = nMCports + 2
      if (nonCycled+2+ncycledPorts /= nports) then
         write(6,'(a)') 'ERROR:: cycled ports not part of MC port list; STOPPING!!'
         stop
      end if
   else
      ! Just include the list of MC ports, and no need to sample feed gas 2
      do i=1,nMCports
         iports(i+1) = MCports(i)
      end do
      nports = nMCports + 1
   end if
         
   return
end subroutine setIports

subroutine initializePorts()
   !This routine sets up communication to equipment attached to RS 232 Serial Ports.
   ! 20 Mar 2010 Ver 2:    Changed setup of valco port.
   ! 03 Aug 2010 Ver 2.1:  Added setup of WTDOT board
   use IFWINTY 
   use IFPORT
   use CommParameters
   implicit None
   integer iresult
   integer baud, parity, dbits, sbits


!  Connect to all ports
   iresult = SPORT_CONNECT (oxigrafPort,DL_TERM_CRLF) !Oxigraf output terminated by crlf
   iresult = SPORT_CONNECT (ch4Port) 
   ! For Valco, toss CR from string on read, adds CR to write string, and end a read when a CR is encountered in buffer.
   iresult = SPORT_CONNECT (valcoPort,(DL_TOSS_CR .or. DL_OUT_CR .or. DL_TERM_CR))  
   iresult = SPORT_CONNECT (mfmPort  ,(DL_OUT_CR .or. DL_OUT_LF .or. DL_TERM_CRLF)) !MKS 660B I/O Terminated with CRLF
   iresult = SPORT_CONNECT (WTDOTport,(DL_TOSS_CR .or. DL_OUT_CR .or. DL_TERM_CR)) 

!  Purge the serial ports
   iresult = SPORT_PURGE (oxigrafPort, (PURGE_TXABORT .or. PURGE_RXABORT .or. PURGE_TXCLEAR .or. PURGE_RXCLEAR)) 
   iresult = SPORT_PURGE (ch4Port,     (PURGE_TXABORT .or. PURGE_RXABORT .or. PURGE_TXCLEAR .or. PURGE_RXCLEAR)) 
   iresult = SPORT_PURGE (valcoPort,   (PURGE_TXABORT .or. PURGE_RXABORT .or. PURGE_TXCLEAR .or. PURGE_RXCLEAR)) 
   iresult = SPORT_PURGE (mfmPort,     (PURGE_TXABORT .or. PURGE_RXABORT .or. PURGE_TXCLEAR .or. PURGE_RXCLEAR)) 
   iresult = SPORT_PURGE (WTDOTport,   (PURGE_TXABORT .or. PURGE_RXABORT .or. PURGE_TXCLEAR .or. PURGE_RXCLEAR)) 

!  First Cancel communications to all ports
   iresult = SPORT_CANCEL_IO (oxigrafPort)
   iresult = SPORT_CANCEL_IO (ch4Port)
   iresult = SPORT_CANCEL_IO (valcoPort)
   iresult = SPORT_CANCEL_IO (mfmPort)
   iresult = SPORT_CANCEL_IO (WTDOTport)

!  Set communcation specifics on all open ports
   baud = 9600; parity = 0; dbits = 8; sbits = 0 !(note, sbits = 0 means stop bits = 1)
   iresult = SPORT_SET_STATE (oxigrafPort , baud, parity, dbits, sbits)
   iresult = SPORT_SET_STATE (ch4Port     , baud, parity, dbits, sbits)
   iresult = SPORT_SET_STATE (valcoPort   , baud, parity, dbits, sbits)
   iresult = SPORT_SET_STATE (mfmPort     , baud, parity, dbits, sbits)
   iresult = SPORT_SET_STATE (WTDOTport   , baud, parity, dbits, sbits)

!  Set the oxigraf to report once
!  Only print info out once (ie. report period)
   iresult = SPORT_WRITE_DATA (oxigrafPort, ESC//'P0;', 0)

return
end subroutine initializePorts

subroutine openUnit (iunit, promptStr)
   ! This opens a file for data output.
   implicit none
   integer iunit, iovar ! which unit to open
   character promptStr*(*) ! String to prompt user with
   ! Local declarations
   character fname*80, yesno*3
   logical fileExists
   
   do
      write(6,'(a,$)') promptStr
      read(5,'(a)') fname
      ! See if the files exist
      inquire(file=trim(fname), exist=fileExists)
      if (fileExists) then
         write(6,'(a,$)') 'File exists, overwrite (y/n): '
         read(5,'(a)') yesno
         if (scan(yesno,'yY') /= 0) then
            open(unit=iunit,file=trim(fname),status='old', iostat=iovar)
            if (iovar /= 0) then
               write(6,'(a)') 'Warning, could not open file, probably invalid filename...'
               cycle
            end if
            exit
         end if
      else
         open(unit=iunit,file=trim(fname),status='new', iostat=iovar)
         if (iovar /= 0) then
            write(6,'(a)') 'Warning, could not open file, probably invalid filename...'
            cycle
         end if
         exit
      end if
   end do
   return
end subroutine openUnit

subroutine readAndStoreAnalzers(iunit)
   ! This routine repetitively reads the gas anlyzers and stores their
   ! value in a user specified file.  Reading stops only when the user
   ! hits the ESC key.
   ! 24 Aug 10:   Ver 2, changed elements of datevec for Ver 3.0 of main changes.
   ! 12 Sep 10:   Ver 3.04, this calls readGasAnalyzersNC instead now.
   use ifwin
   use ifWINTY 
   use ifcore
   use CommParameters, only: ESC, feed2GasPort, mixtankVentPort
   use IFPORT
   implicit none
   integer iunit
   ! Local declarations
   integer readFreq ! Time between analyzer reads (sec).   
   integer iport, dt(8), ifail, j, iread
   character fname*80, lineRead*80, longStr*113, key*1
   real(8) datavec(21), tzero
   logical sampleFeed, keyHit, fileExists, turnPumpOn, turnFeed2On
   ! Data declarations for console cursor access
   integer fhandle
   logical logstat
	Type (T_CONSOLE_SCREEN_BUFFER_INFO) conbuf
   Type (T_COORD)        dwCursorPosition

   ! Get the handle of the console
   fhandle = GetStdHandle(STD_OUTPUT_HANDLE)

   do ! allow the user to rerun this routine   
      call openUnit (iunit, 'Enter name of file to store data: ')   
      write(6,'(a)') 'Enter description for file header and end with a blank line'
      do 
         write(6,'(a,$)') '> '
         read(5,'(a)') lineRead
         if (lineRead == ' ') exit
         write(iunit,'(a)') '# '//lineRead
      end do
      write(iunit,'(a)') ' '
      write(6,'(a,$)') 'Enter delay in seconds between sampling: '
      read(5,*) readFreq
      write(6,'(a,$)') 'Run with sample-loop pump on? (y/n): '
      read(5,'(a)') lineRead
      turnPumpOn = .true. ! use on as default
      if (scan(lineRead,'nN') /= 0) turnPumpOn = .false. 
           
      write(6,'(a,$)') 'Turn feed gas 2 sample valve on? (y/n): '
      read(5,'(a)') lineRead
      turnFeed2On = .false. ! use off as default
      if (scan(lineRead,'yY') /= 0) turnFeed2On = .true. 
      
      write(6,'(a,$)') 'Enter desired Valco port to sample and hit return to start sampling: '
      read(5,*) iport
      if (turnPumpOn .and. iport == mixtankVentPort) then
         write(6,'(a,$)') 'WARNING:: Do you really want to place a vacuum on the mix-tank? (y/n): '
         read(5,'(a)') lineread
         if ( .not.(scan(lineRead,'yN') /= 0) ) cycle     
      end if                 
      
      call valcoControl(iport)
      datavec(14) = dble(iport)
      write(6,'(/,a)') '*** Hit ESC to stop sampling ***'
      ! Set time zero based on current time and date.
      call DATE_AND_TIME (values=dt)
      call julday(dt(1),dt(2),dt(3),dt(5),dt(6),real(dt(7)),tzero,ifail)
      write(iunit,'(a)') 'Variables = "Time (min)" "Flow (sccm)"  "O2 (%)"  "CO2(%)"  "CH4(%)"'
      write(iunit,'(a)') 'Zone'
   
      ! main loop that reads analyzers, reports to screen and stores values      
      if (turnPumpOn) call pumpOn () ! turn pump on first
      if (turnFeed2On) call f2SampleGasOn ! open valve to sample feed 2 gas.
      do
         call cycleFeedGas2 ! This is here to insure that feed gas 2 gets switch over if 
                            ! sample recording goes for a long period.
         call readGasAnalyzers (iport,datavec,tzero,.true.)
         ! Store and display values
         write(iunit,'(f8.2,1x,f6.2,3(1x,f8.4))') 1440.*datavec(13), datavec(15), datavec(17), datavec(19), datavec(21)
         ! get the cursor position
         logstat = GetConsoleScreenBufferInfo(fhandle, conbuf)
         ! write current reading out at the same location
         longStr = '(''Run t:'',f6.1,,'' min'','' Port:'',i2,'' Flow: '',f6.2,'' O2:'',f6.2,'' CO2:'',f7.3,'' CH4:'',f8.4)'
         write(6,longStr)  real(1440.*datavec(13)), int(datavec(14)), real(datavec(15)), (real(datavec(j)),j=17,21,2)
         ! put the cursor back to where it was
         logstat = SetConsoleCursorPosition(fhandle, conbuf.dwCursorPosition)
         ! Go to sleep      
         call SLEEPQQ (1000*readFreq)
         ! See if key was hit
         keyHit = PEEKCHARQQ ( )
         if (keyHit) then
            !See if key hit is ESC, if is so exit.
            key = GETCHARQQ( )
            if (key == ESC) exit
         end if
      end do
      call pumpOff () !turn gas pump off
      call f2SampleGasOff ! close valve to sample feed 2 gas
      
      close(unit=iunit)
      write(6,'(/a,$)') 'Read and store data from a specified port again? (y/n): '
      read(5,'(a)') lineread
      if (scan(lineRead,'nN') /= 0) exit     
   end do
   return
end subroutine readAndStoreAnalzers

subroutine readCH4 (ch4Ave)
   ! This routine reads CH4 concentration from CAI device.
   ! Ver 2: 18 Mar 2010; I have notice problems with reading the CAI device, also
   ! I want to average several reads together.
   ! Once the Measured concentration string is sent (AKON K1), the data string returned will 
   ! look something like:
   !  STX_AKON 0 2.900818 919879468ETXSTX_AKON 0 2.900819 919879469ETX
   ! Where its starts with an STX character and is terminated by an ETX char.
   ! If no errors are returned, then the first 8 characters following the STX are:
   !  "_AKON 0 "
   ! which is then followed by two numbers, the first being the gas concentration and
   ! the second a time stamp in 1/10's of a second.  Also, I've noticed that the 
   ! device seems to send two (or maybe one and a half) return responces.
   ! This current version parses the last completed AKON returned.
   !
   ! Ver 2.1 25 Jun 2010
   ! If the analyzer output exceeds its range (i.e., 0-5%), it will trigger an error
   ! so the response will be
   ! "_AKON # "
   ! where # is any digit between 1 and 9. If there are no errors, then a 0 is returned.
   ! Since I'm now running gas mix at 4.9%, instrument drift can cause it to exceed 5.0,
   ! but the data is still good, so I have removed the check only to look for the first 6
   ! characters as in: 
   ! "_AKON "
   use ifport
   use CommParameters
   implicit none
   real ch4Ave

   ! Local declaratiosn
   integer, parameter:: aveReads = 10 ! Number of samples to average together.
   integer, parameter:: maxTime = 1 ! maximum time to wait for a response (sec).                                 
   integer i, ok2read ,byteCnt, byteCntP, iresult, j, nPts
   integer iSTX, iETX
   real ch4, timeStamp
   real tstart, tend
   character dataStr*1024

   nPts = 0
   do i=1,aveReads
      ! Commands start with STX and are terminated with ETX
      iresult = SPORT_WRITE_DATA (ch4Port, STX//' AKON K1 '//ETX, 11)
      ! Since SPORT_READ_DATA will hold execution until at least on char is in the buffer, make sure
      ! there is one, since if the device dies, we don't want the whole program to hang
      ! here. 
      ! It appears that it can take some time ( ~0.2 sec) for the port to respond
      ! so wait a bit for the response to occur
      call cpu_time(tstart)
      tend = tstart
      do while (tend-tstart <= maxTime)
         iresult = SPORT_PEEK_DATA (ch4Port, ok2read ,byteCntP) 
         if (byteCntP >= 20) exit ! 20 characterts is probably not enough, but it could be if timestamp is only 1 char long.
         call cpu_time(tend)
      end do
      if (.not. ok2read) cycle ! data not in buffer, try again.
      ! Read the entire buffer (or at least dataStr bytes of)
      iresult = SPORT_READ_DATA (ch4Port, dataStr, byteCnt)
      
      ! Starting from the end, find the first complete response between STX and ETX
      iETX = index(trim(dataStr),ETX, back=.true.)
      if (iETX < 20) cycle ! There is not enough info in the string for a complete read, so try again
      iSTX=index(dataStr(1:iETX),STX, back=.true.)
      if (iSTX == 0) cycle ! STX could not be found in sub string, start again
      ! A full string record has been found, make sure response is "_AKON "  
      ! if (index(dataStr(iSTX:iETX),'_AKON 0 ') == 2) then
      if (index(dataStr(iSTX:iETX),'_AKON ') == 2) then
         ! looks good, read values and add to average
         read(dataStr(iSTX+9:iETX-1),*) ch4, timeStamp
         !write(100,'(f8.6,1x,i2,1x,f5.3,1x,i4,1x,a)') ch4, i, tend-tstart, byteCntP,  trim(dataStr) ! remove this 
      else
         ! Error must have been returned in reponse, so try again.  While there may
         ! be good data still in dataStr, just go ahead and read the buffer again.
         cycle
      end if
      
      nPts = nPts + 1
      if (nPts == 1) then ! ch4Ave must be initialized
         ch4Ave = ch4
         cycle
      end if
      ch4Ave = ch4Ave + (ch4 - ch4Ave)/real(nPts) ! running average.
   end do
   
   if (nPts == 0) ch4Ave = badRead
   
   !write(100,'(f8.6,1x,I4)') ch4Ave, nPts ! remove this 
   
   return
end subroutine readCH4

Subroutine readMassFlowMeter (flow)
   ! This routine read the mass flow meter
   ! Ver 2: 20 Mar 2010
   ! Just improving communication
   ! Writes and reads to the 660B PS are terminated by CR/LF.
   use ifport
   use CommParameters

   implicit none
   real flow

   ! Local declarations
   integer, parameter:: maxTime = 1 ! maximum time to wait for a response (sec).  
   integer, parameter:: maxAtmp = 5 ! maximum attempts to read the device.  
   real tstart, tend
   integer iresult, ok2read, byteCnt, byteCntP, i
   character char3*3, dataStr*1024

   flow = badRead
   do i=1,maxAtmp
      ! Command R5<CR/LF> requests the value of the flow.
      ! The responce looks like
      ! P##.###<CR/LF> where the decimal place depends on the setting.
      ! It appears to also take about 0.1-0.2 sec to obtain a response.
      iresult = SPORT_WRITE_LINE (mfmPort, 'R5', 0)
      call cpu_time(tstart)
      tend = tstart
      do while (tend-tstart <= maxTime)
         iresult = SPORT_PEEK_LINE (mfmPort, ok2read, byteCntP)
         if (ok2read) exit
         call cpu_time(tend)
      end do
      if (.not. ok2read) cycle ! data not in buffer, try again.
      iresult = SPORT_READ_LINE (mfmPort, dataStr, byteCnt)
      if (verify(dataStr(2:8),okDigits) /= 0 .OR. len_trim(dataStr(2:8)) == 0) cycle ! Bad read
      read(dataStr(2:8),*) flow
      !write(100,'(f8.4,1x,i2,1x,f5.3,1x,i4,1x,a)') flow, i, tend-tstart, byteCntP,  trim(dataStr) ! remove this 
      exit
   end do
   return
 end Subroutine readMassFlowMeter
    
subroutine readOxigraf(o2Ave, co2Ave)
   ! This routine reads the oxigraf detector for CO2 and O2
   !  oxigraf communcation begins with an ESC and is terminated with a semicolon.
   !  Values that can be read are:
   !     0: System status (16 bit output)
   !     1: Oxygen concentration (0.01%)
   !     2: Sample cell pressure (0.1 mBar)      
   !     3: Sampel cell temperature (0.01 C)
   !     4: Sample flow rate (ml/min)
   !     5: Time stamp counter (9.2 ms?)
   !     6: Alarms (16 bit output)
   !     7: CO2 concentration (0.01%)
   !     8: CO2 cell pressure (0.1 mm Hg)
   !     9: CO2 cell temperature (0.01 C)
   ! Ver 2: 19 Mar 2010
   ! Improving communication and averaging outputs.
   ! The responce from any command sent to the oxigraf is:
   !     C:#####<CR><LF>
   ! Where C is the command sent, #### is any data requested by command
   ! and <CR> and <LF> are the carrage return and line feed (two characters).
   ! Hence, 

   use ifport
   use CommParameters
   real o2Ave, co2Ave

   ! Local declaratiosn
   integer, parameter:: aveReads = 10 ! The number of reads to average together
   integer, parameter:: maxTime = 1 ! maximum time to wait for a response (sec).  
                                  
   integer i, ok2read ,byteCnt, byteCntP, iresult, nPts, lenStr
   real tstart, tend, o2, co2
   character dataStr*100
   

   nPts = 0
   do i=1,aveReads
      ! Get O2, CO2, Cell prssure, and flowrate.  
      !call purgePorts ()    
      iresult = SPORT_WRITE_DATA (oxigrafPort, ESC//'R1,7;', 0)
      call cpu_time(tstart)
      tend = tstart
      do while (tend-tstart <= maxTime)
         iresult = SPORT_PEEK_LINE (oxigrafPort, ok2read ,byteCntP)
         if (ok2read) exit
         call cpu_time(tend)
      end do
      if (.not. ok2read) cycle ! data not in buffer, try again.
      iresult = SPORT_READ_LINE (oxigrafPort, dataStr, byteCnt)
      lenStr = len_trim(dataStr)
      if (dataStr(1:2) == 'R:' .and. dataStr(lenStr-1:lenStr) == crlf) then ! well formed response
         read(dataStr(3:lenStr-2),*) o2, co2
         o2 = o2/100.0
         co2 = co2/100.00
         !write(100,'(2(f8.4,1x),i2,1x,f5.3,1x,i4,1x,a)') o2, co2, i, tend-tstart, byteCntP,  trim(dataStr) ! remove this 
      else
         ! error reading line
         cycle
      end if

      nPts = nPts + 1
      if (nPts == 1) then ! ch4Ave must be initialized
         o2Ave  =  o2
         co2Ave = co2
         cycle
      end if
      o2Ave  =  o2Ave + ( o2 -  o2Ave)/real(nPts) ! running average.
      co2Ave = co2Ave + (co2 - co2Ave)/real(nPts) ! running average.
   end do
   
   if (nPts == 0) then
      o2Ave  = badRead
      co2Ave = badRead
   end if
   
   ! write(100,'(2(f8.4,1x),I4)') o2Ave, co2Ave, nPts ! remove this 
   
   return
end subroutine readOxigraf

Subroutine valcoControl(portNo)
   ! This routine changes the port number of the valco valvue
   ! Ver 2: 20 Mar 2010
   ! Updating serial port read protocal.
   ! The Valco uses only <CR> to end commands and responces
   use ifport
   use CommParameters

   implicit none
   integer portNo

   ! Local declarations
   integer, parameter:: maxTime = 5 ! maximum time to wait for a response (sec). 
                                    ! It takes about 2.5 sec to rotate 1/2 way around valve.  
   integer, parameter:: maxAtmp = 5 ! maximum attempts to read the device.  
   real tstart, tend
   integer iresult, ok2read, byteCnt, portAt, i
   character char4*4, dataStr*1024

   portAt = badRead
   do i=1,maxAtmp
      ! Goto port portNo
      char4 = 'GO'
      write(char4(3:4),'(i2)') portNo
      ! Go to specified port number.
      iresult = SPORT_WRITE_LINE (valcoPort, char4, 0)
      ! Note, an error will be return, such as:
      !     GO 3 = Bad command
      ! if port is already on port 3., so this needs to be flushed from
      ! the buffer.
      
      ! Read port number
      iresult = SPORT_WRITE_LINE (valcoPort, 'CP', 0)
      call cpu_time(tstart)
      tend = tstart
      do while (tend-tstart <= maxTime)
         ! It appears the CP command reponse does not occur until valve
         ! is at its final destination, so this makes the program wait
         ! unitl the valve has stopped moving, unless maxTime is exceeded.
         iresult = SPORT_PEEK_LINE (valcoPort, ok2read ,byteCnt)
         if (ok2read) exit
         call cpu_time(tend)
      end do
      if (.not. ok2read) cycle ! No data in buffer, try again.
      do while (ok2read)
         ! This will continue to read records until buffer is empty
         ! Only the last record in the buffer is used below. 
         iresult = SPORT_READ_LINE (valcoPort, dataStr, byteCnt)
         iresult = SPORT_PEEK_LINE (valcoPort, ok2read ,byteCnt)
      end do
      
      !write(6,'(a,f5.3)') 'Wait time: ', tend-tstart
      !write(6,'(a,i1,a)') 'i=',i,' dataStr: '//trim(dataStr)//'END'
      
      if (dataStr(1:15) == 'Position is  = ' .AND. scan(dataStr(16:17),'0123456789') /= 0) read(dataStr(16:17),*) portAt   
      if (portAt == portNo) exit
   end do

   portNo = portAt
   return
end subroutine valcoControl

Subroutine readWeederT(outStr)
   ! This routine reads outStr from a Weeder Tech module.
   ! Notes:
   !  You must prefix a command with the header character as set by the dipswitch on the board
   !  All command must be terminated by a <CR> character.  Likewise, all data returned is terminated with <CR>
   !  If a command does not return a value, then it appears to echo the command.  Consequently, all input
   !     have some kind of output.
   !  Spaces are not allowed.
   !  Example input: ARB<CR>  !this asks for the state of channel B on board with prefix A
   !  Example responce to above: ABH<CR>
   !  By specifing SPORT_CONNECT and using LINE reads and write, the <CR> character can be easily handled.
   use ifport
   use CommParameters
   implicit none
   character outStr*(*) ! input command string, and output of responce if requrested

   ! Local declarations
   integer, parameter:: maxTime = 1 ! maximum time to wait for a response (sec). 
   integer, parameter:: maxAtmp = 5 ! maximum attempts to read the device.  
   real tstart, tend
   integer iresult, ok2read, byteCnt, i
   character char4*4

   outStr = 'E'  ! set to E in case read fails.
   do i=1,maxAtmp
      ! Read a module, but allow maxTime sectonds for response to be place in buffer.
      call cpu_time(tstart)
      tend = tstart
      do while (tend-tstart <= maxTime)
         iresult = SPORT_PEEK_LINE (WTDOTport, ok2read ,byteCnt)
         if (ok2read==1) exit
         call cpu_time(tend)
      end do
      if (ok2read==0) cycle ! No data in buffer, try again.
      do while (ok2read==1)
         ! This will continue to read records until buffer is empty
         ! Only the last record in the buffer is used below. 
         iresult = SPORT_READ_LINE (WTDOTport, outStr, byteCnt)
         iresult = SPORT_PEEK_LINE (WTDOTport, ok2read ,byteCnt)
      end do
      exit
   end do

   return
end subroutine readWeederT

Subroutine writeWeederT(inStr)
   ! This routine write inStr to a Weeder Tech module.
   ! Only the last record in the buffer is returned.
   use ifport
   use CommParameters
   implicit none
   character inStr*(*) ! input command string, and output of responce if requrested

   ! Local declarations
   integer iresult
   
   iresult = SPORT_WRITE_LINE (WTDOTport, trim(inStr), 0)
   if (iresult /= 0) write(6,*) 'WARNING writing to WTDOT::iresult = ',iresult
   return
end subroutine writeWeederT

subroutine pumpOn ()
   ! This routine turns the gas pump for the gas sampling loop on
   use CommParameters
   implicit none
   integer ierr
   
   character outStr*80
   
   ! write to the WTDOT board
   ! First turn N2 gas on
   call writeWeederT(N2on)
   ! get responce, which would just be the return of command
   call readWeederT(outStr)
   if (outStr /= N2on) then
      write(6,*)      'WARNING:: error turning N2 On'
      write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
   end if         
   
   call writeWeederT(gasPumpOn)
   ! get responce, which would just be the return of command
   call readWeederT(outStr)
   if (outStr /= gasPumpOn) then
      write(6,*)       'WARNING:: error on gas pump On'
      write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
   end if         
return
end subroutine pumpOn

subroutine pumpOff ()
   ! This routine turns the gas pump for the gas sampling loop off
   use CommParameters
   implicit none
   
   character outStr*80
   
   ! write to the WTDOT board   
   ! Turn gas pump off
   call writeWeederT(gasPumpOff)
   ! get responce, which would just be 
   call readWeederT(outStr)
   if (outStr /= gasPumpOff) then
      write(6,*)      'WARNING:: error on gas pump Off'
      write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
   end if         
   
   ! Turn N2 gas off
   call writeWeederT(N2off)
   ! get responce, which would just be the return of command
   call readWeederT(outStr)
   if (outStr /= N2off) then
      write(6,*)      'WARNING:: error turning N2 Off'
      write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
   end if         
   
   return
end subroutine pumpOff

subroutine f2SampleGasOn
   ! This routine turns the Feed 2 gas on for sampling
   use CommParameters, only: f2SmpOn
   implicit none

   character outStr*80
   
   ! Turn ASCO valve on that controls feed 2 gas to sample port
   call writeWeederT(f2SmpOn) 
   ! get responce, which would just be the input string
   call readWeederT(outStr)
   if (outStr /= f2SmpOn) then
      write(6,*)      'WARNING:: error on feed 2 gas sample valve on'
      write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
   end if         
   
   return
end subroutine f2SampleGasOn

subroutine f2SampleGasOff
   ! This routine turns the Feed 2 gas off for sampling
   use CommParameters, only: f2SmpOff
   implicit none

   character outStr*80
   
   ! Turn ASCO valve off that controls feed 2 gas to sample port
   call writeWeederT(f2SmpOff) 
   ! get responce, which would just be the input string
   call readWeederT(outStr)
   if (outStr /= f2SmpOff) then
      write(6,*)      'WARNING:: error on feed 2 gas sample valve off'
      write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
   end if         
   
   return
end subroutine f2SampleGasOff

subroutine cycleFeedGas2 
   ! This routine changes the feed gas based on the current time
   ! the date the experiment started, and the cycle time period.
   ! Note, if gasCycleRunning is false, then state of feed 2 valve is NOT changed.
   ! Modifications
   !  Feed 2 gas is turned ON for even days.
   use CommParameters, only: t0GasCycling, tcycleOnOff, feed2Off, feed2On, co2Feed, o2Feed, ch4Feed, &
                             co2Feed2, o2Feed2, ch4Feed2, co2Feed2cur, o2Feed2cur, ch4Feed2cur, &
                             Feed2GasOn, onOffUnit, tzero, gasCycleRunning
   implicit none

   !local declarations
   character outStr*80
   integer dt(8), ifail, interval
   real(8) eltime, sysTime

   ! if not cycling feed gas 2, then set default values and return
   ! Note, state of feed 2 valve is not changed
   if (.not. gasCycleRunning) then
      co2Feed2cur = co2Feed
      o2Feed2cur  = o2Feed
      ch4Feed2cur = ch4Feed
      return
   end if
   
   ! Get time stamp and julian time
   call DATE_AND_TIME (values=dt)
   call julday(dt(1),dt(2),dt(3),dt(5),dt(6),real(dt(7)),eltime,ifail)
   sysTime = eltime - tzero ! Time relative to start of experiment
   eltime = eltime - t0GasCycling ! elapsed time in days
   interval = int(eltime/tcycleOnOff)
   ! If interval is odd, then main feed gas should be on, if even, then feed gas
   ! 2 should be turned on.
   if (mod(interval,2) /= 0 ) then
      ! In an odd interval, so main gas should be on and feed 2 valve turned off
      call writeWeederT(feed2Off)
      ! get responce, which should just be the return of command
      call readWeederT(outStr)
      if (outStr /= feed2Off) then
         write(6,*)      'WARNING:: error turning Feed gas 2 valve Off'
         write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
      end if         
      ! Update the current feed gas composition
      co2Feed2cur = co2Feed
      o2Feed2cur  = o2Feed
      ch4Feed2cur = ch4Feed
      if (Feed2GasOn) then
         write(6,'(///,a,2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,/)') 'Feed gas 2 turned off at: ', &
                            dt(2),dt(3),dt(1),dt(5),dt(6),dt(7)
         ! Store the time when feed gas 2 was turned off
         write(onOffUnit,'(2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,1x,f9.4,1x,i1)') dt(2),dt(3),dt(1),dt(5),dt(6),dt(7), sysTime, 0                           
      end if
      Feed2GasOn = .false.
   else
      ! In an even interval, so feed gas 2 should be on and valve turned on
      ! Update the current feed gas composition
      call writeWeederT(feed2On)
      ! get responce, which should just be the return of command
      call readWeederT(outStr)
      if (outStr /= feed2On) then
         write(6,*)      'WARNING:: error turning Feed gas 2 valve On'
         write(6,'(2a)') ' WTDOT Output response was: ',trim(outStr)
      end if         
      co2Feed2cur = co2Feed2
      o2Feed2cur  = o2Feed2
      ch4Feed2cur = ch4Feed2
      if (.not. Feed2GasOn) then 
         write(6,'(///,a,2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,/)') 'Feed gas 2 turned on at: ', &
                            dt(2),dt(3),dt(1),dt(5),dt(6),dt(7)
         ! Store the time when feed gas 2 was turned on
         write(onOffUnit,'(2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,1x,f9.4,1x,i1)') dt(2),dt(3),dt(1),dt(5),dt(6),dt(7), sysTime, 1                           
      end if
      Feed2GasOn = .true.
   end if
   return
end subroutine cycleFeedGas2

subroutine forceFeedGas2Valve(offOn) 
   ! This routine changes the feed gas 2 valve on or off based on input.
   ! Note, if gasCycleRunning is true, then valve state will return to its correct
   ! status once control returns to main sample loop.  However, if gasCycleRunning is
   ! false, the valve will retain it's status upon returning to main loop.
   use CommParameters, only: t0GasCycling, tcycleOnOff, feed2Off, feed2On, co2Feed, o2Feed, ch4Feed, &
                             co2Feed2, o2Feed2, ch4Feed2, co2Feed2cur, o2Feed2cur, ch4Feed2cur, &
                             Feed2GasOn, onOffUnit, tzero, gasCycleRunning
   implicit none
   integer offOn ! 0: turn valve off, 1: turn valve on

   !local declarations
   character outStr*80
   integer dt(8), ifail, interval
   real(8) eltime, sysTime
   
   ! Get time stamp and julian time
   call DATE_AND_TIME (values=dt)
   call julday(dt(1),dt(2),dt(3),dt(5),dt(6),real(dt(7)),eltime,ifail)
   sysTime = eltime - tzero ! Time relative to start of experiment
   if (offOn == 0) then
      ! Turn valve off
      call writeWeederT(feed2Off)
      ! get responce, which should just be the return of command
      call readWeederT(outStr)
      if (outStr /= feed2Off) write(6,*) 'WARNING:: error turning Feed gas 2 valve Off'
      ! Update the current feed gas composition
      co2Feed2cur = co2Feed
      o2Feed2cur  = o2Feed
      ch4Feed2cur = ch4Feed
      if (Feed2GasOn) then
         write(6,'(/,a,2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,/)') 'Feed gas 2 turned off at: ', &
                            dt(2),dt(3),dt(1),dt(5),dt(6),dt(7)
         ! Store the time when feed gas 2 was turned off
         write(onOffUnit,'(2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,1x,f9.4,1x,i1)') dt(2),dt(3),dt(1),dt(5),dt(6),dt(7), sysTime, 0                           
      end if
      Feed2GasOn = .false.
   else
      ! Turn valve on
      call writeWeederT(feed2On)
      ! get responce, which should just be the return of command
      call readWeederT(outStr)
      if (outStr /= feed2On) write(6,*) 'WARNING:: error turning Feed gas 2 valve On'
      ! Update the current feed gas composition
      co2Feed2cur = co2Feed2
      o2Feed2cur  = o2Feed2
      ch4Feed2cur = ch4Feed2
      if (.not. Feed2GasOn) then 
         write(6,'(/,a,2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,/)') 'Feed gas 2 turned on at: ', &
                            dt(2),dt(3),dt(1),dt(5),dt(6),dt(7)
         ! Store the time when feed gas 2 was turned on
         write(onOffUnit,'(2(i2.2,''/''),i4,1x,2(i2.2,'':''),i2.2,1x,f9.4,1x,i1)') dt(2),dt(3),dt(1),dt(5),dt(6),dt(7), sysTime, 1                           
      end if
      Feed2GasOn = .true.
   end if
   return
end subroutine forceFeedGas2Valve