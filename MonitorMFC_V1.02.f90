module globalVars 
! This module contains global variables
   implicit none
   integer, parameter:: mfcPort = 8 ! 647C Powersupply/controller is connected to Com Port 8.
   integer, parameter:: noMFC = 8 ! Number of MFC 647C can control (either 4 or 8)
   ! Just handle the first 10 types.
   real(4), parameter:: MFCrangeValues(1:10) = (/1.,2.,5.,10.,20.,50.,100.,200.,500.,1000./)
   integer, parameter:: loopWait = 10 ! number of seconds to wait between reading 647C.
	CHARACTER(LEN=1), PARAMETER :: ESC	= CHAR(27)   
   real(4) MFCranges(noMFC)! Ranges of the MFC recorded in the 647c
   real(4) MFCgcf(noMFC)! Gas correction factors recorded in the 647c
   
   ! Data to map responce of 647C to actual values
   ! Flow ranges:
   ! 0 = 1.000 SCCM, 20 = 1.000 SCFH
   ! 1 = 2.000 SCCM, 21 = 2.000 SCFH
   ! 2 = 5.000 SCCM, 22 = 5.000 SCFH
   ! 3 = 10.00 SCCM, 23 = 10.00 SCFH
   ! 4 = 20.00 SCCM, 24 = 20.00 SCFH
   ! 5 = 50.00 SCCM, 25 = 50.00 SCFH
   ! 6 = 100.0 SCCM, 26 = 100.0 SCFH
   ! 7 = 200.0 SCCM, 27 = 200.0 SCFH
   ! 8 = 500.0 SCCM, 28 = 500.0 SCFH
   ! 9 = 1.000 SLM, 29 = 1.000 SCFM
   ! 10 = 2.000 SLM, 30 = 2.000 SCFM
   ! 11 = 5.000 SLM, 31 = 5.000 SCFM
   ! 12 = 10.00 SLM, 32 = 10.00 SCFM
   ! 13 = 20.00 SLM, 33 = 20.00 SCFM
   ! 14 = 50.00 SLM, 34 = 50.00 SCFM
   ! 15 = 100.0 SLM, 35 = 100.0 SCFM
   ! 16 = 200.0 SLM, 36 = 200.0 SCFM
   ! 17 = 400.0 SLM, 37 = 500.0 SCFM
   ! 18 = 500.0 SLM, 38 = 30.00 SLM
   ! 19 = 1.000 SCMM, 39 = 300.0 SLM
end module globalVars

program MonitorMFC
   ! This program is used to monitor the 647C MFC controller
   ! Ver 1.00 (28 Mar 2010)
   !
   ! Modifications to make (25 Oct 10)
   !  1. Select a gas menu (X, 1 - 4)
   !  2. Turn an MFC on (allow muliple to be specified, as 1,2,3,4,7,8,0)
   !  3. Turn an MFC off (as above)
   !  Make sure to reread MFC setting and factors after a change in the same.
   !  Is it true that no response is returned on a command?
   !  Make sure I'm using the latest algorithm for reading com ports
   !  Package all files into a single file for easier version maintenance.
   !
   ! V 1.02 24 Mar 2011
   !  Finishing implimentation of turning MFCs off/on
   !  Note, selection of gas menu not implemented yet.
    
   use GlobalVars
   use ifport 
   use ifcore
   implicit none
   
   character iostring*1024, wStr*20, longStr*113, key*1
   integer i, j, iValue, dt(8), k, iresult, cntrl
   real(4) MFCflows(noMFC), rValue, lineCnt 
   logical keyHit
   
   write(6,'(a)') 'Version 1.02 (24 Mar 2011)'
   write(6,'(a)') 'Initializing ....'
   
   ! initial port
   call initializePort ()
   
   ! Get the flow ranges for each MFC
   ! Number returned requires table lookup values (see globalVars module).
   wStr = 'RA c R'
   do i=1, noMFC
      write(wStr(4:4),'(i1)') i
      call write647C(trim(wStr))
      call read647c(iostring)
      if ( scan(iostring,'E') /= 0) then
         write(6,*) 'Error getting flow ranges, value = '//trim(iostring)
         stop
      end if
      read(iostring,*) iValue
      MFCranges(i) = MFCrangeValues(iValue+1)
   end do
   
   ! Get the gas correction factors, These are reported at % values.
   wStr = 'GC c R'
   do i=1, noMFC
      write(wStr(4:4),'(i1)') i
      call write647C(trim(wStr))
      call read647c(iostring)
      if ( scan(iostring,'E') /= 0) then
         write(6,*) 'Error getting gas correction factor, value = '//trim(iostring)
         stop
      end if
      read(iostring,*) iValue
      MFCgcf(i) = real(iValue)/100.0
   end do
   
   write(6,'(a)') 'Flow Ranges for MFCs reported by 647C:'
   write(6,'(a)') '   MFC 1   MFC 2   MFC 3   MFC 4   MFC 5   MFC 6   MFC 7   MFC 8'
   write(6,'(8(1x,f7.2)/)')  (MFCranges(j),j=1,noMFC)   
   
   write(6,'(a)') 'Gas Correction Factors for MFCs reported by 647C:'
   write(6,'(a)') '   MFC 1   MFC 2   MFC 3   MFC 4   MFC 5   MFC 6   MFC 7   MFC 8'
   write(6,'(8(1x,f7.2)/)')  (MFCgcf(j),j=1,noMFC)   
   
   ! *******************************************
   ! ********* Begin Sampling Loop *************
   ! *******************************************
   lineCnt = 35
   monitor: do 
      if (lineCnt >= 35) then
         CALL DATE_AND_TIME (values=dt)
         write(6,'(a)') ' *** Hit ESC to bring up menu ***'
         write(6,'(1x,2(i2.2,a),i4,4x,a)') dt(2),'/',dt(3),'/',dt(1),'MFC 1   MFC 2   MFC 3   MFC 4   MFC 5   MFC 6   MFC 7   MFC 8'
         lineCnt = 1
      end if
   
      ! See if ESC key was hit
      keyHit = PEEKCHARQQ ( )
      if (keyHit) then
         !See if key hit is ESC, if is so exit.
         key = GETCHARQQ( )
         if (key == ESC) then
            call menu(cntrl)
            select case (cntrl)
               case(1:8) ! Continue monitoring.
                  write(6,'(/,a)') '               MFC 1   MFC 2   MFC 3   MFC 4   MFC 5   MFC 6   MFC 7   MFC 8'
               case(9) ! exit program
                  exit monitor
            end select ! if neither case is selected, just continue monitoring
         end if
      end if      
      
      ! Get the flow values
      wStr = 'FL #'   
      do i=1, noMFC
         write(wStr(4:4),'(i1)') i
         call write647C(trim(wStr))
         call read647c(iostring)
         if ( scan(iostring,'E') /= 0) then
            write(6,'(a)') 'Error getting flow flows, value = '//trim(iostring)
            MFCflows(i) = -99.99            
            cycle
         end if
         if ( trim(iostring) == '-----' ) then 
            ! assume value is off scale.
            MFCflows(i) = -99.99
         else
            read(iostring,*) iValue
            MFCflows(i) = MFCranges(i)*MFCgcf(i)*real(iValue)/1000.
         end if
      end do
      CALL DATE_AND_TIME (values=dt)

      longStr = '(''At: '',2(i2.2,'':''),i2.2,8(1x,f7.2))'
      write(6,longStr)  (dt(4+j),j=1,3), (MFCflows(j),j=1,noMFC)
      lineCnt = lineCnt + 1
      ! put the cursor back to where it was
      !logstat = SetConsoleCursorPosition(fhandle, conbuf.dwCursorPosition)
      call sleepqq(1000*loopWait)
   end do monitor
      
   !  Release the port
   iresult = SPORT_RELEASE (mfcPort)
   stop
end program MonitorMFC

subroutine changeSetpoint ()
   ! This routine is used to change an MFC set point
   use globalVars
   implicit none
   character ioStr*80, iStr*11
   integer gasMenu, MFC, iValue
   real(4) setPoint, rValue
   
   ! First find out which gas menu is running 
   call write647c('GM R')
   call read647c(ioStr)
   if (scan(ioStr,'E') /= 0 ) then
      write(6,'(a)') 'Error getting MFC gas menu, returning...'
      return
   end if
   read(ioStr,*) gasMenu 
   write(6,'(a,$)') 'Which MFC channel do you want to change SP for: '
   read(5,*) MFC
   
   if (gasMenu == 0) then
      ! 647c is in menu X mode
      ! FS c xxxx
      !  c = 1..8 channel
      !  x = 0..1100 setpoint in 0.1 percent of full scale
      !
      write(6,'(a)') 'Note, gas menu X running'
      ! First get current value:
      iStr = 'FS c R'
      write(iStr(4:4),'(i1)') MFC
      call write647c(iStr)
      call read647c(ioStr)
      if (scan(ioStr,'E') /= 0 ) then
         write(6,'(a)') ' *** Error getting MFC set point, returning...'
         return
      end if
      read(ioStr,*) iValue ! this value is 0.1% of full scale, which depends on GFC also
      rValue = real(iValue)*MFCranges(MFC)*MFCgcf(MFC)/1000.
      
      ! Now set it to user requested value
      write(6,'(3(a,f7.2),a,$)') 'SP currently: ',rValue,' change to (', &
       & MFCranges(MFC)*MFCgcf(MFC)*0.01,'-',MFCranges(MFC)*MFCgcf(MFC)*1.10,'): '
      read(5,*) setPoint
      
      iValue = nint( 1000.*setPoint/(MFCranges(MFC)*MFCgcf(MFC)) )
      write(iStr(6:9),'(i4.4)') iValue
      call write647c(iStr)
      ! Now check it
      iStr(6:9) = 'R   '
      call write647c(iStr)
      call read647c(ioStr)
      if (scan(ioStr,'E') /= 0 ) then
         write(6,'(a)') ' *** Error setting MFC set point, returning...'
         return
      end if
      read(ioStr,*) iValue
      write(6,'(a,i1,a,f8.2)') 'MFC ',MFC,' Set to: ',real(iValue)*MFCranges(MFC)*MFCgcf(MFC)/1000.
   else
      ! a gas menu is being used
      ! GP c s xxxx
      !  c = 1..8 MFC channel
      !  s = 1..5 gas set 1 to 5
      !  x = 0..1100 setpoint in 0.1 percent of full scale
      !
      write(6,'(a,i1,a)') 'Note, gas menu ',gasMenu,' running'
      ! First get current value:
      iStr = 'GP c s R'
      write(iStr(4:4),'(i1)') MFC
      write(iStr(6:6),'(i1)') gasMenu
      call write647c(iStr)
      call read647c(ioStr)
      if (scan(ioStr,'E') /= 0 ) then
         write(6,'(a)') ' *** Error getting MFC set point, returning...'
         return
      end if
      read(ioStr,*) iValue ! this value is 0.1% of full scale, which depends on GFC also
      rValue = real(iValue)*MFCranges(MFC)*MFCgcf(MFC)/1000.

      ! Now set it to user requested value
      write(6,'(3(a,f7.2),a,$)') 'SP currently: ',rValue,' change to (', &
       & MFCranges(MFC)*MFCgcf(MFC)*0.01,'-',MFCranges(MFC)*MFCgcf(MFC)*1.10,'): '
      read(5,*) setPoint
      
      iValue = nint( 1000.*setPoint/(MFCranges(MFC)*MFCgcf(MFC)) )
      write(iStr(8:11),'(i4.4)') iValue
      call write647c(iStr)
      ! Now check it
      iStr(8:11) = 'R   '
      call write647c(iStr)
      call read647c(ioStr)
      if (scan(ioStr,'E') /= 0 ) then
         write(6,'(a)') ' *** Error setting MFC set point, returning...'
         return
      end if
      read(ioStr,*) iValue
      write(6,'(a,i1,a,f8.2)') 'MFC ',MFC,' Set to: ',real(iValue)*MFCranges(MFC)*MFCgcf(MFC)/1000.
   end if
   return
end subroutine changeSetpoint

subroutine MFCsOnOff ()
   ! This routine is used to turn MFCs on or off
   use globalVars
   implicit none
   character MFCstr*80, MFCi*4, MainValveState*7
   integer i, nMFCset, MFCs(9)
   
   ! Turn MFC's on  Note, 0 corresponds to main valve. It must be on
   ! for any flow to occur.
   write(6,'(a)') 'Enter MFC numbers to turn ON. Note, "0" corresponds to main valve,'
   write(6,'(a,$)') 'which must be on for any flow. Enter MFCs (hit return for none): '
   read(5,'(a)') MFCstr
   MFCi = 'ON #'
   nMFCset = 0
   MainValveState = 'UNKNOWN'  ! There is not command to assess the state of the Main Valve
   do i=0,8
      if (scan(trim(MFCstr),char(48+i)) /= 0) then
         ! Turn MFC i on
         write(MFCi(4:4),'(i1)') i
         call write647c(MFCi)
         nMFCset = nMFCset + 1
         MFCs(nMFCset) = i
         if (i == 0) MainValveState = 'ON'
      end if
   end do
   if (nMFCset > 0) then
      write(6,'(/a,9(1x,i1))') 'MFCs set to ON:', (MFCs(i),i=1,nMFCset)
   else
      write(6,'(/a)') 'No MFCs changed to ON'
   end if
      
   ! Turn MFC's OFF  Note, 0 corresponds to main valve. If turned off, no flow occurs
   write(6,'(/,a)') 'Enter MFC numbers to turn OFF. Note, "0" corresponds to main valve.'
   write(6,'(a,$)') 'If set, all flows turned off. Enter MFCs (hit return for none): '
   read(5,'(a)') MFCstr
   MFCi = 'OF #'
   nMFCset = 0
   do i=0,8
      if (scan(trim(MFCstr),char(48+i)) /= 0) then
         ! Turn MFC i off
         write(MFCi(4:4),'(i1)') i
         call write647c(MFCi)
         nMFCset = nMFCset + 1
         MFCs(nMFCset) = i
         if (i == 0) MainValveState = 'OFF'
      end if
   end do
   if (nMFCset > 0) then
      write(6,'(/a,9(1x,i1))') 'MFCs set to OFF:', (MFCs(i),i=1,nMFCset)
   else
      write(6,'(/a)') 'No MFCs changed to OFF'
   end if
   write(6,'(/,a)') 'Note, state of Main Valve: '//trim(MainValveState)

   return
end subroutine MFCsOnOff

subroutine initializePort()
   !This routine sets up communication to equipment attached to RS 232 Serial Ports.
   ! Ver 1: 28 Mar 2010
   use IFWINTY 
   use IFPORT
   use globalVars
   implicit None
   integer iresult
   integer baud, parity, dbits, sbits


   !  Connect to all port
   ! add CR LF to writes, drop CR and LF from reads, and a read with CR/LF is encountered.
   iresult = SPORT_CONNECT (mfcPort,(DL_OUT_CR .or. DL_OUT_LF .or. DL_TERM_CRLF .or. DL_TOSS_CR .or. DL_TOSS_LF)) 

   !  Purge the serial port
   iresult = SPORT_PURGE (mfcPort, (PURGE_TXABORT .or. PURGE_RXABORT .or. PURGE_TXCLEAR .or. PURGE_RXCLEAR)) 
   
   !  First Cancel communications to all ports
   iresult = SPORT_CANCEL_IO (mfcPort)

   !  Set communcation specifics on all open ports
   baud = 9600; parity = 0; dbits = 8; sbits = 0 !(note, sbits = 0 means stop bits = 1)
   iresult = SPORT_SET_STATE (mfcPort , baud, parity, dbits, sbits)

   return
end subroutine initializePort

subroutine menu(optNo)
   ! This routine proves a simple menu for the MFC monitor
   implicit none
   integer optNo
   
   ! Local declarations
   logical getOpt
   
   getOpt = .true.
   opt: do while (getOpt)
      write(6,'(//)') 
      write(6,'(a)') 'Monitoring Paused; Options:'
      write(6,'(a)') '  1. Return to MFC monitoring.'
      write(6,'(a)') '  2. Change the set point of a MFC'
      write(6,'(a)') '  3. Turn ON a MFC or all MFCs'
      write(6,'(a)') '  4. Select gas menu'
      write(6,'(a)') '  9. Exit program'
      read(5,*) optNo
      select case (optNo)
      case(1)
         write(6,*) '  Returning to monitoring. *** Press ESC for menu ***'
         write(6,*) ' '
         getOpt = .false.
      case(2)
         ! Change a MFC set point
         call changeSetpoint ()
         getOpt = .true. ! return to menu.
      case(3)
         ! Turn MFC on.
         call MFCsOnOff ()
         getOpt = .true. ! return to menu.
      case(4)
         ! Select gas menu.
         call selectGasMenu ()
         getOpt = .true. ! return to menu.
      case(9)
         write(6,*) '  Exiting program.' 
         getOpt = .false.
         
      case default
         write(6,*) '  ** Option not recognized, try again **'
      end select
   end do opt
   return
end subroutine menu

Subroutine read647c(outStr)
   ! This routine reads outStr to the 647C MFC controller.
   ! Only the last record in the buffer is returned.
   use ifport
   use globalVars
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
      ! Read the 647C, but allow maxTime sectonds for response to be place in buffer.
      call cpu_time(tstart)
      tend = tstart
      do while (tend-tstart <= maxTime)
         !iresult = SPORT_PEEK_LINE (mfcPort, ok2read ,byteCnt)
         ! It appears if a <CR><LF> is the first thing in the buffer, then
         ! PEEK_LINE returns with ok2read ==0, even if there is data after the CR/LF.
         ! or at least, this appears to be problem in reading 647C after giving a 
         ! command, because the 647C will return just <CR><LF> with no data.
         iresult = SPORT_PEEK_DATA (mfcPort, ok2read ,byteCnt)
         if (ok2read==1) exit
         call cpu_time(tend)
      end do
      if (ok2read==0) cycle ! No data in buffer, try again.
      do while (ok2read==1)
         ! This will continue to read records until buffer is empty
         ! Only the last record in the buffer is used below. 
         iresult = SPORT_READ_LINE (mfcPort, outStr, byteCnt)
         iresult = SPORT_PEEK_LINE (mfcPort, ok2read ,byteCnt)
      end do
      exit
   end do

   return
end subroutine read647c

Subroutine write647c(inStr)
   ! This routine write inStr to the 647C MFC controller.
   use ifport
   use globalVars
   implicit none
   character inStr*(*) ! input command string

   ! Local declarations
   integer iresult
   
   iresult = SPORT_WRITE_LINE (mfcPort, trim(inStr), 0)
   return
end subroutine write647c

subroutine selectGasMenu ()
   ! This routine allows user to set the gas menu
   implicit none
   character ioStr*80, iStr*11, menu*4
   integer gasMenu, MFC, iValue, gasMenuSet   
   
   ! First find out which gas menu is running 
   call write647c('GM R')
   call read647c(ioStr)
   if (scan(ioStr,'E') /= 0 ) then
      write(6,'(a)') 'Error getting MFC gas menu, returning...'
      return
   end if
   read(ioStr,*) gasMenu 
   
   write(6,'(a,i1,a)') 'Note, gas menu ',gasMenu,' running'
   write(6,'(a,$)') 'Enter gas menu to set (1-5 or 0 for X menu): '
   read(5,*) gasMenuSet
   
   menu = 'GM #'
   write(menu(4:4),'(i1)') gasMenuSet
   call write647c(menu)   
   
   ! Check to see if correct
   call write647c('GM R')
   call read647c(ioStr)
   if (scan(ioStr,'E') /= 0 ) then
      write(6,'(a)') 'Error getting MFC gas menu on confirmation check, returning...'
      return
   end if
   read(ioStr,*) gasMenu 
   if (gasMenu /= gasMenuSet) write(6,'(a)') 'WARNING: Gas menu not set as requested. Returning'
   
   return
end subroutine selectGasMenu

