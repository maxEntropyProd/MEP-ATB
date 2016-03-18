      subroutine julday(iye,mon,idy,ihr,min,sec,julian,ifail)
C ******************************************************************************
C *                                                                            *
C *                            FORTRAN SOURCE CODE                             *
C *                                                                            *
C *   PROGRAM SET :  UTILITY                           REF:JRH:23:05:1991      *
C *                                                                            *
C *   REVISION    :  -------------                         JRH:--:--:1991      *
C *                                                                            *
C *   SOURCE      :  FORLIB.FVS                                                *
C *   ROUTINE NAME:  JULDAY                                                    *
C *   TYPE        :  SUBROUTINE                                                *
C *                                                                            *
C *   FUNCTION    :  Converts (year, month, day, hour, minute, second) to      *
C *                  Julian Day (based on Numerical Recipes)                   *
C *                                                                            *
C *                  IYE ........ Year                                         *
C *                  MON ........ Month                                        *
C *                  IDY ........ Day                                          *
C *                  IHR ........ Hour                                         *
C *                  MIN ........ Minute                                       *
C *                  SEC ........ Second                                       *
C *                  JULIAN ..... Julian Day (non-negative, but may be         *
C *                               non-integer)                                 *
C *                  IFAIL ...... 0 for successful execution, otherwise 1      *
C *                                                                            *
C ******************************************************************************
      real*8 julian
      real*4 sec
      integer*4 iye,mon,idy,ihr,min,ifail
      integer*4 iyyy,jy,jm,igreg,ja,ijul
      integer*4 idint2
      parameter (igreg=15+31*(10+12*1582))
C                              ..... Gregorian Calendar was adopted 15 Oct. 1582
      if(iye.eq.0.or.                                    ! There is no year zero
     $   iye.lt.-4713) then                    ! Julian Day must be non-neagtive
        ifail=1
        return
      endif
      if(iye.lt.0) then
        iyyy=iye+1
      else
        iyyy=iye
      endif
      if(mon.gt.2) then
        jy=iyyy
        jm=mon+1
      else
        jy=iyyy-1
        jm=mon+13
      endif
      ijul=idint2(365.25d0*dble(jy))+idint2(30.6001d0*dble(jm))
     $       +idy+1720995
      if(idy+31*(mon+12*iyyy).ge.igreg) then
C                                    ..... Test for change to Gregorian Calendar
        ja=idint(0.01d0*dble(jy))
        ijul=ijul+2-ja+idint(0.25d0*dble(ja))
      endif
      julian=dble(ijul)
     $         +dble(ihr)/24.d0+dble(min)/1440.d0+dble(sec)/86400.d0
     $         -0.5d0                         ! Correction from midnight to noon
      if(julian.lt.0.d0) then                  ! Julian Day must be non-negative
        ifail=1
        return
      endif
      ifail=0
      return
      end

      integer*4 function idint2(x)
C ******************************************************************************
C *                                                                            *
C *                            FORTRAN SOURCE CODE                             *
C *                                                                            *
C *   PROGRAM SET :  UTILITY                           REF:JRH:23:04:1991      *
C *                                                                            *
C *   REVISION    :  -------------                         JRH:--:--:1991      *
C *                                                                            *
C *   SOURCE      :  FORLIB.FVS                                                *
C *   ROUTINE NAME:  IDINT2                                                    *
C *   TYPE        :  FUNCTION                                                  *
C *                                                                            *
C *   FUNCTION    :  As INT, but trunctates towards -(infinity)                *
C                    This is D.P. version of INT2                              *
C *                                                                            *
C ******************************************************************************
      real*8 x
      idint2=idint(x)
      if(dble(idint2).eq.x) return
      if(x.lt.0.d0) idint2=idint2-1
      return
      end
