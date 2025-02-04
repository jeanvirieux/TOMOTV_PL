c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  OLD FORTRAN SUBROUTINES ABOUT UNIVERSAL TIME
c
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_time(year,month,day,hour,min,sec,ier) 
      integer year
      integer month
      integer day
      integer hour
      integer min
      real*4 sec
c
c
c
      ier=0
      if(year.lt.1970) ier=1
      if(month.lt.1.or.month.gt.12) ier=1
      if(day.lt.1.or.day.gt.31) ier=1
      if(hour.lt.0.or.hour.gt.23) ier=1
      if(min.lt.0.or.min.gt.59) ier=1
      if(sec.lt.0. .or. sec.gt.60.000) ier=1
      return
      end
      subroutine str2utime1(string,dbltime)  
      character*(*) string
      real*8 dbltime
      integer daysInMonth(12)
      data daysInMonth/31,28,31,30,31,30,31,31,30,31,30,31/
      integer ts_year
      integer ts_month
      integer ts_day
      integer ts_hour
      integer ts_min
      real*4 ts_sec 
      
      integer sec
      integer msec
      integer year
      integer month
      integer day
      real*8 dbsec
c-------------------------- make sure we got a float
      string(20:20)='.'
c-------------------------- extraction 
c    Input format: yyyy.mm.dd-hh.mm.ss.sss               
c    Example       1995.02.12-11:34:45.123  
c                                     ####    unecessary
      read(string,'(i4)') ts_year
      read(string,'(5x,i2)') ts_month
      read(string,'(8x,i2)') ts_day
      read(string,'(11x,i2)') ts_hour
      read(string,'(14x,i2)') ts_min
      read(string,'(17x,f6.3)') ts_sec
c-------------------------- verification
      call check_time(ts_year,ts_month,ts_day,ts_hour,ts_min,ts_sec,ier)
      if(ier.ne.0) then
        write(*,*) ' error in the date '
        stop 
      endif
c--------------------------
      dbsec=0.0d0
      do year=1970,ts_year-1
      daysInMonth(2)=28
      if(leap_year(year).eq.1) daysInMonth(2)=29
      do month=1,12
      dbsec=dbsec+daysInMonth(month)*86400.0d0
      enddo
      enddo

      daysInMonth(2)=28
      if(leap_year(ts_year).eq.1) daysInMonth(2)=29

      do month=1,ts_month-1
      dbsec=dbsec+daysInMonth(month)*86400.0d0
      enddo
c---------------------- month ts_month
      do day=1,ts_day-1
      dbsec=dbsec+86400.0d0
      enddo
      dbsec=dbsec+ts_hour*3600.0d0
      dbsec=dbsec+ts_min*60.0d0
      dbsec=dbsec+ts_sec
      dbltime=dbsec
      return
      end
c====================================================
      integer function leap_year(year)
      integer year
      leap_year=0
c-------------------- not first year of each century if divided by 4 1700
      if(mod(year,4).eq.0.and.mod(year,100).ne.0) then
        leap_year=1
      endif
c-------------------- but selected millenium as 2000 !
      if(mod(year,400).eq.0) leap_year=1
      return
      end
c=====================================================
      subroutine utime2str1(dbltime,string)
      integer ts_year
      integer ts_month
      integer ts_day
      integer ts_hour
      integer ts_min,iannee
      integer leap_year
      real*4 ts_sec
      real*8 dbltime,db
      character*23 string

      do i=1,23
      string(i:i)=' '
      enddo
      ts_sec=dbltime-int(dbltime/60.d0)*60
      db=int(dbltime/60.d0) !  minute 
      ts_min=db-int(db/60.d0)*60
      db=int(db/60.d0)      !  hour
      ts_hour=db-int(db/24.d0)*24
      db=dbltime/86400.d0    ! day
c--------------------------------------- de 1970 : db est en jours
      ts_year=1969
 1000 continue
      ts_year=ts_year+1
      if(leap_year(ts_year).eq.0) then
        if(int(db).lt.365) goto 2000
        db=db-365.d0
       else
        if(int(db).lt.366) goto 2000
        db=db-366.d0
      endif
      goto 1000
c-------------------------------------- year is deduced compyte day
 2000 continue
      julday=int(db)+1.d0     ! first day start at zero 
      call julien(ts_day,ts_month,iannee,ts_year,julday,-1)
      write(string,'(i4,1h.,i2,1h.,i2,1h-,i2,1h:,i2,1h:f6.3)') 
     &              ts_year,ts_month,ts_day,ts_hour,ts_min,ts_sec
      do i=1,23
      if(string(i:i).eq.' ') string(i:i)='0'
      enddo
      return
      end
C=======================================================================
      SUBROUTINE JULIEN(IJ,IM,IA,IY,ID,ISI)
C=======================================================================
C
C     CONVERTI IJ,IM,IA -> IY,ID       POUR ISI=+1
C     CONVERTI IY,ID    -> IJ,IM,IA    POUR ISI=-1
C
C     EX:  IJ,IM,IA = 21 02 80    IY,ID = 1980 52
C        ou           21 02 1980          1980 52
c=======================================================================
C
      DIMENSION IMO(12)
      DATA IMO/31,28,31,30,31,30,31,31,30,31,30,31/
      integer annee4,leap_year
c
      IF(ISI.EQ.1) then
c
c    calcul jour julien
c
         IF(IA.EQ.0)GO TO 20
         iy=annee4(ia)
         IMO(2)=28
         if(leap_year(iy).eq.1) IMO(2)=29
         ID=0
         DO 1 I=1,12
         IF(I.EQ.IM)GO TO 2
 1       ID=ID+IMO(I)
 2       ID=ID+IJ
         else
c
c    calcul jour, mois
c
         IF(IY.EQ.0)GO TO 20
         ia=annee4(iy)
         IMO(2)=28
         if(leap_year(ia).eq.1) IMO(2)=29
         IDL=ID
         DO 3 I=1,12
         IF(IDL.LE.IMO(I))GO TO 4
 3       IDL=IDL-IMO(I)
 4       IJ=IDL
         IM=I
         endif
      RETURN
c
c    annee nulle: tous mis a zero
c
20    IJ=0
      IM=0
      IA=0
      IY=0
      ID=0
      RETURN
      END
      integer function annee4(ian)
c----------------------------------------------------
c	complete une annee codee sur 2 caracteres
c       ou a partir de 1900 (comme gmtime?)
c----------------------------------------------------
      integer ian
c--------------------------------- if higher than 1900 do nothing
      annee4=ian
      if(annee4.gt.1900) return
c--------------------------------- guess of the year
      if(annee4.gt.70) then
        annee4=annee4+1900
      else
        annee4=annee4+2000
      endif
c
      return
      end
c---------------------------------------------------










