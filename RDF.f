c     This program will calculate g(r) for a system
      PROGRAM raddisfun
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER nmol=864,nbin=200,nframe=1000
      PARAMETER pi=3.14159265
      DIMENSION rat(3,1,nmol*nframe)
      DIMENSION boxl(3)
      DIMENSION rij(3)
      REAL*8 ghist(nbin),rupp,rlow

      sigma=3.4E-10
      boxl(1)=34.7786E-10
      boxl(2)=34.7786E-10
      boxl(3)=34.7786E-10
      rho=nmol/(boxl(1)*boxl(2)*boxl(3))
      dr=(boxl(1)/2.)*(1./nbin)

      OPEN(456,file='gr-864-equil.out',status='unknown')

c Zero the g(r) histogram
      DO i=1,nbin
         ghist(i)=0
      END DO
      
c Read in all the frames
      OPEN(10,file='MCsim.out-frames7',status='unknown')
      DO i=1,nmol*nframe
         READ(10,*)rat(1,1,i),rat(2,1,i),rat(3,1,i)
      END DO
      CLOSE(10)

c Cycle through each frame
      DO k=1,nframe
         istop=nmol*k
         istart=nmol*(k-1)+1

c Calculate the g(r) function
         DO i=istart,istop-1
            DO j=i+1,istop
               xij=rat(1,1,i)-rat(1,1,j)
               xij=xij-NINT(xij/boxl(1))*boxl(1)
               yij=rat(2,1,i)-rat(2,1,j)
               yij=yij-NINT(yij/boxl(2))*boxl(2)
               zij=rat(3,1,i)-rat(3,1,j)
               zij=zij-NINT(zij/boxl(3))*boxl(3)
               gij=((xij)**2.+(yij)**2.+(zij)**2.)**(0.5)
               islab=NINT(gij/dr)+1
               IF(islab.LE.nbin)THEN
                  ghist(islab)=ghist(islab)+2
               END IF
            END DO
         END DO         
      END DO!iframe loop

c Write out the g(r) histogram 
      DO i=1,nbin
         rlow=(4.*rho/3.)*pi*(dr*(i-1))**3
         rupp=(4.*rho/3.)*pi*(dr*(i-1)+dr)**3
         vshell=rupp-rlow
         raxis=(dr*(i-1)+dr/2.)/sigma
      WRITE(456,*)raxis,(ghist(i)/(vshell*nmol*nframe))
      END DO
      
      STOP
      END
      
