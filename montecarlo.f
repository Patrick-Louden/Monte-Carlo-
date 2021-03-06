c This program is a first attempt at a Monte Carlo Simulation
c It will find the lowest energy configuration for an Argon Simulation

      PROGRAM argon
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*4 rand
      PARAMETER boltz=1.38065d-23,ncycle=50000,nmol=864
      DIMENSION rat(3,1,nmol),srat(3,1,nmol)
      DIMENSION boxl(3)
      DIMENSION rij(3)
      DIMENSION ve(ncycle),vesq(ncycle)

c Define constants for sim.
      de=120.*boltz
      epsil=4.*de
      sigma=3.4D-10      
      boxl(1)=34.7786E-10
      boxl(2)=34.7786E-10
      boxl(3)=34.7786E-10
      temp=94.4
      rmax=0.5D-10
      an=6.0221415E23
      itick=0 !keeps track of # of accepted moves

      OPEN(10,file='argon864.pos',status='old')
      DO i=1,nmol
         READ(10,*)rat(1,1,i),rat(2,1,i),rat(3,1,i)
      END DO
      CLOSE(10)

      OPEN(20,file='MCsim.out-energy',status='unknown')
      OPEN(45,file='MCsim.out-frames',status='unknown')

c Clear the potential energy arrays.
      DO i=1,ncycle
         ve(i)=0
         vesq(i)=0
      END DO

c Icycle controls how many loops of moving the whole system.
      DO icycle=1,ncycle
         
c Loop through all Argon.
         DO i=1,nmol
            
            vlj1=0
            vlj2=0

c Store the initial coordinates of the moved Argon.
            srat(1,1,i)=rat(1,1,i)
            srat(2,1,i)=rat(2,1,i)
            srat(3,1,i)=rat(3,1,i)
            
c Calculate the Lennard-Jones Potentials.
            DO j=1,nmol
               IF(j.NE.i)THEN
                  DO idim=1,3
                     rij(idim)=rat(idim,1,i)-rat(idim,1,j)
                     rij(idim)=rij(idim)
     *                    -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
                  END DO
                  rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
                  sr2=(sigma*sigma)/rsq
                  sr6=sr2**3
                  sr12=sr6**2
                  vij=epsil*(sr12-sr6)
                  vlj1=vlj1+vij
               END IF
            END DO
                        
c Assign new random coordinate for Argon.
            rat(1,1,i)=rat(1,1,i)+(2*RAND(,)-1)*rmax
            rat(2,1,i)=rat(2,1,i)+(2*RAND(,)-1)*rmax
            rat(3,1,i)=rat(3,1,i)+(2*RAND(,)-1)*rmax                 

c Periodic Boundary Conditions.
            DO j=1,nmol
               DO idim=1,3
                  IF (rat(idim,1,j).GT.boxl(idim))THEN
                     rat(idim,1,j)=rat(idim,1,j)-boxl(idim)
                  END IF
                  IF (rat(idim,1,j).LT.0)THEN
                     rat(idim,1,j)=rat(idim,1,j)+boxl(idim)
                  END IF
               END DO
            END DO
            
c Calculate the Lennard-Jones Potentials.
            DO j=1,nmol
               IF(j.NE.i)THEN
                  DO idim=1,3
                     rij(idim)=rat(idim,1,i)-rat(idim,1,j)
                     rij(idim)=rij(idim)
     *                    -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
                  END DO
                  rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
                  sr2=(sigma*sigma)/rsq
                  sr6=sr2**3
                  sr12=sr6**2
                  vij=epsil*(sr12-sr6)
                  vlj2=vlj2+vij
               END IF
            END DO       
            
c Checks the energy against random number.
            IF(vlj2.GT.vlj1)THEN
               drand=RAND(,)
               IF(drand.GE.EXP(-(vlj2-vlj1)/(boltz*temp)))THEN
                  rat(1,1,i)=srat(1,1,i)
                  rat(2,1,i)=srat(2,1,i)
                  rat(3,1,i)=srat(3,1,i)
                  itick=itick+1
               END IF
            END IF
                        
         END DO !Argon Loop

c Calculate the Lennard-Jones Potentials to be written out.
         vlj=0
         DO i=1,nmol
            DO j=i,nmol
               IF(j.NE.i)THEN
               DO idim=1,3
                  rij(idim)=rat(idim,1,i)-rat(idim,1,j)
                  rij(idim)=rij(idim)
     *                 -IDNINT(rij(idim)/boxl(idim))*boxl(idim)
               END DO
               rsq=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
               sr2=(sigma*sigma)/rsq
               sr6=sr2**3
               sr12=sr6**2
               vij=epsil*(sr12-sr6)
               vlj=vlj+vij
               END IF
            END DO
         END DO
         WRITE(20,*)icycle,vlj

         ve(icycle)=vlj
         vesq(icycle)=vlj*vlj

         IF(MOD(icycle,50).EQ.0)THEN
            DO iat=1,nmol
               WRITE(45,*)rat(1,1,iat),rat(2,1,iat),rat(3,1,iat)
            END DO
            IF(MOD(icycle,100).EQ.0)THEN
               PRINT*,icycle
            END IF
         END IF
         
      END DO !cycle Loop

c Compute Cv.
      vljsum=0
      vljsqsum=0
      DO i=1,ncycle
         vljsum=vljsum+ve(i)
         vljsqsum=vljsqsum+vesq(i)
      END DO
      vljav=vljsum/ncycle
      vljsqav=vljsqsum/ncycle
      delV=vljsqav-(vljav*vljav)
      Cv=delV/(boltz*temp*temp)+(3.*nmol*boltz/2.)
      Cv=Cv/nmol*an
      OPEN(33,file='MCsim.out-Cv',status='unknown')
      WRITE(33,*)Cv
      CLOSE(33)

      OPEN(30,file='MCsim.out-pos',status='unknown')
      DO i=1,nmol
         WRITE(30,*)rat(1,1,i),rat(2,1,i),rat(3,1,i)
      END DO

      OPEN(40,file='MCsim.out-effic',status='unknown')
      WRITE(40,*)'There were',itick,'rejected moves out
     *     of',ncycle*nmol,'moves'

      PRINT*,itick,ncycle*nmol

      OPEN(44,file='argon.xyz',status='unknown')
      WRITE(44,*) nmol
      WRITE(44,*)'argon.xyz          0.000'
      DO i=1,nmol
      WRITE(44,*)'Ar',rat(1,1,i)*1D10,rat(2,1,i)*1D10,rat(3,1,i)*1D10
      END DO
      
      CLOSE(44)
      CLOSE(40)
      CLOSE(30)      
      CLOSE(20)
      CLOSE(45)

      STOP
      END
