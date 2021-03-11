       program main 
       include   "green_com.inc"

!c  (0). Basic constants:
!c  -----------------------------
        aj = cmplx(0d0, 1d0)
        pi = 4d0*atan(1d0)
        pi2= 2d0*pi

!c  (1). Reading & Checking input:
!c  -------------------------------   
        call green_input
!c  (2). Basic parameters: 
!c  ----------------------
        call green_basic
!c      
        print*,'ndist',ndist
        do i = 1,ndist
            print*,'Cal the ',i,'th record'
            print*,dist(i),t0(i),Outname(i)
            call green(i)
        enddo
        end
