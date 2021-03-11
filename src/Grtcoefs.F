!c
!c "GRT_COEFS.F"
!c
        SUBROUTINE GRT_COEFS(KN,O)

        include   "green_com.inc"

        COMPLEX*16   a(4,4),b(4,4),ab(4,4),ed(4,4),eu(4,4),cu,cd
        COMPLEX*16   tem1, tem2, tem3
        REAL*8       hu, hd, count

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   1. Calculate modified R/T coefficients(j=1,2,...,N)  ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c---------------------------------------   
!c        (1). SH case(cf.(31) and (32)).
!c---------------------------------------

        do lay=1,nly-1
            cu=csn(lay)*mu(lay)
            cd=csn(lay+1)*mu(lay+1)
            tdd0(lay)=2.*cu/(cu+cd)
            rdu0(lay)=(cu-cd)/(cu+cd)
            rud0(lay)=-rdu0(lay)
            tuu0(lay)=2.*cd/(cu+cd)
        end do

!c---------------------------------------   
!c        (2). PSV case(cf.(41) and (42)).
!c---------------------------------------
        do lay=1,nly-1
          call mtxe(lay,kn,o,eu)
          call mtxe(lay+1,kn,o,ed)
          do j = 1,2
            do i = 1,2
              a(i,j)     = ed(i,j)
              a(2+i,j)   = ed(2+i,j)
              a(i,j+2)   =-eu(i,j+2)
              a(i+2,j+2) =-eu(i+2,j+2)
              b(i,j)     = eu(i,j)
              b(2+i,j)   = eu(2+i,j)
              b(i,j+2)   =-ed(i,j+2)
              b(i+2,j+2) =-ed(i+2,j+2)
            end do
          end do

          call mtxinv_4(a)
          call mtxmtp(4, a, b, ab)

          do j=1,2
            do i=1,2
                tdd(lay,i,j)=ab(i,j)
                rdu(lay,i,j)=ab(i+2,j)
                rud(lay,i,j)=ab(i,j+2)
                tuu(lay,i,j)=ab(i+2,j+2)
            end do
          end do
        end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   2. Calculate generalized R/T coefficients(j=1,2,...,N)  ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c.......GRUD AT THE FREE-SURFACE (LAY=0):
!c       
!c       SH WAVES
        grud0(0) = 1.0

!c       PSV-WAVES:
        call mtxe(1, kn, O, eu)
        do j=1,2
           do i=1,2
              a(i,j) = eu(i+2,j)
              b(i,j) = eu(i+2,j+2)
           end do
        end do

        call mtxinv_2(a)

        do j=1,2
           do i=1,2
             grud(0,i,j) =-(a(i,1)*b(1,j) + a(i,2)*b(2,j))
           end do
        end do

!c.......GRUD AND GTUU FOR LAY=1,2,......,(LS-1)-TH INTERFACES:

        do 1000 lay=1,ls-1,1
!c       --------------------------------------------
        hu = z(lay)-z(lay-1)
        exu(1) = cdexp(-cpn(lay)*hu)
        exu(2) = cdexp(-csn(lay)*hu)
        ex(lay,1) = exu(1)
        ex(lay,2) = exu(2) 

!c    SH-WAVE:
        tem1 = grud0(lay-1)*exu(2)
        tem2 = rdu0(lay)*exu(2)
        tem3 = tdd0(lay)*exu(2)
        gtu0(lay)  = tuu0(lay)/(1-tem1*tem2)
        grud0(lay) = rud0(lay) + tem3*tem1*gtu0(lay)

!c       PSV-WAVES:
        do j=1,2
           do i=1,2
              a(i,j) = - (rdu(lay,i,1)*exu(1)*grud(lay-1,1,j)
     &                 +  rdu(lay,i,2)*exu(2)*grud(lay-1,2,j))*exu(j)
              b(i,j) =   (tdd(lay,i,1)*exu(1)*grud(lay-1,1,j)
     &                 +  tdd(lay,i,2)*exu(2)*grud(lay-1,2,j))*exu(j)
           end do
           a(j,j)   = 1.0 + a(j,j)
        end do

        call mtxinv_2(a)

        do j=1,2
           do i=1,2
              gtu(lay,i,j) = a(i,1)*tuu(lay,1,j)+a(i,2)*tuu(lay,2,j)     
           end do
        end do
          do j=1,2
             do i=1,2
               grud(lay,i,j)  = rud(lay,i,j) + b(i,1)*gtu(lay,1,j)
     &                                     + b(i,2)*gtu(lay,2,j)
             end do
        end do
!c       -------------------------------------
 1000   CONTINUE


!c.......GRDU AND GTDD FOR LAY=(NLY-1)-TH INTERFACE:

!c       call mtxe(nly, kn, o, ed)
        lay = nly-1
        hd  = z(lay+1)-z(lay)
        exd(1) = cdexp(-cpn(lay+1)*hd)
        exd(2) = cdexp(-csn(lay+1)*hd)
        ex(lay,1) = exd(1)
        ex(lay,2) = exd(2)

!c       (1).SH-WAVES:
        grdu0(lay) = rdu0(lay)
        gtd0 (lay) = tdd0(lay)        
      
!c    (2). PSV-WAVES:
        do j = 1,2
           do i = 1,2
              gtd (lay,i,j) = tdd(lay,i,j)
              grdu(lay,i,j) = rdu(lay,i,j)
           end do
        end do
 
!c        do j=1,4
!c            do i=1,4
!c                ed(i,j) = eu(i,j)
!c            end do
!c        end do

!c.......GRDU AND GTDD FOR LAY=NLY-2,....,LS-TH INTERFACES:

        do 2000 lay = nly-2,ls,-1
!c       --------------------------------------------
        hd     = z(lay+1) - z(lay)
        exd(1) = cdexp(-cpn(lay+1)*hd)
        exd(2) = cdexp(-csn(lay+1)*hd)
        ex(lay,1) = exd(1)
        ex(lay,2) = exd(2)
 
!c       SH-WAVES:
        tem1      = rud0(lay)*exd(2)
        tem2      = grdu0(lay+1)*exd(2)
        tem3      = tuu0(lay)*exd(2)
        gtd0(lay) = tdd0(lay)/(1-tem1*tem2)
        grdu0(lay)= rdu0(lay) + tem3*tem2*gtd0(lay)
 
!c       PSV-WAVES:
        do j=1,2
           do i=1,2
              a(i,j) = -(rud(lay,i,1)*exd(1)*grdu(lay+1,1,j)
     &                 + rud(lay,i,2)*exd(2)*grdu(lay+1,2,j))*exd(j)
              b(i,j) =  (tuu(lay,i,1)*exd(1)*grdu(lay+1,1,j)
     &                 + tuu(lay,i,2)*exd(2)*grdu(lay+1,2,j))*exd(j)
           end do
           a(j,j)    = 1d0 + a(j,j)
        end do
        call mtxinv_2(a)
        do j=1,2
           do i=1,2
              gtd(lay,i,j) = a(i,1)*tdd(lay,1,j) + a(i,2)*tdd(lay,2,j)
           end do
        end do
        do i=1,2
           do j=1,2
              grdu(lay,i,j) = rdu(lay,i,j) + b(i,1)*gtd(lay,1,j)
     &                                   + b(i,2)*gtd(lay,2,j)
           end do
        end do
       
!c       -------------------------------------
 2000   CONTINUE

        ex(ls, 1) = cdexp( -cpn(ls)*(z(ls)-z(ls-1)) )
        ex(ls, 2) = cdexp( -csn(ls)*(z(ls)-z(ls-1)) )

        RETURN
        END
        



