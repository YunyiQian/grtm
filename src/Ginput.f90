subroutine green_input
    include "green_com.inc"
    integer layernumb
    write(0,'(a)') 'Input nlay src_depth rcv_depth'
    read(*,*) nly,zs,z0
    print *, nly,zs,z0
    do i = 1,nly
        write(0,'(a)') 'Input depth rho Vs Vp Qb Qa for layer'
        read(*,*) z(i-1),rho(i),vs0(i),vp0(i),qs(i),qp(i)
        print *, z(i-1),rho(i),vs0(i),vp0(i),qs(i),qp(i)
    enddo
    print*,'znly',z(nly)
    write(0,'(a)') 'Input m dt taper f1 f2'
    read(*,*) m,dt,taper,f1,f2
    print*, m,dt,taper,f1,f2
    write(0,'(a)') 'Input number of distance to compute'
    read(*,*) ndist
    print*,'ndist',ndist
    do i = 1,ndist
        write(0,'(a)') 'Input dist t0 output_name'
        read(*,'(2f10.3,1x,a)') dist(i),t0(i),Outname(i)
        print *, dist(i),t0(i),Outname(i)
    enddo

    ls = layernumb(zs,nly,z)
    lo = layernumb(z0,nly,z)
    print*,'ls lo',ls,lo
    if(nly.lt.ls.or.ls.lt.1)then
        print*,'Error: layer number'
    endif
    if (nly.eq.ls)then
        if(zs.gt.z(nly-1))then
            print*,'Add one frictious interface'
            z(nly) = zs + 5.0
            vs0(nly+1) = vs0(nly)
               vp0(nly+1) = vp0(nly)
               qs(nly+1)  = qs(nly)
               qp(nly+1)  = qp(nly)
               rho(nly+1) = rho(nly)
               nly=nly+1
               print*,'nly=',nly,'     z(nly-1)=',z(nly-1)
        else
        endif
    endif

end
