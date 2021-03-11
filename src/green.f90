subroutine green(indx)

include   "green_com.inc"
integer indx
integer   fff1, fff2
real*8    f_1, f_2, amp0 ,win ,freal, window, t, ifftcoef
real*8    random1,random2,random3
complex*16   st0, sw, rrr(n), fff(n), zzz(n)
complex*16   integf(n/2,10), integt(n,10)
real*8    err, tmp(n+5), tmpmax
integer   imax, ifmax
integer wc,iblank,idx0
character*80 sacout
real dtout,shiftout,r0out,dfac,tmp1
real,dimension(:),allocatable::outz,outr,outt
wc = m1*(1-taper)
if (wc.LT.1) wc = 1

nf2=m1
nf1=1
         
ifmax = m1
fff1 = 1
fff2 = ifmax
shift = t0(indx)
r0 = dist(indx)
call green_spectra_r(fff1, fff2, integf)


do i = nf1,nf2
    freal = df*(i-1)         
    o     = pi2*freal - aj*oi
    st0   = 1
    !if(i.GT.wc)then
    !    st0 = cos((i-wc)*pi/(m1-wc)/2)
    !endif
    if(i.GT.wc)then
        st0 = 0.5*(1.+cos((i-wc)*pi/(m1-wc+1)))
    endif
    ifftcoef = df*mt/pi2
    integt(i,:) = st0*integf(i,:)*ifftcoef
end do

do i=nf2+1, mt
    integt(i,:) = conjg(integt(mt+2-i,:))
end do

! ifft
do i=1,10
    call fft(integt(:,i),m,-1)
end do

        
iblank = index(Outname(indx),' ')
Outname(indx)(iblank+1:iblank+1) = char(0)
dtout = dt
shiftout=shift
r0out = r0
allocate(outz(mt))
allocate(outr(mt))
allocate(outt(mt))
!!!!!!!!!!!!!!!! G2x       
zzz = -integt(:,5)
rrr = -integt(:,2) 
fff = 0
do i = 1,mt
    t=(i-1)*dt
    !t = i*dt
    outz(i) = real(zzz(i))*exp(oi*t)*exp(oi*shift)
    outr(i) = real(rrr(i))*exp(oi*t)*exp(oi*shift)
    outt(i) = real(fff(i))*exp(oi*t)*exp(oi*shift)
enddo
idx0=47
Outname(indx)(iblank:iblank) = char(idx0+1)
call wrtsac0(Outname(indx),dtout,mt,shiftout,r0out,outz)
Outname(indx)(iblank:iblank) = char(idx0+2)
call wrtsac0(Outname(indx),dtout,mt,shiftout,r0out,outr)
Outname(indx)(iblank:iblank) = char(idx0+3)
call wrtsac0(Outname(indx),dtout,mt,shiftout,r0out,outt)
!!!!!!!!!!!!!!!!!!
        
        idx0 = 47
!!!!!!!!!!!!!!!!  G1x   
        zzz = integt(:,4)
        rrr = integt(:,1)
        fff = -integt(:,3)
        do i = 1,mt
            t=(i-1)*dt
            !t = i*dt
            outz(i) = real(zzz(i))*exp(oi*t)*exp(oi*shift)
            outr(i) = real(rrr(i))*exp(oi*t)*exp(oi*shift)
            outt(i) = real(fff(i))*exp(oi*t)*exp(oi*shift)
        enddo
        Outname(indx)(iblank:iblank) = char(idx0+4)
        call wrtsac0(Outname(indx),dtout,mt,shiftout,r0out,outz)
        Outname(indx)(iblank:iblank) = char(idx0+5)
        call wrtsac0(Outname(indx),dtout,mt,shiftout,r0out,outr)
        Outname(indx)(iblank:iblank) = char(idx0+6)
        call wrtsac0(Outname(indx),dtout,mt,shiftout,r0out,outt)
!!!!!!!!!!!!!!!!!!

end subroutine
