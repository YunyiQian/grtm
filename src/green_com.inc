

!c  "GREENMAIN_COM.INC"
!c       

        implicit  none
        integer   NM, N, MAX, i, j, lay, nf1, nf2
        parameter ( NM=25, N=2**14+2, MAX=5 )

!c  R/T coefficients:
        complex*16 rdu(NM,2,2), rud(NM,2,2), tdd(NM,2,2), tuu(NM,2,2)
        complex*16 grdu(NM,2,2),grud(0:NM,2,2),gtu(NM,2,2),gtd(NM,2,2)
        complex*16 rdu0(NM), rud0(NM), tdd0(NM), tuu0(NM)
        complex*16 exu(2), exd(2), ex(NM,2)
        complex*16 grdu0(NM), grud0(0:NM), gtu0(NM), gtd0(NM)
        common     /grt_r/grud, grdu, gtu, gtd, ex
        common     /grt_l/grud0, grdu0, gtd0, gtu0

!c  Basic quantities:
        complex*16 cpn(NM), csn(NM)
        complex*16 aj, o, unit(2,2)
        real*8     kn, pi, pi2
        common     /cn_ps/cpn, csn
        integer    m, m1, mt, nt
        common     /m_n/m1, m, mt, nt
        real*8      L,Twin,oi,df,dt,dk
        common     /basic_1/L,Twin,oi,dt,df,dk
        common     /basic_2/aj, pi, pi2, unit

!c  Quantities related to wave files:
        complex*16 yd(2,2), yu(2,2), yd0, yu0
        common     /ydu/yd, yu, yd0, yu0
        complex*16 ur2,ur1,ur0,ur01,ur02,uf2,uf1,uz2,uz1,uz01,uz02
        common     /urfz/ur2,ur1,ur0,ur01,ur02,uf2,uf1,uz2,uz1,uz01,uz02
        complex*16 su_0(2), sd_0(2), su_1(2), sd_1(2), su_2(2), sd_2(2)
        complex*16 su0_1, su0_2, sd0_1, sd0_2
        common     /SouVec_PSV/su_0, su_1, su_2, sd_0, sd_1, sd_2
        common     /SouVec_SH /su0_1, su0_2, sd0_1, sd0_2
        complex*16 cigma0(2),cigma1(2), cigma2(2), cigma1_sh, cigma2_sh
        common     /cigma/cigma0, cigma1, cigma2, cigma1_sh, cigma2_sh
!c        complex*16 ur_z(nvsp,n/2), uf_z(nvsp,n/2), uz_z(nvsp,n/2)
        complex*16 ur(n/2), uf(n/2), uz(n/2)	
        integer    n_z
        real*8     dz0, z0_vsp
!c        common     /u_vsp/n_z, dz0, z0_vsp, ur_z, uf_z, uz_z
        common     /u_rfz/ur, uf, uz
        
!c  Source & receiver parameters:
        real*8     fc, tou, zs, fai0, z0,r0, taper,shift
		real*8     dist(100),t0(100)
        integer    lo, ls, nly,ndist
        common     /sou_rec/nly, lo, ls,ndist
		common     /source_parameter/zs, z0, r0,dist, fai0, fc, tou, shift,t0 ,taper
        real*8     s0, nx, ny, nz
        common     /single_force/s0, nx, ny, nz
        real*8     D0s, A0s, dip, rake, str, Me, M0
	real*8     mxx, myy, mzz, mxy, mxz, myz
        common     /double_couple /M0, D0s, A0s, dip, rake, str
	common	   /cmt/ mxx, myy, mzz, mxy, mxz, myz
        common     /explosive_source/Me
        complex*16 fpsv01, fpsv02, fpsv1, fpsv2, fsh1, fsh2
        common    /radiation_pattern/fpsv01,fpsv02,fpsv1,fpsv2,fsh1,fsh2
       
!c  Media parameters:
        real*8     mu(NM),z(0:NM),vs0(NM),vp0(NM),qs(NM),rho(NM),qp(NM)
        complex*16 vs(NM), vp(NM)
        common     /c_velocity/vs, vp
        common     /media/z, vp0, vs0, rho, mu, qs, qp

!c  Others:
        real*8      coef, criterion, vmin,vmax,fmax, f1,f2,f3,f4
        common      /others/coef,criterion,vmin,vmax,fmax
        common      /f1234/f1,f2,f3,f4
        character*80 SourceType, Case, Type_STF, Outname(100), OutFormat
        character*80 WinSwitch,WinType
        common     /strings/SourceType,Case,Type_STF,Outname,OutFormat
        common      /filter/ WinSwitch, WinType
