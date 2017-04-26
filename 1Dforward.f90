program main
    use coefficient
    dimension V(200),m(15)
    double precision V,m
    double precision test
    m=(/500,500,500,500,500,500,500,800,800,800,800,800,800,800,800/)
    
    test=21.0
    print*,cccos(1),W(1)
    call coeffready
!    V=forward(m)
    print*,cccos(1),W(2)
    print*,t(42)
    print*,t(test*2),t(test*2-1)
    print*,tline
    end
!读取参数文件：
subroutine coeffready()
    use coefficient
    tline=0
    !读取余弦变换滤波系数
    open(11,file='C:\study\code\G_S.txt')
    do i=1,12
    read(11,100) cccos(i)
    enddo
100 format(e22.15)    
    close(11,status='keep')
    !读取汉克尔变换滤波系数
    open(12,file='C:\study\code\J1_140.txt')
    read(12,*)          !skip the first line
    do i=1,140
        read(12,101) W(i)
    enddo
101 format(e18.12)
    close(12,status='keep')
    !读取装置参数
    open(13,file='C:\study\code\parameter.txt')
    read(13,*)
    read(13,102) L
    read(13,*)
    read(13,102) Ie
    read(13,*)
    read(13,102) q
    read(13,*)
    read(13,102) N
    read(13,*)
    read(13,102) Tmax  
    read(13,102) Tmin
    read(13,*)
    read(13,*) Layer 
102 format(f5.1)    
    close(13,status='keep')
    !计算时间序列
    do i=1,200
        t(i)=0.0
    enddo
    
    do i=1,Tmin*10.0-Tmax*10.0+1
        tline=1+tline
        t(i)=10.0**Tmax
        Tmax=Tmax+0.10
    enddo
    
    !计算每一层的厚度
    h(1)=20.0
    do i=2,Layer
        h(i)=h(i-1)**1.05!层厚增长因子
    enddo
    end
    
function forward(m)
    use coefficient
    dimension lambda(140),zz(12,140),z0(12,140),u(300,140),z(300,140),forward(200)
    double precision forward,lambda,zz,z0,u,z
    double precision u1,z1,sumCos,sumHankel,ii,jj,kk
    integer cs
    double precision ws,wa,a,pi,mu0
    pi=3.141592653
    mu0=4*pi*10**(-7)
    sumCos=0
    sumHankel=0
    a=l/sqrt(pi)!发射线圈等效半径
    wa=-7.91001919000d-00
    ws=8.79671439570d-02
    do i=1,140
        lambda(i)=(10**(wa+(i-1)*ws))/a
    enddo
    do ii=1,tline
        do cs=Layer,1,-1
            do jj=1,12
                do kk=1,140
                    z0(jj,kk)=-log(2.0)*jj/t(ii)*mu0/lambda(kk)
                    u((cs-1)*12+jj,kk)=sqrt((lambda(kk)**2)+log(2.0)*jj/t(ii)*mu0/m(cd))
                    z((cs-1)*12+jj,kk)=-log(2.0)*jj*mu0/u((cs-1)*12+jj,kk)/t(ii)
                    if (cs==Layer) then
                        zz(jj,kk)=z((cs-1)*12+jj,kk)
                    endif
                enddo
            enddo
        enddo
        do cs=Layer-1,1,-1
            do jj=1,12
                do kk=1,140
                    u1=u((cs-1)*12+jj,kk)
                    z1=z((cs-1)*12+jj,kk)
                    zz(jj,kk)=z1*(zz(jj,kk)+z1*tanh(u1*h(cs)))/(z1+zz(jj,kk)*tanh(u1*h(cs)))
                enddo
            enddo
        enddo
        do jj=1,12
            sumHankel=0.0
            do kk=1,140
                sumHankel=sumHankel+lambda(kk)*W(kk)*(zz(jj,kk)/(zz(jj,kk)+z0(jj,kk)))
            enddo
            sumCos=sumCos+sumHankel*cccos(jj)
        enddo
        sumCos=sumCos/t(ii)*log(2.0)
        forward(ii)=sumCos*Ie*q*mu0*N
    enddo
end function forward