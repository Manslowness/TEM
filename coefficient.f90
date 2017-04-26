module coefficient
    implicit none
    dimension cccos(12),W(1:140),t(200),h(20)
    double precision cccos,W
    double precision Ie,L,q,N!电流值,发射线圈边长，接受线圈面积，线圈匝数
    double precision t!时间序列
    double precision Tmax,Tmin,h
    integer Layer,tline
    end module