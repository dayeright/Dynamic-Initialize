begin

info=asciiread("./di_info.txt", -1, "float")
lat_c=info(0)
lon_c=info(1)
di_time=asciiread("./di_time.txt", -1, "string")
year=di_time(0)
month=di_time(1)
day=di_time(2)
hour=di_time(3)
minute=di_time(4)
second=di_time(5)

f=addfile("./wrfout_d03_"+str_get_field(year,1," ")+"-"+str_get_field(month,1," ")+"-"+str_get_field(day,1," ")+"_"+str_get_field(hour,1," ")+":"+str_get_field(minute,1," ")+":"+str_get_field(second,1," "),"r")
;p=f->P
s_p=wrf_user_getvar(f,"slp",-1)

xlon=f->XLONG
xlat=f->XLAT

ds=dimsizes(s_p)
ntim=ds(0)
nlat=ds(1)
mlon=ds(2)

minp=min(s_p(1,:,:))
m_i=0
m_j=0
    do i=0,nlat-1
        do j=0,mlon-1
           if(s_p(1,i,j).eq.minp) then
                m_i=i
                m_j=j
                end if 
        end do
    end do 
lat_co=xlat(1,m_i,m_j)
lon_co=xlon(1,m_i,m_j)

delete(f)
delete(ds)
delete(m_i)
delete(m_j)
delete(xlon)
delete(xlat)
;;;
;;;
f=addfile("./wrfout_d01_"+str_get_field(year,1," ")+"-"+str_get_field(month,1," ")+"-"+str_get_field(day,1," ")+"_"+str_get_field(hour,1," ")+":"+str_get_field(minute,1," ")+":"+str_get_field(second,1," "),"r")

xlon=f->XLONG
xlat=f->XLAT

ds=dimsizes(xlon)
ntim=ds(0)
nlat=ds(1)
mlon=ds(2)

time = 1    ; first time step
u_a      = wrf_user_getvar(f, "U",time)          ; grid point height
v_a      = wrf_user_getvar(f, "V",time)          ; grid point height
u_aa     = u_a(:,:,0:mlon-1)
v_aa     = v_a(:,0:nlat-1,:)
  p      = wrf_user_getvar(f, "pressure",time)   ; total pressure
  u  = wrf_user_interp_level(u_aa,p,500,False)
  v  = wrf_user_interp_level(v_aa,p,500,False)
printVarSummary(u)

    do i=0,nlat-1
        do j=0,mlon-1
           if((xlat(1,i,j)-lat_co).lt.0.1.and.(xlon(1,i,j)-lon_co).lt.0.1) then
                mm_i=i
                mm_j=j
                end if 
        end do
    end do 

D=min((/mm_i,mm_j,abs(mm_i-nlat),abs(mm_j-mlon)/))-5
print(D)
vb=new(D,"float")
do i=1,D
vb(i-1)=(abs(u(mm_i+i,mm_j))+abs(u(mm_i-i,mm_j))+abs(v(mm_i,mm_j+i))+abs(v(mm_i,mm_j-i)))/4.
end do
i=0
do while ((vb(i).ge.6.or.abs(vb(i)-vb(i+1)).ge.0.072))
;do while (abs(vb(i)).ge.6.or.abs(vb(i)-vb(i+1)).ge.0.072)
;do while (vb(i).ge.3)
i=i+1
if (i.gt.(D-2)) then
print("i overloaded")
break
end if
end do
;R1=floattointeger(z(4))
R1=i




lat_co=xlat(1,mm_i,mm_j)
lon_co=xlon(1,mm_i,mm_j)
;print(lat_co)
;print(lon_co)
;print(mm_i)
;print(mm_j)
asciiwrite("./bogusvortex.txt", (/lat_c,lon_c,lat_co,lon_co,R1/))
print(R1)
print("Then substitute")
end
