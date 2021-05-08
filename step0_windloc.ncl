;step0_windloc
begin
info=asciiread("./di_info.txt", -1, "float")
uv10_ini_max=info(4)

di_time=asciiread("./di_time.txt", -1, "string")
year=di_time(0)
month=di_time(1)
day=di_time(2)
hour=di_time(3)
minute=di_time(4)
second=di_time(5)

f=addfile("./wrfout_d03_"+str_get_field(year,1," ")+"-"+str_get_field(month,1," ")+"-"+str_get_field(day,1," ")+"_"+str_get_field(hour,1," ")+":"+str_get_field(minute,1," ")+":"+str_get_field(second,1," "),"r")
;f=addfile("./wrfinput_d03","w")
p=f->P
ds=dimsizes(p)
ntim=ds(0)
nlev=ds(1)
nlat=ds(2)
mlon=ds(3)


u10=f->U10(1,:,:) 
v10=f->V10(1,:,:)

s_p=p(1,0,:,:)
;printVarSummary(s_p)


xlon=f->XLONG(1,:,:)
xlat=f->XLAT(1,:,:)
 


 minp=min(s_p(:,:))
 ;print(minp)
 uv10=new((/nlat,mlon/),"float")
 ii=0
 jj=0
 do i=0,nlat-1
        do j=0,mlon-1
   if(s_p(i,j).eq.minp) then
        ii=i
        jj=j
         end if
         uv10(i,j)=sqrt(u10(i,j)^2+v10(i,j)^2)
 end do
 end do
xp_wrf=xlon(ii,jj)
yp_wrf=xlat(ii,jj)

uv10_max=max(uv10)

diff_uv10=(uv10_max-uv10_ini_max)
asciiwrite("./di_judge.txt", diff_uv10)
;print(xp_wrf)
;print(yp_wrf)
print(uv10_max)
print(diff_uv10)



end