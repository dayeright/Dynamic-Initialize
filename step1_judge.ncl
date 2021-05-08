;判断bogus vortex
begin

info=asciiread("./di_info.txt", -1, "float")
lat_c=info(0)
lon_c=info(1)

;第一步，确定圆心
f=addfile("./wrfinput_d03","r")

s_p=wrf_user_getvar(f,"slp",-1)
;printVarSummary(s_p)

xlon=f->XLONG
xlat=f->XLAT

ds=dimsizes(s_p)
ntim=ds(0)
nlat=ds(1)
mlon=ds(2)

minp=min(s_p(0,:,:))
m_i=0
m_j=0
    do i=0,nlat-1
        do j=0,mlon-1
           if(s_p(0,i,j).eq.minp) then
                m_i=i
                m_j=j
                end if 
        end do
    end do 
lat_c1=xlat(0,m_i,m_j)
lon_c1=xlon(0,m_i,m_j)
;print(lat_c1)
;print(lon_c1)
;print(m_i)
;print(m_j)

delete(f)
delete(xlon)
delete(xlat)
delete(ds)
delete(nlat)
delete(mlon)
delete(m_i)
delete(m_j)

;;;
;---------改变初始场的中心位置
;;;

f=addfile("./wrfinput_d01","w")

v_p=f->P 
v_u=f->U 
v_v=f->V 
v_t=f->T 
v_ph=f->PH  
v_vapor=f->QVAPOR
xlon=f->XLONG
xlat=f->XLAT

delta=18
b_t=40
b_t_s=41
pi=3.1415926
ds=dimsizes(v_p)
ntim=ds(0)
nlev=ds(1)
nlat=ds(2)
mlon=ds(3)


M=round(((111*cos(lat_c))/delta),3)
;print(M)

;定位输出中的坐标,并确定半径

m_i=0
m_j=0
    do i=M,nlat-M-1 
        do j=M,mlon-M-1
           if((xlat(0,i,j)-lat_c1).lt.0.1.and.(xlon(0,i,j)-lon_c1).lt.0.1) then
                m_i=i
                m_j=j
                end if 
        end do
    end do 
 ;   print(m_i)
 ;  print(m_j)

m_ii=0
m_jj=0
    do i=M,nlat-M-1 
        do j=M,mlon-M-1
           if((xlat(0,i,j)-lat_c).lt.0.1.and.(xlon(0,i,j)-lon_c).lt.0.1) then
                m_ii=i
                m_jj=j
                end if 
        end do
    end do 
 ;   print(m_ii)
  ;  print(m_jj)

time = 0    ; first time step
u_a      = wrf_user_getvar(f, "U",time)          ; grid point height
v_a      = wrf_user_getvar(f, "V",time)          ; grid point height
u_aa     = u_a(:,:,0:mlon-1)
v_aa     = v_a(:,0:nlat-1,:)
  p      = wrf_user_getvar(f, "pressure",time)   ; total pressure
  u  = wrf_user_interp_level(u_aa,p,500,False)
  v  = wrf_user_interp_level(v_aa,p,500,False)
;printVarSummary(u)

D=min((/m_i,m_j,abs(m_i-nlat),abs(m_j-mlon)/))-5
;print(D)
vb=new(D,"float")

do i=1,D
vb(i-1)=(abs(u(m_i+i,m_j))+abs(u(m_i-i,m_j))+abs(v(m_i,m_j+i))+abs(v(m_i,m_j-i)))/4.
end do

i=0
do while (abs(vb(i)).ge.6.or.abs(vb(i)-vb(i+1)).ge.0.072)
;do while (vb(i).ge.3)
i=i+1
if (i.gt.(D-2)) then
;print("i overloaded")
end if
end do
R1=i
;print(R1)
R=R1
;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;

hM=floattoint(M^(1./2))
;print(hM)
q=new(hM+1,"integer")

kkm=new(11,"float")
q(0)=M 
do i=1,hM 
q(i)=floattoint(M^(1./(2*i)))
end do 
;print(q)
km=(/2,3,4,2,5,6,7,2,8,9,2/)


do i=0,10
kkm(i)=(0.5*(1-cos((2./km(i))*pi))^(-1))
end do 
;print(kkm)

uB=v_u(0,:,:,:)
pB=v_p(0,:,:,:)
vB=v_v(0,:,:,:)
tB=v_t(0,:,:,:)
phB=v_ph(0,:,:,:)
vaporB=v_vapor(0,:,:,:)

ubar=v_u(0,:,:,:)
pbar=v_p(0,:,:,:)
vbar=v_v(0,:,:,:)
tbar=v_t(0,:,:,:)
phbar=v_ph(0,:,:,:)
vaporbar=v_vapor(0,:,:,:)

do ii=0,hM
do jj=0,10
 M=q(ii)
    do i=M,nlat-M-1 
        do j=M,mlon-M-1
pbar(:,i,j)=pB(:,i,j)+kkm(jj)*(pB(:,i,j-M)+pB(:,i,j+M)-2*pB(:,i,j))
ubar(:,i,j)=uB(:,i,j)+kkm(jj)*(uB(:,i,j-M)+uB(:,i,j+M)-2*uB(:,i,j))
vbar(:,i,j)=vB(:,i,j)+kkm(jj)*(vB(:,i,j-M)+vB(:,i,j+M)-2*vB(:,i,j))
tbar(:,i,j)=tB(:,i,j)+kkm(jj)*(tB(:,i,j-M)+tB(:,i,j+M)-2*tB(:,i,j))
phbar(:,i,j)=phB(:,i,j)+kkm(jj)*(phB(:,i,j-M)+phB(:,i,j+M)-2*phB(:,i,j))
vaporbar(:,i,j)=vaporB(:,i,j)+kkm(jj)*(vaporB(:,i,j-M)+vaporB(:,i,j+M)-2*vaporB(:,i,j))
end do
end do

    do i=M,nlat-M-1 
        do j=M,mlon-M-1
pB(:,i,j)=pbar(:,i,j)+kkm(jj)*(pbar(:,i-M,j)+pbar(:,i+M,j)-2*pbar(:,i,j))
uB(:,i,j)=ubar(:,i,j)+kkm(jj)*(ubar(:,i-M,j)+ubar(:,i+M,j)-2*ubar(:,i,j))
vB(:,i,j)=vbar(:,i,j)+kkm(jj)*(vbar(:,i-M,j)+vbar(:,i+M,j)-2*vbar(:,i,j))
tB(:,i,j)=tbar(:,i,j)+kkm(jj)*(tbar(:,i-M,j)+tbar(:,i+M,j)-2*tbar(:,i,j))
phB(:,i,j)=phbar(:,i,j)+kkm(jj)*(phbar(:,i-M,j)+phbar(:,i+M,j)-2*phbar(:,i,j))
vaporB(:,i,j)=vaporbar(:,i,j)+kkm(jj)*(vaporbar(:,i-M,j)+vaporbar(:,i+M,j)-2*vaporbar(:,i,j))
end do 
end do
end do 
end do



pD=v_p(0,:,:,:)
uD=v_u(0,:,:,:)
vD=v_v(0,:,:,:)
tD=v_t(0,:,:,:)
phD=v_ph(0,:,:,:)
vaporD=v_vapor(0,:,:,:)

    do i=0,nlat-1 
        do j=0,mlon-1
        pD(:,i,j)=pD(:,i,j)-pB(:,i,j)        
        uD(:,i,j)=uD(:,i,j)-uB(:,i,j)
        vD(:,i,j)=vD(:,i,j)-vB(:,i,j)
        tD(:,i,j)=tD(:,i,j)-tB(:,i,j)
        phD(:,i,j)=phD(:,i,j)-phB(:,i,j)
        vaporD(:,i,j)=vaporD(:,i,j)-vaporB(:,i,j) 
        end do 
    end do 


;=========================================================
;以上计算出了FD，全部算是准备过程，接下来的操作是对半径R内的操作
;=========================================================

;----------------------------------------------------------
E=new((R+1),"float")  
do i=0,R
    E(i)=(exp(-((R-i)^2)/((0.2*R)^2))-exp(-(R^2)/((0.2*R)^2)))/(1-exp(-(R^2)/((0.2*R)^2)))       
end do 

pDS=new(nlev,"float")
pDbar=new(nlev,"float")
do i=0,nlev-1
pDS(i)=0.
end do

do i=-R,R 
pDS(:)=pDS(:)+pD(:,m_i+R,m_j+i)
end do
do i=-R,R 
pDS(:)=pDS(:)+pD(:,m_i-R,m_j+i)
end do
do i=-R,R 
pDS(:)=pDS(:)+pD(:,m_i+i,m_j+R)
end do
do i=-R,R 
pDS(:)=pDS(:)+pD(:,m_i+i,m_j-R)
end do 
pDbar(:)=pDS(:)/(8.*R)

uDS=new(nlev,"float")
uDbar=new(nlev,"float")
do i=0,nlev-1
uDS(i)=0.
end do

do i=-R,R 
uDS(:)=uDS(:)+abs(uD(:,m_i+R,m_j+i))
end do
do i=-R,R 
uDS(:)=uDS(:)+abs(uD(:,m_i-R,m_j+i))
end do
do i=-R,R 
uDS(:)=uDS(:)+abs(uD(:,m_i+i,m_j+R))
end do
do i=-R,R 
uDS(:)=uDS(:)+abs(uD(:,m_i+i,m_j-R))
end do 
uDbar(:)=uDS(:)/(8.*R)

vDS=new(nlev,"float")
vDbar=new(nlev,"float")
do i=0,nlev-1
vDS(i)=0.
end do

do i=-R,R 
vDS(:)=vDS(:)+abs(vD(:,m_i+R,m_j+i))
end do
do i=-R,R 
vDS(:)=vDS(:)+abs(vD(:,m_i-R,m_j+i))
end do
do i=-R,R 
vDS(:)=vDS(:)+abs(vD(:,m_i+i,m_j+R))
end do
do i=-R,R 
vDS(:)=vDS(:)+abs(vD(:,m_i+i,m_j-R))
end do 
vDbar(:)=vDS(:)/(8.*R)

tDS=new(nlev,"float")
tDbar=new(nlev,"float")
do i=0,nlev-1
tDS(i)=0.
end do

do i=-R,R 
tDS(:)=tDS(:)+tD(:,m_i+R,m_j+i)
end do
do i=-R,R 
tDS(:)=tDS(:)+tD(:,m_i-R,m_j+i)
end do
do i=-R,R 
tDS(:)=tDS(:)+tD(:,m_i+i,m_j+R)
end do
do i=-R,R 
tDS(:)=tDS(:)+tD(:,m_i+i,m_j-R)
end do 
tDbar(:)=tDS(:)/(8.*R)
;print(tDbar)

phDS=new(b_t_s,"float")
phDbar=new(b_t_s,"float")
do i=0,b_t_s-1
phDS(i)=0.
end do

do i=-R,R 
phDS(:)=phDS(:)+phD(:,m_i+R,m_j+i)
end do
do i=-R,R 
phDS(:)=phDS(:)+phD(:,m_i-R,m_j+i)
end do
do i=-R,R 
phDS(:)=phDS(:)+phD(:,m_i+i,m_j+R)
end do
do i=-R,R 
phDS(:)=phDS(:)+phD(:,m_i+i,m_j-R)
end do 
phDbar(:)=phDS(:)/(8.*R)

vaporDS=new(nlev,"float")
vaporDbar=new(nlev,"float")
do i=0,nlev-1
vaporDS(i)=0.
end do

do i=-R,R 
vaporDS(:)=vaporDS(:)+vaporD(:,m_i+R,m_j+i)
end do
do i=-R,R 
vaporDS(:)=vaporDS(:)+vaporD(:,m_i-R,m_j+i)
end do
do i=-R,R 
vaporDS(:)=vaporDS(:)+vaporD(:,m_i+i,m_j+R)
end do
do i=-R,R 
vaporDS(:)=vaporDS(:)+vaporD(:,m_i+i,m_j-R)
end do 
vaporDbar(:)=vaporDS(:)/(8.*R)



;计算FV
pV=pD 
uV=uD
vV=vD
tV=tD
phV=phD 
vaporV=vaporD 

do i=-R,R
        do j=-R,R 
        r=max((/abs(i),abs(j)/))    
        if (abs(i).ge.abs(j)) then
      pV(:,m_i+i,m_j+j)=pD(:,m_i+i,m_j+j)-(pD(:,m_i+sign_matlab(i)*R,m_j+j)*E(r)+pDbar(:)*(1-E(r)))
      uV(:,m_i+i,m_j+j)=uD(:,m_i+i,m_j+j)-(abs(uD(:,m_i+sign_matlab(i)*R,m_j+j))*E(r)+uDbar(:)*(1-E(r)))*sign_matlab(uD(:,m_i+i,m_j+j))
      vV(:,m_i+i,m_j+j)=vD(:,m_i+i,m_j+j)-(abs(vD(:,m_i+sign_matlab(i)*R,m_j+j))*E(r)+vDbar(:)*(1-E(r)))*sign_matlab(vD(:,m_i+i,m_j+j))
      tV(:,m_i+i,m_j+j)=tD(:,m_i+i,m_j+j)-(tD(:,m_i+sign_matlab(i)*R,m_j+j)*E(r)+tDbar(:)*(1-E(r)))
      phV(:,m_i+i,m_j+j)=phD(:,m_i+i,m_j+j)-(phD(:,m_i+sign_matlab(i)*R,m_j+j)*E(r)+phDbar(:)*(1-E(r)))
      vaporV(:,m_i+i,m_j+j)=vaporD(:,m_i+i,m_j+j)-(vaporD(:,m_i+sign_matlab(i)*R,m_j+j)*E(r)+vaporDbar(:)*(1-E(r)))
        else
      pV(:,m_i+i,m_j+j)=pD(:,m_i+i,m_j+j)-(pD(:,m_i+i,m_j+sign_matlab(j)*R)*E(r)+pDbar(:)*(1-E(r)))
      uV(:,m_i+i,m_j+j)=uD(:,m_i+i,m_j+j)-(abs(uD(:,m_i+i,m_j+sign_matlab(j)*R))*E(r)+uDbar(:)*(1-E(r)))*sign_matlab(uD(:,m_i+i,m_j+j))
      vV(:,m_i+i,m_j+j)=vD(:,m_i+i,m_j+j)-(abs(vD(:,m_i+i,m_j+sign_matlab(j)*R))*E(r)+vDbar(:)*(1-E(r)))*sign_matlab(vD(:,m_i+i,m_j+j))
      tV(:,m_i+i,m_j+j)=tD(:,m_i+i,m_j+j)-(tD(:,m_i+i,m_j+sign_matlab(j)*R)*E(r)+tDbar(:)*(1-E(r)))
      phV(:,m_i+i,m_j+j)=phD(:,m_i+i,m_j+j)-(phD(:,m_i+i,m_j+sign_matlab(j)*R)*E(r)+phDbar(:)*(1-E(r)))
      vaporV(:,m_i+i,m_j+j)=vaporD(:,m_i+i,m_j+j)-(vaporD(:,m_i+i,m_j+sign_matlab(j)*R)*E(r)+vaporDbar(:)*(1-E(r)))
        end if
	end do       
    end do 

;计算FVax

;substitute

v_p_n=v_p
v_u_n=v_u  
v_v_n=v_v  
v_t_n=v_t 
v_ph_n=v_ph 
v_vapor_n=v_vapor




do i=-R,R
    do j=-R,R 
        r=max((/abs(i),abs(j)/)) 
      v_p_n(0,:,m_i+i,m_j+j)=v_p(0,:,m_i+i,m_j+j)-pV(:,m_i+i,m_j+j)
      v_u_n(0,:,m_i+i,m_j+j)=v_u(0,:,m_i+i,m_j+j)-uV(:,m_i+i,m_j+j)
      v_v_n(0,:,m_i+i,m_j+j)=v_v(0,:,m_i+i,m_j+j)-vV(:,m_i+i,m_j+j)
      v_t_n(0,:,m_i+i,m_j+j)=v_t(0,:,m_i+i,m_j+j)-tV(:,m_i+i,m_j+j)
      v_ph_n(0,:,m_i+i,m_j+j)=v_ph(0,:,m_i+i,m_j+j)-phV(:,m_i+i,m_j+j)
      v_vapor_n(0,:,m_i+i,m_j+j)=v_vapor(0,:,m_i+i,m_j+j)-vaporV(:,m_i+i,m_j+j)
    end do       
end do 

do i=-R,R
    do j=-R,R 
        r=max((/abs(i),abs(j)/)) 
      v_p_n(0,:,m_ii+i,m_jj+j)=v_p_n(0,:,m_ii+i,m_jj+j)+pV(:,m_i+i,m_j+j)
      v_u_n(0,:,m_ii+i,m_jj+j)=v_u_n(0,:,m_ii+i,m_jj+j)+uV(:,m_i+i,m_j+j)
      v_v_n(0,:,m_ii+i,m_jj+j)=v_v_n(0,:,m_ii+i,m_jj+j)+vV(:,m_i+i,m_j+j)
      v_t_n(0,:,m_ii+i,m_jj+j)=v_t_n(0,:,m_ii+i,m_jj+j)+tV(:,m_i+i,m_j+j)
      v_ph_n(0,:,m_ii+i,m_jj+j)=v_ph_n(0,:,m_ii+i,m_jj+j)+phV(:,m_i+i,m_j+j)
      v_vapor_n(0,:,m_ii+i,m_jj+j)=v_vapor_n(0,:,m_ii+i,m_jj+j)+vaporV(:,m_i+i,m_j+j)
    end do       
 end do 
;print(v_vapor_n)
;printVarSummary(v_ph)
f->P=v_p_n   
f->U=v_u_n   
f->V=v_v_n   
f->T=v_t_n
f->PH=v_ph_n  
f->QVAPOR=v_vapor_n
delete(f)
delete(v_p)
delete(v_u)
delete(v_v)
delete(v_t)
delete(v_ph)
delete(v_vapor)
delete(s_p)
delete(minp)
delete(ds)
delete(m_i)
delete(m_j)
delete(m_ii)
delete(m_jj)

print("First relocation has been finished.")

end

