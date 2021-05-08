;step2_substitute
begin
z=asciiread("./bogusvortex.txt", -1, "float")

di_time=asciiread("./di_time.txt", -1, "string")
year=di_time(0)
month=di_time(1)
day=di_time(2)
hour=di_time(3)
minute=di_time(4)
second=di_time(5)

f=addfile("./wrfout_d01_"+str_get_field(year,1," ")+"-"+str_get_field(month,1," ")+"-"+str_get_field(day,1," ")+"_"+str_get_field(hour,1," ")+":"+str_get_field(minute,1," ")+":"+str_get_field(second,1," "),"r")
R=floattointeger(z(4))
lat_c=z(0)
lon_c=z(1)
delta=18
b_t=40
b_t_s=41
pi=3.1415926
;print(R)

v_p=f->P 
v_u=f->U 
v_v=f->V 
v_t=f->T 
v_ph=f->PH  
v_vapor=f->QVAPOR 

;printVarSummary(v_vapor)
ds=dimsizes(v_p)
ntim=ds(0)
nlev=ds(1)
nlat=ds(2)
mlon=ds(3)
;print(nlat)
;calculate FB 
M=round(((111*cos(lat_c))/delta),3)
;print(M)
hM=floattoint(M^(1./2))
q=new(hM+1,"integer")

kkm=new(11,"float")
q(0)=M 
do i=1,hM 
q(i)=floattoint(M^(1./(2*i)))
end do 

km=(/2,3,4,2,5,6,7,2,8,9,2/)

do i=0,10
kkm(i)=(0.5*(1-cos((2./km(i))*pi))^(-1))
end do 
;print(kkm)

pbar=v_p(1,:,:,:)
pB=v_p(1,:,:,:)
ubar=v_u(1,:,:,:)
uB=v_u(1,:,:,:)
vbar=v_v(1,:,:,:)
vB=v_v(1,:,:,:)
tbar=v_t(1,:,:,:)
tB=v_t(1,:,:,:)
phbar=v_ph(1,:,:,:)
phB=v_ph(1,:,:,:)
vaporbar=v_vapor(1,:,:,:)
vaporB=v_vapor(1,:,:,:)

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

pD=v_p(1,:,:,:)
uD=v_u(1,:,:,:)
vD=v_v(1,:,:,:)
tD=v_t(1,:,:,:)
phD=v_ph(1,:,:,:)
vaporD=v_vapor(1,:,:,:)

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

;定位输出中的坐标
xlon=f->XLONG
xlat=f->XLAT

m_i=0
m_j=0
    do i=M,nlat-M-1 
        do j=M,mlon-M-1
           if((xlat(1,i,j)-z(2)).lt.0.1.and.(xlon(1,i,j)-z(3)).lt.0.1) then
                m_i=i
                m_j=j
                end if 
        end do
    end do 
 ;   print(m_i)
 ;   print(m_j)

;计算常量ER和FDbar
;p=new((R+1),"float") 
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

;pE=v_p(1,:,:,:)
;uE=v_u(1,:,:,:)
;vE=v_v(1,:,:,:)
;tE=v_t(1,:,:,:)
;phE=v_ph(1,:,:,:)
;vaporE=v_vapor(1,:,:,:)

;     do i=-R,R
;        do j=-R,R 
;pE(:,m_i+i,m_j+j)=v_p(1,:,m_i+i,m_j+j)-pV(:,m_i+i,m_j+j)
;uE(:,m_i+i,m_j+j)=v_u(1,:,m_i+i,m_j+j)-uV(:,m_i+i,m_j+j)
;vE(:,m_i+i,m_j+j)=v_v(1,:,m_i+i,m_j+j)-vV(:,m_i+i,m_j+j)
;tE(:,m_i+i,m_j+j)=v_t(1,:,m_i+i,m_j+j)-tV(:,m_i+i,m_j+j)
;phE(:,m_i+i,m_j+j)=v_ph(1,:,m_i+i,m_j+j)-phV(:,m_i+i,m_j+j)
;vaporE(:,m_i+i,m_j+j)=v_vapor(1,:,m_i+i,m_j+j)-vaporV(:,m_i+i,m_j+j)
;end do
;end do

;计算FVax

pVS=new((/nlev,R+1/),"float")
pVax=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
pVS(i,k)=0.
pVax(i,k)=0.
end do 
end do 

pVax(:,0)=pV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
pVS(:,k)=pVS(:,k)+pV(:,m_i+k,m_j+i)
end do
do i=-k,k 
pVS(:,k)=pVS(:,k)+pV(:,m_i-k,m_j+i)
end do
do i=-k,k 
pVS(:,k)=pVS(:,k)+pV(:,m_i+i,m_j+k)
end do
do i=-k,k 
pVS(:,k)=pVS(:,k)+pV(:,m_i+i,m_j-k)
end do  
pVax(:,k)=pVS(:,k)/(8.*k)
end do 

uVS=new((/nlev,R+1/),"float")
uVax=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
uVS(i,k)=0.
uVax(i,k)=0.
end do 
end do 

uVax(:,0)=abs(uV(:,m_i,m_j))
do k=1,R 
do i=-k,k 
uVS(:,k)=uVS(:,k)+abs(uV(:,m_i+k,m_j+i))
end do
do i=-k,k 
uVS(:,k)=uVS(:,k)+abs(uV(:,m_i-k,m_j+i))
end do
do i=-k,k 
uVS(:,k)=uVS(:,k)+abs(uV(:,m_i+i,m_j+k))
end do
do i=-k,k 
uVS(:,k)=uVS(:,k)+abs(uV(:,m_i+i,m_j-k))
end do  
uVax(:,k)=uVS(:,k)/(8.*k)
end do 

vVS=new((/nlev,R+1/),"float")
vVax=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
vVS(i,k)=0.
vVax(i,k)=0.
end do 
end do 

vVax(:,0)=abs(vV(:,m_i,m_j))
do k=1,R 
do i=-k,k 
vVS(:,k)=vVS(:,k)+abs(vV(:,m_i+k,m_j+i))
end do
do i=-k,k 
vVS(:,k)=vVS(:,k)+abs(vV(:,m_i-k,m_j+i))
end do
do i=-k,k 
vVS(:,k)=vVS(:,k)+abs(vV(:,m_i+i,m_j+k))
end do
do i=-k,k 
vVS(:,k)=vVS(:,k)+abs(vV(:,m_i+i,m_j-k))
end do  
vVax(:,k)=vVS(:,k)/(8.*k)
end do 

tVS=new((/nlev,R+1/),"float")
tVax=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
tVS(i,k)=0.
tVax(i,k)=0.
end do 
end do 

tVax(:,0)=tV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
tVS(:,k)=tVS(:,k)+tV(:,m_i+k,m_j+i)
end do
do i=-k,k 
tVS(:,k)=tVS(:,k)+tV(:,m_i-k,m_j+i)
end do
do i=-k,k 
tVS(:,k)=tVS(:,k)+tV(:,m_i+i,m_j+k)
end do
do i=-k,k 
tVS(:,k)=tVS(:,k)+tV(:,m_i+i,m_j-k)
end do  
tVax(:,k)=tVS(:,k)/(8.*k)
end do 


phVS=new((/b_t_s,R+1/),"float")
phVax=new((/b_t_s,R+1/),"float")
do i=0,b_t_s-1
do k=0,R 
phVS(i,k)=0.
phVax(i,k)=0.
end do 
end do 

phVax(:,0)=phV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
phVS(:,k)=phVS(:,k)+phV(:,m_i+k,m_j+i)
end do
do i=-k,k 
phVS(:,k)=phVS(:,k)+phV(:,m_i-k,m_j+i)
end do
do i=-k,k 
phVS(:,k)=phVS(:,k)+phV(:,m_i+i,m_j+k)
end do 
do i=-k,k 
phVS(:,k)=phVS(:,k)+phV(:,m_i+i,m_j-k)
end do  
phVax(:,k)=phVS(:,k)/(8.*k)
end do 

vaporVS=new((/nlev,R+1/),"float")
vaporVax=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
vaporVS(i,k)=0.
vaporVax(i,k)=0.
end do 
end do 

vaporVax(:,0)=vaporV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
vaporVS(:,k)=vaporVS(:,k)+vaporV(:,m_i+k,m_j+i)
end do
do i=-k,k 
vaporVS(:,k)=vaporVS(:,k)+vaporV(:,m_i-k,m_j+i)
end do
do i=-k,k 
vaporVS(:,k)=vaporVS(:,k)+vaporV(:,m_i+i,m_j+k)
end do
do i=-k,k 
vaporVS(:,k)=vaporVS(:,k)+vaporV(:,m_i+i,m_j-k)
end do  
vaporVax(:,k)=vaporVS(:,k)/(8.*k)
end do 

;;;
;;;
;;;
;=============================================================
;接下来是替换过程，这个时候会在wrfinput文件里重复上述过程，然后替换
;=============================================================
;;;
;;;
;;;

delete(f)
delete(v_p)
delete(v_u)
delete(v_v)
delete(v_t)
delete(v_ph)
delete(v_vapor)
delete(ds)
delete(m_i)
delete(m_j)

f=addfile("./wrfinput_d01","w")

v_p=f->P 
v_u=f->U 
v_v=f->V 
v_t=f->T
v_ph=f->PH  
v_vapor=f->QVAPOR


ds=dimsizes(v_p)
ntim=ds(0)
nlev=ds(1)
nlat=ds(2)
mlon=ds(3)

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

M=round(((111*cos(lat_c))/delta),3)
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

;定位输出中的坐标
m_i=0
m_j=0
    do i=M,nlat-M-1 
        do j=M,mlon-M-1
           if((xlat(0,i,j)-z(0)).lt.0.1.and.(xlon(0,i,j)-z(1)).lt.0.1) then
                m_i=i
                m_j=j
                end if 
        end do
    end do 
 ;   print(m_i)
 ;   print(m_j)



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

pVS_in=new((/nlev,R+1/),"float")
pVax_in=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
pVS_in(i,k)=0.
pVax_in(i,k)=0.
end do 
end do 

pVax_in(:,0)=pV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
pVS_in(:,k)=pVS_in(:,k)+pV(:,m_i+k,m_j+i)
end do
do i=-k,k 
pVS_in(:,k)=pVS_in(:,k)+pV(:,m_i-k,m_j+i)
end do
do i=-k,k 
pVS_in(:,k)=pVS_in(:,k)+pV(:,m_i+i,m_j+k)
end do
do i=-k,k 
pVS_in(:,k)=pVS_in(:,k)+pV(:,m_i+i,m_j-k)
end do  
pVax_in(:,k)=pVS_in(:,k)/(8.*k)
end do 

uVS_in=new((/nlev,R+1/),"float")
uVax_in=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
uVS_in(i,k)=0.
uVax_in(i,k)=0.
end do 
end do 

uVax_in(:,0)=abs(uV(:,m_i,m_j))
do k=1,R 
do i=-k,k 
uVS_in(:,k)=uVS_in(:,k)+abs(uV(:,m_i+k,m_j+i))
end do
do i=-k,k 
uVS_in(:,k)=uVS_in(:,k)+abs(uV(:,m_i-k,m_j+i))
end do
do i=-k,k 
uVS_in(:,k)=uVS_in(:,k)+abs(uV(:,m_i+i,m_j+k))
end do
do i=-k,k 
uVS_in(:,k)=uVS_in(:,k)+abs(uV(:,m_i+i,m_j-k))
end do  
uVax_in(:,k)=uVS_in(:,k)/(8.*k)
end do 

vVS_in=new((/nlev,R+1/),"float")
vVax_in=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
vVS_in(i,k)=0.
vVax_in(i,k)=0.
end do 
end do 

vVax_in(:,0)=abs(vV(:,m_i,m_j))
do k=1,R 
do i=-k,k 
vVS_in(:,k)=vVS_in(:,k)+abs(vV(:,m_i+k,m_j+i))
end do
do i=-k,k 
vVS_in(:,k)=vVS_in(:,k)+abs(vV(:,m_i-k,m_j+i))
end do
do i=-k,k 
vVS_in(:,k)=vVS_in(:,k)+abs(vV(:,m_i+i,m_j+k))
end do
do i=-k,k 
vVS_in(:,k)=vVS_in(:,k)+abs(vV(:,m_i+i,m_j-k))
end do  
vVax_in(:,k)=vVS_in(:,k)/(8.*k)
end do 

tVS_in=new((/nlev,R+1/),"float")
tVax_in=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
tVS_in(i,k)=0.
tVax_in(i,k)=0.
end do 
end do 

tVax_in(:,0)=tV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
tVS_in(:,k)=tVS_in(:,k)+tV(:,m_i+k,m_j+i)
end do
do i=-k,k 
tVS_in(:,k)=tVS_in(:,k)+tV(:,m_i-k,m_j+i)
end do
do i=-k,k 
tVS_in(:,k)=tVS_in(:,k)+tV(:,m_i+i,m_j+k)
end do
do i=-k,k 
tVS_in(:,k)=tVS_in(:,k)+tV(:,m_i+i,m_j-k)
end do  
tVax_in(:,k)=tVS_in(:,k)/(8.*k)
end do 


phVS_in=new((/b_t_s,R+1/),"float")
phVax_in=new((/b_t_s,R+1/),"float")
do i=0,b_t_s-1
do k=0,R 
phVS_in(i,k)=0.
phVax_in(i,k)=0.
end do 
end do 

phVax_in(:,0)=phV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
phVS_in(:,k)=phVS_in(:,k)+phV(:,m_i+k,m_j+i)
end do
do i=-k,k 
phVS_in(:,k)=phVS_in(:,k)+phV(:,m_i-k,m_j+i)
end do
do i=-k,k 
phVS_in(:,k)=phVS_in(:,k)+phV(:,m_i+i,m_j+k)
end do 
do i=-k,k 
phVS_in(:,k)=phVS_in(:,k)+phV(:,m_i+i,m_j-k)
end do  
phVax_in(:,k)=phVS_in(:,k)/(8.*k)
end do 

vaporVS_in=new((/nlev,R+1/),"float")
vaporVax_in=new((/nlev,R+1/),"float")
do i=0,nlev-1
do k=0,R 
vaporVS_in(i,k)=0.
vaporVax_in(i,k)=0.
end do 
end do 

vaporVax_in(:,0)=vaporV(:,m_i,m_j)
do k=1,R 
do i=-k,k 
vaporVS_in(:,k)=vaporVS_in(:,k)+vaporV(:,m_i+k,m_j+i)
end do
do i=-k,k 
vaporVS_in(:,k)=vaporVS_in(:,k)+vaporV(:,m_i-k,m_j+i)
end do
do i=-k,k 
vaporVS_in(:,k)=vaporVS_in(:,k)+vaporV(:,m_i+i,m_j+k)
end do
do i=-k,k 
vaporVS_in(:,k)=vaporVS_in(:,k)+vaporV(:,m_i+i,m_j-k)
end do  
vaporVax_in(:,k)=vaporVS_in(:,k)/(8.*k)
end do 

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
        if (r.gt.(R/2.)) then 
        wr=exp(-(((2.*(r-R))/(0.4*R))^2))
        else 
        wr=0
        end if  
      v_p_n(0,:,m_i+i,m_j+j)=v_p(0,:,m_i+i,m_j+j)-pVax_in(:,r)+wr*pVax_in(:,r)+(1-wr)*pVax(:,r)
      v_u_n(0,:,m_i+i,m_j+j)=v_u(0,:,m_i+i,m_j+j)+(-uVax_in(:,r)+wr*uVax_in(:,r)+(1-wr)*uVax(:,r))*sign_matlab(v_u(0,:,m_i+i,m_j+j))
      v_v_n(0,:,m_i+i,m_j+j)=v_v(0,:,m_i+i,m_j+j)+(-vVax_in(:,r)+wr*vVax_in(:,r)+(1-wr)*vVax(:,r))*sign_matlab(v_v(0,:,m_i+i,m_j+j))
      v_t_n(0,:,m_i+i,m_j+j)=v_t(0,:,m_i+i,m_j+j)-tVax_in(:,r)+wr*tVax_in(:,r)+(1-wr)*tVax(:,r)
      v_ph_n(0,:,m_i+i,m_j+j)=v_ph(0,:,m_i+i,m_j+j)-phVax_in(:,r)+wr*phVax_in(:,r)+(1-wr)*phVax(:,r)
      v_vapor_n(0,:,m_i+i,m_j+j)=v_vapor(0,:,m_i+i,m_j+j)-vaporVax_in(:,r)+wr*vaporVax_in(:,r)+(1-wr)*vaporVax(:,r)
    end do       
end do 
;print(uVax_in(0,:))
;print(uVax(0,:))
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
delete(ds)
delete(m_i)
delete(m_j)

print("Substitution has been completed!")
end