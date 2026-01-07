program Mtable
implicit none

integer          :: num(2),j
real, parameter  :: delr=0.002
double precision :: r,f(3),g(3),Rc,la,ep,R0,al,be,gm

open(unit=1,file='table.xvg',status='unknown',action='write')
open(unit=2,file='tabled.xvg',status='unknown',action='write')
write(*,*) 'Input parameters: D(kBT),R0,alpha,beta,lamda,rcut'
read(*,*) ep,R0,al,be,la,Rc
write(*,'(8f10.5)') ep,R0,al,be,la,Rc
if ((ep .le. 0) .or. (R0 .le. 0)) g = 0
write(1,'(7e15.6)') 0, 0, 0, 0, 0, 0, 0
write(2,'(7e15.6)') 0, 0, 0, 0, 0, 0, 0
num(1)=int(0.5*R0/delr)
num(2)=int((Rc+1)/delr)+1
ep=ep*2.477
do j=1,num(1)
    r=delr*j
    gm=al*(3+2*cos(2*3.14159265358979324*r/R0))
    f(1) = 1/((1-la)**2+r**2)**0.5+((1-la)**2+r**2)**0.5/Rc**2-2/Rc
    f(2) = r/((1-la)**2+r**2)**1.5-r/((Rc**2)*((1-la)**2+r**2)**0.5)
    f(3) = (1-la)/((1-la)**2+r**2)**1.5-(1-la)/((Rc**2)*((1-la)**2+r**2)**0.5)
    if ((ep .gt. 0) .and. (R0 .gt. 0)) then
        if (r.le.R0*la) then
            g(1)=la*ep*(exp(gm*(la-r/R0))-2.*exp(0.5*gm*(la-r/R0)))
            g(2)=la*ep*gm/R0*(exp(gm*(la-r/R0))-1.*exp(0.5*gm*(la-r/R0)))
            g(3)=(ep+gm*la*ep)*exp(gm*(la-r/R0))-(2*ep+gm*la*ep)*exp(0.5*gm*(la-r/R0))
        else if (r.gt.R0*la) then
            g(1)=la*ep*(exp( be*(la-r/R0))-2.*exp(0.5* be*(la-r/R0)))
            g(2)=la*ep* be/R0*(exp( be*(la-r/R0))-1.*exp(0.5* be*(la-r/R0)))
            g(3)=(ep+be*la*ep)*exp(be*(la-r/R0))-(2*ep+be*la*ep)*exp(0.5*be*(la-r/R0))
        end if
    end if
    write(1,'(7e15.6)') r,f(1),f(2),g(1),g(2),g(1),g(2)
    write(2,'(7e15.6)') r,f(3),f(2),g(3),g(2),g(3),g(2)
end do
do j=num(1)+1,num(2)
    r=delr*j
    f(1) = 1/((1-la)**2+r**2)**0.5+((1-la)**2+r**2)**0.5/Rc**2-2/Rc
    f(2) = r/((1-la)**2+r**2)**1.5-r/((Rc**2)*((1-la)**2+r**2)**0.5)
    f(3) = (1-la)/((1-la)**2+r**2)**1.5-(1-la)/((Rc**2)*((1-la)**2+r**2)**0.5)
    if ((ep .gt. 0) .and. (R0 .gt. 0)) then
        if (r.le.R0*la) then
            g(1)=la*ep*(exp(al*(la-r/R0))-2.*exp(0.5*al*(la-r/R0)))
            g(2)=la*ep*al/R0*(exp(al*(la-r/R0))-1.*exp(0.5*al*(la-r/R0)))
            g(3)=(ep+al*la*ep)*exp(al*(la-r/R0))-(2*ep+al*la*ep)*exp(0.5*al*(la-r/R0))
        else if (r.gt.R0*la) then
            g(1)=la*ep*(exp( be*(la-r/R0))-2.*exp(0.5* be*(la-r/R0)))
            g(2)=la*ep* be/R0*(exp( be*(la-r/R0))-1.*exp(0.5* be*(la-r/R0)))
            g(3)=(ep+be*la*ep)*exp(be*(la-r/R0))-(2*ep+be*la*ep)*exp(0.5*be*(la-r/R0))
        end if
    end if
    write(1,'(7e15.6)') r,f(1),f(2),g(1),g(2),g(1),g(2)
    write(2,'(7e15.6)') r,f(3),f(2),g(3),g(2),g(3),g(2)
end do
close(1)
close(2)
end program Mtable
