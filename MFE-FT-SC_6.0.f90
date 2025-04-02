![A].b calculater with sparse efficient algorithm (used in BiCGSTAB subroutine):


subroutine mesh_generator(nx,ny,nt,nb,block_id,ghost_cell_size,&
                         x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                         u_gs,v_gs,p_gs,&
                         uvbc,pbc,winds,&
                         u_w,u_e,u_s,u_n,&
                         v_w,v_e,v_s,v_n,&
                         p_w,p_e,p_s,p_n,&
                         ubc_e,ubc_w,&
                         vbc_n,vbc_s,&
                         pbc_e,pbc_w,pbc_n,pbc_s)
    
integer :: i,j,k,k_start
integer :: nx,ny,nt,nb,block_id,ghost_cell_size

real(kind=8) dx,dy

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter




real(kind=8),dimension(nb)::u_w,u_e,u_s,u_n
real(kind=8),dimension(nb)::v_w,v_e,v_s,v_n
real(kind=8),dimension(nb)::p_w,p_e,p_s,p_n



character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds

!local variable
integer :: step_i_start = 1!i
integer :: step_j_start = 1!j
integer :: step_i_end = 10!i
integer :: step_j_end = 10!j



end subroutine mesh_generator






subroutine boundary_mesh(nx,ny,nt,&
                         x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                         u,u_alter,v,v_alter,&
                         p,p_alter,p_correct,&
                         f_former,f,f_alter,&
                         a_former,a,a_alter,&
                         b_former,b,b_alter,&
                         u_gs,v_gs,p_gs,&
                         uvbc,pbc,winds,&
                         u_w,u_e,u_s,u_n,&
                         v_w,v_e,v_s,v_n,&
                         p_w,p_e,p_s,p_n,&
                         ubc_e,ubc_w,&
                         vbc_n,vbc_s,&
                         pbc_e,pbc_w,pbc_n,pbc_s,&
                         u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                         u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                         p_bound_w,p_bound_e,&
                         p_bound_s,p_bound_n,&
                         f_bound_w,f_bound_e,&
                         f_bound_s,f_bound_n,&
                         a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                         a_bound_s,a_bound_n,b_bound_s,b_bound_n)
    
integer :: i,j
integer :: nx,ny,nt

real(kind=8) dx,dy

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter



real(kind=8)::u_w,u_e,u_s,u_n
real(kind=8)::v_w,v_e,v_s,v_n
real(kind=8)::p_w,p_e,p_s,p_n


real(kind=8),dimension(0:ny+2)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1)::a_bound_s,a_bound_n,b_bound_s,b_bound_n



character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds


dx = (domain_x_end-domain_x_start)/nx
dy = (domain_y_end-domain_y_start)/ny


x(0,:)    = domain_x_start
x(1,:)    = domain_x_start + dx/2.0d0
x(nx,:)   = domain_x_end   - dx/2.0d0
x(nx+1,:) = domain_x_end
    
y(:,0)    = domain_y_start 
y(:,1)    = domain_y_start + dy/2.0d0
y(:,ny)   = domain_y_end   - dy/2.0d0
y(:,ny+1) = domain_y_end

do i=2,nx-1
    x(i,:)=domain_x_start + dx/2.0d0 + 1.0d0*(i-1)*dx
end do

do j=2,ny-1
    y(:,j)=domain_y_start + dy/2.0d0 + 1.0d0*(j-1)*dy
end do
    

!init uv/initial condition
    do i=0,nx+2
        do j=0,ny+2
            u(i,j)=1e-20
            u_alter(i,j)=1e-20
            u_gs(i,j,:)="c"
        end do
    end do
    
    do i=0,nx+2
        do j=0,ny+2
            v(i,j)=1e-20
            v_alter(i,j)=1e-20
            v_gs(i,j,:)="c"
        end do
    end do
!init pf/initial condition
    do i=0,nx+1
        do j=0,ny+1
            p(i,j)=0.0d0
            p_alter(i,j)=0.0d0
            p_gs(i,j,:)="c"
        end do
    end do
    


    do i=0,nx+1
        do j=0,ny+1
            f(i,j)=0.0d0
            f_former(i,j)=0.0d0
            f_alter(i,j)=0.0d0
            
            a(i,j)=0.0d0
            a_former(i,j)=0.0d0
            a_alter(i,j)=0.0d0

            b(i,j)=0.0d0
            b_former(i,j)=0.0d0
            b_alter(i,j)=0.0d0
        end do
    end do
    
print*, "Array initialized!"


end subroutine boundary_mesh



subroutine bicg_stab_solver
  implicit none

  integer, parameter :: n = 5  ! size of the matrix A
  real(8) :: A(n,n), B(n), X(n), R(n), P(n), V(n), S(n), T(n), alpha, beta, rho, rho_old, omega
  integer :: i, j, max_iter
  real(8) :: tol, norm_r

  ! Initialize guess for X
  X = 0.0d0

  ! Set tolerance and maximum number of iterations
  tol = 1.0d-6
  max_iter = 1000

  ! Initial residual
  call mat_vec(A, X, R)
  R = B - R

  ! Initial rho
  rho = dot_product(R, R)
  
  P = R
  V = 0.0d0
  S = 0.0d0
  T = 0.0d0
  omega = 1.0d0
  norm_r = sqrt(rho)

  print *, "Initial residual norm: ", norm_r

  ! BICGSTAB iterations
  do i = 1, max_iter
     rho_old = rho
     rho = dot_product(R, R)
     
     beta = (rho / rho_old) * (alpha / omega)
     P = R + beta * (P - omega * V)
     
     call mat_vec(A, P, V)
     
     alpha = rho / dot_product(R, V)
     S = R - alpha * V
     
     call mat_vec(A, S, T)
     
     omega = dot_product(T, S) / dot_product(T, T)
     X = X + alpha * P + omega * S
     
     R = S - omega * T
     
     norm_r = sqrt(dot_product(R, R))
     print *, "Iteration: ", i, " Residual norm: ", norm_r
     
     if (norm_r < tol) then
        print *, "Convergence reached after ", i, " iterations"
        exit
     end if
  end do

  if (norm_r >= tol) then
     print *, "Solution did not converge after ", max_iter, " iterations"
  end if

!print *, "Solution X: ", X

contains

  subroutine mat_vec(A, X, R)
    real(8), dimension(n,n), intent(in) :: A
    real(8), dimension(n), intent(in) :: X
    real(8), dimension(n), intent(out) :: R
    integer :: i, j

    R = 0.0d0
    do i = 1, n
       R(i) = 0.0d0
       do j = 1, n
          R(i) = R(i) + A(i,j) * X(j)
       end do
    end do
  end subroutine mat_vec

end subroutine bicg_stab_solver






subroutine Iteration_solver(n,A,X,B)
    integer :: i,j,iter,max_iter=1000
    integer :: n
    real(kind=8),dimension(n,n)::A
    real(kind=8),dimension(n)::X,X_new
    real(kind=8),dimension(n)::B
    real(kind=8) :: tolerance=1E-15,error
    X=0.0d0
    X_new=0.0d0

    ! Iterate using the Jacobi method
    do iter = 1, max_iter
        ! Compute the new values of X using the Jacobi formula
        !$OMP PARALLEL DO PRIVATE(i) SHARED(A, X, B, X_new, n)
        do i = 1, n
            X_new(i) = (B(i) - sum(A(i, 1:i-1) * X(1:i-1)) - sum(A(i, i+1:n) * X(i+1:n))) / A(i, i)
        end do
        !$OMP END PARALLEL DO
        
        ! Calculate the error (the difference between the old and new X values)
        error = maxval(abs(X_new - X))

        ! Update the solution vector X
        X = X_new

        ! Check if the solution has converged
        if (error < tolerance) then
            print *, "Converged after ", iter, " iterations."
            !print *, "Solution vector X:"
            !print *, X
            exit
        end if
        
        if ( mod(iter,100).eq. 0) then
            print *, "solver iteration process is: ", iter, " iterations."
        end if
        
        if (iter .eq. max_iter .and. error > tolerance) then
        print *, "WARNING!!! Did not converge after: ", max_iter, " iterations, WARNING!!!"
        !print *, "Last computed solution vector X:"
        !print *, X
        end if

    end do

end subroutine Iteration_solver
    
subroutine matvec(a,ai,aj,nnzero,b,n,c)
	implicit none
	INTEGER, INTENT(inout) :: n !n*n
	INTEGER, INTENT(inout) :: nnzero !none zero
	REAL*8,  INTENT(inout), DIMENSION(1:nnzero) :: a  !A matrix
	INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: ai  !A matrix i
	INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: aj  !A matrix j
	REAL*8,  INTENT(inout), DIMENSION(1:n) :: b !b vector (right hand side)
	REAL*8,  INTENT(inout), DIMENSION(1:n) :: c !Answer vector
	integer i
	c=0.0d0


	do i=1,nnzero
		c(aj(i))=c(aj(i))+a(i)*b(ai(i))
       
	enddo
        
end subroutine

!To solve [A].[x]=[b] system of equations by BiCGSTAB method:
subroutine Bicgstab_solver(a,ai,aj,nnzero,b,n,x,bires,bimaxit)
    implicit none
    INTEGER, INTENT(inout) :: n  !Matrix size (n*n)
    INTEGER, INTENT(inout) :: nnzero !No. of none zero elements in [A]
    real*8,  INTENT(inout) :: bires !Max. Residual to converge
    INTEGER, INTENT(inout) :: bimaxit !Max. iterations
    REAL*8,  INTENT(inout), DIMENSION(1:nnzero) :: a  !A matrix
    INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: ai  !A matrix i
    INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: aj  !A matrix j
    REAL*8,  INTENT(inout), DIMENSION(1:n) :: b !Right hand side
    REAL*8,  INTENT(inout), DIMENSION(1:n) :: x !Answer

    allocatable temp(:),r0(:),r(:),p(:),s(:),rr(:)
    integer i
    real(8) temp,r0,r,p,alpha,s,w,rr,beta,res
    allocate (temp(1:n),r0(1:n),r(1:n),p(1:n),s(1:n),rr(1:n))
    i=0

	!BiCGSTAB Algorithm:
	!1:
    x=0.0 !First Guess
    call matvec(a,ai,aj,nnzero,x,n,temp)
    r0=b-temp

	!2:
    p=r0
    r=r0
    
	!3:
	do
		i=i+1

		!4:
		call matvec(a,ai,aj,nnzero,p,n,temp)
		alpha=(dot_product(r,r0))/(dot_product(temp,r0))

		!5:
		s=r-alpha*temp

		!6:
		call matvec(a,ai,aj,nnzero,s,n,temp)
		w=(dot_product(temp,s))/(dot_product(temp,temp))

		!7:
		x=x+alpha*p+w*s

		!8:
		rr=s-w*temp

		!9:
		beta=((dot_product(rr,r0))/(dot_product(r,r0)))*(alpha/w)
	
		!10:
		call matvec(a,ai,aj,nnzero,p,n,temp)
		p=rr+beta*(p-w*temp)

		!11:
		r=rr

		!12:
		res=dsqrt(dot_product(r,r))/n


		if ((res<bires).or.(i>bimaxit)) then !Convergance or Max iterations Checking
			!print*,"BiCGSTAB Res.=",res ,"by",i,"iters"
			exit
		endif
	enddo

end subroutine

subroutine Gaussian_solver (n, A, x, coe)
!The aim of this sub is to solve Ax=coe

real(kind=8),dimension(n,n) ::A
real(kind=8),dimension(n) ::x
real(kind=8),dimension(n) ::coe

real(kind=8),dimension(n,n) ::A_temp
real(kind=8),dimension(n) ::coe_temp
real(kind=8),dimension(n) ::A_diagonalized

real(kind=8)::temp=0
integer::i,j,k
!duijiaohua

do i=1,n
coe_temp(i) = coe(i)
  do j=1,n
  A_temp(i,j) = A(i,j)
  end do
end do


do k=1,n-1

  do j=k+1,n
  coe_temp(j) = coe_temp(j)-A_temp(k,j)/A_temp(k,k)*coe_temp(k)
  end do
  
  
  do j=k+1,n
  temp=A_temp(k,j)
  
    do i=1,n
    A_temp(i,j)=A_temp(i,j)-temp/A_temp(k,k)*A_temp(i,k)
    end do

  end do
  
end do

  do i=1,n
    A_diagonalized(i)=A_temp(i,i)
    coe_temp(i)=coe_temp(i) / A_diagonalized(i)
  end do
 
  do j=1,n
  do i=j,n
    A_temp(i,j)=A_temp(i,j) / A_diagonalized(j)
  end do
  end do

temp=0
do j=n,1,-1
  do i=n,j+1,-1
  temp=temp+A_temp(i,j)*x(i)
  end do
x(j)=(coe_temp(j)-temp)
temp=0
end do


  !do j=1,n
  !print *, A(:,j),"x",x(j),"=", coe(j)
  !end do

end subroutine Gaussian_solver

subroutine Upwind_calculator(delta,wind,down,center,up,term)
    
real(kind=8)::delta,down,center,up,wind,term
if (wind .ge. 0.0) then 
    term = wind*( center - down )/delta
else
    term = wind*( up - center )/delta
end if

end subroutine Upwind_calculator

subroutine Central_calculator(delta,multiplier,down,center,up,term)
real(kind=8)::delta,down,center,up,multiplier,term
    term = multiplier*( up - down )/2.0/delta
end subroutine Central_calculator

subroutine Upwind_time_calculator(delta,time,wind,down1,center1,up1,down2,center2,up2,term)
    
real(kind=8)::delta,time,down1,center1,up1,down2,center2,up2,wind,term
if (wind .ge. 0.0) then 
    term = wind*( center2 - center1 - down2 + down1 )/delta/time
else
    term = wind*( up2 - up1 - center2 + center1 )/delta/time
end if

end subroutine Upwind_time_calculator

subroutine Central_first_time_calculator(delta,time,multiplier,down1,center1,up1,down2,center2,up2,term)
    
real(kind=8)::delta,time,down1,center1,up1,down2,center2,up2,multiplier,term

    term = multiplier*( ( up2 - down2 ) - ( up1 - down1 ) ) / time / 2.0d0 / delta
    !term = wind*( center2 - center1 ) / time 
end subroutine Central_first_time_calculator

subroutine Central_second_calculator(delta,multiplier,down,center,up,term)
    
real(kind=8)::delta,down,center,up,multiplier,term
    term = multiplier*( up + down - 2.0 * center )/delta/delta
    
end subroutine Central_second_calculator

subroutine Upwind_2D_calculator(dx,dy,order,wind,x11,x12,x13,x21,x22,x23,x31,x32,x33,term)
    
real(kind=8)::dx,dy,wind,term,x11,x12,x13,x21,x22,x23,x31,x32,x33
integer :: order

if (order .eq. 1) then !dx first, dy second
    if ( wind .ge. 0.0 ) then 
        if ( x12 .ge. 0.0 ) then 
            term = wind*( x22 - x21 - (x12 - x11) )/dx/dy
        end if
        if ( x12 .lt. 0.0 ) then 
            term = wind*( x22 - x21 - (x13 - x12) )/dx/dy
        end if
    end if
    
    if ( wind .lt. 0.0 ) then 
        if ( x32 .ge. 0.0 ) then 
            term = wind*( x32 - x31 - (x23 - x22) )/dx/dy
        end if
        if ( x32 .lt. 0.0 ) then 
            term = wind*( x33 - x32 - (x23 - x22) )/dx/dy
        end if
    end if
else if (order .eq. 2) then !dy first, dx second
    if ( wind .ge. 0.0 ) then 
        if ( x21 .ge. 0.0 ) then 
            term = wind*( x22 - x12 - (x21 - x11) )/dx/dy
        end if
        if ( x21 .lt. 0.0 ) then 
            term = wind*( x22 - x12 - (x31 - x21) )/dx/dy
        end if
    end if
    
    if ( wind .lt. 0.0 ) then 
        if ( x23 .ge. 0.0 ) then 
            term = wind*( x23 - x13 - (x32 - x22) )/dx/dy
        end if
        if ( x23 .lt. 0.0 ) then 
            term = wind*( x33 - x23 - (x32 - x22) )/dx/dy
        end if
    end if
else 
    print *, "error in Upwind_2D_calculator"
end if 


end subroutine Upwind_2D_calculator

subroutine Central_2D_calculator(dx,dy,order,multiplier,x11,x12,x13,x21,x22,x23,x31,x32,x33,term)
real(kind=8)::dx,dy,multiplier,term,x11,x12,x13,x21,x22,x23,x31,x32,x33
integer :: order
if (order .eq. 1) then
    term = multiplier * ( (x33-x31)/2.0/dy - (x13-x11)/2.0/dy )/2.0/dx
else if (order .eq. 2) then
    term = multiplier * ( (x33-x13)/2.0/dx - (x31-x11)/2.0/dx )/2.0/dy
else 
print *, "wrong in central 2d"
end if


end subroutine Central_2D_calculator

subroutine Grid(nx,ny,domain_x_start,domain_x_end,domain_y_start,domain_y_end,min_cell_size,cell_size_inc,layer_count,&
dx,dy)

    integer i,j
    integer :: layer_count
    real(kind=8) :: domain_x_start,domain_x_end,domain_y_start,domain_y_end
    real(kind=8) :: min_cell_size,cell_size_inc
    
    real(kind=8),dimension(0:nx+1)  :: cell_size_x
    real(kind=8),dimension(0:ny+1)  :: cell_size_y
    
    real(kind=8),dimension(0:nx)  :: cell_size_x_u
    real(kind=8),dimension(0:ny)  :: cell_size_y_v
    
    real(kind=8),dimension(0:nx+1,0:ny+1)  :: dx
    real(kind=8),dimension(0:nx+1,0:ny+1)  :: dy
    
    real(kind=8) :: cell_size_x_total=0,cell_size_y_total=0
    real(kind=8),dimension(0:nx+2,0:ny+2)  :: grid_x
    real(kind=8),dimension(0:nx+2,0:ny+2)  :: grid_y
    
    real(kind=8),dimension(0:nx+1,0:ny+2)  :: grid_u_x
    real(kind=8),dimension(0:nx+1,0:ny+2)  :: grid_u_y
    
    real(kind=8),dimension(0:nx+2,0:ny+1)  :: grid_v_x
    real(kind=8),dimension(0:nx+2,0:ny+1)  :: grid_v_y
    
    do i=0,layer_count

        cell_size_x(i)=min_cell_size+i*cell_size_inc
    
        cell_size_x(nx+1-i)=cell_size_x(i)
    
    end do 
    
    do i=layer_count,nx+1-layer_count
    
        cell_size_x(i)=cell_size_x(layer_count)
    
    end do
    
    do i=0,nx+1
        
        cell_size_x_total=cell_size_x_total+cell_size_x(i)
        
    end do
    
    do i=0,nx+1
    
        cell_size_x(i)=cell_size_x(i)/cell_size_x_total*(domain_x_end-domain_x_start)
        dx(i,:)=cell_size_x(i)
    
    end do
    











    
    do j=0,layer_count

        cell_size_y(j)=min_cell_size+j*cell_size_inc
    
        cell_size_y(ny+1-j)=cell_size_y(j)
    
    end do 
    
    do j=layer_count,ny+1-layer_count
    
        cell_size_y(j)=cell_size_y(layer_count)
    
    end do
    
    do j=0,ny+1
        
        cell_size_y_total=cell_size_y_total+cell_size_y(j)
        
    end do
    
    do j=0,ny+1
    
        cell_size_y(j)=cell_size_y(j)/cell_size_y_total*(domain_y_end-domain_y_start)
        dy(:,j)=cell_size_y(j)
    end do
    
    
    
    grid_x(0,:)=domain_x_start
    do i=0,nx+1
        grid_x(i+1,:)=domain_x_start+cell_size_x(i)
    end do
    
    grid_y(:,0)=domain_y_start
    do j=0,ny+1
        grid_y(:,j+1)=domain_y_start+cell_size_y(j)
    end do
    
    
    


    
end subroutine 







subroutine initial_array(nx,ny,nt,timestep_id,&
                         x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                         u,u_alter,v,v_alter,&
                         p,p_alter,p_correct,&
                         f_former,f,f_alter,&
                         a_former,a,a_alter,&
                         b_former,b,b_alter,&
                         u_gs,v_gs,p_gs,&
                         uvbc,pbc,winds,&
                         u_w,u_e,u_s,u_n,&
                         v_w,v_e,v_s,v_n,&
                         p_w,p_e,p_s,p_n,&
                         ubc_e,ubc_w,&
                         vbc_n,vbc_s,&
                         pbc_e,pbc_w,pbc_n,pbc_s,&
                         u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                         u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                         p_bound_w,p_bound_e,&
                         p_bound_s,p_bound_n,&
                         f_bound_w,f_bound_e,&
                         f_bound_s,f_bound_n,&
                         a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                         a_bound_s,a_bound_n,b_bound_s,b_bound_n)
    
integer :: i,j
integer :: nx,ny,nt
integer :: timestep_id

real(kind=8) dx,dy

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter



real(kind=8)::u_w,u_e,u_s,u_n
real(kind=8)::v_w,v_e,v_s,v_n
real(kind=8)::p_w,p_e,p_s,p_n


real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n



character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds

!init domain
!x1=0 x101=100

!DOMAIN SHAPE IS
!1,ny+1     nx+1,ny+1
!
!
!
!1,1        1,nx+1

dx = (domain_x_end-domain_x_start)/nx
dy = (domain_y_end-domain_y_start)/ny


x(0,:)    = domain_x_start
x(1,:)    = domain_x_start + dx/2.0d0
x(nx,:)   = domain_x_end   - dx/2.0d0
x(nx+1,:) = domain_x_end
    
y(:,0)    = domain_y_start 
y(:,1)    = domain_y_start + dy/2.0d0
y(:,ny)   = domain_y_end   - dy/2.0d0
y(:,ny+1) = domain_y_end

do i=2,nx-1
    x(i,:)=domain_x_start + dx/2.0d0 + 1.0d0*(i-1)*dx
end do

do j=2,ny-1
    y(:,j)=domain_y_start + dy/2.0d0 + 1.0d0*(j-1)*dy
end do
    

!init uv/initial condition
    do i=0,nx+2
        do j=0,ny+2
            u(i,j)=1e-20
            u_alter(i,j)=1e-20
            u_gs(i,j,:)="c"
        end do
    end do
    
    do i=0,nx+2
        do j=0,ny+2
            v(i,j)=1e-20
            v_alter(i,j)=1e-20
            v_gs(i,j,:)="c"
        end do
    end do
!init pf/initial condition
    do i=0,nx+1
        do j=0,ny+1
            p(i,j)=0.0d0
            p_alter(i,j)=0.0d0
            p_gs(i,j,:)="c"
        end do
    end do
    


    do i=0,nx+1
        do j=0,ny+1
            f(i,j)=0.0d0
            f_former(i,j)=0.0d0
            f_alter(i,j)=0.0d0
            
            a(i,j)=0.0d0
            a_former(i,j)=0.0d0
            a_alter(i,j)=0.0d0

            b(i,j)=0.0d0
            b_former(i,j)=0.0d0
            b_alter(i,j)=0.0d0
        end do
    end do
    
print*, "Array initialized!"


end subroutine initial_array

subroutine initial_boundary(nx,ny,nt,nb,block_id,timestep_id,ghost_cell_size,&
                         x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                         u,u_alter,v,v_alter,&
                         p,p_alter,p_correct,&
                         f_former,f,f_alter,&
                         a_former,a,a_alter,&
                         b_former,b,b_alter,&
                         u_gs,v_gs,p_gs,&
                         uvbc,pbc,winds,&
                         u_w,u_e,u_s,u_n,&
                         v_w,v_e,v_s,v_n,&
                         p_w,p_e,p_s,p_n,&
                         ubc_e,ubc_w,&
                         vbc_n,vbc_s,&
                         pbc_e,pbc_w,pbc_n,pbc_s,&
                         u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                         u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                         p_bound_w,p_bound_e,&
                         p_bound_s,p_bound_n,&
                         f_bound_w,f_bound_e,&
                         f_bound_s,f_bound_n,&
                         a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                         a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                         u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                         u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                         p_ghost_cell_w,p_ghost_cell_e,&
                         p_ghost_cell_s,p_ghost_cell_n)
    
integer :: i,j,k,k_start
integer :: nx,ny,nt,nb,block_id,ghost_cell_size,timestep_id

real(kind=8) dx,dy

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter




real(kind=8),dimension(nb)::u_w,u_e,u_s,u_n
real(kind=8),dimension(nb)::v_w,v_e,v_s,v_n
real(kind=8),dimension(nb)::p_w,p_e,p_s,p_n


real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


    !The output ghost cell information for block_id
    real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
    

character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds

integer :: step_i_start = 1
integer :: step_i_end = 100
integer :: step_j_start = 1
integer :: step_j_end = 25
!init uv/initial condition
    do i=0,nx+2
        do j=0,ny+2
            u_gs(i,j,:)="c"
            v_gs(i,j,:)="c"
        end do
    end do
    
!init pf/initial condition
    do i=0,nx+1
        do j=0,ny+1
            p_gs(i,j,:)="c"
        end do
    end do
    
    
!init uv/boundary condition
do j=0,ny+2
        u_bound_e(j,:)=u_e(block_id)!0.0d0! inlet velocity
        u_bound_w(j,:)=u_w(block_id)!1.0d0! inlet velocity
        v_bound_e(j,:)=v_e(block_id)!0.0d0! wall velocity
        v_bound_w(j,:)=v_w(block_id)!0.0d0! wall velocity
end do

do i=0,nx+2
        u_bound_n(i,:)=u_n(block_id)!0.0d0! wall velocity
        u_bound_s(i,:)=u_s(block_id)!0.0d0! wall velocity
        v_bound_n(i,:)=v_n(block_id)!0.0d0! inlet velocity
        v_bound_s(i,:)=v_s(block_id)!0.0d0! inlet velocity
end do    

do j=0,ny+1
        p_bound_e(j,:)=p_e(block_id)!0.0d0
        p_bound_w(j,:)=p_w(block_id)!0.0d0
end do

do i=0,nx+1
        p_bound_n(i,:)=p_n(block_id)!0.0d0
        p_bound_s(i,:)=p_s(block_id)!0.0d0
end do    



do j=0,ny+1
    do k=0,nt
        f_bound_e(j,k)=0.0d0!-0.5! u east
        f_bound_w(j,k)=0.0d0!0.5! u west
        
        a_bound_e(j,k)=0.0d0
        a_bound_w(j,k)=0.0d0
        b_bound_e(j,k)=0.0d0
        b_bound_w(j,k)=0.0d0
    end do
end do
do i=0,nx+1
    do k=0,nt
        f_bound_n(i,k)=0.0d0! u north
        f_bound_s(i,k)=0.0d0! u south
        
        a_bound_n(i,k)=0.0d0
        a_bound_s(i,k)=0.0d0
        b_bound_n(i,k)=1.0E-8*sin(1.0*k)!0.0d0
        b_bound_s(i,k)=0.0d0
    end do
end do

!init /boundary condition
do j=0,ny+1
        f_former(0,j)=f_bound_w(j,timestep_id-1)
        f_former(nx+1,j)=f_bound_e(j,timestep_id-1)
        
        a_former(0,j)=a_bound_w(j,timestep_id-1)
        a_former(nx+1,j)=a_bound_e(j,timestep_id-1)
        
        b_former(0,j)=b_bound_w(j,timestep_id-1)
        b_former(nx+1,j)=b_bound_e(j,timestep_id-1)
        
        
        
        
        f(0,j)=f_bound_w(j,timestep_id)
        f(nx+1,j)=f_bound_e(j,timestep_id)
        
        a(0,j)=a_bound_w(j,timestep_id)
        a(nx+1,j)=a_bound_e(j,timestep_id)
        
        b(0,j)=b_bound_w(j,timestep_id)
        b(nx+1,j)=b_bound_e(j,timestep_id)
end do

do i=0,nx+1
        f_former(i,0)=f_bound_s(i,timestep_id-1)
        f_former(i,ny+1)=f_bound_n(i,timestep_id-1)
        
        a_former(i,0)=a_bound_s(i,timestep_id-1)
        a_former(i,ny+1)=a_bound_n(i,timestep_id-1)
        
        b_former(i,0)=b_bound_s(i,timestep_id-1)
        b_former(i,ny+1)=b_bound_n(i,timestep_id-1)
        
        
        
        f(i,0)=f_bound_s(i,timestep_id)
        f(i,ny+1)=f_bound_n(i,timestep_id)
        
        a(i,0)=a_bound_s(i,timestep_id)
        a(i,ny+1)=a_bound_n(i,timestep_id)
        
        b(i,0)=b_bound_s(i,timestep_id)
        b(i,ny+1)=b_bound_n(i,timestep_id)
end do    





!transfer boundary conditions to calculation
!1=VELOCITY INLET
!ubc_e=0
!ubc_w=1
!vbc_n=0
!vbc_s=0

!1=CONSTANT PRESSURE
!pbc_e=1
!pbc_w=0
!pbc_n=0
!pbc_s=0

if (uvbc(1,1,3).eq."I") then 
    ubc_w="I"
else if (uvbc(1,1,3).eq."O") then 
    ubc_w="O"
else 
    ubc_w="X"
end if
if (uvbc(1,2,3).eq."I") then
    ubc_e="I"
else if (uvbc(1,2,3).eq."O") then
    ubc_e="O"
else 
    ubc_e="X"
end if
if (uvbc(1,3,3).eq."I") then
    vbc_s="I"
else if (uvbc(1,3,3).eq."O") then
    vbc_s="O"
else 
    vbc_s="X"
end if
if (uvbc(1,4,3).eq."I") then
    vbc_n="I"
else if (uvbc(1,4,3).eq."O") then
    vbc_n="O"
else 
    vbc_n="X"
end if

if (pbc(1,1,3).eq."C") then
    pbc_w="C"
else 
    pbc_w="X"
end if
if (pbc(1,2,3).eq."C") then
    pbc_e="C"
else 
    pbc_e="X"
end if
if (pbc(1,3,3).eq."C") then
    pbc_s="C"
else 
    pbc_s="X"
end if
if (pbc(1,4,3).eq."C") then
    pbc_n="C"
else 
    pbc_n="X"
end if

!INLET LOCATED ON THE NEXT OF BOUNDARY
if (ubc_w.eq."I") then
    !u(1,:)=u_bound_w(:) ; u_alter(1,:)=u_bound_w(:)
    do j=step_j_end+1,ny!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP HEIGHT!!!!!!!!!!!!!!!!!
       u(1,j)=u_bound_w(j,timestep_id) ; u_alter(1,j)=u_bound_w(j,timestep_id) 
    end do
    
end if

if (ubc_e.eq."I") then
    u(nx+1,:)=u_bound_e(:,timestep_id) ; u_alter(nx+1,:)=u_bound_e(:,timestep_id)
end if
if (vbc_s.eq."I") then
    v(:,1)=v_bound_s(:,timestep_id)    ; v_alter(:,1)=v_bound_s(:,timestep_id)
end if
if (vbc_n.eq."I") then
    v(:,ny+1)=v_bound_n(:,timestep_id) ; v_alter(:,ny+1)=v_bound_n(:,timestep_id)
end if

if (pbc_e.eq."C") then
    p(nx,:)=p_bound_e(:,timestep_id)
end if
if (pbc_w.eq."C") then
    p(1,:)=p_bound_w(:,timestep_id)
end if
if (pbc_n.eq."C") then
    p(:,ny)=p_bound_n(:,timestep_id)
end if
if (pbc_s.eq."C") then
    p(:,1)=p_bound_s(:,timestep_id)
end if




!ASSIGN GRID STATUS grid may equal to w e s n or c
!Inner
do j=1,ny
    u_gs(2,  j, 1)="w" ; u_gs(nx, j, 2)="e" ; u_gs(2, j, 5)="B"  ; u_gs(nx, j, 5)="B"
end do
do i=2,nx
    u_gs(i, 1,  3)="s" ; u_gs(i, ny, 4)="n" ; u_gs(i, 1,  5)="B" ; u_gs(i, ny, 5)="B"
end do
do j=2,ny
    v_gs(1,  j, 1)="w" ; v_gs(nx, j, 2)="e" ; v_gs(1, j, 5)="B"  ; v_gs(nx, j, 5)="B"
end do
do i=1,nx
    v_gs(i, 2,  3)="s" ; v_gs(i, ny, 4)="n" ; v_gs(i, 2,  5)="B" ; v_gs(i, ny, 5)="B"
end do
do j=1,ny
    p_gs(1,  j, 1)="w" ; p_gs(nx, j, 2)="e" ; p_gs(1, j, 5)="B"  ; p_gs(nx, j, 5)="B"
end do
do i=1,nx
    p_gs(i, 1,  3)="s" ; p_gs(i, ny, 4)="n" ; p_gs(i, 1,  5)="B" ; p_gs(i, ny, 5)="B"
end do
print*, "Boundary initialized!"
!Ghost cell






do i=step_i_start,step_i_end
    do j=step_j_start,step_j_end
        !p_gs(i,j,1)="c"
        !p_gs(i,j,2)="c"
        !p_gs(i,j,3)="c"
        !p_gs(i,j,4)="c"
        !p_gs(i,j,5)="X"!As the inner block
    end do
end do

do i=step_i_end+1,step_i_end+1
    do j=step_j_start,step_j_end+1
        !p_gs(i, j, 1)="w" ; p_gs(i, j, 5)="B" ; p_gs(i, j, 6)="T"
    end do
end do
do i=step_i_start,step_i_end+1
    do j=step_j_end+1,step_j_end+1
        !p_gs(i, j, 3)="s" ; p_gs(i, j, 5)="B" ; p_gs(i, j, 6)="T"
    end do
end do

do i=1,nx
    do j=1,ny
        if ( p_gs(i,j,5) .eq. "X" ) then 
            !u_gs(i+1,j,5)="X"
            !v_gs(i,j+1,5)="X"
        else if ( p_gs(i, j, 5) .eq. "B" .and. p_gs(i, j, 6) .eq. "T") then
            !u_gs(i+1, j, 1)=p_gs(i, j, 1)
            !u_gs(i+1, j, 2)=p_gs(i, j, 2)
            !u_gs(i+1, j, 3)=p_gs(i, j, 3)
            !u_gs(i+1, j, 4)=p_gs(i, j, 4)
            !u_gs(i+1, j, 5)="B"
            !u_gs(i+1, j, 6)="T"
            
            !v_gs(i, j+1, 1)=p_gs(i, j, 1)
            !v_gs(i, j+1, 2)=p_gs(i, j, 2)
            !v_gs(i, j+1, 3)=p_gs(i, j, 3)
            !v_gs(i, j+1, 4)=p_gs(i, j, 4)
            !v_gs(i, j+1, 5)="B"
            !v_gs(i, j+1, 6)="T"
        end if
    end do
end do

!To avoid unidentified grid properties
do i=step_i_start,step_i_end
    do j=step_j_start,step_j_end
        !p_gs(i,j,1)="c"
        !p_gs(i,j,2)="c"
        !p_gs(i,j,3)="c"
        !p_gs(i,j,4)="c"
        !p_gs(i,j,5)="X"!As the inner block
         
        !u_gs(i+1, j, 1)=p_gs(i, j, 1)
        !u_gs(i+1, j, 2)=p_gs(i, j, 2)
        !u_gs(i+1, j, 3)=p_gs(i, j, 3)
        !u_gs(i+1, j, 4)=p_gs(i, j, 4)
        
        !v_gs(i, j+1, 1)=p_gs(i, j, 1)
        !v_gs(i, j+1, 2)=p_gs(i, j, 2)
        !v_gs(i, j+1, 3)=p_gs(i, j, 3)
        !v_gs(i, j+1, 4)=p_gs(i, j, 4)    
    end do
end do













if (uvbc(1,1,3).eq."L") then
    do j=1,ny
        u_gs(2,  j, 1)="w" ; u_gs(2,  j, 5)="L"
    end do
    do j=2,ny
        v_gs(1,  j, 1)="w" ; v_gs(1,  j, 5)="L"
    end do
else if (uvbc(1,1,3).eq."O") then
    do j=1,ny
        u_gs(2,  j, 1)="w" ; u_gs(2,  j, 5)="O"
    end do
    do j=2,ny
        v_gs(1,  j, 1)="w" ; v_gs(1,  j, 5)="O"
    end do
else
end if


if (uvbc(1,2,3).eq."L") then
    do j=1,ny
        u_gs(nx, j, 2)="e" ; u_gs(nx, j, 5)="L"
    end do
    do j=2,ny
        v_gs(nx, j, 2)="e" ; v_gs(nx, j, 5)="L"
    end do
else if (uvbc(1,2,3).eq."O") then
    do j=1,ny
        u_gs(nx, j, 2)="e" ; u_gs(nx, j, 5)="O"
    end do
    do j=2,ny
        v_gs(nx, j, 2)="e" ; v_gs(nx, j, 5)="O"
    end do
else
end if


if (uvbc(1,3,3).eq."L") then
    do i=2,nx
        u_gs(i, 1,  3)="s" ; u_gs(i, 1,  5)="L"
    end do
    do i=1,nx
        v_gs(i, 2,  3)="s" ; v_gs(i, 2,  5)="L"
    end do
else if (uvbc(1,3,3).eq."O") then
    do i=2,nx
        u_gs(i, 1,  3)="s" ; u_gs(i, 1,  5)="O"
    end do
    do i=1,nx
        v_gs(i, 2,  3)="s" ; v_gs(i, 2,  5)="O"
    end do
else 
end if


if (uvbc(1,4,3).eq."L") then
    do i=2,nx
        u_gs(i, ny, 4)="n" ; u_gs(i, ny, 5)="L"
    end do
    do i=1,nx
        v_gs(i, ny, 4)="n" ; v_gs(i, ny, 5)="L"
    end do
else if (uvbc(1,4,3).eq."O") then
    do i=2,nx
        u_gs(i, ny, 4)="n" ; u_gs(i, ny, 5)="O"
    end do
    do i=1,nx
        v_gs(i, ny, 4)="n" ; v_gs(i, ny, 5)="O"
    end do
else 
end if


            



if (pbc(1,1,3).eq."L") then
    do j=1,ny
        p_gs(1,  j, 1)="w" ; p_gs(1,  j, 5)="L"
    end do
end if

if (pbc(1,2,3).eq."L") then
    do j=1,ny
        p_gs(nx, j, 2)="e" ; p_gs(nx, j, 5)="L"
    end do
end if
if (pbc(1,3,3).eq."L") then
    do i=1,nx
        p_gs(i, 1,  3)="s" ; p_gs(i, 1,  5)="L"
    end do
end if
if (pbc(1,4,3).eq."L") then
    do i=1,nx
        p_gs(i, ny, 4)="n" ; p_gs(i, ny, 5)="L"
    end do
end if


end subroutine initial_boundary

subroutine ghost_cell_finder(block_id,nb,nx,ny,gx,gy,ghost_cell_size,&
                             block_grid_count1,block_grid_count2,&
                             block_grid_id1,block_grid_id2,&
                             linkbc,&
                             x_block,y_block,u_block,v_block,p_block,f_block,a_block,b_block,&
                             u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                             u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                             p_ghost_cell_w,p_ghost_cell_e,&
                             p_ghost_cell_s,p_ghost_cell_n)
    integer :: i,j
    integer :: nb,nx,ny,ghost_cell_size
    integer,dimension(1:nb)::gx,gy
    integer ::block_id
    integer ::block_grid_count1,block_grid_count2
    integer,dimension(1:nb)::block_grid_id1,block_grid_id2

    !For ghost cell
    integer :: ghost_cell_w_block_id,ghost_cell_e_block_id,ghost_cell_s_block_id,ghost_cell_n_block_id
    character(len=1),dimension(4)::ghost_cell_boundary="X"

    !Local array to storage the data that grabbed out
    !real(kind=8),dimension(0:nx+1, 0:ny+1)::x_local,y_local
    !real(kind=8),dimension(0:nx+2, 0:ny+2)::u_local,v_local
    !real(kind=8),dimension(0:nx+1, 0:ny+1)::p_local,f_local,a_local,b_local
    
    !All data list
    real(kind=8),dimension(block_grid_count1)::x_block,y_block,p_block,f_block,a_block,b_block
    real(kind=8),dimension(block_grid_count2)::u_block,v_block


    !Local array to storage boundary topology information
    integer,dimension(4)::boundary_id_this_block
    character(len=2),dimension(4)::boundary_link_this_block
    integer,dimension(4*nb)::boundary_id
    character(len=2),dimension(4*nb)::boundary_link
    !Boundary link passed to this subroutine
    character(len=2),dimension(nb,4) ::linkbc
    
    !The output ghost cell information for block_id
    real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
!ONLY REFLESH GHOST CELL ONCE A TIME
    u_ghost_cell_w=0.0d0 ; v_ghost_cell_w=0.0d0 ; p_ghost_cell_w=0.0d0 ; 
    u_ghost_cell_e=0.0d0 ; v_ghost_cell_e=0.0d0 ; p_ghost_cell_e=0.0d0
    u_ghost_cell_s=0.0d0 ; v_ghost_cell_s=0.0d0 ; p_ghost_cell_s=0.0d0
    u_ghost_cell_n=0.0d0 ; v_ghost_cell_n=0.0d0 ; p_ghost_cell_n=0.0d0
            
    !ID in this block
    boundary_id_this_block(1)=1+4*(block_id-1)
    boundary_id_this_block(2)=2+4*(block_id-1)
    boundary_id_this_block(3)=3+4*(block_id-1)
    boundary_id_this_block(4)=4+4*(block_id-1)
    !Unique ID links to topology in this block
    boundary_link_this_block(1)=linkbc(block_id,1)!W
    boundary_link_this_block(2)=linkbc(block_id,2)!E
    boundary_link_this_block(3)=linkbc(block_id,3)!S
    boundary_link_this_block(4)=linkbc(block_id,4)!N

    
    do k=1,nb
    !ID list for all blocks
        boundary_id(1+4*(k-1))=1+4*(k-1)
        boundary_id(2+4*(k-1))=2+4*(k-1)
        boundary_id(3+4*(k-1))=3+4*(k-1)
        boundary_id(4+4*(k-1))=4+4*(k-1)
        boundary_link(1+4*(k-1))=linkbc(k,1)!W
        boundary_link(2+4*(k-1))=linkbc(k,2)!E
        boundary_link(3+4*(k-1))=linkbc(k,3)!S
        boundary_link(4+4*(k-1))=linkbc(k,4)!N
    end do
    !Compare ID in this block with the list
    !W only links to E
    !S only links to N
    do k=1,nb
        if ( ( boundary_link_this_block(1) .eq. boundary_link(2+4*(k-1)) ) .and. &
             ( boundary_id_this_block(1)   .ne. boundary_id(2+4*(k-1)) )   .and. &
               boundary_link_this_block(1) .ne. "xx" ) then !W as ghost cell connects to E
            ghost_cell_e_block_id=k
            ghost_cell_boundary(2)="e"
            print *, "Find ghost cell at block",ghost_cell_e_block_id,"boundary  ",ghost_cell_boundary(2),"  Now fill ghost cell!"
            call ghost_cell_filler(ghost_cell_e_block_id,block_id,nb,nx,ny,gx,gy,ghost_cell_size,ghost_cell_boundary,&
                                   block_grid_count1,block_grid_count2,&
                                   block_grid_id1,block_grid_id2,&
                                   x_block,y_block,u_block,v_block,p_block,f_block,a_block,b_block,&
                                   u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                                   u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                                   p_ghost_cell_w,p_ghost_cell_e,&
                                   p_ghost_cell_s,p_ghost_cell_n)
        end if
               
        if ( ( boundary_link_this_block(2) .eq. boundary_link(1+4*(k-1)) ) .and. &
             ( boundary_id_this_block(2)   .ne. boundary_id(1+4*(k-1)) )   .and. &
               boundary_link_this_block(2) .ne. "xx" ) then !E as ghost cell connects to W
            ghost_cell_w_block_id=k
            ghost_cell_boundary(1)="w"
            print *, "Find ghost cell at block",ghost_cell_w_block_id,"boundary  ",ghost_cell_boundary(1),"  Now fill ghost cell!"
            call ghost_cell_filler(ghost_cell_w_block_id,block_id,nb,nx,ny,gx,gy,ghost_cell_size,ghost_cell_boundary,&
                                   block_grid_count1,block_grid_count2,&
                                   block_grid_id1,block_grid_id2,&
                                   x_block,y_block,u_block,v_block,p_block,f_block,a_block,b_block,&
                                   u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                                   u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                                   p_ghost_cell_w,p_ghost_cell_e,&
                                   p_ghost_cell_s,p_ghost_cell_n)
        end if 

        if ( ( boundary_link_this_block(3) .eq. boundary_link(4+4*(k-1)) ) .and. &
             ( boundary_id_this_block(3)   .ne. boundary_id(4+4*(k-1)) )   .and. &
               boundary_link_this_block(3) .ne. "xx" ) then !W as ghost cell connects to E
            ghost_cell_n_block_id=k
            ghost_cell_boundary(4)="n"
            print *, "Find ghost cell at block",ghost_cell_n_block_id,"boundary  ",ghost_cell_boundary(4),"  Now fill ghost cell!"
            call ghost_cell_filler(ghost_cell_n_block_id,block_id,nb,nx,ny,gx,gy,ghost_cell_size,ghost_cell_boundary,&
                                   block_grid_count1,block_grid_count2,&
                                   block_grid_id1,block_grid_id2,&
                                   x_block,y_block,u_block,v_block,p_block,f_block,a_block,b_block,&
                                   u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                                   u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                                   p_ghost_cell_w,p_ghost_cell_e,&
                                   p_ghost_cell_s,p_ghost_cell_n)
        end if
        
        if ( ( boundary_link_this_block(4) .eq. boundary_link(3+4*(k-1)) ) .and. &
             ( boundary_id_this_block(4)   .ne. boundary_id(3+4*(k-1)) )   .and. &
               boundary_link_this_block(4) .ne. "xx" ) then !E as ghost cell connects to W
            ghost_cell_s_block_id=k
            ghost_cell_boundary(3)="s"
            print *, "Find ghost cell at block",ghost_cell_s_block_id,"boundary  ",ghost_cell_boundary(3),"  Now fill ghost cell!"
            call ghost_cell_filler(ghost_cell_s_block_id,block_id,nb,nx,ny,gx,gy,ghost_cell_size,ghost_cell_boundary,&
                                   block_grid_count1,block_grid_count2,&
                                   block_grid_id1,block_grid_id2,&
                                   x_block,y_block,u_block,v_block,p_block,f_block,a_block,b_block,&
                                   u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                                   u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                                   p_ghost_cell_w,p_ghost_cell_e,&
                                   p_ghost_cell_s,p_ghost_cell_n)
        end if
 
    end do
    
    ghost_cell_boundary(:)="X"

end subroutine ghost_cell_finder

subroutine ghost_cell_filler(ghost_cell_block_id,block_id,nb,nx,ny,gx,gy,ghost_cell_size,ghost_cell_boundary,&
                             block_grid_count1,block_grid_count2,&
                             block_grid_id1,block_grid_id2,&
                             x_block,y_block,u_block,v_block,p_block,f_block,a_block,b_block,&
                             u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                             u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                             p_ghost_cell_w,p_ghost_cell_e,&
                             p_ghost_cell_s,p_ghost_cell_n)
    
    !LOCAL DATA

    
    integer :: i,j,k,k_start
    integer :: nb,nx,ny,ngx,ngy,ghost_cell_size
    integer,dimension(1:nb)::gx,gy
    integer ::ghost_cell_block_id,block_id
    character(len=1),dimension(4)::ghost_cell_boundary
    integer ::block_grid_count1,block_grid_count2
    integer,dimension(1:nb)::block_grid_id1,block_grid_id2
    !All data list
    real(kind=8),dimension(block_grid_count1)::x_block,y_block,p_block,f_block,a_block,b_block
    real(kind=8),dimension(block_grid_count2)::u_block,v_block

    !The output ghost cell information for block_id
    real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
    
    real(kind=8),allocatable,dimension(:,:)::x_local,y_local
    real(kind=8),allocatable,dimension(:,:)::u_local,v_local
    real(kind=8),allocatable,dimension(:,:)::p_local,f_local,a_local,b_local
    
    ngx=gx(ghost_cell_block_id)
    ngy=gy(ghost_cell_block_id)
    print *,"this block id",block_id,"size", nx,ny,"ghost block id",ghost_cell_block_id,"size", ngx,ngy
    print *, block_grid_id1(:),block_grid_id2(:)
    
    allocate(x_local(0:ngx+1, 0:ngy+1),y_local(0:ngx+1, 0:ngy+1))
    allocate(u_local(0:ngx+2, 0:ngy+2),v_local(0:ngx+2, 0:ngy+2))
    allocate(p_local(0:ngx+1, 0:ngy+1),f_local(0:ngx+1, 0:ngy+1),a_local(0:ngx+1, 0:ngy+1),b_local(0:ngx+1, 0:ngy+1))

    !u_local = 0.0d0 ; v_local = 0.0d0 
    !print *, ghost_cell_block_id,ghost_cell_boundary(:)
    !Grab data out from all 

    do i=1,ngx+2
        do j=1,ngy+2
            x_local(i-1,j-1)=x_block( j+(ngy+2)*(i-1)+block_grid_id1(ghost_cell_block_id) )!i+nx*(j-1)
            y_local(i-1,j-1)=y_block( j+(ngy+2)*(i-1)+block_grid_id1(ghost_cell_block_id) )!i+nx*(j-1)
            p_local(i-1,j-1)=p_block( j+(ngy+2)*(i-1)+block_grid_id1(ghost_cell_block_id) )!i+nx*(j-1)
            f_local(i-1,j-1)=f_block( j+(ngy+2)*(i-1)+block_grid_id1(ghost_cell_block_id) )!i+nx*(j-1)
            a_local(i-1,j-1)=a_block( j+(ngy+2)*(i-1)+block_grid_id1(ghost_cell_block_id) )!i+nx*(j-1)
            b_local(i-1,j-1)=b_block( j+(ngy+2)*(i-1)+block_grid_id1(ghost_cell_block_id) )!i+nx*(j-1)
        end do
    end do
    
        
    do i=1,ngx+3
        do j=1,ngy+3
            u_local(i-1,j-1)=u_block( j+(ngy+3)*(i-1)+block_grid_id2(ghost_cell_block_id) )!i+nx*(j-1)
            v_local(i-1,j-1)=v_block( j+(ngy+3)*(i-1)+block_grid_id2(ghost_cell_block_id) )!i+nx*(j-1)
            !print *, i,j,j+(ngy+3)*(i-1)+block_grid_id2(ghost_cell_block_id),u_local(i-1,j-1)
            !print *,i-1,j-1,u_local(i-1,j-1),"GS U"
        end do
    end do
    
    
    !CHECK compaction of the ghost cell
    if (ghost_cell_boundary(1) .eq. "w") then
        if (ny .ne. ngy) then
            print *, "critical mistake, please check ghost cell finder"
        end if 
    end if
    if (ghost_cell_boundary(2) .eq. "e") then
        if (ny .ne. ngy) then
            print *, "critical mistake, please check ghost cell finder"
        end if 
    end if
    if (ghost_cell_boundary(3) .eq. "s") then
        if (nx .ne. ngx) then
            print *, "critical mistake, please check ghost cell finder"
        end if 
    end if
    if (ghost_cell_boundary(4) .eq. "n") then
        if (nx .ne. ngx) then
            print *, "critical mistake, please check ghost cell finder"
        end if 
    end if

    
    !init uv/boundary condition
    do j=0+1,ngy+2+1
        if (ghost_cell_boundary(1) .eq. "w") then
            !u_ghost_cell_e(1,j-1)=u_local(2,j-1)
            !u_ghost_cell_e(2,j-1)=u_local(3,j-1)
            !u_ghost_cell_e(3,j-1)=u_local(4,j-1)
            !u_ghost_cell_e(4,j-1)=u_local(5,j-1)
            do k = 1, ghost_cell_size
                k_start = k + 1
                u_ghost_cell_e(k,j-1) = u_local(k_start,j-1)
            end do
            
            !v_ghost_cell_e(1,j-1)=v_local(1,j-1)
            !v_ghost_cell_e(2,j-1)=v_local(2,j-1)
            !v_ghost_cell_e(3,j-1)=v_local(3,j-1)
            !v_ghost_cell_e(4,j-1)=v_local(4,j-1)
            do k = 1, ghost_cell_size
                k_start = k 
                v_ghost_cell_e(k,j-1)=v_local(k_start,j-1)
            end do
        end if
        if (ghost_cell_boundary(2) .eq. "e") then
            !u_ghost_cell_w(7,j-1)=u_local(ngx-3,j-1)
            !u_ghost_cell_w(8,j-1)=u_local(ngx-2,j-1)
            !u_ghost_cell_w(9,j-1)=u_local(ngx-1,j-1)
            !u_ghost_cell_w(10,j-1)=u_local(ngx,j-1)
            
            do k = 1, ghost_cell_size
                k_start = ngx - ghost_cell_size + k
                u_ghost_cell_w(k,j-1) = u_local(k_start,j-1)
            end do

            !v_ghost_cell_w(7,j-1)=v_local(ngx-3,j-1)
            !v_ghost_cell_w(8,j-1)=v_local(ngx-2,j-1)
            !v_ghost_cell_w(9,j-1)=v_local(ngx-1,j-1)
            !v_ghost_cell_w(10,j-1)=v_local(ngx,j-1)
            
            do k = 1, ghost_cell_size
                k_start = ngx - ghost_cell_size + k
                v_ghost_cell_w(k,j-1) = v_local(k_start,j-1)
            end do
        end if
    end do
    do i=0+1,ngx+2+1
        if (ghost_cell_boundary(3) .eq. "s") then
            !u_ghost_cell_n(1,i-1)=u_local(i-1,1)
            !u_ghost_cell_n(2,i-1)=u_local(i-1,2)
            !u_ghost_cell_n(3,i-1)=u_local(i-1,3)
            !u_ghost_cell_n(4,i-1)=u_local(i-1,4)
            
            do k = 1, ghost_cell_size
                k_start = k 
                u_ghost_cell_n(k,i-1)=u_local(i-1,k_start)
            end do
        
            !v_ghost_cell_n(1,i-1)=v_local(i-1,2)
            !v_ghost_cell_n(2,i-1)=v_local(i-1,3)
            !v_ghost_cell_n(3,i-1)=v_local(i-1,4)
            !v_ghost_cell_n(4,i-1)=v_local(i-1,5)
            do k = 1, ghost_cell_size
                k_start = k + 1
                v_ghost_cell_n(k,i-1)=v_local(i-1,k_start)
            end do
        end if
        if (ghost_cell_boundary(4) .eq. "n") then
            !u_ghost_cell_s(7,i-1)=u_local(i-1,ngy-3)
            !u_ghost_cell_s(8,i-1)=u_local(i-1,ngy-2)
            !u_ghost_cell_s(9,i-1)=u_local(i-1,ngy-1)
            !u_ghost_cell_s(10,i-1)=u_local(i-1,ngy)
            
            do k = 1, ghost_cell_size
                k_start = ngy - ghost_cell_size + k
                u_ghost_cell_s(k,i-1)=u_local(i-1,k_start)
            end do


            !v_ghost_cell_s(7,i-1)=v_local(i-1,ngy-3)
            !v_ghost_cell_s(8,i-1)=v_local(i-1,ngy-2)
            !v_ghost_cell_s(9,i-1)=v_local(i-1,ngy-1)
            !v_ghost_cell_s(10,i-1)=v_local(i-1,ngy)

            do k = 1, ghost_cell_size
                k_start = ngy - ghost_cell_size + k
                v_ghost_cell_s(k,i-1)=v_local(i-1,k_start)
            end do
        end if
    end do    


    do j=0,ngy+1
        if (ghost_cell_boundary(1) .eq. "w") then
            !p_ghost_cell_e(1,j)=p_local(0,j)!p_w!0.0d0
            !p_ghost_cell_e(2,j)=p_local(1,j)!p_w!0.0d0
            !p_ghost_cell_e(3,j)=p_local(2,j)!p_w!0.0d0
            
            do k = 1, ghost_cell_size
                k_start = k - 1 
                p_ghost_cell_e(k,j) = p_local(k_start,j)
            end do
        end if
        if (ghost_cell_boundary(2) .eq. "e") then
            !p_ghost_cell_w(1,j)=p_local(ngx-1,j)!p_e!0.0d0
            !p_ghost_cell_w(2,j)=p_local(ngx  ,j)!p_e!0.0d0
            !p_ghost_cell_w(3,j)=p_local(ngx+1,j)!p_e!0.0d0
            
            do k = 1, ghost_cell_size
                k_start = ngx + 1 - ghost_cell_size + k
                p_ghost_cell_w(k,j) = p_local(k_start,j)
            end do
        end if
    end do
    do i=0,ngx+1
        if (ghost_cell_boundary(3) .eq. "s") then
            !p_ghost_cell_n(1,i)=p_local(i,0)!p_s!0.0d0
            !p_ghost_cell_n(2,i)=p_local(i,1)
            !p_ghost_cell_n(3,i)=p_local(i,2)
            
            do k = 1, ghost_cell_size
                k_start = k - 1 
                p_ghost_cell_n(k,i) = p_local(k_start,i)
            end do
        end if
        if (ghost_cell_boundary(4) .eq. "n") then
            !p_ghost_cell_s(1,i)=p_local(i,ngy-1)
            !p_ghost_cell_s(2,i)=p_local(i,ngy  )!p_n!0.0d0
            !p_ghost_cell_s(3,i)=p_local(i,ngy+1)!p_n!0.0d0
            
            do k = 1, ghost_cell_size
                k_start = ngy + 1 - ghost_cell_size + k
                p_ghost_cell_s(k,i) = p_local(k_start,i)
            end do
        end if
    end do
    
    
    deallocate(x_local,y_local)
    deallocate(u_local,v_local)
    deallocate(p_local,f_local,a_local,b_local)

end subroutine ghost_cell_filler


subroutine coefficient_filler_u(nx,ny,nn,A_coe,B_bound,AA_coe,BB_bound)

integer nx,ny,nn
real(kind=8),dimension(2:nx,1:ny,5):: AA_coe
real(kind=8),dimension(2:nx,1:ny):: BB_bound

real(kind=8),dimension(1:nx-1,1:ny,5):: AA_coe_temp
real(kind=8),dimension(1:nx-1,1:ny):: BB_bound_temp

real(kind=8),dimension(nn,nn):: A_coe
real(kind=8),dimension(nn):: B_bound

real(kind=8),dimension(nn):: temp

integer counter 
A_coe =0.0d0

do i=2,nx
    do j=1,ny
        AA_coe_temp(i-1,j,1)=AA_coe(i,j,1)
        AA_coe_temp(i-1,j,2)=AA_coe(i,j,2)
        AA_coe_temp(i-1,j,3)=AA_coe(i,j,3)
        AA_coe_temp(i-1,j,4)=AA_coe(i,j,4)
        AA_coe_temp(i-1,j,5)=AA_coe(i,j,5)
        BB_bound_temp(i-1,j)=BB_bound(i,j)
    end do
end do


counter = 0

    do j=1,ny
        do i=1,nx-1
            counter = counter + 1
            temp(:)=0.0d0
            temp(i+(j-1)*(nx-1))=AA_coe_temp(i,j,5)
            
            if (i+(j-1)*(nx-1)-1 .ge. 1) then
            temp(i+(j-1)*(nx-1)-1)=AA_coe_temp(i,j,1)
            end if
                        
            if (i+(j-1)*(nx-1)+1 .le. nn) then
            temp(i+(j-1)*(nx-1)+1)=AA_coe_temp(i,j,2)
            end if
            
            if (i+(j-1)*(nx-1)-(nx-1) .ge. 1) then
            temp(i+(j-1)*(nx-1)-(nx-1))=AA_coe_temp(i,j,3)
            end if

            if (i+(j-1)*(nx-1)+(nx-1) .le. nn) then
            temp(i+(j-1)*(nx-1)+nx-1)=AA_coe_temp(i,j,4)
            end if
            
            A_coe(counter,:)=temp(:)
            B_bound(i+(j-1)*(nx-1))=BB_bound_temp(i,j)
        end do
    end do
    
    do i=1,nn
            !print *,A_coe(i,:), B_bound(i)
    end do

end subroutine coefficient_filler_u

subroutine u_solver(nx,ny,nt,timestep_id,density,viscosity,scale_factor,ghost_cell_size,&
                    x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                    u_gs,v_gs,p_gs,&
                    u,u_alter,v,v_alter,&
                    p,p_alter,p_correct,&
                    f_former,f,f_alter,&
                    a_former,a,a_alter,&
                    b_former,b,b_alter,&
                    ap_u,ap_v,ap_u1,ap_v1,&
                    uvbc,pbc,winds,&
                    ubc_e,ubc_w,&
                    vbc_n,vbc_s,&
                    pbc_e,pbc_w,pbc_n,pbc_s,&
                    u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                    u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                    p_bound_w,p_bound_e,&
                    p_bound_s,p_bound_n,&
                    f_bound_w,f_bound_e,&
                    f_bound_s,f_bound_n,&
                    a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                    a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                    u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                    u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                    p_ghost_cell_w,p_ghost_cell_e,&
                    p_ghost_cell_s,p_ghost_cell_n)
    

integer :: i,j,k,counter,ghost_cell_size
integer :: nx,ny,nt,timestep_id
integer :: bimaxit
real(kind=8)::bires

real(kind=8)::dx,dy
real(kind=8)::density
real(kind=8)::viscosity
real(kind=8)::scale_factor

real(kind=8)::temp_new1,temp_new2,temp_new3,temp_new4
real(kind=8)::temp_square,temp_square_derive

real(kind=8)::a_n,a_s,a_w,a_e
real(kind=8)::d_n,d_s,d_w,d_e
real(kind=8)::source_term,bound
real(kind=8)::at_n,at_s,at_w,at_e,at_c,at_f
real(kind=8)::bt_n,bt_s,bt_w,bt_e,bt_c

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end


character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter

real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v

real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s

character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds
    !The output ghost cell information for block_id
    real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
    
!USED FOR SOLVER
integer :: un,vn,pn
integer :: unnzero,vnnzero,pnnzero
real(kind=8),allocatable,dimension(:):: umatv,umrhs,vmatv,vmrhs,pmatv,pmrhs
integer,allocatable,dimension(:):: umati,umatj,vmati,vmatj,pmati,pmatj
real(kind=8),allocatable,dimension(:)::solution_u,solution_v,solution_p

!real(kind=8),dimension(2:nx,1:ny,5)::AA_coe
!real(kind=8),dimension(2:nx,1:ny)::BB_bound

!real(kind=8),dimension((nx-1)*ny,(nx-1)*ny)::A_coe
!real(kind=8),dimension((nx-1)*ny)::B_bound

    real(kind=8),allocatable,dimension(:,:,:)::AA_coe
    real(kind=8),allocatable,dimension(:,:)::BB_bound
    real(kind=8),allocatable,dimension(:,:)::A_coe
    real(kind=8),allocatable,dimension(:)::B_bound

un=(nx-1)*ny ; vn=nx*(ny-1) ; pn=nx*ny

unnzero=5*(nx-3)*(ny-2)+4*(nx-3)*2+4*(ny-2)*2+3*4
vnnzero=5*(nx-2)*(ny-3)+4*(nx-2)*2+4*(ny-3)*2+3*4
pnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4

bires=1.0e-20
bimaxit=300

allocate (AA_coe(2:nx,1:ny,5),BB_bound(2:nx,1:ny),A_coe((nx-1)*ny,(nx-1)*ny),B_bound((nx-1)*ny))

allocate (umatv(1:unnzero),umati(1:unnzero),umatj(1:unnzero),umrhs(1:(nx-1)*(ny)))
allocate (vmatv(1:vnnzero),vmati(1:vnnzero),vmatj(1:vnnzero),vmrhs(1:nx*(ny-1)))
allocate (pmatv(1:pnnzero),pmati(1:pnnzero),pmatj(1:pnnzero),pmrhs(1:nx*ny))

allocate(solution_u( (nx-1)*(ny) ) ) ; allocate(solution_v( nx*(ny-1) ) ) ; allocate(solution_p( nx*ny ) )
    umatv=0 ; umati=0 ; umatj=0 
    vmatv=0 ; vmati=0 ; vmatj=0 
    pmatv=0 ; pmati=0 ; pmatj=0 
!USED FOR SOLVER
    solution_u=0.0 ; solution_v=0.0 ; solution_p=0.0
    AA_coe = 0.0 ; BB_bound=0.0 ; A_coe = 0.0 ; B_bound = 0.0
!START CALCULATION
!for main body
    counter = 0
    do i=2,nx
        do j=1,ny
            dx=x(i,j)-x(i-1,j) ; dy=y(i,j)-y(i,j-1)
            dx=(domain_x_end-domain_x_start)/nx ; dy=(domain_y_end-domain_y_start)/ny
            !Compute coefficients
            a_e=0.0d0 ; a_w=0.0d0 ; a_n=0.0d0 ; a_s=0.0d0
            d_e=0.0d0 ; d_w=0.0d0 ; d_n=0.0d0 ; d_s=0.0d0
            at_e=0.0d0 ; at_w=0.0d0 ; at_n=0.0d0 ;at_s=0.0d0 ; at_c=0.0d0
            bt_e=0.0d0 ; bt_w=0.0d0 ; bt_n=0.0d0 ;bt_s=0.0d0 ; bt_c=0.0d0

            a_e = density*( u(i,j)     + u(i+1,j) )*dy/2.0d0
            a_w = density*( u(i-1,j)   + u(i,j)   )*dy/2.0d0
            a_n = density*( v(i-1,j+1) + v(i,j+1) )*dx/2.0d0
            a_s = density*( v(i-1,j)   + v(i,j)   )*dx/2.0d0
            
            d_e=viscosity*dy/dx ; d_w=viscosity*dy/dx ; d_n=viscosity*dx/dy ; d_s=viscosity*dx/dy

            !hybrid format
            at_e=max(-a_e,(d_e-a_e/2.0d0),0.0d0)
            at_w=max(+a_w,(d_w+a_w/2.0d0),0.0d0)
            at_n=max(-a_n,(d_n-a_n/2.0d0),0.0d0)
            at_s=max(+a_s,(d_s+a_s/2.0d0),0.0d0)
            at_c=at_e+at_w+at_n+at_s+(a_e-a_w+a_n-a_s)
            at_f=1.0d0
            temp_square=viscosity*( a(i,j)**2+b(i,j)**2 )*( u(i,j)+u(i+1,j) )/2.0d0
            
            temp_square_derive=viscosity*scale_factor*( ( b(i+1,j)-b(i-1,j) )/2.0d0/dx - &
                                                        ( a(i,j+1)-a(i,j-1) )/2.0d0/dy )**2*( u(i,j)+u(i+1,j) )/2.0d0

            temp_new1=viscosity*2*a(i,j)*(v(i,j)+v(i+1,j)+v(i,j+1)+v(i+1,j+1)-v(i-1,j+1)-v(i,j+1)-v(i-1,j)-v(i,j))/4.0d0/dx
            temp_new2=viscosity*2*b(i,j)*(v(i,j+1)-v(i,j))/dy
            temp_new3=viscosity*(v(i,j)+v(i,j+1))/2.0d0*(a(i+1,j)-a(i-1,j))/2.0d0/dx
            temp_new4=viscosity*(v(i,j)+v(i,j+1))/2.0d0*(b(i,j+1)-b(i,j-1))/2.0d0/dy
        
            !Compute source term
            source_term=0.0d0
            source_term = - temp_square - temp_square_derive - temp_new1 - temp_new2 - temp_new3 - temp_new4
            bound=0.0d0
            

            if (u_gs(i,j,1) .ne. "w") then !W
                counter=counter+1
                umatv(counter)=-at_w
                umati(counter)=(j-1)*(nx-1)+(i-1-1)
                umatj(counter)=(j-1)*(nx-1)+(i-1)
                AA_coe(i,j,1)=-at_w
                !print *, i,j,counter,umatv(counter)
            else if (u_gs(i,j,1) .eq. "w" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "c") then
                bound = bound + at_w * u(i-1,j)
            else if (u_gs(i,j,1) .eq. "w" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "T") then
                at_w = 0.0 
                AA_coe(i,j,1) = -at_w
                bound = bound + at_w * u(i-1,j)
            else if (u_gs(i,j,1) .eq. "w" .and. u_gs(i,j,5) .eq. "O") then
                bound = bound
                at_w = 0.0
                
            else if (u_gs(i,j,1) .eq. "w" .and. u_gs(i,j,5) .eq. "L"  .and. u_gs(i,j,6) .eq. "U") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0
                !bound = bound + at_w * u(i-1,j)
            else if (u_gs(i,j,1) .eq. "w" .and. u_gs(i,j,5) .eq. "L"  .and. u_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = bound + at_w * u(i-1,j)
            else
                print *, "unidentified grid property at",i,j
            end if
            

            
            if (u_gs(i,j,2) .ne. "e") then !E
                counter=counter+1
                umatv(counter)=-at_e
                umati(counter)=(j-1)*(nx-1)+(i-1+1)
                umatj(counter)=(j-1)*(nx-1)+(i-1)
                AA_coe(i,j,2)=-at_e
                !print *, i,j,counter,umatv(counter)
            else if (u_gs(i,j,2) .eq. "e" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "c") then
                bound = bound + at_e * u(i+1,j)
            else if (u_gs(i,j,2) .eq. "e" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,5) .eq. "B") then
                at_e = 0.0
                AA_coe(i,j,2) = -at_e
                bound = bound + at_e * u(i+1,j)
            else if (u_gs(i,j,2) .eq. "e" .and. u_gs(i,j,5) .eq. "O") then
                bound = bound 
                at_e = 0.0
                                
            else if (u_gs(i,j,2) .eq. "e" .and. u_gs(i,j,5) .eq. "L" .and. u_gs(i,j,6) .eq. "U") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0
                !bound = bound + at_e * u(i+1,j)
            else if (u_gs(i,j,2) .eq. "e" .and. u_gs(i,j,5) .eq. "L" .and. u_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = bound + at_e * u(i+1,j)
            else 
                print *, "unidentified grid property at",i,j
            end if

            
            if (u_gs(i,j,3) .ne. "s") then !S
                counter=counter+1
                umatv(counter)=-at_s
                umati(counter)=(j-1-1)*(nx-1)+(i-1)
                umatj(counter)=(j-1)*(nx-1)+(i-1)
                AA_coe(i,j,3)=-at_s
                !print *, i,j,counter,umatv(counter)
            else if (u_gs(i,j,3) .eq. "s" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "c") then
                at_c = at_c + viscosity * dx /(0.5*dy)
                bound = viscosity * u_bound_s(i,timestep_id) * dx / (0.5*dy)
            else if (u_gs(i,j,3) .eq. "s" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "T") then
                at_s = 0.0
                AA_coe(i,j,3) = -at_s
                at_c = at_c + viscosity * dx /(0.5*dy)
                bound = viscosity * u_bound_s(i,timestep_id) * dx / (0.5*dy)
            else if (u_gs(i,j,3) .eq. "s" .and. u_gs(i,j,5) .eq. "O") then
                bound = bound
                at_s = 0.0
                
            else if (u_gs(i,j,3) .eq. "s" .and. u_gs(i,j,5) .eq. "L" .and. u_gs(i,j,6) .eq. "U") then
                !bound = bound + at_s * u(i,j-1)
            else if (u_gs(i,j,3) .eq. "s" .and. u_gs(i,j,5) .eq. "L" .and. u_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = bound + at_s * u(i,j-1) !viscosity * u_bound_s(i) * dx / (0.5*dy) !
            else
                print *, "unidentified grid property at",i,j
            end if
            

            
            if (u_gs(i,j,4) .ne. "n") then !N
                counter=counter+1
                umatv(counter)=-at_n
                umati(counter)=(j-1+1)*(nx-1)+(i-1)
                umatj(counter)=(j-1)*(nx-1)+(i-1)
                AA_coe(i,j,4)=-at_n
                !print *, i,j,counter,umatv(counter)
            else if (u_gs(i,j,4) .eq. "n" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "c") then
                at_c = at_c + viscosity * dx /(0.5*dy)
                bound = viscosity * u_bound_n(i,timestep_id) * dx / (0.5*dy)
             else if (u_gs(i,j,4) .eq. "n" .and. u_gs(i,j,5) .eq. "B" .and. u_gs(i,j,6) .eq. "T") then
                at_n = 0.0
                AA_coe(i,j,4) = -at_n
                at_c = at_c + viscosity * dx /(0.5*dy)
                bound = viscosity * u_bound_n(i,timestep_id) * dx / (0.5*dy)
            else if (u_gs(i,j,4) .eq. "n" .and. u_gs(i,j,5) .eq. "O") then
                bound = bound
                at_n = 0.0
                
            else if (u_gs(i,j,4) .eq. "n" .and. u_gs(i,j,5) .eq. "L" .and. u_gs(i,j,6) .eq. "U") then
                !bound = bound + at_n * u(i,j+1)
            else if (u_gs(i,j,4) .eq. "n" .and. u_gs(i,j,5) .eq. "L" .and. u_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = bound + at_n * u(i,j+1) !viscosity * u_bound_n(i) * dx / (0.5*dy)!
            else
                print *, "unidentified grid property at",i,j
            end if
            


            AA_coe(i,j,5)=at_c/0.7
            
            counter=counter+1
            
            umatv(counter)=at_c/0.7!-sp
            umati(counter)=(j-1)*(nx-1)+(i-1)
            umatj(counter)=(j-1)*(nx-1)+(i-1)
            ap_u(i,j)=at_c
            !print *, i,j,counter,umatv(counter)
            !print *, i,j,counter,-at_w,umatv(counter-4),-at_e,umatv(counter-3),-at_s,umatv(counter-2),-at_n,umatv(counter-1)

            BB_bound(i,j)=(p(i-1,j)-p(i,j))*dy+((1-0.7)*at_f*at_c/0.7)*u(i,j)+bound+source_term
            
            umrhs((j-1)*(nx-1)+(i-1))=(p(i-1,j)-p(i,j))*dy+((1-0.7)*at_f*at_c/0.7)*u(i,j)+bound+source_term
            
            if (u_gs(i,j,5) .eq. "X") then
                AA_coe(i,j,1)=0.0
                AA_coe(i,j,2)=0.0
                AA_coe(i,j,3)=0.0
                AA_coe(i,j,4)=0.0
                AA_coe(i,j,5)=1.0
                BB_bound(i,j) = 1E-20
            end if
            
            
            !print *,i,j,u_gs(i,j,:),AA_coe(i,j,:),BB_bound(i,j)
            
            
            !print *, BB_bound(i,j)
!PRINT *, i,j,at_w,at_e,at_s,at_n,at_c,umatv(counter),umrhs((j-1)*(nx-1)+(i-1))!(p(i-1,j)-p(i,j))*dy+((1-0.7)*at_c/0.7)*u(i,j)+bound+source_term

            !For ghost cell
            !Assign u ghost cell w to u(1,1:ny)
            !Assign u ghost cell e to u(nx+1,1:ny)
            !u_bound_s(2:nx)=u ghost cell s
            !u_bound_n(2:nx)=u ghost cell n
            !Assign v ghost cell s to v(1:nx,1)
            !Assign v ghost cell n to v(1:nx,ny+1)
            
            !Assign p ghost cell w to p(1,1:ny)
            !Assign p ghost cell e to p(nx,1:ny)
        end do     
    end do


    do j=1,ny
        do i=1,nx+1
            !Recalculate APU
            a_e=0.0d0 ; a_w=0.0d0 ; a_n=0.0d0 ; a_s=0.0d0
            d_e=0.0d0 ; d_w=0.0d0 ; d_n=0.0d0 ; d_s=0.0d0
            at_e=0.0d0 ; at_w=0.0d0 ; at_n=0.0d0 ;at_s=0.0d0 ; at_c=0.0d0
            bt_e=0.0d0 ; bt_w=0.0d0 ; bt_n=0.0d0 ;bt_s=0.0d0 ; bt_c=0.0d0

            a_e = density*( u(i,j)     + u(i+1,j) )*dy/2.0d0
            a_w = density*( u(i-1,j)   + u(i,j)   )*dy/2.0d0
            a_n = density*( v(i-1,j+1) + v(i,j+1) )*dx/2.0d0
            a_s = density*( v(i-1,j)   + v(i,j)   )*dx/2.0d0
            
            d_e=viscosity*dy/dx
            d_w=viscosity*dy/dx
            d_n=viscosity*dx/dy
            d_s=viscosity*dx/dy

            !hybrid format
            at_e=max(-a_e,(d_e-a_e/2.0d0),0.0d0)
            at_w=max(+a_w,(d_w+a_w/2.0d0),0.0d0)
            at_n=max(-a_n,(d_n-a_n/2.0d0),0.0d0)
            at_s=max(+a_s,(d_s+a_s/2.0d0),0.0d0)
            at_c=at_e+at_w+at_n+at_s+(a_e-a_w+a_n-a_s)
            ap_u1(i,j)=at_c
            
        end do     
    end do 

    call coefficient_filler_u(nx,ny,(nx-1)*ny,A_coe,B_bound,AA_coe,BB_bound)

    call Iteration_solver((nx-1)*ny,A_coe,solution_u,B_bound)

    !call Bicgstab_solver(umatv,umati,umatj,unnzero,umrhs,un,solution_u,bires,bimaxit)

!To solve [A].[x]=[b] system of equations by BiCGSTAB method:

    !GET SOLUTION
    do i=2,nx
        do j=1,ny
            u_alter(i,j)=solution_u( i-1+(nx-1)*(j-1) )
            !print *, i,j,u(i,j),u_alter(i,j),"U"
        end do
    end do
    
    !do i=0,nx+2
        !do j=0,ny+2
            !print *, i,j,u(i,j),u_alter(i,j),"U"
        !end do
    !end do
print *, "Solve U complete"

    deallocate(AA_coe)
    print *, "Solve U complete"
    deallocate(BB_bound)
    print *, "Solve U complete"
    deallocate(A_coe)
        print *, "Solve U complete"
    deallocate(B_bound)

    
    deallocate(umatv,umati,umatj,umrhs)
    deallocate(vmatv,vmati,vmatj,vmrhs)
    deallocate(pmatv,pmati,pmatj,pmrhs)

    deallocate(solution_u)
    deallocate(solution_v)
    deallocate(solution_p)
print *, "Solve U complete"
    
        !print *, "Solve U complete"
end subroutine u_solver

subroutine coefficient_filler_v(nx,ny,nn,A_coe,B_bound,AA_coe,BB_bound)

integer nx,ny,nn
real(kind=8),dimension(1:nx,2:ny,5):: AA_coe
real(kind=8),dimension(1:nx,2:ny):: BB_bound

real(kind=8),dimension(1:nx,1:ny-1,5):: AA_coe_temp
real(kind=8),dimension(1:nx,1:ny-1):: BB_bound_temp

real(kind=8),dimension(nn,nn):: A_coe
real(kind=8),dimension(nn):: B_bound

real(kind=8),dimension(nn):: temp

integer counter 
A_coe =0.0d0

do i=1,nx
    do j=2,ny
        AA_coe_temp(i,j-1,1)=AA_coe(i,j,1)
        AA_coe_temp(i,j-1,2)=AA_coe(i,j,2)
        AA_coe_temp(i,j-1,3)=AA_coe(i,j,3)
        AA_coe_temp(i,j-1,4)=AA_coe(i,j,4)
        AA_coe_temp(i,j-1,5)=AA_coe(i,j,5)
        BB_bound_temp(i,j-1)=BB_bound(i,j)
    end do
end do


counter = 0

        do j=1,ny-1
            do i=1,nx
            counter = counter + 1
            temp(:)=0.0d0
            temp(i+(j-1)*(nx))=AA_coe_temp(i,j,5)
            
            if (i+(j-1)*(nx)-1 .ge. 1) then
            temp(i+(j-1)*(nx)-1)=AA_coe_temp(i,j,1)
            end if
                        
            if (i+(j-1)*(nx)+1 .le. nn) then
            temp(i+(j-1)*(nx)+1)=AA_coe_temp(i,j,2)
            end if
            
            if (i+(j-1)*(nx)-nx .ge. 1) then
            temp(i+(j-1)*(nx)-nx)=AA_coe_temp(i,j,3)
            end if

            if (i+(j-1)*(nx)+nx .le. nn) then
            temp(i+(j-1)*(nx)+nx)=AA_coe_temp(i,j,4)
            end if
            
            A_coe(counter,:)=temp(:)
            B_bound(i+(j-1)*(nx))=BB_bound_temp(i,j)
        end do
    end do
    
    do i=1,nn
            !print *,A_coe(i,:), B_bound(i)
    end do

end subroutine coefficient_filler_v

subroutine v_solver(nx,ny,nt,timestep_id,density,viscosity,scale_factor,ghost_cell_size,&
                    x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                    u_gs,v_gs,p_gs,&
                    u,u_alter,v,v_alter,&
                    p,p_alter,p_correct,&
                    f_former,f,f_alter,&
                    a_former,a,a_alter,&
                    b_former,b,b_alter,&
                    ap_u,ap_v,ap_u1,ap_v1,&
                    uvbc,pbc,winds,&
                    ubc_e,ubc_w,&
                    vbc_n,vbc_s,&
                    pbc_e,pbc_w,pbc_n,pbc_s,&
                    u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                    u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                    p_bound_w,p_bound_e,&
                    p_bound_s,p_bound_n,&
                    f_bound_w,f_bound_e,&
                    f_bound_s,f_bound_n,&
                    a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                    a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                    u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                    u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                    p_ghost_cell_w,p_ghost_cell_e,&
                    p_ghost_cell_s,p_ghost_cell_n)
    


integer :: i,j,k,counter,ghost_cell_size
integer :: nx,ny,nt,timestep_id
integer::bimaxit
real(kind=8)::bires

real(kind=8)::dx,dy
real(kind=8)::density
real(kind=8)::viscosity
real(kind=8)::scale_factor

real(kind=8)::temp_new1,temp_new2,temp_new3,temp_new4
real(kind=8)::temp_square,temp_square_derive

real(kind=8)::a_n,a_s,a_w,a_e
real(kind=8)::d_n,d_s,d_w,d_e
real(kind=8)::source_term,bound
real(kind=8)::at_n,at_s,at_w,at_e,at_c
real(kind=8)::bt_n,bt_s,bt_w,bt_e,bt_c

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter

real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v

real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s

character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds
    !The output ghost cell information for block_id
    real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
    
!USED FOR SOLVER
integer :: un,vn,pn
integer :: unnzero,vnnzero,pnnzero
real(kind=8),allocatable,dimension(:):: umatv,umrhs,vmatv,vmrhs,pmatv,pmrhs
integer,allocatable,dimension(:):: umati,umatj,vmati,vmatj,pmati,pmatj
real(kind=8),allocatable,dimension(:)::solution_u,solution_v,solution_p


real(kind=8),allocatable,dimension(:,:,:)::AA_coe
real(kind=8),allocatable,dimension(:,:)::BB_bound
real(kind=8),allocatable,dimension(:,:)::A_coe
real(kind=8),allocatable,dimension(:)::B_bound

!real(kind=8),dimension(1:nx,2:ny,5)::AA_coe
!real(kind=8),dimension(1:nx,2:ny)::BB_bound

!real(kind=8),dimension(nx*(ny-1),nx*(ny-1))::A_coe
!real(kind=8),dimension(nx*(ny-1))::B_bound

un=(nx-1)*ny ; vn=nx*(ny-1) ; pn=nx*ny

unnzero=5*(nx-3)*(ny-2)+4*(nx-3)*2+4*(ny-2)*2+3*4
vnnzero=5*(nx-2)*(ny-3)+4*(nx-2)*2+4*(ny-3)*2+3*4
pnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4

bires=1.0e-20
bimaxit=300

allocate (AA_coe(1:nx,2:ny,5),BB_bound(1:nx,2:ny),A_coe(nx*(ny-1),nx*(ny-1)),B_bound(nx*(ny-1)))

allocate (umatv(1:unnzero),umati(1:unnzero),umatj(1:unnzero),umrhs(1:(nx-1)*ny))
allocate (vmatv(1:vnnzero),vmati(1:vnnzero),vmatj(1:vnnzero),vmrhs(1:nx*(ny-1)))
allocate (pmatv(1:pnnzero),pmati(1:pnnzero),pmatj(1:pnnzero),pmrhs(1:nx*ny))

allocate(solution_u( (nx-1)*ny ) ) ; allocate(solution_v( nx*(ny-1) ) ) ; allocate(solution_p( nx*ny ) )

    umatv=0 ; umati=0 ; umatj=0
    vmatv=0 ; vmati=0 ; vmatj=0
    pmatv=0 ; pmati=0 ; pmatj=0
!USED FOR SOLVER
    solution_u=0.0 ; solution_v=0.0 ; solution_p=0.0
    
        AA_coe = 0.0 ; BB_bound=0.0 ; A_coe = 0.0 ; B_bound = 0.0
!START CALCULATION
!for main body

    counter = 0
    do i=1,nx
        do j=2,ny

            dx=x(i,j)-x(i-1,j) ; dy=y(i,j)-y(i,j-1)
            dx=(domain_x_end-domain_x_start)/nx ; dy=(domain_y_end-domain_y_start)/ny
            !Compute coefficients
            a_e=0.0d0 ; a_w=0.0d0 ; a_n=0.0d0 ; a_s=0.0d0
            d_e=0.0d0 ; d_w=0.0d0 ; d_n=0.0d0 ; d_s=0.0d0
            at_e=0.0d0 ; at_w=0.0d0 ; at_n=0.0d0 ;at_s=0.0d0 ; at_c=0.0d0
            bt_e=0.0d0 ; bt_w=0.0d0 ; bt_n=0.0d0 ;bt_s=0.0d0 ; bt_c=0.0d0
            
            a_e = density*( u(i+1,j-1) + u(i+1,j) )*dy/2.0d0
            a_w = density*( u(i,j-1)   + u(i,j)   )*dy/2.0d0
            a_n = density*( v(i,j)     + v(i,j+1) )*dx/2.0d0
            a_s = density*( v(i,j-1)   + v(i,j)   )*dx/2.0d0
        
            d_e=viscosity*dy/dx ; d_w=viscosity*dy/dx ; d_n=viscosity*dx/dy ; d_s=viscosity*dx/dy

            !hybrid format
            at_e=max(-a_e,(d_e-a_e/2.0),0.0d0)
            at_w=max(+a_w,(d_w+a_w/2.0),0.0d0)
            at_n=max(-a_n,(d_n-a_n/2.0),0.0d0)
            at_s=max(+a_s,(d_s+a_s/2.0),0.0d0)
            at_c=at_e+at_w+at_n+at_s+(a_e-a_w+a_n-a_s)
            at_f=1.0d0
            temp_square=viscosity*( a(i,j)**2+b(i,j)**2 )*( v(i,j)+v(i,j+1) )/2.0d0
            
            temp_square_derive=viscosity*scale_factor*( ( b(i+1,j)-b(i-1,j) )/2.0d0/dx - &
                                                        ( a(i,j+1)-a(i,j-1) )/2.0d0/dy )**2*( v(i,j)+v(i+1,j) )/2.0d0

            temp_new1=viscosity*2*a(i,j)*(u(i+1,j)-u(i,j))/dx
            temp_new2=viscosity*2*b(i,j)*(u(i,j)+u(i,j+1)+u(i+1,j)+u(i+1,j+1)-u(i,j-1)-u(i,j)-u(i+1,j-1)-u(i+1,j))/4.0d0/dy
            temp_new3=viscosity*(u(i,j)+u(i+1,j))/2.0d0*(a(i+1,j)-a(i-1,j))/2.0d0/dx
            temp_new4=viscosity*(u(i,j)+u(i+1,j))/2.0d0*(b(i,j+1)-b(i,j-1))/2.0d0/dy
        
            !Compute source term
            source_term=0.0d0
            source_term = - temp_square - temp_square_derive - temp_new1 - temp_new2 - temp_new3 - temp_new4
            bound=0.0d0
            
            if (v_gs(i,j,1) .ne. "w") then !W
                counter=counter+1
                vmatv(counter)=-at_w
                vmati(counter)=(j-2)*(nx)+(i-1)
                vmatj(counter)=(j-2)*(nx)+(i)
                AA_coe(i,j,1)=-at_w
            else if (v_gs(i,j,1) .eq. "w" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "c") then
                at_c = at_c + viscosity * dy / (0.5*dx)
                bound = viscosity * v_bound_w(j,timestep_id) * dy / (0.5*dx)
            else if (v_gs(i,j,1) .eq. "w" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "T") then
                at_w = 0.0
                AA_coe(i,j,1) = -at_w
                at_c = at_c + viscosity * dy / (0.5*dx)
                bound = viscosity * v_bound_w(j,timestep_id) * dy / (0.5*dx)
            else if (v_gs(i,j,1) .eq. "w" .and. v_gs(i,j,5) .eq. "O") then
                !at_c = at_c + viscosity * dy / (0.5*dx)
                bound = bound
                at_w = 0.0
            
            else if (v_gs(i,j,1) .eq. "w" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "U") then
                !bound = bound + at_w * v(i-1,j) !viscosity * v_bound_w(j) * dy / (0.5*dx) !
            else if (v_gs(i,j,1) .eq. "w" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = 0.0 !bound + at_w * v(i-1,j) !viscosity * v_bound_w(j) * dy / (0.5*dx) !
            else
                print *, "unidentified grid property at",i,j
            end if
            
    
            if (v_gs(i,j,2) .ne. "e") then !E
                counter=counter+1
                vmatv(counter)=-at_e
                vmati(counter)=(j-2)*(nx)+(i+1)
                vmatj(counter)=(j-2)*(nx)+(i)
                AA_coe(i,j,2)=-at_e
            else if (v_gs(i,j,2) .eq. "e" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "c") then
                at_c = at_c + viscosity * dy / (0.5*dx)
                bound = viscosity * v_bound_e(j,timestep_id) * dy / (0.5*dx)
            else if (v_gs(i,j,2) .eq. "e" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "T") then
                at_e = 0.0
                AA_coe(i,j,2) = -at_e
                at_c = at_c + viscosity * dy / (0.5*dx)
                bound = viscosity * v_bound_e(j,timestep_id) * dy / (0.5*dx)
            else if (v_gs(i,j,2) .eq. "e" .and. v_gs(i,j,5) .eq. "O") then
                !at_c = at_c + viscosity * dy / (0.5*dx)
                bound = bound
                at_e = 0.0
                
            else if (v_gs(i,j,2) .eq. "e" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "U") then
                !bound = bound + at_e * v(i+1,j)
            else if (v_gs(i,j,2) .eq. "e" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = 0.0 !bound !+ at_e * v(i+1,j) !viscosity * v_bound_e(j) * dy / (0.5*dx)!
            else 
                print *, "unidentified grid property at",i,j
            end if


            if (v_gs(i,j,3) .ne. "s") then !S
                counter=counter+1
                vmatv(counter)=-at_s
                vmati(counter)=(j-2-1)*(nx)+(i)
                vmatj(counter)=(j-2)*(nx)+(i)
                AA_coe(i,j,3)=-at_s
            else if (v_gs(i,j,3) .eq. "s" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "c") then
                bound = bound + at_s * v(i,j-1)
            else if (v_gs(i,j,3) .eq. "s" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "T") then
                at_s = 0.0 
                AA_coe(i,j,3) = -at_s
                bound = bound + at_s * v(i,j-1)
            else if (v_gs(i,j,3) .eq. "s" .and. v_gs(i,j,5) .eq. "O") then
                bound = bound   
                at_s = 0.0
                
            else if (v_gs(i,j,3) .eq. "s" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "U") then
                !bound = bound + at_s * v(i,j-1)
            else if (v_gs(i,j,3) .eq. "s" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = 0.0 
            else
                print *, "unidentified grid property at",i,j
            end if

            if (v_gs(i,j,4) .ne. "n") then !N
                counter=counter+1
                vmatv(counter)=-at_n
                vmati(counter)=(j-2+1)*(nx)+(i)
                vmatj(counter)=(j-2)*(nx)+(i)
                AA_coe(i,j,4)=-at_n
            else if (v_gs(i,j,4) .eq. "n" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "c") then
                bound = bound + at_n * v(i,j+1)
            else if (v_gs(i,j,4) .eq. "n" .and. v_gs(i,j,5) .eq. "B" .and. v_gs(i,j,6) .eq. "T") then
                at_n = 0.0 
                AA_coe(i,j,4) = -at_n
                bound = bound + at_n * v(i,j+1)
            else if (v_gs(i,j,4) .eq. "n" .and. v_gs(i,j,5) .eq. "O") then
                bound = bound
                at_n = 0.0

            else if (v_gs(i,j,4) .eq. "n" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "U") then
                !bound = bound + at_n * v(i,j+1)
            else if (v_gs(i,j,4) .eq. "n" .and. v_gs(i,j,5) .eq. "L" .and. v_gs(i,j,6) .eq. "D") then
                !at_e = 0.0; at_w = 0.0; at_n = 0.0; at_s = 0.0; at_c = 1.0*0.7
                !at_f=1.0/0.3
                !bound = 0.0 
            else
                print *, "unidentified grid property at",i,j
            end if
            
            AA_coe(i,j,5)=at_c/0.7
            
            counter=counter+1

            vmatv(counter)=at_c/0.7!-sp
            vmati(counter)=(j-2)*(nx)+(i)
            vmatj(counter)=(j-2)*(nx)+(i)
            
            ap_v(i,j)=at_c
            
            !PRINT *, i,j,at_w,at_e,at_s,at_n,at_c,(p(i,j-1)-p(i,j))*dx+((1-0.7)*at_c/0.7)*v(i,j)+bound+source_term
            BB_bound(i,j)=(p(i,j-1)-p(i,j))*dx+((1-0.7)*at_f*at_c/0.7)*v(i,j)+bound+source_term
            
            vmrhs((j-2)*(nx)+(i))=(p(i,j-1)-p(i,j))*dx+((1-0.7)*at_f*at_c/0.7)*v(i,j)+bound+source_term
            
            if (v_gs(i,j,5) .eq. "X") then
                AA_coe(i,j,1)=0.0
                AA_coe(i,j,2)=0.0
                AA_coe(i,j,3)=0.0
                AA_coe(i,j,4)=0.0
                AA_coe(i,j,5)=1.0
                BB_bound(i,j) = 1E-20
            end if
            
            !print *,i,j,v_gs(i,j,:),AA_coe(i,j,:),BB_bound(i,j)
            !For ghost cell
            !Assign u ghost cell w to u(1,1:ny)
            !Assign u ghost cell e to u(nx+1,1:ny)
            !Assign v ghost cell s to v(1:nx,1)
            !Assign v ghost cell n to v(1:nx,ny+1)
            !v_bound_w(2:ny)=v ghost cell w
            !v_bound_e(2:ny)=v ghost cell w
            
            !Assign p ghost cell s to p(1:nx,1)
            !Assign p ghost cell n to p(1:nx,ny)
        end do     
    end do 

        do j=1,ny+1
            do i=1,nx
            dx=x(i,j)-x(i-1,j)
            dy=y(i,j)-y(i,j-1)
            dx=(domain_x_end-domain_x_start)/nx
            dy=(domain_y_end-domain_y_start)/ny
            !Compute coefficients
            a_e=0.0d0 ; a_w=0.0d0 ; a_n=0.0d0 ; a_s=0.0d0
            d_e=0.0d0 ; d_w=0.0d0 ; d_n=0.0d0 ; d_s=0.0d0
            at_e=0.0d0 ; at_w=0.0d0 ; at_n=0.0d0 ;at_s=0.0d0 ; at_c=0.0d0
            bt_e=0.0d0 ; bt_w=0.0d0 ; bt_n=0.0d0 ;bt_s=0.0d0 ; bt_c=0.0d0
            
            a_e = density*( u(i+1,j-1) + u(i+1,j) )*dy/2.0d0
            a_w = density*( u(i,j-1)   + u(i,j)   )*dy/2.0d0
            a_n = density*( v(i,j)     + v(i,j+1) )*dx/2.0d0
            a_s = density*( v(i,j-1)   + v(i,j)   )*dx/2.0d0
        
            d_e=viscosity*dy/dx
            d_w=viscosity*dy/dx
            d_n=viscosity*dx/dy
            d_s=viscosity*dx/dy

            !hybrid format
            at_e=max(-a_e,(d_e-a_e/2.0),0.0d0)
            at_w=max(+a_w,(d_w+a_w/2.0),0.0d0)
            at_n=max(-a_n,(d_n-a_n/2.0),0.0d0)
            at_s=max(+a_s,(d_s+a_s/2.0),0.0d0)
            at_c=at_e+at_w+at_n+at_s+(a_e-a_w+a_n-a_s)
            
            ap_v1(i,j)=at_c
        end do     
    end do 
    
    call coefficient_filler_v(nx,ny,nx*(ny-1),A_coe,B_bound,AA_coe,BB_bound)

    call Iteration_solver(nx*(ny-1),A_coe,solution_v,B_bound)
    
    !call Bicgstab_solver(vmatv,vmati,vmatj,vnnzero,vmrhs,vn,solution_v,bires,bimaxit)
    !get solution
    do i=1,nx
        do j=2,ny
            v_alter(i,j)=solution_v( i+nx*(j-1-1) )
            !print *, i,j,v(i,j),v_alter(i,j),"V"
        end do
    end do
    

    do i=0,nx+2
        do j=0,ny+2
          !print *, i,j,v(i,j),v_alter(i,j),"V"
        end do
    end do

    deallocate (umatv,umati,umatj,umrhs)
    deallocate (vmatv,vmati,vmatj,vmrhs)
    deallocate (pmatv,pmati,pmatj,pmrhs)

    deallocate(solution_u)
    deallocate(solution_v)
    deallocate(solution_p)

    deallocate (AA_coe,BB_bound,A_coe,B_bound)
    
    print *, "Solve V complete"
    
end subroutine v_solver

subroutine uv_correction(nx,ny,u_alter,v_alter,pbc_e,pbc_w,pbc_n,pbc_s)
    
integer ::i,j,nx,ny
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s

    !Boundary Velocity Correction for Peressure BCs:
    do j=1,ny
        if (pbc_w.eq."C") then !W
            u_alter(1,j)    = u_alter(2,j)+v_alter(1,j+1)-v_alter(1,j)
        end if
        if (pbc_e.eq."C") then !E
            u_alter(nx+1,j) = u_alter(nx,j)+v_alter(nx,j)-v_alter(nx,j+1)
        end if
    end do
    
    do i=1,nx
        if (pbc_n.eq."C") then !N
            v_alter(i,ny+1) = v_alter(i,ny)+u_alter(i,ny)-u_alter(i+1,ny)
        end if
        if (pbc_s.eq."C") then !S
            v_alter(i,1)    = v_alter(i,2)+u_alter(i,2)-u_alter(i+1,2)
        end if
    end do

end subroutine uv_correction

subroutine coefficient_filler_p(nx,ny,nn,A_coe,B_bound,AA_coe,BB_bound)

integer nx,ny,nn
real(kind=8),dimension(1:nx,1:ny,5):: AA_coe
real(kind=8),dimension(1:nx,1:ny):: BB_bound

real(kind=8),dimension(1:nx,1:ny,5):: AA_coe_temp
real(kind=8),dimension(1:nx,1:ny):: BB_bound_temp

real(kind=8),dimension(nn,nn):: A_coe
real(kind=8),dimension(nn):: B_bound

real(kind=8),dimension(nn):: temp

integer counter 
A_coe =0.0d0

do i=1,nx
    do j=1,ny
        AA_coe_temp(i,j,1)=AA_coe(i,j,1)
        AA_coe_temp(i,j,2)=AA_coe(i,j,2)
        AA_coe_temp(i,j,3)=AA_coe(i,j,3)
        AA_coe_temp(i,j,4)=AA_coe(i,j,4)
        AA_coe_temp(i,j,5)=AA_coe(i,j,5)
        BB_bound_temp(i,j)=BB_bound(i,j)
    end do
end do


counter = 0
    do j=1,ny
        do i=1,nx
            counter = counter + 1
            temp(:)=0.0d0
            temp(i+(j-1)*(nx))=AA_coe_temp(i,j,5)
            
            if (i+(j-1)*(nx)-1 .ge. 1) then
            temp(i+(j-1)*(nx)-1)=AA_coe_temp(i,j,1)
            end if
                        
            if (i+(j-1)*(nx)+1 .le. nn) then
            temp(i+(j-1)*(nx)+1)=AA_coe_temp(i,j,2)
            end if
            
            if (i+(j-1)*(nx)-nx .ge. 1) then
            temp(i+(j-1)*(nx)-nx)=AA_coe_temp(i,j,3)
            end if

            if (i+(j-1)*(nx)+nx .le. nn) then
            temp(i+(j-1)*(nx)+nx)=AA_coe_temp(i,j,4)
            end if
            
            A_coe(counter,:)=temp(:)
            B_bound(i+(j-1)*(nx))=BB_bound_temp(i,j)
        end do
    end do
    
    do i=1,nn
            !print *,A_coe(i,:), B_bound(i)
    end do

end subroutine coefficient_filler_p

subroutine p_solver(nx,ny,nt,timestep_id,density,viscosity,scale_factor,ghost_cell_size,&
                    x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                    u_gs,v_gs,p_gs,&
                    u,u_alter,v,v_alter,&
                    p,p_alter,p_correct,&
                    f_former,f,f_alter,&
                    a_former,a,a_alter,&
                    b_former,b,b_alter,&
                    ap_u,ap_v,ap_u1,ap_v1,&
                    uvbc,pbc,winds,&
                    ubc_e,ubc_w,&
                    vbc_n,vbc_s,&
                    pbc_e,pbc_w,pbc_n,pbc_s,&
                    u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                    u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                    p_bound_w,p_bound_e,&
                    p_bound_s,p_bound_n,&
                    f_bound_w,f_bound_e,&
                    f_bound_s,f_bound_n,&
                    a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                    a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                    u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                    u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                    p_ghost_cell_w,p_ghost_cell_e,&
                    p_ghost_cell_s,p_ghost_cell_n)



integer :: i,j,k,counter,ghost_cell_size
integer :: nx,ny,nt,timestep_id
integer::bimaxit
real(kind=8)::bires,con_res

real(kind=8)::dx,dy
real(kind=8)::density
real(kind=8)::viscosity
real(kind=8)::scale_factor

real(kind=8)::temp_new1,temp_new2,temp_new3,temp_new4
real(kind=8)::temp_square,temp_square_derive

real(kind=8)::a_n,a_s,a_w,a_e
real(kind=8)::d_n,d_s,d_w,d_e
real(kind=8)::source_term,bound
real(kind=8)::at_n,at_s,at_w,at_e,at_c
real(kind=8)::bt_n,bt_s,bt_w,bt_e,bt_c

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter

real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v

real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds

    !The output ghost cell information for block_id
    real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
    
!USED FOR SOLVER
integer :: un,vn,pn
integer :: unnzero,vnnzero,pnnzero
real(kind=8),allocatable,dimension(:):: umatv,umrhs,vmatv,vmrhs,pmatv,pmrhs
integer,allocatable,dimension(:):: umati,umatj,vmati,vmatj,pmati,pmatj
real(kind=8),allocatable,dimension(:)::solution_u,solution_v,solution_p

real(kind=8),allocatable,dimension(:,:,:)::AA_coe
real(kind=8),allocatable,dimension(:,:)::BB_bound
real(kind=8),allocatable,dimension(:,:)::A_coe
real(kind=8),allocatable,dimension(:)::B_bound
allocate (AA_coe(1:nx,1:ny,5),BB_bound(1:nx,1:ny),A_coe(nx*ny,nx*ny),B_bound(nx*ny))

AA_coe = 0.0
BB_bound = 0.0 

!real(kind=8),dimension(1:nx,1:ny,5)::AA_coe
!real(kind=8),dimension(1:nx,1:ny)::BB_bound

!real(kind=8),dimension(nx*ny,nx*ny)::A_coe
!real(kind=8),dimension(nx*ny)::B_bound

un=(nx-1)*ny ; vn=nx*(ny-1); pn=nx*ny

unnzero=5*(nx-3)*(ny-2)+4*(nx-3)*2+4*(ny-2)*2+3*4
vnnzero=5*(nx-2)*(ny-3)+4*(nx-2)*2+4*(ny-3)*2+3*4
pnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4

bires=1.0e-20
bimaxit=300

allocate (umatv(1:unnzero),umati(1:unnzero),umatj(1:unnzero),umrhs(1:(nx-1)*ny))
allocate (vmatv(1:vnnzero),vmati(1:vnnzero),vmatj(1:vnnzero),vmrhs(1:nx*(ny-1)))
allocate (pmatv(1:pnnzero),pmati(1:pnnzero),pmatj(1:pnnzero),pmrhs(1:nx*ny))

allocate(solution_u( (nx-1)*ny ) ) ; allocate(solution_v( nx*(ny-1) ) ) ; allocate(solution_p( nx*ny ) )
    solution_u=0.0 ; solution_v=0.0 ; solution_p=0.0

    pmatv = 0 ; pmati = 0 ; pmatj = 0
!START CALCULATION
!for main body

    counter=0
    do i=1,nx 
        do j=1,ny
                       
            dx=x(i+1,j)-x(i,j) ; dy=y(i,j+1)-y(i,j)
            dx=(domain_x_end-domain_x_start)/nx
            dy=(domain_y_end-domain_y_start)/ny
            !SET ALL TO ZERO
            at_e=0.0d0 ; at_w=0.0d0 ; at_n=0.0d0 ; at_s=0.0d0 ; at_c=0.0d0
            !SET NON-BOUNDARY
            !if (i/=nx) at_e = density*dy/ap_u(i+1,j) *dy*0.7d0
            !if (i/=1 ) at_w = density*dy/ap_u(i,j)   *dy*0.7d0
            !if (j/=ny) at_n = density*dx/ap_v(i,j+1) *dx*0.7d0
            !if (j/=1 ) at_s = density*dx/ap_v(i,j)   *dx*0.7d0
            bound = 0.0d0
            if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "B") then
            at_w = 0.0d0   
            else
            at_w = density*dy/ap_u1(i,j)   *dy*0.7d0
            end if
            
            if (p_gs(i,j,2) .eq. "e" .and. p_gs(i,j,5) .eq. "B") then
            at_e = 0.0d0
            else 
            at_e = density*dy/ap_u1(i+1,j) *dy*0.7d0
            end if
            
            if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "B") then
            at_s = 0.0d0
            else
            at_s = density*dx/ap_v1(i,j)   *dx*0.7d0
            end if
            
            if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "B") then
            at_n = 0.0d0
            else 
            at_n = density*dx/ap_v1(i,j+1) *dx*0.7d0
            end if
            !Compute coefficient
            
            !print *, i,j,ap_u1(i,j),ap_u1(i+1,j),ap_v1(i,j),ap_v1(i,j+1)
            !print *,i,j,AA_coe(i,j,1)
            if ((p_gs(i,j,1) .eq. "w").and.(p_gs(i,j,5) .eq. "B").and.(pbc_w.eq."C")) then !P BC on W
!PRINT *, i,j,at_e,at_w,at_n,at_s,at_c,density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) )
                counter=counter+1
                at_c=0.01d0
                pmatv(counter)=at_c
                pmati(counter)=(j-1)*nx+(i)
                pmatj(counter)=(j-1)*nx+(i)
                pmrhs((j-1)*nx+(i))=0.0d0
                AA_coe(i,j,5)=at_c
                !print *, i,j,"W"
                
                cycle
            end if
            
            if ((p_gs(i,j,2) .eq. "e").and.(p_gs(i,j,5) .eq. "B").and.(pbc_e.eq."C")) then !P BC on E
!PRINT *, i,j,at_e,at_w,at_n,at_s,at_c,density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) )
                counter=counter+1
                at_c=0.01d0
                pmatv(counter)=at_c
                pmati(counter)=(j-1)*nx+(i)
                pmatj(counter)=(j-1)*nx+(i)
                pmrhs((j-1)*nx+(i))=0.0d0
                AA_coe(i,j,5)=at_c
               ! print *, i,j,"E"
                cycle
            end if

            if ((p_gs(i,j,3) .eq. "s").and.(p_gs(i,j,5) .eq. "B").and.(pbc_s.eq."C")) then !P BC on S
!PRINT *, i,j,at_e,at_w,at_n,at_s,at_c,density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) )
                counter=counter+1
                at_c=0.01d0
                pmatv(counter)=at_c
                pmati(counter)=(j-1)*nx+(i)
                pmatj(counter)=(j-1)*nx+(i)
                pmrhs((j-1)*nx+(i))=0.0d0
                AA_coe(i,j,5)=at_c
                !print *, i,j,"S"
                cycle
            end if

            if ((p_gs(i,j,4) .eq. "n").and.(p_gs(i,j,5) .eq. "B").and.(pbc_n.eq."C")) then !P BC on N
!PRINT *, i,j,at_e,at_w,at_n,at_s,at_c,density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) )
                counter=counter+1
                at_c=0.01d0
                pmatv(counter)=at_c
                pmati(counter)=(j-1)*nx+(i)
                pmatj(counter)=(j-1)*nx+(i)
                pmrhs((j-1)*nx+(i))=0.0d0
                AA_coe(i,j,5)=at_c
                !print *, i,j,"N"
                cycle
            end if
            !Pressure free condition

            

            
            if (p_gs(i,j,1) .ne. "w") then
                counter=counter+1
                pmatv(counter)=-at_w
                pmati(counter)=(j-1)*nx+(i-1)
                pmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,1)=-at_w
                
            else if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "L") then
                !bound = bound + at_w * p_ghost_cell_w(2,j)
            else
                !print *,"Error at",i,j,p_gs(i,j,1),p_gs(i,j,5)
            end if

            if (p_gs(i,j,2) .ne. "e") then
                counter=counter+1
                pmatv(counter)=-at_e
                pmati(counter)=(j-1)*nx+(i+1)
                pmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,2)=-at_e
            else if (p_gs(i,j,2) .eq. "e") then
                !bound = bound + at_e * p_ghost_cell_e(2,j)
            else
                !print *,"Error at",i,j,p_gs(i,j,2),p_gs(i,j,5)
            end if

            if (p_gs(i,j,3) .ne. "s") then
                counter=counter+1
                pmatv(counter)=-at_s
                pmati(counter)=(j-1-1)*nx+(i)
                pmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,3)=-at_s
            else if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "L") then
                !bound = bound + at_s * p_ghost_cell_s(2,i)
            else
                !print *,"Error at",i,j,p_gs(i,j,3),p_gs(i,j,5)
            end if
        
            if (p_gs(i,j,4) .ne. "n") then
                counter=counter+1
                pmatv(counter)=-at_n
                pmati(counter)=(j-1+1)*nx+(i)
                pmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,4)=-at_n
            else if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "L") then
                !bound = bound + at_n * p_ghost_cell_n(2,i)
            else
                !print *,"Error at",i,j,p_gs(i,j,4),p_gs(i,j,5)
            end if
            
            !print *,i,j,AA_coe(i,j,1)
                    
            at_c = at_c + at_e + at_w + at_n + at_s
            counter=counter+1
            pmatv(counter)=at_c
            pmati(counter)=(j-1)*nx+(i)
            pmatj(counter)=(j-1)*nx+(i)
            
            AA_coe(i,j,5)=at_c/0.7
            BB_bound(i,j)=density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) ) + bound
!PRINT *, i,j,ap_u1(i,j),ap_v1(i,j)
!PRINT *, i,j,at_e,at_w,at_n,at_s,at_c,density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) )
!PRINT *, i,j,u_alter(i,j),u_alter(i+1,j),v_alter(i,j),v_alter(i,j+1)
            !print *, pmatj(counter)
            pmrhs((j-1)*nx+(i))=density*dy*( u_alter(i,j) - u_alter(i+1,j) ) + density*dx*( v_alter(i,j) - v_alter(i,j+1) ) + bound
            !print *, i,j, p_gs(i,j,5) 
            
            if (p_gs(i,j,5) .eq. "X") then
                AA_coe(i,j,:)=0.0
                AA_coe(i,j,5)=1.0
                BB_bound(i,j) = p(i,j)
            end if
            
            !print *, i,j,u_alter(i,j),v_alter(i,j)
            !print *,i,j,p_gs(i,j,:),AA_coe(i,j,:),BB_bound(i,j)
            
        end do     
    end do 
    
    
    

            

    !REFLESH P CORRECT
    p_correct=0.0

    !call Bicgstab_solver(pmatv,pmati,pmatj,pnnzero,pmrhs,pn,solution_p,bires,bimaxit)
    
    call coefficient_filler_p(nx,ny,nx*ny,A_coe,B_bound,AA_coe,BB_bound)

    call Iteration_solver(nx*ny,A_coe,solution_p,B_bound)
    
    do i=1,nx
        do j=1,ny
            !correct value
            p_correct(i,j)=solution_p( i+nx*(j-1) )
            !p_correct(i,j)=p_correct(i,j)-p_correct(50,20)
            !print *, i,j,p_correct(i,j),"P"
            !p_alter(i,j)=p(i,j)+p_correct(i,j)*0.3d0
        end do
    end do
    

     do i=1,nx
        do j=1,ny
            p_alter(i,j)=p(i,j)+p_correct(i,j)*0.3d0
            !print *, i,j,p_alter(i,j)
        end do
    end do
    
    print *, "Solve P complete"
    
    con_res=sqrt(dot_product(pmrhs,pmrhs))/pn
    print*,"SIMPLE Iter=",timestep_id,"Continuity Res.=",con_res
    !write(88,*),"SIMPLE Iter=",nt,"Continuity Res.=",con_res
    write(88,*) timestep_id,con_res
    

    flush(unit=88)

    
    deallocate (umatv,umati,umatj,umrhs)
    deallocate (vmatv,vmati,vmatj,vmrhs)
    deallocate (pmatv,pmati,pmatj,pmrhs)

    deallocate(solution_u)
    deallocate(solution_v)
    deallocate(solution_p)
    
    deallocate (AA_coe,BB_bound,A_coe,B_bound)
end subroutine p_solver

subroutine f_solver(nx,ny,nt,dt,timestep_id,density,viscosity,scale_factor,ghost_cell_size,&
                    x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                    u_gs,v_gs,p_gs,&
                    u,u_alter,v,v_alter,&
                    p,p_alter,p_correct,&
                    f_former,f,f_alter,&
                    a_former,a,a_alter,&
                    b_former,b,b_alter,&
                    ap_u,ap_v,ap_u1,ap_v1,&
                    uvbc,pbc,winds,&
                    ubc_e,ubc_w,&
                    vbc_n,vbc_s,&
                    pbc_e,pbc_w,pbc_n,pbc_s,&
                    u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                    u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                    p_bound_w,p_bound_e,&
                    p_bound_s,p_bound_n,&
                    f_bound_w,f_bound_e,&
                    f_bound_s,f_bound_n,&
                    a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                    a_bound_s,a_bound_n,b_bound_s,b_bound_n)

integer :: i,j,k,counter,ghost_cell_size
integer :: nx,ny,nt,timestep_id
integer::bimaxit
real(kind=8)::bires,con_res

real(kind=8)::dx,dy,dt
real(kind=8)::density,viscosity,scale_factor

real(kind=8)::temp_f1,temp_f2
real(kind=8)::temp_square,temp_square_derive

real(kind=8)::a_n,a_s,a_w,a_e,d_n,d_s,d_w,d_e
real(kind=8)::at_n,at_s,at_w,at_e,at_c,bt_n,bt_s,bt_w,bt_e,bt_c
real(kind=8)::source_term,bound

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter

real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v

real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds


    
!USED FOR SOLVER
integer :: un,vn,pn,fn
integer :: unnzero,vnnzero,pnnzero,fnnzero
real(kind=8),allocatable,dimension(:):: umatv,umrhs,vmatv,vmrhs,pmatv,pmrhs,fmatv,fmrhs
integer,allocatable,dimension(:):: umati,umatj,vmati,vmatj,pmati,pmatj,fmati,fmatj
real(kind=8),allocatable,dimension(:)::solution_u,solution_v,solution_p,solution_f

real(kind=8),allocatable,dimension(:,:,:)::AA_coe
real(kind=8),allocatable,dimension(:,:)::BB_bound
real(kind=8),allocatable,dimension(:,:)::A_coe
real(kind=8),allocatable,dimension(:)::B_bound
allocate (AA_coe(1:nx,1:ny,5),BB_bound(1:nx,1:ny),A_coe(nx*ny,nx*ny),B_bound(nx*ny))

AA_coe = 0.0
BB_bound = 0.0 

un=(nx-1)*ny ; vn=nx*(ny-1); pn=nx*ny ; fn=nx*ny

unnzero=5*(nx-3)*(ny-2)+4*(nx-3)*2+4*(ny-2)*2+3*4
vnnzero=5*(nx-2)*(ny-3)+4*(nx-2)*2+4*(ny-3)*2+3*4
pnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4
fnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4

bires=1.0e-20
bimaxit=300

allocate (umatv(1:unnzero),umati(1:unnzero),umatj(1:unnzero),umrhs(1:(nx-1)*ny))
allocate (vmatv(1:vnnzero),vmati(1:vnnzero),vmatj(1:vnnzero),vmrhs(1:nx*(ny-1)))
allocate (pmatv(1:pnnzero),pmati(1:pnnzero),pmatj(1:pnnzero),pmrhs(1:nx*ny))
allocate (fmatv(1:pnnzero),fmati(1:pnnzero),fmatj(1:pnnzero),fmrhs(1:nx*ny))

allocate(solution_u( (nx-1)*ny ) );allocate(solution_v( nx*(ny-1) ) );allocate(solution_p( nx*ny ) );allocate(solution_f( nx*ny ) )

    solution_u=0.0 ; solution_v=0.0 ; solution_p=0.0 ; solution_f=0.0
    fmatv = 0 ; fmati = 0 ; fmatj = 0

!PRINT *, nt,"NT"
!init /boundary condition

!START CALCULATION
!for main body
    !PRINT *, "FFF,",b(:,ny+1)   

write (99,*)  "timestep",nt
    counter=0
    do i=1,nx 
        do j=1,ny
            dx=x(i+1,j)-x(i,j)
            dy=y(i,j+1)-y(i,j)
            dx=(domain_x_end-domain_x_start)/nx
            dy=(domain_y_end-domain_y_start)/ny
            
            at_e=scale_factor/dx/dx
            at_w=scale_factor/dx/dx
            at_n=scale_factor/dy/dy
            at_s=scale_factor/dy/dy
            at_c= 1.0 + 2.0d0*scale_factor/dx/dx + 2.0d0*scale_factor/dy/dy
            
            bound = 0.0d0
            temp_f1 = 0.0d0 ; temp_f2 = 0.0d0 
            
            if (p_gs(i,j,1) .ne. "w") then !W
                counter=counter+1
                fmatv(counter)=-at_w
                fmati(counter)=(j-1)*nx+(i-1)
                fmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,1)=-at_w
            else if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "c") then
                bound = bound + at_w * f(i-1,j)!f_bound_w(j)
            else if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "T") then
                at_w = 0.0 
                AA_coe(i,j,1) = -at_w
                bound = bound + at_w * 0.0!f_bound_w(j)
            else if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "O") then
                at_w = 0.0
                bound = bound
            else if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "U") then

            else if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "D") then

            else
                print *, "unidentified grid property at",i,j
            end if
            
            
            if (p_gs(i,j,2) .ne. "e") then !E
                counter=counter+1
                fmatv(counter)=-at_e
                fmati(counter)=(j-1)*nx+(i+1)
                fmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,2)=-at_e
            else if (p_gs(i,j,2) .eq. "e" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "c") then
                bound = bound + at_e * f(i+1,j)!f_bound_e(j)
            else if (p_gs(i,j,2) .eq. "e" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "T") then
                at_e = 0.0 
                AA_coe(i,j,2) = -at_e
                bound = bound + at_e * 0.0!f_bound_e(j)
            else if (p_gs(i,j,2) .eq. "e" .and. p_gs(i,j,5) .eq. "O") then
                at_e = 0.0
                bound = bound
            else if (p_gs(i,j,2) .eq. "e" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "U") then

            else if (p_gs(i,j,2) .eq. "e" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "D") then

            else
                print *, "unidentified grid property at",i,j
            end if
            
            
            if (p_gs(i,j,3) .ne. "s") then !s
                counter=counter+1
                fmatv(counter)=-at_s
                fmati(counter)=(j-1-1)*nx+(i)
                fmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,3)=-at_s
            else if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "c") then
                bound = bound + at_s * f(i,j-1)!f_bound_s(i)
            else if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "T") then
                at_s = 0.0 
                AA_coe(i,j,3) = -at_s
                bound = bound + at_s * 0.0!f_bound_s(i)
            else if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "O") then
                at_s = 0.0
                bound = bound
            else if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "U") then

            else if (p_gs(i,j,3) .eq. "s" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "D") then

            else
                print *, "unidentified grid property at",i,j
            end if


            if (p_gs(i,j,4) .ne. "n") then !s
                counter=counter+1
                fmatv(counter)=-at_n
                fmati(counter)=(j-1+1)*nx+(i)
                fmatj(counter)=(j-1)*nx+(i)
                AA_coe(i,j,4)=-at_n
            else if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "c") then
                bound = bound + at_n * f(i,j+1)!f_bound_n(i)
            else if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "B" .and. p_gs(i,j,6) .eq. "T") then
                at_n = 0.0 
                AA_coe(i,j,4) = -at_n
                bound = bound + at_n * 0.0!f_bound_n(i)
            else if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "O") then
                at_n = 0.0
                bound = bound
            else if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "U") then

            else if (p_gs(i,j,4) .eq. "n" .and. p_gs(i,j,5) .eq. "L"  .and. p_gs(i,j,6) .eq. "D") then

            else
                print *, "unidentified grid property at",i,j
            end if
            
            
            AA_coe(i,j,5)=at_c/0.7
            
            counter=counter+1
            fmatv(counter)=at_c*0.7
            fmati(counter)=(j-1)*nx+(i)
            fmatj(counter)=(j-1)*nx+(i)
            
            temp_f1=(scale_factor*(a(i+1,j)-a(i-1,j))/2.0d0/dx-scale_factor*(a_former(i+1,j)-a_former(i-1,j))/2.0d0/dx)/dt
            temp_f2=(scale_factor*(b(i,j+1)-b(i,j-1))/2.0d0/dy-scale_factor*(b_former(i,j+1)-b_former(i,j-1))/2.0d0/dy)/dt
            
            BB_bound(i,j)=-1.0*temp_f1 - 1.0*temp_f2 + bound

            fmrhs((j-1)*nx+(i))=-1.0*temp_f1 - 1.0*temp_f2 + bound
            
            !print *, i,j,fmrhs((j-1)*nx+(i))
            if (p_gs(i,j,5) .eq. "X") then
                AA_coe(i,j,:)=0.0
                AA_coe(i,j,5)=1.0
                BB_bound(i,j) = f(i,j)
            end if
            
            !print *, i,j,a(i,j),b(i,j),a_former(i,j),b_former(i,j)
            !print *, i,j,temp_f1,temp_f2,b(i,j+1),b(i,j-1),b_former(i,j+1),b_former(i,j-1)
            !print *, i,j,AA_coe(i,j,:),BB_bound(i,j)
            !write (99,*) i,j,AA_coe(i,j,:),BB_bound(i,j),a_former(i-1,j),a_former(i+1,j),a(i-1,j),a(i+1,j),&
            !b_former(i,j-1),b_former(i,j+1),b(i,j-1),b(i,j+1)
        end do     
    end do 

    do i =0,nx+1
        do j=0,ny+1
            !print *, i,j,temp_f1,temp_f2,b(i,j+1),b(i,j-1),b_former(i,j+1),b_former(i,j-1)
        end do
    end do
                
    call coefficient_filler_p(nx,ny,nx*ny,A_coe,B_bound,AA_coe,BB_bound)

    call Iteration_solver(nx*ny,A_coe,solution_f,B_bound)

    !call Bicgstab_solver(fmatv,fmati,fmatj,fnnzero,fmrhs,fn,solution_f,bires,bimaxit)

    
    
        do i=1,nx
            do j=1,ny
                f_alter(i,j)=solution_f( i+nx*(j-1) )
                !print *, i,j,f_alter(i,j),"F"
                !print *, i,j,f(i,j)
            end do
        end do
    
    !f_alter(0,:)    = f_bound_w(:)
    !f_alter(nx+1,:) = f_bound_e(:)
    !f_alter(:,0)    = f_bound_s(:)
    !f_alter(:,ny+1) = f_bound_n(:)


    !f_alter(0,0)=       1.0*(f_alter(0,1)+f_alter(1,0))         -0.5*(f_alter(0,2)+f_alter(2,0))
    !f_alter(0,ny+1)=    1.0*(f_alter(0,ny)+f_alter(1,ny+1))     -0.5*(f_alter(0,ny-1)+f_alter(2,ny+1))
    !f_alter(nx+1,0)=    1.0*(f_alter(nx+1,1)+f_alter(nx,0))     -0.5*(f_alter(nx+1,2)+f_alter(nx-1,0))
    !f_alter(nx+1,ny+1)= 1.0*(f_alter(nx+1,ny)+f_alter(nx,ny+1)) -0.5*(f_alter(nx+1,ny-1)+f_alter(nx-1,ny+1))
    
    print *, "Solve F complete"
    
    deallocate (umatv,umati,umatj,umrhs)
    deallocate (vmatv,vmati,vmatj,vmrhs)
    deallocate (pmatv,pmati,pmatj,pmrhs)
    deallocate (fmatv,fmati,fmatj,fmrhs)
    
    deallocate(solution_u)
    deallocate(solution_v)
    deallocate(solution_p)
    deallocate(solution_f)
    
    deallocate (AA_coe,BB_bound,A_coe,B_bound)
    
end subroutine f_solver

subroutine ab_solver(nx,ny,nt,dt,timestep_id,density,viscosity,scale_factor,ghost_cell_size,&
                    x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                    u_gs,v_gs,p_gs,&
                    u,u_alter,v,v_alter,&
                    p,p_alter,p_correct,&
                    f_former,f,f_alter,&
                    a_former,a,a_alter,&
                    b_former,b,b_alter,&
                    ap_u,ap_v,ap_u1,ap_v1,&
                    uvbc,pbc,winds,&
                    ubc_e,ubc_w,&
                    vbc_n,vbc_s,&
                    pbc_e,pbc_w,pbc_n,pbc_s,&
                    u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                    u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                    p_bound_w,p_bound_e,&
                    p_bound_s,p_bound_n,&
                    f_bound_w,f_bound_e,&
                    f_bound_s,f_bound_n,&
                    a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                    a_bound_s,a_bound_n,b_bound_s,b_bound_n)

integer :: i,j,k,counter,ghost_cell_size
integer :: nx,ny,nt,timestep_id
integer::bimaxit
real(kind=8)::bires,con_res

real(kind=8)::dx,dy,dt
real(kind=8)::density,viscosity,scale_factor

real(kind=8)::temp_square
real(kind=8)::temp_con1,temp_con2,temp_ft,temp_uu,temp_vv


real(kind=8)::a_n,a_s,a_w,a_e,d_n,d_s,d_w,d_e
real(kind=8)::at_n,at_s,at_w,at_e,at_c,bt_n,bt_s,bt_w,bt_e,bt_c
real(kind=8)::source_term,bound

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8)::domain_x_start,domain_x_end,domain_y_start,domain_y_end

character(len=1),dimension(0:nx+2, 0:ny+2, 6) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 6) ::p_gs
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct
real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter

real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v

real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s,pbc_e,pbc_w,pbc_n,pbc_s
character(len=1),dimension(1,4,3)::uvbc,pbc
character(len=1),dimension(1,4,1)::winds


    
!USED FOR SOLVER
integer :: un,vn,pn,fn
integer :: unnzero,vnnzero,pnnzero,fnnzero
real(kind=8),allocatable,dimension(:):: umatv,umrhs,vmatv,vmrhs,pmatv,pmrhs,fmatv,fmrhs
integer,allocatable,dimension(:):: umati,umatj,vmati,vmatj,pmati,pmatj,fmati,fmatj
real(kind=8),allocatable,dimension(:)::solution_u,solution_v,solution_p,solution_f

real(kind=8),allocatable,dimension(:,:,:)::AA_coe
real(kind=8),allocatable,dimension(:,:)::BB_bound
real(kind=8),allocatable,dimension(:,:)::A_coe
real(kind=8),allocatable,dimension(:)::B_bound
allocate (AA_coe(1:nx,1:ny,5),BB_bound(1:nx,1:ny),A_coe(nx*ny,nx*ny),B_bound(nx*ny))

AA_coe = 0.0
BB_bound = 0.0 

un=(nx-1)*ny ; vn=nx*(ny-1); pn=nx*ny ; fn=nx*ny

unnzero=5*(nx-3)*(ny-2)+4*(nx-3)*2+4*(ny-2)*2+3*4
vnnzero=5*(nx-2)*(ny-3)+4*(nx-2)*2+4*(ny-3)*2+3*4
pnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4
fnnzero=5*(nx-2)*(ny-2)+4*(nx-2)*2+4*(ny-2)*2+3*4

bires=1.0e-20
bimaxit=300

allocate (umatv(1:unnzero),umati(1:unnzero),umatj(1:unnzero),umrhs(1:(nx-1)*ny))
allocate (vmatv(1:vnnzero),vmati(1:vnnzero),vmatj(1:vnnzero),vmrhs(1:nx*(ny-1)))
allocate (pmatv(1:pnnzero),pmati(1:pnnzero),pmatj(1:pnnzero),pmrhs(1:nx*ny))
allocate (fmatv(1:pnnzero),fmati(1:pnnzero),fmatj(1:pnnzero),fmrhs(1:nx*ny))


allocate(solution_u( (nx-1)*ny ) );allocate(solution_v( nx*(ny-1) ) );allocate(solution_p( nx*ny ) );allocate(solution_f( nx*ny ) )

    solution_u=0.0 ; solution_v=0.0 ; solution_p=0.0 ; solution_f=0.0
    fmatv = 0 ; fmati = 0 ; fmatj = 0
!START CALCULATION
!for main body
    !SOLVE A

    
    
    counter=0
    do i=1,nx 
        do j=1,ny
            dx=x(i+1,j)-x(i,j)
            dy=y(i,j+1)-y(i,j)
            dx=(domain_x_end-domain_x_start)/nx
            dy=(domain_y_end-domain_y_start)/ny
            
            bound = 0.0d0

            temp_con1=0.0d0
            temp_con2=0.0d0
            temp_square=0.0d0
            temp_ft=0.0d0
            temp_uu=0.0d0
            temp_vv=0.0d0

            temp_con1=( u(i,j)+u(i+1,j) )/2.0d0*( v(i,j)+v(i+1,j)+v(i,j+1)+v(i+1,j+1)-v(i-1,j+1)-v(i,j+1)-v(i-1,j)-v(i,j) )/4.0d0/dx
            temp_con2=( v(i,j)+v(i,j+1) )/2.0d0*( u(i+1,j)-u(i,j) )/dx
            
            temp_square=( ( u(i,j)/2.0d0+u(i+1,j)/2.0d0 )**2+( v(i,j)/2.0d0+v(i,j+1)/2.0d0 )**2 )*a(i,j)
            
            temp_ft=scale_factor*( ( f_alter(i+1,j)-f_alter(i-1,j) )/2.0d0/dx-( f(i+1,j)-f(i-1,j) )/2.0d0/dx )/dt
            
            temp_uu=scale_factor*( ( u(i,j)/2.0d0+u(i+1,j)/2.0d0 )**2+( v(i,j)/2.0d0+v(i,j+1)/2.0d0 )**2 )*&
                                 ((b(i+1,j+1)-b(i+1,j-1))/2.0d0/dy-(b(i-1,j+1)-b(i-1,j-1))/2.0d0/dy)/2.0d0/dx
                                 
            temp_vv=scale_factor*( ( u(i,j)/2.0d0+u(i+1,j)/2.0d0 )**2+( v(i,j)/2.0d0+v(i,j+1)/2.0d0 )**2 )*&
                                 (a(i,j-1)+a(i,j+1)-2.0d0*a(i,j))/dy/dy
            !print *, i,j,u(i,j),v(i,j),f(i,j),f_alter(i,j),"UVF"
            !print *, i,j,temp_con1,temp_con2,temp_square,temp_ft,temp_uu,temp_vv
            a_alter(i,j)=(-temp_con1+temp_con2-temp_square+temp_ft-temp_uu+temp_vv)*dt*dt+2.0d0*a(i,j)-a_former(i,j)
            if (p_gs(i,j,5) .eq. "X") then
                a_alter(i,j)=0.0
            end if
            
        end do
    end do
    !SOLVE B
    counter=0
    do i=1,nx 
        do j=1,ny
            dx=x(i+1,j)-x(i,j)
            dy=y(i,j+1)-y(i,j)
            dx=(domain_x_end-domain_x_start)/nx
            dy=(domain_y_end-domain_y_start)/ny
            
            bound = 0.0d0
            
            temp_con1=0.0d0
            temp_con2=0.0d0
            temp_square=0.0d0
            temp_ft=0.0d0
            temp_uu=0.0d0
            temp_vv=0.0d0
            
            temp_con1=( u(i,j)+u(i+1,j) )/2.0d0*( v(i,j+1)-v(i,j) )/dy
            temp_con2=( v(i,j)+v(i,j+1) )/2.0d0*( u(i,j)+u(i,j+1)+u(i+1,j)+u(i+1,j+1)-u(i,j-1)-u(i,j)-u(i+1,j-1)-u(i+1,j) )/4.0d0/dy
            
            temp_square=( ( u(i,j)/2.0d0+u(i+1,j)/2.0d0 )**2+( v(i,j)/2.0d0+v(i,j+1)/2.0d0 )**2 )*b(i,j)
            
            temp_ft=scale_factor*( ( f_alter(i,j+1)-f_alter(i,j-1) )/2.0d0/dy-( f(i,j+1)-f(i,j-1) )/2.0d0/dy )/dt
            
            temp_uu=scale_factor*( ( u(i,j)/2.0d0+u(i+1,j)/2.0d0 )**2+( v(i,j)/2.0d0+v(i,j+1)/2.0d0 )**2 )*&
                                 (b(i-1,j)+b(i+1,j)-2.0d0*b(i,j))/dy/dy
            
            temp_vv=scale_factor*( ( u(i,j)/2.0d0+u(i+1,j)/2.0d0 )**2+( v(i,j)/2.0d0+v(i,j+1)/2.0d0 )**2 )*&
                                 ((a(i+1,j+1)-a(i-1,j+1))/2.0d0/dx-(a(i+1,j-1)-a(i-1,j-1))/2.0d0/dx)/2.0d0/dy


            b_alter(i,j)=(-temp_con1+temp_con2-temp_square+temp_ft+temp_uu-temp_vv)*dt*dt+2.0d0*b(i,j)-b_former(i,j)
            
            if (p_gs(i,j,5) .eq. "X") then
                b_alter(i,j)=0.0
            end if
            
        end do
    end do


    !a_alter(0,:)    = a_bound_w(:)
    !a_alter(nx+1,:) = a_bound_e(:)
    !a_alter(:,0)    = a_bound_s(:)
    !a_alter(:,ny+1) = a_bound_n(:)
    
    !b_alter(0,:)    = b_bound_w(:)
    !b_alter(nx+1,:) = b_bound_e(:)
    !b_alter(:,0)    = b_bound_s(:)
    !b_alter(:,ny+1) = b_bound_n(:)
    
    do i=1,nx
        do j=1,ny
            !print *, i,j,a_alter(i,j),b_alter(i,j),"AB"
        end do
    end do
    
    !a_alter(0,0)=       1.0*(a_alter(0,1)+a_alter(1,0))         -0.5*(a_alter(0,2)+a_alter(2,0))
    !a_alter(0,ny+1)=    1.0*(a_alter(0,ny)+a_alter(1,ny+1))     -0.5*(a_alter(0,ny-1)+a_alter(2,ny+1))
    !a_alter(nx+1,0)=    1.0*(a_alter(nx+1,1)+a_alter(nx,0))     -0.5*(a_alter(nx+1,2)+a_alter(nx-1,0))
    !a_alter(nx+1,ny+1)= 1.0*(a_alter(nx+1,ny)+a_alter(nx,ny+1)) -0.5*(a_alter(nx+1,ny-1)+a_alter(nx-1,ny+1))
    
    !b_alter(0,0)=       1.0*(b_alter(0,1)+b_alter(1,0))         -0.5*(b_alter(0,2)+b_alter(2,0))
    !b_alter(0,ny+1)=    1.0*(b_alter(0,ny)+b_alter(1,ny+1))     -0.5*(b_alter(0,ny-1)+b_alter(2,ny+1))
    !b_alter(nx+1,0)=    1.0*(b_alter(nx+1,1)+b_alter(nx,0))     -0.5*(b_alter(nx+1,2)+b_alter(nx-1,0))
    !b_alter(nx+1,ny+1)= 1.0*(b_alter(nx+1,ny)+b_alter(nx,ny+1)) -0.5*(b_alter(nx+1,ny-1)+b_alter(nx-1,ny+1))
            
            
    print *, "Solve A & B complete"
    
    deallocate (umatv,umati,umatj,umrhs)
    deallocate (vmatv,vmati,vmatj,vmrhs)
    deallocate (pmatv,pmati,pmatj,pmrhs)
    deallocate (fmatv,fmati,fmatj,fmrhs)
    
    deallocate(solution_u)
    deallocate(solution_v)
    deallocate(solution_p)
    deallocate(solution_f)
    
    deallocate (AA_coe,BB_bound,A_coe,B_bound)
end subroutine ab_solver

subroutine uvp_process(nx,ny,nt,ghost_cell_size,&
                            ap_u,ap_v,ap_u1,ap_v1,&
                            u,u_alter,v,v_alter,&
                            p,p_alter,p_correct,&
                            u_gs,v_gs,p_gs,&
                            u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                            u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                            p_bound_w,p_bound_e,&
                            p_bound_s,p_bound_n)
integer :: i,j
integer :: nx,ny,nt,ghost_cell_size



real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct

real(kind=8),dimension(0:ny+2)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1)::p_bound_s,p_bound_n

real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v
real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

character(len=1),dimension(0:nx+2, 0:ny+2, 5) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 5) ::p_gs

    character(len=1),dimension(1,4,3)::uvbc,pbc
    
    !real(kind=8),dimension(ghost_cell_size,0:ny+2)::u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e
    !real(kind=8),dimension(ghost_cell_size,0:nx+2)::u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n
    !real(kind=8),dimension(ghost_cell_size,0:ny+1)::p_ghost_cell_w,p_ghost_cell_e
    !real(kind=8),dimension(ghost_cell_size,0:nx+1)::p_ghost_cell_s,p_ghost_cell_n
    
    
    
    do i=1,1
        do j=2,2
            !print *, "...",i,j,k,v(i,j),v_alter(i,j),"v3"
        END do
    END do 
    
    !correct u and v
    do i=2,nx
        do j=1,ny
            if (u_gs(i,j,5) .ne. "X") then
            u_alter(i,j)=u_alter(i,j) + dy*0.7/ap_u1(i,j)*( p_correct(i-1,j) - p_correct(i,j) )
            !print *,i,j,u_alter(i,j),"us"
            end if
        end do
    end do    
    
    do i=1,nx
        do j=2,ny
            if (v_gs(i,j,5) .ne. "X") then
            v_alter(i,j)=v_alter(i,j) + dx*0.7/ap_v1(i,j)*( p_correct(i,j-1) - p_correct(i,j) )
            !print *,i,j,v_alter(i,j),"vs"
            end if
        end do
    end do


    !INTERPOLATE
    !2p(1,:)-2p(2,:)=p(2,:)-p(3,:)
    
    do i=1,nx
        do j=1,ny
            if (p_gs(i,j,1) .eq. "w") then
                p_alter(i-1,j)=  1.5*p_alter(i,j)-0.5*p_alter(i+1,j)
            end if 
            
            if (p_gs(i,j,2) .eq. "e") then
                p_alter(i+1,j)=  1.5*p_alter(i,j)-0.5*p_alter(i-1,j)
            end if 
            
            if (p_gs(i,j,3) .eq. "s") then
                p_alter(i,j-1)=  1.5*p_alter(i,j)-0.5*p_alter(i,j+1)
            end if 
            
            if (p_gs(i,j,4) .eq. "n") then
                p_alter(i,j+1)=  1.5*p_alter(i,j)-0.5*p_alter(i,j-1)
            end if 
            
            
            if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,3) .eq. "s") then
                p_alter(i-1,j-1) = 1.0*(p_alter(i-1,j)+p_alter(i,j-1)) - 0.5*(p_alter(i-1,j+1)+p_alter(i+1,j-1))
            end if 
            
            if (p_gs(i,j,1) .eq. "w" .and. p_gs(i,j,3) .eq. "n") then
                p_alter(i-1,j+1) = 1.0*(p_alter(i-1,j)+p_alter(i,j+1)) - 0.5*(p_alter(i-1,j-1)+p_alter(i+1,j+1))
            end if 
            
            if (p_gs(i,j,1) .eq. "e" .and. p_gs(i,j,3) .eq. "s") then
                p_alter(i+1,j-1) = 1.0*(p_alter(i+1,j)+p_alter(i,j-1)) - 0.5*(p_alter(i+1,j+1)+p_alter(i-1,j-1))
            end if 
            
            if (p_gs(i,j,1) .eq. "e" .and. p_gs(i,j,3) .eq. "n") then
                p_alter(i+1,j+1) = 1.0*(p_alter(i+1,j)+p_alter(i,j+1)) - 0.5*(p_alter(i+1,j-1)+p_alter(i-1,j+1))
            end if 
            
            
        end do
    end do
    
    
    


    
    u_alter(0,:)    = u_alter(1,:) !W
    u_alter(nx+2,:) = u_alter(nx+1,:) !E
    
    u_alter(:,0)    = u_bound_s(:)
    u_alter(:,ny+1) = u_bound_n(:)
    
    v_alter(:,0)    = v_alter(:,1) !S
    v_alter(:,ny+2) = v_alter(:,ny+1) !N
    
    v_alter(0,:)    = v_bound_w(:)
    v_alter(nx+1,:) = v_bound_e(:)
    

    !Transfer matrices for time iternation
    
    do i=0,nx+2
        do j=0,ny+2
        u(i,j)=u_alter(i,j)
        v(i,j)=v_alter(i,j)
        !print *, i+1,j+1,u_alter(i,j),v_alter(i,j)
        end do
    end do

    do i=0,nx+1
        do j=0,ny+1
        p(i,j)=p_alter(i,j)
        end do
    end do


end subroutine uvp_process




subroutine fab_process(nx,ny,nt,timestep_id,ghost_cell_size,&
                            ap_u,ap_v,ap_u1,ap_v1,&
                            u,u_alter,v,v_alter,&
                            p,p_alter,p_correct,&
                            f_former,f,f_alter,&
                            a_former,a,a_alter,&
                            b_former,b,b_alter,&
                            u_gs,v_gs,p_gs,&
                            u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                            u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                            p_bound_w,p_bound_e,&
                            p_bound_s,p_bound_n,&
                            f_bound_w,f_bound_e,&
                            f_bound_s,f_bound_n,&
                            a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                            a_bound_s,a_bound_n,b_bound_s,b_bound_n)
integer :: i,j
integer :: nx,ny,nt,ghost_cell_size,timestep_id



real(kind=8),dimension(0:nx+2, 0:ny+2)::u,u_alter,v,v_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,p_alter,p_correct

real(kind=8),dimension(0:nx+1, 0:ny+1)::f_former,f,f_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::a_former,a,a_alter
real(kind=8),dimension(0:nx+1, 0:ny+1)::b_former,b,b_alter

real(kind=8),dimension(0:ny+2,nt)::u_bound_w,u_bound_e,v_bound_w,v_bound_e
real(kind=8),dimension(0:nx+2,nt)::u_bound_s,u_bound_n,v_bound_s,v_bound_n
real(kind=8),dimension(0:ny+1,nt)::p_bound_w,p_bound_e
real(kind=8),dimension(0:nx+1,nt)::p_bound_s,p_bound_n

real(kind=8),dimension(0:ny+1,0:nt)::f_bound_w,f_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::f_bound_s,f_bound_n
real(kind=8),dimension(0:ny+1,0:nt)::a_bound_w,a_bound_e,b_bound_w,b_bound_e
real(kind=8),dimension(0:nx+1,0:nt)::a_bound_s,a_bound_n,b_bound_s,b_bound_n


real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v
real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1

character(len=1),dimension(0:nx+2, 0:ny+2, 5) ::u_gs,v_gs
character(len=1),dimension(0:nx+1, 0:ny+1, 5) ::p_gs


    f_former(0,:)    = f_bound_w(:,timestep_id-1)
    f_former(nx+1,:) = f_bound_e(:,timestep_id-1)
    
    f_former(:,0)    = f_bound_s(:,timestep_id-1)
    f_former(:,ny+1) = f_bound_n(:,timestep_id-1)
    
    f(0,:)    = f_bound_w(:,timestep_id)
    f(nx+1,:) = f_bound_e(:,timestep_id)
    
    f(:,0)    = f_bound_s(:,timestep_id)
    f(:,ny+1) = f_bound_n(:,timestep_id)

    f_alter(0,:)    = f_bound_w(:,timestep_id)
    f_alter(nx+1,:) = f_bound_e(:,timestep_id)
    
    f_alter(:,0)    = f_bound_s(:,timestep_id)
    f_alter(:,ny+1) = f_bound_n(:,timestep_id)
    
    
    

    a_former(0,:)    = a_bound_w(:,timestep_id-1)
    a_former(nx+1,:) = a_bound_e(:,timestep_id-1)
    
    a_former(:,0)    = a_bound_s(:,timestep_id-1)
    a_former(:,ny+1) = a_bound_n(:,timestep_id-1)
    
    a(0,:)    = a_bound_w(:,timestep_id)
    a(nx+1,:) = a_bound_e(:,timestep_id)
    
    a(:,0)    = a_bound_s(:,timestep_id)
    a(:,ny+1) = a_bound_n(:,timestep_id)

    a_alter(0,:)    = a_bound_w(:,timestep_id)
    a_alter(nx+1,:) = a_bound_e(:,timestep_id)
    
    a_alter(:,0)    = a_bound_s(:,timestep_id)
    a_alter(:,ny+1) = a_bound_n(:,timestep_id)
    
    
    
    
    
    b_former(0,:)    = b_bound_w(:,timestep_id-1)
    b_former(nx+1,:) = b_bound_e(:,timestep_id-1)
    
    b_former(:,0)    = b_bound_s(:,timestep_id-1)
    b_former(:,ny+1) = b_bound_n(:,timestep_id-1)
    
    b(0,:)    = b_bound_w(:,timestep_id)
    b(nx+1,:) = b_bound_e(:,timestep_id)

    b(:,0)    = b_bound_s(:,timestep_id)
    b(:,ny+1) = b_bound_n(:,timestep_id)

    b_alter(0,:)    = b_bound_w(:,timestep_id)
    b_alter(nx+1,:) = b_bound_e(:,timestep_id)
    
    b_alter(:,0)    = b_bound_s(:,timestep_id)
    b_alter(:,ny+1) = b_bound_n(:,timestep_id)
    
    
    

    !Transfer matrices for time iternation


    do i=0,nx+1
        do j=0,ny+1
            f_former(i,j)=f(i,j)
            f(i,j)=f_alter(i,j)
            a_former(i,j)=a(i,j)
            a(i,j)=a_alter(i,j)
            b_former(i,j)=b(i,j)
            b(i,j)=b_alter(i,j)
        end do
    end do


end subroutine fab_process
!this subroutine storage or grab field data of block_id to all 
subroutine data_collect(data_collect_status,data_collect_property,block_id,nb,nx,ny,&
                              block_grid_count1,block_grid_count2,&
                              block_grid_id1,block_grid_id2,&
                              x,y,u,v,p,f,a_former,b_former,a,b,&
                              x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                              ap_u,ap_v,ap_u1,ap_v1,&
                              ap_u_block,ap_v_block)
                              
character(len=7)::data_collect_status
character(len=3)::data_collect_property
integer ::i,j
integer ::nb,nx,ny
integer ::block_id
integer ::block_grid_count1,block_grid_count2

integer,dimension(1:nb)::block_grid_id1,block_grid_id2

real(kind=8),dimension(0:nx+1, 0:ny+1)::x,y
real(kind=8),dimension(0:nx+2, 0:ny+2)::u,v
real(kind=8),dimension(0:nx+1, 0:ny+1)::p,f,a_former,b_former,a,b
real(kind=8),dimension(2:nx,1:ny)::ap_u
real(kind=8),dimension(1:nx,2:ny)::ap_v
real(kind=8),dimension(1:nx+1,1:ny)::ap_u1
real(kind=8),dimension(1:nx,1:ny+1)::ap_v1
real(kind=8),dimension(0:nx+1, 0:ny+1)::ap_u_temp,ap_v_temp
real(kind=8),dimension(block_grid_count1)::x_block,y_block,p_block,f_block,a_block,b_block,af_block,bf_block
real(kind=8),dimension(block_grid_count2)::u_block,v_block
real(kind=8),dimension(block_grid_count1)::ap_u_block,ap_v_block

ap_u_temp=0.0d0
ap_v_temp=0.0d0



    !Put all block data to storage
if (data_collect_status .eq. "Storage") then
    
    do i=1,nx+2
        do j=1,ny+2
             if (data_collect_property .eq. "X&Y") then
             x_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=x(i-1,j-1)!i+nx*(j-1)
             y_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=y(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "UVP") then
             p_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=p(i-1,j-1)!i+nx*(j-1)
             f_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=f(i-1,j-1)!i+nx*(j-1)
             
             af_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=a_former(i-1,j-1)!i+nx*(j-1)
             bf_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=b_former(i-1,j-1)!i+nx*(j-1)
             
             a_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=a(i-1,j-1)!i+nx*(j-1)
             b_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=b(i-1,j-1)!i+nx*(j-1)
             
             end if
             if (data_collect_property .eq. "-P-") then
             p_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=p(i-1,j-1)!i+nx*(j-1)
             end if
             
             if (data_collect_property .eq. "-F-") then
             f_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=f(i-1,j-1)!i+nx*(j-1)
             end if
             
             if (data_collect_property .eq. "A&B") then
             af_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=a_former(i-1,j-1)!i+nx*(j-1)
             bf_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=b_former(i-1,j-1)!i+nx*(j-1)
             
             a_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=a(i-1,j-1)!i+nx*(j-1)
             b_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=b(i-1,j-1)!i+nx*(j-1)
             end if
             
             if (data_collect_property .eq. "ALL") then
             x_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=x(i-1,j-1)!i+nx*(j-1)
             y_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=y(i-1,j-1)!i+nx*(j-1)
             p_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=p(i-1,j-1)!i+nx*(j-1)
             f_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=f(i-1,j-1)!i+nx*(j-1)
             a_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=a(i-1,j-1)!i+nx*(j-1)
             b_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=b(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "APU") then
             ap_u_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=ap_u_temp(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "APV") then
             ap_v_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=ap_v_temp(i-1,j-1)!i+nx*(j-1)
             end if
        end do
    end do
    
    do i=1,nx+3
        do j=1,ny+3
             if (data_collect_property .eq. "UVP") then
             u_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )=u(i-1,j-1)!i+nx*(j-1)
             v_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )=v(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "-U-") then
             u_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )=u(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "-V-") then
             v_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )=v(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "ALL") then
             u_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )=u(i-1,j-1)!i+nx*(j-1)
             v_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )=v(i-1,j-1)!i+nx*(j-1)
             end if
        end do
    end do
    
    
    do i=1,nx+1
        do j=1,ny
        if (data_collect_property .eq. "APU") then
        ap_u_temp(i,j)=ap_u1(i,j)
        end if
        end do
    end do

    do i=1,nx
        do j=1,ny+1
        if (data_collect_property .eq. "APV") then
        ap_v_temp(i,j)=ap_v1(i,j)
        end if
        end do
    end do
    
    do i=1,nx+2
        do j=1,ny+2
             if (data_collect_property .eq. "APU") then
             ap_u_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=ap_u_temp(i-1,j-1)!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "APV") then
             ap_v_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )=ap_v_temp(i-1,j-1)!i+nx*(j-1)
             end if
        end do
    end do

    
else if (data_collect_status .eq. "GrabOut") then

    do i=1,nx+2
        do j=1,ny+2
             if (data_collect_property .eq. "X&Y") then
             x(i-1,j-1)=x_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             y(i-1,j-1)=y_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "UVP") then
             p(i-1,j-1)=p_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             f(i-1,j-1)=f_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             
             a_former(i-1,j-1)=af_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             b_former(i-1,j-1)=bf_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             
             a(i-1,j-1)=a_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             b(i-1,j-1)=b_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "-P-") then
             p(i-1,j-1)=p_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
             
             if (data_collect_property .eq. "-F-") then
             f(i-1,j-1)=f_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
             
             if (data_collect_property .eq. "A&B") then
             a_former(i-1,j-1)=af_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             b_former(i-1,j-1)=bf_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             
             a(i-1,j-1)=a_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             b(i-1,j-1)=b_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
             
             
             if (data_collect_property .eq. "ALL") then
             x(i-1,j-1)=x_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             y(i-1,j-1)=y_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             p(i-1,j-1)=p_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             f(i-1,j-1)=f_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             a_former(i-1,j-1)=af_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             b_former(i-1,j-1)=bf_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             a(i-1,j-1)=a_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             b(i-1,j-1)=b_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
        end do
    end do
    
    do i=1,nx+3
        do j=1,ny+3
             if (data_collect_property .eq. "UVP") then
             u(i-1,j-1)=u_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )!i+nx*(j-1)
             v(i-1,j-1)=v_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "-U-") then
             u(i-1,j-1)=u_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "-V-") then
             v(i-1,j-1)=v_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "ALL") then
             u(i-1,j-1)=u_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )!i+nx*(j-1)
             v(i-1,j-1)=v_block( j+(ny+3)*(i-1)+block_grid_id2(block_id) )!i+nx*(j-1)
             end if
        end do
    end do
    
    
    do i=1,nx+2
        do j=1,ny+2
             if (data_collect_property .eq. "APU") then
             ap_u_temp(i-1,j-1)=ap_u_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
             if (data_collect_property .eq. "APV") then
             ap_v_temp(i-1,j-1)=ap_v_block( j+(ny+2)*(i-1)+block_grid_id1(block_id) )!i+nx*(j-1)
             end if
        end do
    end do
    
    do i=1,nx+1
        do j=1,ny
        if (data_collect_property .eq. "APU") then
        ap_u1(i,j)=ap_u_temp(i,j)
        end if
        end do
    end do

    do i=1,nx
        do j=1,ny+1
        if (data_collect_property .eq. "APV") then
        ap_v1(i,j)=ap_v_temp(i,j)
        end if
        end do
    end do
else
    print *, "status unclear, no action will be executed!"
end if

end subroutine data_collect



program New_2D
    
character(len=10) ::print_iteration
character(len=100) ::output_filename

integer :: i,j,k,number_block,t,counter
integer :: ghost_cell_size

integer ::nx,ny,nt,nb
real(kind=8)::dt,dx,dy
real(kind=8)::density=1000.0d0
real(kind=8)::viscosity=0.01d0
real(kind=8)::scale_factor=1.0d0
integer ::dump_every,skip_new_theory

integer,allocatable,dimension(:):: gx,gy
real(kind=8),allocatable,dimension(:)::gx_start,gx_end,gy_start,gy_end
!real(kind=8),allocatable,dimension(:,:,:)::gu,gv,gp
character(len=1),allocatable,dimension(:,:,:) ::uvbc,pbc
character(len=1),allocatable,dimension(:,:,:) ::winds
character(len=2),allocatable,dimension(:,:) ::linkbc


real(kind=8)::domain_x_start,domain_x_end,domain_y_start, domain_y_end

character(len=1)::ubc_e,ubc_w,vbc_n,vbc_s
character(len=1)::pbc_e,pbc_w,pbc_n,pbc_s

real(kind=8),allocatable,dimension(:)::u_w,u_e,u_s,u_n
real(kind=8),allocatable,dimension(:)::v_w,v_e,v_s,v_n
real(kind=8),allocatable,dimension(:)::p_w,p_e,p_s,p_n

!NEED TO ALLOCATE AFTER READING nx ny
character(len=1),allocatable,dimension(:,:,:) ::u_gs,v_gs,p_gs

real(kind=8),allocatable,dimension(:,:)::x,y
real(kind=8),allocatable,dimension(:,:)::u,u_alter,v,v_alter,p,p_alter,p_correct

real(kind=8),allocatable,dimension(:,:)::u_bound_w,u_bound_e,u_bound_n,u_bound_s
real(kind=8),allocatable,dimension(:,:)::v_bound_w,v_bound_e,v_bound_n,v_bound_s
real(kind=8),allocatable,dimension(:,:)::p_bound_w,p_bound_e,p_bound_n,p_bound_s

real(kind=8),allocatable,dimension(:,:)::ap_u,ap_v
real(kind=8),allocatable,dimension(:,:)::ap_u1,ap_v1

real(kind=8),allocatable,dimension(:,:)::f_former,f,f_alter
real(kind=8),allocatable,dimension(:,:)::f_bound_w,f_bound_e,f_bound_n,f_bound_s

real(kind=8),allocatable,dimension(:,:)::a_former,a,a_alter,b_former,b,b_alter
real(kind=8),allocatable,dimension(:,:)::a_bound_w,a_bound_e,a_bound_n,a_bound_s
real(kind=8),allocatable,dimension(:,:)::b_bound_w,b_bound_e,b_bound_n,b_bound_s
!NEED TO ALLOCATE AFTER READING nx ny

!combined_block
real(kind=8),allocatable,dimension(:)::x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block

real(kind=8),allocatable,dimension(:)::ap_u_block,ap_v_block

real(kind=8),allocatable,dimension(:)::x0_block,y0_block,u0_block,v0_block,p0_block,f0_block,af0_block,bf0_block,a0_block,b0_block!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


integer::block_grid_count1=0
integer,allocatable,dimension(:)::block_grid_id1
integer::block_grid_count2=0
integer,allocatable,dimension(:)::block_grid_id2
!
real(kind=8),allocatable,dimension(:,:)::u_ghost_cell_w,u_ghost_cell_e,u_ghost_cell_s,u_ghost_cell_n
real(kind=8),allocatable,dimension(:,:)::v_ghost_cell_w,v_ghost_cell_e,v_ghost_cell_s,v_ghost_cell_n
real(kind=8),allocatable,dimension(:,:)::p_ghost_cell_w,p_ghost_cell_e,p_ghost_cell_s,p_ghost_cell_n



open (1,file='control2.txt',status='old')
open (2,file='block2.txt',status='old')
open (88,file='residual2.txt',status='replace')
open (99,file='coefficient2.txt',status='replace')


read (1,*)
read (1,*) nt
read (1,*) 
read (1,*) dt
read (1,*) 
read (1,*) dump_every
read (1,*) 
read (1,*) density,viscosity
read (1,*) 
read (1,*) scale_factor
read (1,*)
read (1,*) skip_new_theory
read (1,*)


write (88,*) "VARIABLES= T,R"
write (88,*) "ZONE I=",nt,",F=POINT"


read (2,*)
read (2,*) nb, ghost_cell_size

allocate (gx(1:nb),gy(1:nb))
allocate (gx_start(1:nb),gx_end(1:nb))
allocate (gy_start(1:nb),gy_end(1:nb))
allocate (uvbc(nb,4,3),pbc(nb,4,3))
allocate (linkbc(nb,4))
allocate (winds(nb,4,2))

allocate (u_w(1:nb),u_e(1:nb),u_s(1:nb),u_n(1:nb))
allocate (v_w(1:nb),v_e(1:nb),v_s(1:nb),v_n(1:nb))
allocate (p_w(1:nb),p_e(1:nb),p_s(1:nb),p_n(1:nb))

!allocate (gu(nb,nt,4),gv(nb,nt,4),gp(nb,nt,4))
do k=1,nb
read (2,*)
read (2,*) gx(k),gy(k)
read (2,*)
read (2,*) gx_start(k),gx_end(k)
read (2,*)
read (2,*) gy_start(k),gy_end(k)
read (2,*)
read (2,*) uvbc(k,1,1),uvbc(k,1,2),uvbc(k,1,3),linkbc(k,1),winds(k,1,1),u_w(k),v_w(k)
read (2,*) uvbc(k,2,1),uvbc(k,2,2),uvbc(k,2,3),linkbc(k,2),winds(k,2,1),u_e(k),v_e(k)
read (2,*) uvbc(k,3,1),uvbc(k,3,2),uvbc(k,3,3),linkbc(k,3),winds(k,3,1),u_s(k),v_s(k)
read (2,*) uvbc(k,4,1),uvbc(k,4,2),uvbc(k,4,3),linkbc(k,4),winds(k,4,1),u_n(k),v_n(k)
read (2,*)
read (2,*) pbc(k,1,1),pbc(k,1,2),pbc(k,1,3),linkbc(k,1),winds(k,1,2),p_w(k)
read (2,*) pbc(k,2,1),pbc(k,2,2),pbc(k,2,3),linkbc(k,2),winds(k,2,2),p_e(k)
read (2,*) pbc(k,3,1),pbc(k,3,2),pbc(k,3,3),linkbc(k,3),winds(k,3,2),p_s(k)
read (2,*) pbc(k,4,1),pbc(k,4,2),pbc(k,4,3),linkbc(k,4),winds(k,4,2),p_n(k)
end do

allocate (block_grid_id1(1:nb))
allocate (block_grid_id2(1:nb))
block_grid_count1 = 0
block_grid_count2 = 0
do k=1,nb
!ONE BLOCK SIZE LENGTH
    block_grid_id1(k)=block_grid_count1
    block_grid_id2(k)=block_grid_count2
    print*,"Block grid P id starts at ",block_grid_id1(k),"for block",k
    print*,"Block grid UV id starts at ",block_grid_id2(k),"for block",k
    !TOTAL BLOCK GRID SIZE
    block_grid_count1=block_grid_count1+(gx(k)+2)*(gy(k)+2)
    block_grid_count2=block_grid_count2+(gx(k)+3)*(gy(k)+3)
end do

print*,"Total Block grid",block_grid_count

allocate(u_block(block_grid_count2),v_block(block_grid_count2))
allocate(x_block(block_grid_count1),y_block(block_grid_count1))
allocate(p_block(block_grid_count1),f_block(block_grid_count1))
allocate(a_block(block_grid_count1),b_block(block_grid_count1))
allocate(af_block(block_grid_count1),bf_block(block_grid_count1))

allocate(ap_u_block(block_grid_count1),ap_v_block(block_grid_count1))

x_block=0.0d0
y_block=0.0d0
u_block=0.0d0
v_block=0.0d0
p_block=0.0d0
f_block=0.0d0
a_block=0.0d0
b_block=0.0d0
af_block=0.0d0
bf_block=0.0d0

ap_u_block=0.0d0
ap_v_block=0.0d0

allocate(u0_block(block_grid_count2),v0_block(block_grid_count2))
allocate(x0_block(block_grid_count1),y0_block(block_grid_count1))
allocate(p0_block(block_grid_count1),f0_block(block_grid_count1))
allocate(a0_block(block_grid_count1),b0_block(block_grid_count1))
allocate(af0_block(block_grid_count1),bf0_block(block_grid_count1))

x0_block=0.0d0
y0_block=0.0d0
u0_block=0.0d0
v0_block=0.0d0
p0_block=0.0d0
f0_block=0.0d0
a0_block=0.0d0
b0_block=0.0d0
af0_block=0.0d0
bf0_block=0.0d0

print *, "Start!"







do k=1,nb
    nx=gx(k)
    ny=gy(k)
    domain_x_start = gx_start(k)
    domain_x_end   = gx_end(k)
    domain_y_start = gy_start(k)
    domain_y_end   = gy_end(k)
    allocate(x(0:nx+1,0:ny+1)) ; allocate(y(0:nx+1,0:ny+1))

    allocate(      u_gs(0:nx+2, 0:ny+2, 6)) ; allocate(      v_gs(0:nx+2, 0:ny+2, 6)) ; allocate(      p_gs(0:nx+1, 0:ny+1, 6))!0, 1,..,nx,nx+1

    allocate(      u(0:nx+2,  0:ny+2)) ; allocate(u_alter(0:nx+2,  0:ny+2))
    allocate(      v(0:nx+2,  0:ny+2)) ; allocate(v_alter(0:nx+2,  0:ny+2))
    allocate(      p(0:nx+1,0:ny+1)  ) ; allocate(p_alter(0:nx+1,0:ny+1)) ; allocate(p_correct(0:nx+1,0:ny+1))

    allocate( f_former(0:nx+1,0:ny+1)) ; allocate(        f(0:nx+1,0:ny+1)) ; allocate(  f_alter(0:nx+1,0:ny+1))
    allocate( a_former(0:nx+1,0:ny+1)) ; allocate(        a(0:nx+1,0:ny+1)) ; allocate(  a_alter(0:nx+1,0:ny+1))
    allocate( b_former(0:nx+1,0:ny+1)) ; allocate(        b(0:nx+1,0:ny+1)) ; allocate(  b_alter(0:nx+1,0:ny+1))

allocate(u_bound_e(0:ny+2,nt)) ; allocate(u_bound_w(0:ny+2,nt)) ; allocate(u_bound_n(0:nx+2,nt)) ; allocate(u_bound_s(0:nx+2,nt))
allocate(v_bound_e(0:ny+2,nt)) ; allocate(v_bound_w(0:ny+2,nt)) ; allocate(v_bound_n(0:nx+2,nt)) ; allocate(v_bound_s(0:nx+2,nt))
allocate(p_bound_e(0:ny+1,nt)) ; allocate(p_bound_w(0:ny+1,nt)) ; allocate(p_bound_n(0:nx+1,nt)) ; allocate(p_bound_s(0:nx+1,nt))

allocate(f_bound_e(0:ny+1,nt)) ; allocate(f_bound_w(0:ny+1,nt)) ; allocate(f_bound_n(0:nx+1,nt)) ; allocate(f_bound_s(0:nx+1,nt))
allocate(a_bound_e(0:ny+1,nt)) ; allocate(a_bound_w(0:ny+1,nt)) ; allocate(a_bound_n(0:nx+1,nt)) ; allocate(a_bound_s(0:nx+1,nt))
allocate(b_bound_e(0:ny+1,nt)) ; allocate(b_bound_w(0:ny+1,nt)) ; allocate(b_bound_n(0:nx+1,nt)) ; allocate(b_bound_s(0:nx+1,nt))

    allocate(ap_u(2:nx,1:ny)) ; allocate(ap_v(1:nx,2:ny)) ; allocate(ap_u1(1:nx+1,1:ny)) ; allocate(ap_v1(1:nx,1:ny+1))

    !print *, "Allocated!"
    !allocate ghost cell
    allocate(u_ghost_cell_e(ghost_cell_size,0:ny+2))
    allocate(u_ghost_cell_w(ghost_cell_size,0:ny+2))
    allocate(u_ghost_cell_n(ghost_cell_size,0:nx+2))
    allocate(u_ghost_cell_s(ghost_cell_size,0:nx+2))

    allocate(v_ghost_cell_e(ghost_cell_size,0:ny+2))
    allocate(v_ghost_cell_w(ghost_cell_size,0:ny+2))
    allocate(v_ghost_cell_n(ghost_cell_size,0:nx+2))
    allocate(v_ghost_cell_s(ghost_cell_size,0:nx+2))

    allocate(p_ghost_cell_e(ghost_cell_size,0:ny+1))
    allocate(p_ghost_cell_w(ghost_cell_size,0:ny+1))
    allocate(p_ghost_cell_n(ghost_cell_size,0:nx+1))
    allocate(p_ghost_cell_s(ghost_cell_size,0:nx+1))
    

!ghost cell allocated

    call initial_array(nx,ny,nt,1,&
                         x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                         u,u_alter,v,v_alter,&
                         p,p_alter,p_correct,&
                         f_former,f,f_alter,&
                         a_former,a,a_alter,&
                         b_former,b,b_alter,&
                         u_gs,v_gs,p_gs,&
                         uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                         u_w(k),u_e(k),u_s(k),u_n(k),&
                         v_w(k),v_e(k),v_s(k),v_n(k),&
                         p_w(k),p_e(k),p_s(k),p_n(k),&
                         ubc_e,ubc_w,&
                         vbc_n,vbc_s,&
                         pbc_e,pbc_w,pbc_n,pbc_s,&
                         u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                         u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                         p_bound_w,p_bound_e,&
                         p_bound_s,p_bound_n,&
                         f_bound_w,f_bound_e,&
                         f_bound_s,f_bound_n,&
                         a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                         a_bound_s,a_bound_n,b_bound_s,b_bound_n)
    !DIFFERENT
    !Storage all to block
    call data_collect("Storage","ALL",k,nb,nx,ny,&
                         block_grid_count1,block_grid_count2,&
                         block_grid_id1,block_grid_id2,&
                         x,y,u,v,p,f,a_former,b_former,a,b,&
                         x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                         ap_u,ap_v,ap_u1,ap_v1,&
                         ap_u_block,ap_v_block)
                         
    call data_collect("Storage","ALL",k,nb,nx,ny,&
                         block_grid_count1,block_grid_count2,&
                         block_grid_id1,block_grid_id2,&
                         x,y,u,v,p,f,a_former,b_former,a,b,&
                         x0_block,y0_block,u0_block,v0_block,p0_block,f0_block,af0_block,bf0_block,a0_block,b0_block,&
                         ap_u,ap_v,ap_u1,ap_v1,&
                         ap_u_block,ap_v_block)
                         
    
    PRINT *,"Initialize all arrays!"   
    
    
    
    !PRINT *,"Deallocated!"
    
    deallocate(ap_u,ap_v,ap_u1,ap_v1)
    deallocate(u_bound_w,u_bound_e,u_bound_n,u_bound_s) 
    deallocate(v_bound_w,v_bound_e,v_bound_n,v_bound_s)
    deallocate(p_bound_w,p_bound_e,p_bound_n,p_bound_s)
    deallocate(f_bound_w,f_bound_e,f_bound_n,f_bound_s)
    deallocate(a_bound_w,a_bound_e,a_bound_n,a_bound_s)
    deallocate(b_bound_w,b_bound_e,b_bound_n,b_bound_s)
    deallocate(f_former,f,f_alter,a_former,a,a_alter,b_former,b,b_alter)
    deallocate(u_alter,v_alter,p_alter,p_correct)
    deallocate(u,v,p)
    deallocate(u_gs,v_gs,p_gs)
    deallocate(x,y)
    deallocate(u_ghost_cell_w,u_ghost_cell_e,u_ghost_cell_n,u_ghost_cell_s) 
    deallocate(v_ghost_cell_w,v_ghost_cell_e,v_ghost_cell_n,v_ghost_cell_s)
    deallocate(p_ghost_cell_w,p_ghost_cell_e,p_ghost_cell_n,p_ghost_cell_s)

end do

print *,"Start time advance!"










do k=1,nb
!ghost cell allocated
end do



        k=1
        nx=gx(k)
        ny=gy(k)

        domain_x_start = gx_start(k)
        domain_x_end   = gx_end(k)
        domain_y_start = gy_start(k)
        domain_y_end   = gy_end(k)
    
        allocate(x(0:nx+1,0:ny+1)) ; allocate(y(0:nx+1,0:ny+1))
        allocate(     u_gs(0:nx+2, 0:ny+2, 6)) ; allocate( v_gs(0:nx+2, 0:ny+2, 6)) ; allocate(   p_gs(0:nx+1,0:ny+1, 6))!0, 1,..,nx,nx+1
        allocate(        u(0:nx+2, 0:ny+2) ) ; allocate(u_alter(0:nx+2, 0:ny+2))
        allocate(        v(0:nx+2, 0:ny+2) ) ; allocate(v_alter(0:nx+2, 0:ny+2))
        allocate(        p(0:nx+1, 0:ny+1) ) ; allocate(p_alter(0:nx+1, 0:ny+1)) ; allocate( p_correct(0:nx+1,0:ny+1))
        allocate( f_former(0:nx+1,0:ny+1)) ; allocate(        f(0:nx+1, 0:ny+1)) ; allocate( f_alter(0:nx+1,0:ny+1))
        allocate( a_former(0:nx+1,0:ny+1)) ; allocate(        a(0:nx+1, 0:ny+1)) ; allocate( a_alter(0:nx+1,0:ny+1))
        allocate( b_former(0:nx+1,0:ny+1)) ; allocate(        b(0:nx+1, 0:ny+1)) ; allocate( b_alter(0:nx+1,0:ny+1))
allocate(u_bound_e(0:ny+2,nt)) ; allocate(u_bound_w(0:ny+2,nt)) ; allocate(u_bound_n(0:nx+2,nt)) ; allocate(u_bound_s(0:nx+2,nt))
allocate(v_bound_e(0:ny+2,nt)) ; allocate(v_bound_w(0:ny+2,nt)) ; allocate(v_bound_n(0:nx+2,nt)) ; allocate(v_bound_s(0:nx+2,nt))
allocate(p_bound_e(0:ny+1,nt)) ; allocate(p_bound_w(0:ny+1,nt)) ; allocate(p_bound_n(0:nx+1,nt)) ; allocate(p_bound_s(0:nx+1,nt))
allocate(f_bound_e(0:ny+1,nt)) ; allocate(f_bound_w(0:ny+1,nt)) ; allocate(f_bound_n(0:nx+1,nt)) ; allocate(f_bound_s(0:nx+1,nt))
allocate(a_bound_e(0:ny+1,nt)) ; allocate(a_bound_w(0:ny+1,nt)) ; allocate(a_bound_n(0:nx+1,nt)) ; allocate(a_bound_s(0:nx+1,nt))
allocate(b_bound_e(0:ny+1,nt)) ; allocate(b_bound_w(0:ny+1,nt)) ; allocate(b_bound_n(0:nx+1,nt)) ; allocate(b_bound_s(0:nx+1,nt))
        allocate(ap_u(2:nx,1:ny)) ; allocate(ap_v(1:nx,2:ny)) ; allocate(ap_u1(1:nx+1,1:ny)) ; allocate(ap_v1(1:nx,1:ny+1))


!allocate ghost cell
        allocate(u_ghost_cell_e(ghost_cell_size,0:ny+2))
        allocate(u_ghost_cell_w(ghost_cell_size,0:ny+2))
        allocate(u_ghost_cell_n(ghost_cell_size,0:nx+2))
        allocate(u_ghost_cell_s(ghost_cell_size,0:nx+2))

        allocate(v_ghost_cell_e(ghost_cell_size,0:ny+2))
        allocate(v_ghost_cell_w(ghost_cell_size,0:ny+2))
        allocate(v_ghost_cell_n(ghost_cell_size,0:nx+2))
        allocate(v_ghost_cell_s(ghost_cell_size,0:nx+2))

        allocate(p_ghost_cell_e(ghost_cell_size,0:ny+1))
        allocate(p_ghost_cell_w(ghost_cell_size,0:ny+1))
        allocate(p_ghost_cell_n(ghost_cell_size,0:nx+1))
        allocate(p_ghost_cell_s(ghost_cell_size,0:nx+1))


    call initial_array(nx,ny,nt,1,&
                         x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                         u,u_alter,v,v_alter,&
                         p,p_alter,p_correct,&
                         f_former,f,f_alter,&
                         a_former,a,a_alter,&
                         b_former,b,b_alter,&
                         u_gs,v_gs,p_gs,&
                         uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                         u_w(k),u_e(k),u_s(k),u_n(k),&
                         v_w(k),v_e(k),v_s(k),v_n(k),&
                         p_w(k),p_e(k),p_s(k),p_n(k),&
                         ubc_e,ubc_w,&
                         vbc_n,vbc_s,&
                         pbc_e,pbc_w,pbc_n,pbc_s,&
                         u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                         u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                         p_bound_w,p_bound_e,&
                         p_bound_s,p_bound_n,&
                         f_bound_w,f_bound_e,&
                         f_bound_s,f_bound_n,&
                         a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                         a_bound_s,a_bound_n,b_bound_s,b_bound_n)
    !DIFFERENT
    !Storage all to block
    call data_collect("Storage","ALL",k,nb,nx,ny,&
                         block_grid_count1,block_grid_count2,&
                         block_grid_id1,block_grid_id2,&
                         x,y,u,v,p,f,a_former,b_former,a,b,&
                         x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                         ap_u,ap_v,ap_u1,ap_v1,&
                         ap_u_block,ap_v_block)
                         
    call data_collect("Storage","ALL",k,nb,nx,ny,&
                         block_grid_count1,block_grid_count2,&
                         block_grid_id1,block_grid_id2,&
                         x,y,u,v,p,f,a_former,b_former,a,b,&
                         x0_block,y0_block,u0_block,v0_block,p0_block,f0_block,af0_block,bf0_block,a0_block,b0_block,&
                         ap_u,ap_v,ap_u1,ap_v1,&
                         ap_u_block,ap_v_block)
                         
    
    PRINT *,"Initialize all arrays!"   



!time advance
!INTERNAL LOOP
do t=1,nt
    
    

                
    
!extract data from ghost cell
        call data_collect("GrabOut","ALL",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        call data_collect("GrabOut","-U-",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        print *, "Extracted data from all! Now solve U"
                         
        call initial_boundary(nx,ny,nt,nb,k,t,ghost_cell_size,&
                              x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                              u,u_alter,v,v_alter,&
                              p,p_alter,p_correct,&
                              f_former,f,f_alter,&
                              a_former,a,a_alter,&
                              b_former,b,b_alter,&
                              u_gs,v_gs,p_gs,&
                              uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                              u_w,u_e,u_s,u_n,&
                              v_w,v_e,v_s,v_n,&
                              p_w,p_e,p_s,p_n,&
                              ubc_e,ubc_w,&
                              vbc_n,vbc_s,&
                              pbc_e,pbc_w,pbc_n,pbc_s,&
                              u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                              u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                              p_bound_w,p_bound_e,&
                              p_bound_s,p_bound_n,&
                              f_bound_w,f_bound_e,&
                              f_bound_s,f_bound_n,&
                              a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                              a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                              u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                              u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                              p_ghost_cell_w,p_ghost_cell_e,&
                              p_ghost_cell_s,p_ghost_cell_n)

        call u_solver(nx,ny,nt,t,density,viscosity,scale_factor,ghost_cell_size,&
                      x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                      u_gs,v_gs,p_gs,&
                      u,u_alter,v,v_alter,&
                      p,p_alter,p_correct,&
                      f_former,f,f_alter,&
                      a_former,a,a_alter,&
                      b_former,b,b_alter,&
                      ap_u,ap_v,ap_u1,ap_v1,&
                      uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                      ubc_e,ubc_w,&
                      vbc_n,vbc_s,&
                      pbc_e,pbc_w,pbc_n,pbc_s,&
                      u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                      u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                      p_bound_w,p_bound_e,&
                      p_bound_s,p_bound_n,&
                      f_bound_w,f_bound_e,&
                      f_bound_s,f_bound_n,&
                      a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                      a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                      u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                      u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                      p_ghost_cell_w,p_ghost_cell_e,&
                      p_ghost_cell_s,p_ghost_cell_n)
                    
        print *, "Finish U iternation for block",k,"timestep",t
!STORE DATA OF ALL BLOCKS TO STORGE
        call data_collect("Storage","APU",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p_alter,f_alter,a,b,a_alter,b_alter,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
        print *, "Finish U iternation for block",k,"timestep",t
        call data_collect("Storage","-U-",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
         
        print *, "Finish U iternation for block",k,"timestep",t
    
x0_block = x_block ; y0_block = y_block
u0_block = u_block ; v0_block = v_block
p0_block = p_block ; f0_block = f_block
a0_block = a_block ; b0_block = b_block
af0_block = af_block ; bf0_block = bf_block

!extract data from ghost cell
        call data_collect("GrabOut","ALL",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        call data_collect("GrabOut","-V-",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u,v_alter,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)

        print *, "Extracted data from all! Now solve V"

        call initial_boundary(nx,ny,nt,nb,k,t,ghost_cell_size,&
                              x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                              u,u_alter,v,v_alter,&
                              p,p_alter,p_correct,&
                              f_former,f,f_alter,&
                              a_former,a,a_alter,&
                              b_former,b,b_alter,&
                              u_gs,v_gs,p_gs,&
                              uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                              u_w,u_e,u_s,u_n,&
                              v_w,v_e,v_s,v_n,&
                              p_w,p_e,p_s,p_n,&
                              ubc_e,ubc_w,&
                              vbc_n,vbc_s,&
                              pbc_e,pbc_w,pbc_n,pbc_s,&
                              u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                              u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                              p_bound_w,p_bound_e,&
                              p_bound_s,p_bound_n,&
                              f_bound_w,f_bound_e,&
                              f_bound_s,f_bound_n,&
                              a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                              a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                              u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                              u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                              p_ghost_cell_w,p_ghost_cell_e,&
                              p_ghost_cell_s,p_ghost_cell_n)
    
        call v_solver(nx,ny,nt,t,density,viscosity,scale_factor,ghost_cell_size,&
                      x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                      u_gs,v_gs,p_gs,&
                      u,u_alter,v,v_alter,&
                      p,p_alter,p_correct,&
                      f_former,f,f_alter,&
                      a_former,a,a_alter,&
                      b_former,b,b_alter,&
                      ap_u,ap_v,ap_u1,ap_v1,&
                      uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                      ubc_e,ubc_w,&
                      vbc_n,vbc_s,&
                      pbc_e,pbc_w,pbc_n,pbc_s,&
                      u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                      u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                      p_bound_w,p_bound_e,&
                      p_bound_s,p_bound_n,&
                      f_bound_w,f_bound_e,&
                      f_bound_s,f_bound_n,&
                      a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                      a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                      u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                      u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                      p_ghost_cell_w,p_ghost_cell_e,&
                      p_ghost_cell_s,p_ghost_cell_n)
                    

!STORE DATA OF ALL BLOCKS TO STORGE
        call data_collect("Storage","APV",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p_alter,f_alter,a,b,a_alter,b_alter,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)

        call data_collect("Storage","-V-",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p_alter,f_alter,a,b,a_alter,b_alter,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                              
        print *, "Finish V iternation for block",k,"timestep",t

x0_block = x_block ; y0_block = y_block
u0_block = u_block ; v0_block = v_block
p0_block = p_block ; f0_block = f_block
a0_block = a_block ; b0_block = b_block
af0_block = af_block ; bf0_block = bf_block
                         
!extract data from ghost cell
         !PRINT *, b(:,ny+1)

                          
        call data_collect("GrabOut","ALL",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        print *, "Extracted data from all! Now solve P"
                        
        !PRINT *, b(:,ny+1)
                        
                        
        call initial_boundary(nx,ny,nt,nb,k,t,ghost_cell_size,&
                              x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                              u,u_alter,v,v_alter,&
                              p,p_alter,p_correct,&
                              f_former,f,f_alter,&
                              a_former,a,a_alter,&
                              b_former,b,b_alter,&
                              u_gs,v_gs,p_gs,&
                              uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                              u_w,u_e,u_s,u_n,&
                              v_w,v_e,v_s,v_n,&
                              p_w,p_e,p_s,p_n,&
                              ubc_e,ubc_w,&
                              vbc_n,vbc_s,&
                              pbc_e,pbc_w,pbc_n,pbc_s,&
                              u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                              u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                              p_bound_w,p_bound_e,&
                              p_bound_s,p_bound_n,&
                              f_bound_w,f_bound_e,&
                              f_bound_s,f_bound_n,&
                              a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                              a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                              u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                              u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                              p_ghost_cell_w,p_ghost_cell_e,&
                              p_ghost_cell_s,p_ghost_cell_n)
                         
        call data_collect("GrabOut","UVP",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                          

                        
        call data_collect("GrabOut","APU",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
        
        call data_collect("GrabOut","APV",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                 
        call uv_correction(nx,ny,u_alter,v_alter,pbc_e,pbc_w,pbc_n,pbc_s)
        
        print *, "Corrected velocity u & v"

        call p_solver(nx,ny,nt,t,density,viscosity,scale_factor,ghost_cell_size,&
                      x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                      u_gs,v_gs,p_gs,&
                      u,u_alter,v,v_alter,&
                      p,p_alter,p_correct,&
                      f_former,f,f_alter,&
                      a_former,a,a_alter,&
                      b_former,b,b_alter,&
                      ap_u,ap_v,ap_u1,ap_v1,&
                      uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                      ubc_e,ubc_w,&
                      vbc_n,vbc_s,&
                      pbc_e,pbc_w,pbc_n,pbc_s,&
                      u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                      u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                      p_bound_w,p_bound_e,&
                      p_bound_s,p_bound_n,&
                      f_bound_w,f_bound_e,&
                      f_bound_s,f_bound_n,&
                      a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                      a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                      u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                      u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                      p_ghost_cell_w,p_ghost_cell_e,&
                      p_ghost_cell_s,p_ghost_cell_n)
                
        print *, "Now processing the data"

        call uvp_process(nx,ny,nt,ghost_cell_size,&
                              ap_u,ap_v,ap_u1,ap_v1,&
                              u,u_alter,v,v_alter,&
                              p,p_alter,p_correct,&
                              u_gs,v_gs,p_gs,&
                              u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                              u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                              p_bound_w,p_bound_e,&
                              p_bound_s,p_bound_n)
                            
!STORE DATA OF ALL BLOCKS TO STORGE

        call data_collect("Storage","UVP",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p_alter,f_alter,a,b,a_alter,b_alter,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        print *, "Finish P iternation for block",k,"timestep",t
    
x0_block = x_block ; y0_block = y_block
u0_block = u_block ; v0_block = v_block
p0_block = p_block ; f0_block = f_block
a0_block = a_block ; b0_block = b_block
af0_block = af_block ; bf0_block = bf_block

if (skip_new_theory .eq. 1) then
f0_block = 0
f_block = 0
a0_block = 0
a_block = 0
b0_block = 0
b_block = 0
af0_block = 0
af_block = 0
bf0_block = 0
bf_block = 0

f_former = 0
a_former = 0
b_former = 0

f = 0
a = 0
b = 0

f_alter = 0
a_alter = 0
b_alter = 0

goto 11111



end if


!solve finite element based flow theory

!extract data from ghost cell

        call data_collect("GrabOut","ALL",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        print *, "Extracted data from all! Now solve F"
                         
        call initial_boundary(nx,ny,nt,nb,k,t,ghost_cell_size,&
                              x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                              u,u_alter,v,v_alter,&
                              p,p_alter,p_correct,&
                              f_former,f,f_alter,&
                              a_former,a,a_alter,&
                              b_former,b,b_alter,&
                              u_gs,v_gs,p_gs,&
                              uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                              u_w,u_e,u_s,u_n,&
                              v_w,v_e,v_s,v_n,&
                              p_w,p_e,p_s,p_n,&
                              ubc_e,ubc_w,&
                              vbc_n,vbc_s,&
                              pbc_e,pbc_w,pbc_n,pbc_s,&
                              u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                              u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                              p_bound_w,p_bound_e,&
                              p_bound_s,p_bound_n,&
                              f_bound_w,f_bound_e,&
                              f_bound_s,f_bound_n,&
                              a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                              a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                              u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                              u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                              p_ghost_cell_w,p_ghost_cell_e,&
                              p_ghost_cell_s,p_ghost_cell_n)
                  !PRINT *, "F,",b(:,ny+1)                

                          
                  !PRINT *, "F2,",b(:,ny+1)   
 
        call f_solver(nx,ny,nt,dt,t,density,viscosity,scale_factor,ghost_cell_size,&
                      x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                      u_gs,v_gs,p_gs,&
                      u,u_alter,v,v_alter,&
                      p,p_alter,p_correct,&
                      f_former,f,f_alter,&
                      a_former,a,a_alter,&
                      b_former,b,b_alter,&
                      ap_u,ap_v,ap_u1,ap_v1,&
                      uvbc,pbc,winds,&
                      ubc_e,ubc_w,&
                      vbc_n,vbc_s,&
                      pbc_e,pbc_w,pbc_n,pbc_s,&
                      u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                      u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                      p_bound_w,p_bound_e,&
                      p_bound_s,p_bound_n,&
                      f_bound_w,f_bound_e,&
                      f_bound_s,f_bound_n,&
                      a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                      a_bound_s,a_bound_n,b_bound_s,b_bound_n)
        
        call data_collect("Storage","-F-",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p_alter,f_alter,a,b,a_alter,b_alter,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        print *, "Finish F iternation for block",k,"timestep",t
    
x0_block = x_block ; y0_block = y_block
u0_block = u_block ; v0_block = v_block
p0_block = p_block ; f0_block = f_block
a0_block = a_block ; b0_block = b_block
af0_block = af_block ; bf0_block = bf_block
                         
!extract data from ghost cell

        call data_collect("GrabOut","ALL",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u,v,p,f,a_former,b_former,a,b,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
        print *, "Extracted data from all! Now solve A&B"
                         
        call initial_boundary(nx,ny,nt,nb,k,t,ghost_cell_size,&
                              x,y,dx,dy,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                              u,u_alter,v,v_alter,&
                              p,p_alter,p_correct,&
                              f_former,f,f_alter,&
                              a_former,a,a_alter,&
                              b_former,b,b_alter,&
                              u_gs,v_gs,p_gs,&
                              uvbc(k,:,:),pbc(k,:,:),winds(k,:,1),&
                              u_w,u_e,u_s,u_n,&
                              v_w,v_e,v_s,v_n,&
                              p_w,p_e,p_s,p_n,&
                              ubc_e,ubc_w,&
                              vbc_n,vbc_s,&
                              pbc_e,pbc_w,pbc_n,pbc_s,&
                              u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                              u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                              p_bound_w,p_bound_e,&
                              p_bound_s,p_bound_n,&
                              f_bound_w,f_bound_e,&
                              f_bound_s,f_bound_n,&
                              a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                              a_bound_s,a_bound_n,b_bound_s,b_bound_n,&
                              u_ghost_cell_w,u_ghost_cell_e,v_ghost_cell_w,v_ghost_cell_e,&
                              u_ghost_cell_s,u_ghost_cell_n,v_ghost_cell_s,v_ghost_cell_n,&
                              p_ghost_cell_w,p_ghost_cell_e,&
                              p_ghost_cell_s,p_ghost_cell_n)
                         
                        

        call ab_solver(nx,ny,nt,dt,t,density,viscosity,scale_factor,ghost_cell_size,&
                       x,y,domain_x_start,domain_x_end,domain_y_start,domain_y_end,&
                       u_gs,v_gs,p_gs,&
                       u,u_alter,v,v_alter,&
                       p,p_alter,p_correct,&
                       f_former,f,f_alter,&
                       a_former,a,a_alter,&
                       b_former,b,b_alter,&
                       ap_u,ap_v,ap_u1,ap_v1,&
                       uvbc,pbc,winds,&
                       ubc_e,ubc_w,&
                       vbc_n,vbc_s,&
                       pbc_e,pbc_w,pbc_n,pbc_s,&
                       u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                       u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                       p_bound_w,p_bound_e,&
                       p_bound_s,p_bound_n,&
                       f_bound_w,f_bound_e,&
                       f_bound_s,f_bound_n,&
                       a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                       a_bound_s,a_bound_n,b_bound_s,b_bound_n)
        
        call fab_process(nx,ny,nt,t,ghost_cell_size,&
                            ap_u,ap_v,ap_u1,ap_v1,&
                            u,u_alter,v,v_alter,&
                            p,p_alter,p_correct,&
                            f_former,f,f_alter,&
                            a_former,a,a_alter,&
                            b_former,b,b_alter,&
                            u_gs,v_gs,p_gs,&
                            u_bound_w,u_bound_e,v_bound_w,v_bound_e,&
                            u_bound_s,u_bound_n,v_bound_s,v_bound_n,&
                            p_bound_w,p_bound_e,&
                            p_bound_s,p_bound_n,&
                            f_bound_w,f_bound_e,&
                            f_bound_s,f_bound_n,&
                            a_bound_w,a_bound_e,b_bound_w,b_bound_e,&
                            a_bound_s,a_bound_n,b_bound_s,b_bound_n)
        
        
        call data_collect("Storage","A&B",k,nb,nx,ny,&
                          block_grid_count1,block_grid_count2,&
                          block_grid_id1,block_grid_id2,&
                          x,y,u_alter,v_alter,p_alter,f_alter,a,b,a_alter,b_alter,&
                          x_block,y_block,u_block,v_block,p_block,f_block,af_block,bf_block,a_block,b_block,&
                          ap_u,ap_v,ap_u1,ap_v1,&
                          ap_u_block,ap_v_block)
                        
    print *, "Finish A&B iternation for block",k,"timestep",t
    
11111 continue






    !OUTPUT 
    if (mod(t,dump_every).eq.0) then
                          
                                     
        write (print_iteration,'(i5)') t
        output_filename="All_tecplot_iter"//trim(adjustl(print_iteration))//".plt"
        output_filename=trim(adjustl(output_filename))
        
        open(100,file=output_filename,status='replace')
        write (100,*) "VARIABLES= X,Y,U,V,P,F,A,B"
        write (100,*) "ZONE I=",nx+2,",J=",ny+2,",F=POINT"
       
            do number_block=1,nb
                do j=1,gy(number_block)+2
                    do i=1,gx(number_block)+2
                        write(100,*) x_block( j+(ny+2)*(i-1)+block_grid_id1(number_block) ),&
                                     y_block( j+(ny+2)*(i-1)+block_grid_id1(number_block) ),&
        (u_block( j+(ny+3)*(i-1)+block_grid_id2(number_block) )+u_block( j+(ny+3)*(i)+block_grid_id2(number_block) ))/2.0d0,&
        (v_block( j+(ny+3)*(i-1)+block_grid_id2(number_block) )+v_block( j+1+(ny+3)*(i-1)+block_grid_id2(number_block) ))/2.0d0,&
                                     p_block( j+(ny+2)*(i-1)+block_grid_id1(number_block) ),&
                                     f_block( j+(ny+2)*(i-1)+block_grid_id1(number_block) ),&
                                     a_block( j+(ny+2)*(i-1)+block_grid_id1(number_block) ),&
                                     b_block( j+(ny+2)*(i-1)+block_grid_id1(number_block) )
                    end do
                end do
            end do

        close(100)
    end if


end do

    deallocate(ap_u,ap_v,ap_u1,ap_v1)
    deallocate(u_bound_w,u_bound_e,u_bound_n,u_bound_s) 
    deallocate(v_bound_w,v_bound_e,v_bound_n,v_bound_s)
    deallocate(p_bound_w,p_bound_e,p_bound_n,p_bound_s)
    deallocate(f_bound_w,f_bound_e,f_bound_n,f_bound_s)
    deallocate(a_bound_w,a_bound_e,a_bound_n,a_bound_s)
    deallocate(b_bound_w,b_bound_e,b_bound_n,b_bound_s)
    deallocate(f_former,f,f_alter,a_former,a,a_alter,b_former,b,b_alter)
    deallocate(u_alter,v_alter,p_alter,p_correct)
    deallocate(u,v,p)
    deallocate(u_gs,v_gs,p_gs)
    deallocate(x,y)
    deallocate(u_ghost_cell_w,u_ghost_cell_e,u_ghost_cell_n,u_ghost_cell_s) 
    deallocate(v_ghost_cell_w,v_ghost_cell_e,v_ghost_cell_n,v_ghost_cell_s)
    deallocate(p_ghost_cell_w,p_ghost_cell_e,p_ghost_cell_n,p_ghost_cell_s)


print *,"end of calculation"




close(88)
close (2)
close (1)

end program New_2D


