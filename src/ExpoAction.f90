module Expo
	
	use Nullarbor

	
	implicit none
		
	contains
		
	subroutine Read_Options(m_min, m_max, p_min, p_max, l, theta, & 
			& ms_and_ps)
			
		integer, intent(out) :: m_min, m_max, p_min, p_max, l
		real(8), dimension(:), allocatable, intent(out) :: theta
		integer, dimension(:,:), allocatable, intent(out) :: ms_and_ps
		character(128) :: matrix_name, vector_name
		
		integer :: ms_and_ps_length
		integer :: i, j, indx

		open(unit = 10, file = "Options/ExpoOptions.txt", status = 'old', &
				& form = 'formatted')
						
		read(10,*)
		read(10,*)
		read(10,*) m_min
		read(10,*)
		read(10,*) m_max
		read(10,*)
		read(10,*) p_min
		read(10,*)
		read(10,*) p_max
		read(10,*)
		read(10,*) l
		
		close(10)
		
		ms_and_ps_length = 0
		indx = 1
				
		do i = 2, p_max
		 do j = (i*(i-1)-1)+2, m_max
	
				ms_and_ps_length = ms_and_ps_length + 1
	
		 end do
		end do
		
		allocate(ms_and_ps(ms_and_ps_length,2))
				
		do i = 2, p_max
		 do j = (i*(i-1)-1)+2, m_max
	
			ms_and_ps(indx, 1) = j
			ms_and_ps(indx,2) = i 
			indx = indx + 1
			
		 end do
		end do
	
		open(unit = 10, file = "Options/data/theta.dat", status = 'old', &
			& form = 'formatted')
			
		allocate(theta(m_min:m_max))
	
		do i = m_min, m_max
		
			read(10,*) theta(i) 

		end do
		
		close(10)
		
	end subroutine Read_Options
	
	subroutine Parameters(A, m_min, m_max, p_min, p_max, l, theta, ms_and_ps, m, s)
	
		type(CSR), intent(in) :: A
		integer, intent(in) :: m_min, m_max, p_min, p_max, l
		real(8), dimension(m_min:m_max), intent(in) :: theta
		integer, dimension(:,:), intent(in) :: ms_and_ps
		integer, intent(out) :: m, s
		
		type(CSR) :: A_Temp, A_Dagger
		real(8) :: one_norm	
		real(8), dimension(p_min:p_max+1) :: one_norm_array
		real(8), dimension(p_min:p_max) :: alpha_array
		integer, dimension(m_min:m_max) :: m_list
		integer(8), dimension(size(ms_and_ps)) :: c_array
		integer :: min_c
		integer :: i, j, indx
		
		indx = 1
		A_Temp = A
		A_Dagger = Sparse_Conjugate(Sparse_Transpose(A))
		one_norm = Sparse_OneNormEst(A, A_Dagger, 1, l)
		
		if(one_norm.eq.0) then	
			
			m = 0
			s = 1
			
		else if(Test(one_norm, theta, m_max, p_max, l)) then
			do i = m_min, m_max
				m_list(i) = i*ceiling(one_norm/theta(i))
			end do
					
			m = minloc(m_list,1) + m_min - 1
			s = ceiling(one_norm/theta(m))
			
		else
	
			do i = p_min, p_max + 1 
				one_norm_array(i) = Sparse_OneNormEst(A, A_Dagger, i, l)**(1/real(i,8))
			end do
			
			do i = p_min, p_max
				if(one_norm_array(i).gt.one_norm_array(i+1)) then
					alpha_array(i) = one_norm_array(i)
				else
					alpha_array(i) = one_norm_array(i+1)
				end if
			end do

			do i = 1, size(ms_and_ps,1)
				c_array(i) = ms_and_ps(i,1)*ceiling(alpha_array(ms_and_ps(i,2)) &
						& /theta(ms_and_ps(i,1)))
			end do
			
			min_c = minval(c_array,1)
			m = ms_and_ps(minloc(c_array,1),1)
			s = max(min_c/m, 1)
	
		end if

	end subroutine Parameters

	function Test(one_norm, theta, m_max, p_max, l)
		
		logical :: test
		real(8), intent(in) :: one_norm
		real(8), dimension(:), intent(in) :: theta(3:m_max)
		
		integer :: m_max, p_max, l
		If(one_norm.le.(2*l*p_max*(p_max+3)*theta(m_max)/m_max))then
			
			test = .true.
		
		else
		
			test = .false.
		
		end if
		
	end function Test

	function Matrix_Exponential_Action(A, b_vec, m, s)
	
		complex(8), dimension(:), allocatable :: Matrix_Exponential_Action
		type(CSR), intent(in) :: A
		complex(8), dimension(:), intent(in) :: b_vec
		integer, intent(in) :: m, s
		
		complex(8), dimension(:), allocatable :: b_vec_temp
		real(8) :: c1, c2 = 0
		integer :: i, j, k, length
		real(8) :: TOLERENCE = 2.d0**(-53)

		length = size(b_vec)
		
		allocate(Matrix_Exponential_Action(size(b_vec)))
		
		Matrix_Exponential_Action = 0
		
		allocate(b_vec_temp(size(b_vec)))
		
		!$omp parallel do schedule(guided) 
		do i = 1, length
			b_vec_temp(i) = b_vec(i)
		end do
		!$omp end parallel do
		
		!$omp parallel do schedule(guided) 
		do i = 1, length
			Matrix_Exponential_Action(i) = b_vec(i)
		end do
		!$omp end parallel do
			
		
		do i = 1, s
				c1 = Vector_InfNorm(b_vec_temp)
		do j = 1, m
		 
				b_vec_temp = (A.dot.b_vec_temp)/(real(s*j,8))
				c2 = Vector_InfNorm(b_vec_temp)
				Matrix_Exponential_Action = Matrix_Exponential_Action + b_vec_temp
				if((c1+c2).lt.(TOLERENCE*Vector_InfNorm(Matrix_Exponential_Action))) &
						then
					exit 
				end if
				
				c1 = c2
			
			end do
			
			!$omp parallel do schedule(guided) 
			do k = 1, length
				b_vec_temp(k) = Matrix_Exponential_Action(k)
			end do
			!$omp end parallel do
				
		end do
		
	end function Matrix_Exponential_Action

end module Expo


program ExpoAction
		
	use Expo
	use omp_lib
	
	implicit none
	
	type(CSR) :: A		
	complex(8), dimension(:), allocatable :: ExpAction	
	complex(8), dimension(:,:), allocatable :: dense
	complex(8), dimension(:), allocatable :: v
	integer :: m_min, m_max, p_min, p_max, l, i, m, s, input_number
	integer, dimension(:,:), allocatable :: ms_and_ps
	real(8) :: t, tstart, tend
	real(8), dimension(:), allocatable :: theta
	character(128) :: matrix_name, vector_name, matrix_type, input_file, log_name
	character(128) :: time
	character(8)  :: date
  character(6) :: current_time
  
	call get_command_argument(1,input_file)
		
	!read in options
	call Read_Options(m_min, m_max, p_min, p_max, l, theta, & 
			& ms_and_ps)

	!set up logfile
	call date_and_time(date,current_time)

  log_name = "Output/LogFiles/"//date//"_"//current_time
  
	open(30,file=log_name,status='replace',form='formatted')
	write(30,*) "ExpoAction LogFile: " // log_name
	write(30,*)

	!opening input file list
	open(20,file= trim(input_file), status = 'old', form = 'formatted')
	read(20,*)
	read(20,*)
	read(20,*)
	read(20,*) input_number
	
	do i = 1, input_number
		  
		read(20,*)	
		read(20,*)
		read(20,*) matrix_name
		read(20,*)
		read(20,*) vector_name
		read(20,*)
		read(20,*) t
		read(20,*)
		read(20,*) matrix_type
		
		!writing to logfile
		write(30,*) "Matrix: ", trim(matrix_name)
		write(30,*) "Vector: ", trim(vector_name)
		write(30,*) "Time: ", t
		write(30,*) "Matrix Type: ", trim(matrix_type)
		
		!for result output file name
		write(time,*) t
		
		!load matrix and vector data
		call Nullarbor_Import_Vector(v, "Input/Vectors/" // trim(vector_name))
		
		if(trim(matrix_type).eq."dense")then
			
			call Nullarbor_Import_Dense_Matrix(dense, &
					& "Input/Matrices/" // trim(matrix_name))
					
			A = Sparse(dense)
			
			deallocate(dense)
		
		else if(trim(matrix_type).eq."Harwell-Boeing")then
		
			call Nullarbor_Import_Harwell_Boeing(A, &
					& "Input/Matrices/" // trim(matrix_name))
			
		else
		
			call Nullarbor_Import_Sparse_Matrix(A, "Input/Matrices/" // trim(matrix_name))
			
			A = A
			
		end if
		
		!TIMING START
		tstart = omp_get_wtime()
		
		A = t*A

		!determine scaling and squaring parameters
		call Parameters(A, m_min, m_max, p_min, p_max, l, theta, ms_and_ps, m, s)

		ExpAction = Matrix_Exponential_Action(A, v, m, s)

		tend = omp_get_wtime()
		!TIMING END
		
		!more logfile output
		write(30,*) "Compute Time: ", tend - tstart
		write(30,*) "Threads: ", omp_get_max_threads()
	
		!output result
		call Nullarbor_Export_Vector(ExpAction,"Output/"//trim(matrix_name)&
			//"_"//trim(vector_name)//trim(time(1:9)))
	
		call Sparse_Deallocate(A)
		deallocate(v)
	
	end do
	
	close(20)
	close(30)
	
end program ExpoAction
