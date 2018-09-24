program ExpoAction_Test

	implicit none
	
	complex(16), dimension(:), allocatable :: expo, mma
	real(16), dimension(:), allocatable :: delta_real, delta_complex
	real(16) :: mean_real, mean_complex, allowed_diff
	complex(16) :: average_difference
	real(16) :: maximum_difference
	integer :: vector_length, io, i
	
	open(10,file='../ExpoAction/Output/out.rsa_v.dat   1.0000')
	open(20,file='Mvexp.dat')
	
	vector_length = 0
	
	allowed_diff = 7.99361Q-13
		
	do
		read(10,*,iostat=io)
		if (io.ne.0) exit
		vector_length = vector_length + 1
	end do

	allocate(expo(vector_length), mma(vector_length), &
			& delta_real(vector_length))
	
	rewind(10)
	
	read(10,*) (expo(i),i=1,vector_length)
	read(20,*) (mma(i),i=1,vector_length)

	delta_real = abs((real(expo) - real(mma))/real(expo))*100
	delta_complex = abs((aimag(expo) - aimag(mma))/aimag(expo))*100
	
	mean_real = sum(delta_real)/vector_length
	mean_complex = sum(delta_complex)/vector_length
	
	write(*,'(A13,E15.5,A17,E15.5)') "Real % Diff: ", mean_real, &
			& " Complex % Diff: ", mean_complex
			
	write(*,*)
	
	write(*,'(A17,E15.5,A21,E15.5)') "Max Real % Diff: ", maxval(delta_real), &
			& " Max Complex % Diff: ", maxval(delta_complex)
			
	if((mean_real.le.allowed_diff).and.(mean_complex.le.allowed_diff))then
		write(*,*) "OK"
	end if

end program ExpoAction_Test
