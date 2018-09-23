program test
	
	use Nullarbor
	
	implicit none
	
	type(CSR) :: A

	call Nullarbor_Import_Harwell_Boeing ('out.rsa', A)
	
	call Sparse_Write(A,1,1,10,10)

end program test
