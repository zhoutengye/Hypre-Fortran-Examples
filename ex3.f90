!---------------------------------------------------------------------
! Notes:
! The function 
!    Call HYPRE_StructPCGSetPrecond(solver, precond_id, precond, ier)
! is not consistent with the C version interface, it should be written
! as the above form.
! Ite seems that the SMG pre_conditioner would be activated when
!     predond_id = 0
!---------------------------------------------------------------------
program main
  Implicit None

  ! MPI variables
  Integer :: i, j, k
  Integer :: myid, num_procs, MPI_COMM_world
  Integer :: ier

  Integer :: n, NN, pi, pj
  Real(8) :: h, h2
  Real(8) :: PHI(33,33)
  Real(8) :: PHI_FULL(33,33)
  Integer :: ilower(2), iupper(2)

  Integer :: solver_id
  Integer :: n_pre, n_post

  Integer*8:: grid
  Integer*8 :: stencil
  Integer*8 :: A
  Integer*8 :: b
  Integer*8 :: x
  Integer*8 :: solver
  Integer*8 :: precond
  Integer*8 :: precond_id

  Integer :: num_iterations
  Real(8) :: final_res_norm

  Call MPI_INIT(MPI_COMM_WORLD, myid, ier)
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ier)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ier)

  n = 33
  solver_id = 0
  n_pre  = 1
  n_post = 1

  ! Figure out the processor grid (N x N).  The local problem
  ! size for the interior nodes is indicated by n (n x n).
  ! pi and pj indicate position in the processor grid.
  NN = sqrt(float(num_procs))
  h  = 1.0 / float(NN*n+1) ! note that when calculating h we must remember to count the boundary nodes 
  h2 = h*h
  pj = myid / float(NN)
  pi = myid - pj*NN

  ! Figure out the extents of each processor's piece of the grid.
  ilower(1) = pi*n;
  ilower(2) = pj*n;

  iupper(1) = ilower(1) + n-1
  iupper(2) = ilower(2) + n-1



  !-1. Set up a grid. 

  ! Create an empty 2D grid object */
  Call HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, grid, ier)
  ! Add a new box to the grid */
  Call HYPRE_StructGridSetExtents(grid, ilower, iupper, ier)
  ! This is a collective call finalizing the grid assembly.
  ! The grid is now ``ready to be used'' */
  Call HYPRE_StructGridAssemble(grid, ier)

  !--2. Define the discretization stencil------------------
  ! Create an empty 2D, 5-pt stencil object
  Block
    Integer :: entry
    Integer :: offsets(5,2)
    Integer :: offset(2)
    Call HYPRE_StructStencilCreate(2, 5, stencil, ier)

    ! Define the geometry of the stencil. Each represents a
    ! relative offset (in the index space)
    offsets(1,:) = [0,0]
    offsets(2,:) = [-1,0]
    offsets(3,:) = [1,0]
    offsets(4,:) = [0,-1]
    offsets(5,:) = [0,1]

    ! Assign each of the 5 stencil entries
    Do entry = 1,5
      offset = offsets(entry,:)
      Call HYPRE_StructStencilSetElement(stencil, entry-1, offset, ier);
    End Do
  End Block


  !-----3. Set up a Struct Matrix--------------------------
  Block
    Integer :: nentries = 5
    Integer, Allocatable :: stencil_indices(:)
    Integer :: nvalues
    Real(8), Allocatable :: values(:)
    nvalues = nentries*n*n
    Allocate(stencil_indices(nentries))
    Allocate(values(nvalues))
    ! Create an empty matrix object 
    Call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, A, ier);

    Call HYPRE_StructMatrixInitialize(A, ier)

    ! labels for the stencil entries these correspond to the offsets
    ! defined above 
    stencil_indices = [0,1,2,3,4]
    nentries = 5


    ! We have 6 grid points, each with 5 stencil entries 
    Do i = 1, nvalues, nentries
      values(i) = 4.0
      Do j = 2, nentries
        values(i+j-1) = -1.0
      End Do
    End Do

    Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, &
        stencil_indices, values, ier)

    ! Deallocate(stencil_indices)
    ! Deallocate(values)

  End Block

  ! 4. Incorporate the zero boundary conditions: go along each edge of
  ! the domain and set the stencil entry that reaches to the boundary to
  ! zero.
  Block
    Integer bc_ilower(2)
    Integer bc_iupper(2)
    Integer :: nentries = 1
    Integer, Allocatable :: stencil_indices(:)
    Integer :: nvalues
    Real(8), Allocatable :: values(:)
    nvalues = nentries*n
    Allocate(stencil_indices(nentries))
    Allocate(values(nvalues))
    do i = 1, nvalues
      values(i) = 0.0
    enddo
    if (pj .eq. 0) Then
      ! Bottom row of grid points 
      Block
        bc_ilower(1) = pi*n
        bc_ilower(2) = pj*n
        bc_iupper(1) = bc_ilower(1) + n - 1
        bc_iupper(2) = bc_ilower(2)
        stencil_indices = 3
        Call HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1, &
            stencil_indices, values, ier)
      End Block
    Else if (pj .eq. NN-1) Then
      ! Upper row of grid points 
      Block
        bc_ilower(1) = pi*n
        bc_ilower(2) = pj*n + n - 1
        bc_iupper(1) = bc_ilower(1) + n - 1
        bc_iupper(2) = bc_ilower(2)
        stencil_indices = 4
        Call HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1, &
            stencil_indices, values, ier)
      End Block
    Else if (pi .eq. 0) Then
      ! Left row of grid points 
      Block
        bc_ilower(1) = pi*n 
        bc_ilower(2) = pj*n 
        bc_iupper(1) = bc_ilower(1)
        bc_iupper(2) = bc_ilower(2) + n - 1
        stencil_indices = 1
        Call HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1, &
            stencil_indices, values, ier)
      End Block
    Else if (pi .eq. NN-1) Then
      ! Right row of grid points 
      Block
        bc_ilower(1) = pi*n + n - 1
        bc_ilower(2) = pj*n 
        bc_iupper(1) = bc_ilower(1)
        bc_iupper(2) = bc_ilower(2) + n - 1
        stencil_indices = 2
        Call HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1, &
            stencil_indices, values, ier)
      End Block
    end if

    ! Deallocate(stencil_indices)
    ! Deallocate(values)

    ! This is a collective call finalizing the matrix assembly.
    ! The matrix is now ``ready to be used'' 
    Call HYPRE_StructMatrixAssemble(A, ier)

  End Block

  ! 5. Set up Struct Vectors for b and x 
  Block
    Integer :: nvalues 
    Real(8), Allocatable :: values(:)
    nvalues = n*n
    Allocate(values(nvalues))
    Call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, b, ier)
    Call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, x, ier)

    Call HYPRE_StructVectorInitialize(b, ier)
    Call HYPRE_StructVectorInitialize(x, ier)

    Do i = 1, nvalues
      values(i) = h2
    End Do
    Call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values, ier)

    Do i = 1, nvalues
      values(i) = 0.0
    end Do
    Call HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values, ier)

    ! Deallocate(values)

    Call HYPRE_StructVectorAssemble(b, ier)
    Call HYPRE_StructVectorAssemble(x, ier)

  end Block

  !---6. Set up and use a solver---------------------------------
  ! (See the Reference Manual for descriptions of all of the options.)

  ! Create an empty PCG Struct solver 
  Call HYPRE_StructPCGCreate(MPI_COMM_WORLD, solver, ier)

  ! Set some parameters 
  Call HYPRE_StructPCGSetMaxIter(solver, 50, ier)
  Call HYPRE_StructPCGSetTol(solver, 1.0d-06, ier)
  Call HYPRE_StructPCGSetTwoNorm(solver, 1, ier)
  Call HYPRE_StructPCGSetRelChange(solver, 0, ier)
  Call HYPRE_StructPCGSetPrintLevel(solver, 2, ier)
  Call HYPRE_StructPCGSetLogging(solver, 1, ier)


  ! Use symmetric SMG as preconditioner 
  Call HYPRE_StructSMGCreate(MPI_COMM_WORLD, precond, ier)
  Call HYPRE_StructSMGSetMemoryUse(precond, 0, ier)
  Call HYPRE_StructSMGSetMaxIter(precond, 1, ier)
  Call HYPRE_StructSMGSetTol(precond, 0.0, ier)
  Call HYPRE_StructSMGSetZeroGuess(precond, ier)
  Call HYPRE_StructSMGSetNumPreRelax(precond, 1, ier)
  Call HYPRE_StructSMGSetNumPostRelax(precond, 1, ier)

  

  ! Setup and solve 
  precond_id = 0
  Call HYPRE_StructPCGSetPrecond(solver, precond_id, precond, ier)
  Call HYPRE_StructPCGSetup(solver, A, b, x, ier)
  Call HYPRE_StructPCGSolve(solver, A, b, x, ier)

  ! Get some info on the run 
  Call HYPRE_StructPCGGetNumIterations(solver, num_iterations, ier)
  ! Thic function does not seem exist
  ! Call HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, final_res_norm, ier)

  ! Free memory 
  Call HYPRE_StructGridDestroy(grid, ier)
  Call HYPRE_StructStencilDestroy(stencil, ier)
  Call HYPRE_StructMatrixDestroy(A, ier)
  Call HYPRE_StructVectorDestroy(b, ier)
  Call HYPRE_StructVectorDestroy(x, ier)
  Call HYPRE_StructPCGDestroy(solver, ier)
  Call HYPRE_StructSMGDestroy(precond, ier)


  If (myid .eq. 0) Then
    print *, "Iterations = ", num_iterations
    print *, "Final Relative Residual Norm = ", final_res_norm
  End If
  
  Call MPI_Finalize(ier)

end program main

