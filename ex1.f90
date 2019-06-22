program main
  Implicit None

  ! MPI variables
  Integer :: i, j, myid, num_procs, MPI_COMM_world
  Integer :: ier

  Integer :: vis = 0

  Integer*8 :: grid
  Integer*8 :: stencil
  Integer*8 :: A
  Integer*8 :: b
  Integer*8 :: x
  Integer*8 :: solver

  Call MPI_INIT(MPI_COMM_WORLD, myid, ier)
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ier)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ier)

  If (num_procs .ne. 2) Then
    If (myid .eq.0) print *, ("Must run with 2 processors!")
    Call MPI_Finalize(ier)
    stop
  End If

  !-1. Set up a grid. Each processor describes the piece-----
  ! of the grid that it owns. 

  ! Create an empty 2D grid object */
  Block
    Integer :: ilower(2)
    Integer :: iupper(2)
    Call HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, grid, ier)
    ! Add boxes to the grid
    If (myid .eq. 0) Then
      ilower = [-3,1]
      iupper = [-1,2]
      Call HYPRE_StructGridSetExtents(grid, ilower, iupper, ier)
    Else if (myid .eq. 1) Then
      ilower = [0,1]
      iupper = [2,4]
      Call HYPRE_StructGridSetExtents(grid, ilower, iupper, ier)
    End If

    Call HYPRE_StructGridAssemble(grid, ier)
  End Block

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
    ! Create an empty matrix object 
    Call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, A, ier);

    ! Indicate that the matrix coefficients are ready to be set 
    Call HYPRE_StructMatrixInitialize(A, ier)

    ! Set the matrix coefficients.  Each processor assigns coefficients
    ! for the boxes in the grid that it owns. Note that the coefficients
    ! associated with each stencil entry may vary from grid point to grid
    ! point if desired.  Here, we first set the same stencil entries for
    ! each grid point.  Then we make modifications to grid points near
    ! the boundary.

    if (myid .eq. 0) Then
      Block
        Integer :: ilower(2)
        Integer :: iupper(2)
        Integer :: stencil_indices(5)
        Integer :: nentries
        Integer :: nvalues
        Real(8) :: values(30)
        ! labels for the stencil entries these correspond to the offsets
        ! defined above 
        ilower = [-3,1]
        iupper = [-1,2]
        stencil_indices = [0,1,2,3,4]
        nentries = 5
        nvalues = 30 ! 6 grid points, each with 5 stencil entries 

        ! We have 6 grid points, each with 5 stencil entries 
        Do i = 1, nvalues, nentries
          values(i) = 4.0
          Do j = 2, nentries
            values(i+j-1) = -1.0
          End Do
        End Do

        Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, &
            stencil_indices, values, ier)
      End Block
    else if (myid .eq. 1) Then
      Block
        Integer :: ilower(2)
        Integer :: iupper(2)
        Integer :: stencil_indices(5)
        Integer :: nentries
        Integer :: nvalues
        Real(8) :: values(60)

        ! labels for the stencil entries these correspond to the offsets
        ! defined above
        ilower = [0,1]
        iupper = [2,4]
        stencil_indices = [0,1,2,3,4]
        nentries = 5
        nvalues = 60 ! 12 grid points, each with 5 stencil entries 

        ! We have 12 grid points, each with 5 stencil entries 
        Do i = 1, nvalues, nentries
          values(i) = 4.0
          Do j = 2, nentries
            values(i+j-1) = -1.0
          End Do
        End Do



        Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, &
            stencil_indices, values, ier)
      End Block
    end if

    ! Set the coefficients reaching outside of the boundary to 0
    if (myid .eq. 0) Then
      Block
        real(8) :: values(3)
        Integer :: ilower(2)
        Integer :: iupper(2)
        Integer :: stencil_indices(1)
        do i = 1,3
          values(i) = 0.0
        End Do
        Block
          ilower = [-3,1]
          iupper = [-1,1]
          stencil_indices = 3
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
        Block
          ilower = [-3,1]
          iupper = [-3,2]
          stencil_indices = 1
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
        Block
          ilower = [-3,2]
          iupper = [-1,2]
          stencil_indices = 4
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
      End Block
    else if (myid .eq. 1) Then
      Block
        real(8) :: values(4)
        Integer :: ilower(2)
        Integer :: iupper(2)
        Integer :: stencil_indices(1)
        do i = 1,4
          values(i) = 0.0
        End Do
        Block
          ilower = [0,1]
          iupper = [2,1]
          stencil_indices = 3
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
        Block
          ilower = [2,1]
          iupper = [2,4]
          stencil_indices = 2
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
        Block
          ilower = [0,4]
          iupper = [2,4]
          stencil_indices = 4
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
        Block
          ilower = [0,3]
          iupper = [0,4]
          stencil_indices = 1
          Call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, &
              stencil_indices, values, ier)
        End Block
      End Block
    end if

    ! This is a collective call finalizing the matrix assembly.
    ! The matrix is now ``ready to be used'' 
    Call HYPRE_StructMatrixAssemble(A, ier)
  End Block

  !--4. Set up Struct Vectors for b and x.
  ! Each processor sets the vectors corresponding to its boxes.
  Block
    Call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, b, ier)
    Call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, x, ier)

    Call HYPRE_StructVectorInitialize(b, ier)
    Call HYPRE_StructVectorInitialize(x, ier)
    If (myid .eq. 0) Then
      Block
        Integer :: ilower(2)
        Integer :: iupper(2)
        Real(8) :: values(6)
        ilower = [-3,1]
        iupper = [-1,2]
        Do i = 1,6
          values(i) = 1.0
        End Do

        Call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values, ier)

        Do i = 1,6
          values(i) = 0.0
        End Do

        Call HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values, ier)
      End Block
    Else If (myid .eq. 1) Then
      Block
        Integer :: ilower(2)
        Integer :: iupper(2)
        Real(8) :: values(12)
        ilower = [0,1]
        iupper = [2,4]
        Do i = 1,12
          values(i) = 1.0
        End Do

        Call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values, ier)

        Do i = 1,12
          values(i) = 0.0
        End Do

        Call HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values, ier)
      End Block
    End If

    Call HYPRE_StructVectorAssemble(b, ier)
    Call HYPRE_StructVectorAssemble(x, ier)

  end Block

  !---5. Set up and use a solver---------------------------------
  ! (See the Reference Manual for descriptions of all of the options.)

  ! Create an empty PCG Struct solver 
  Call HYPRE_StructPCGCreate(MPI_COMM_WORLD, solver, ier)

  ! Set some parameters 
  Call HYPRE_StructPCGSetTol(solver, 1.0d-06, ier)
  Call HYPRE_StructPCGSetPrintLevel(solver, 2, ier)

  ! Setup and solve 
  Call HYPRE_StructPCGSetup(solver, A, b, x, ier)
  Call HYPRE_StructPCGSolve(solver, A, b, x, ier)


  ! Free memory 
  Call HYPRE_StructGridDestroy(grid, ier)
  Call HYPRE_StructStencilDestroy(stencil, ier)
  Call HYPRE_StructMatrixDestroy(A, ier)
  Call HYPRE_StructVectorDestroy(b, ier)
  Call HYPRE_StructVectorDestroy(x, ier)
  Call HYPRE_StructPCGDestroy(solver, ier)

  If (myid .eq. 0) print *, 'yes'
  Call MPI_Finalize(ier)

end program main

