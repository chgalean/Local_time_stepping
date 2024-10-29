function solve_linear_system(Aglo, Bglo, Fglo)
    ############################## MUMPS ##########################################
    MPI.Init()
    root = 0
    comm = MPI.COMM_WORLD

    icntl = default_icntl[:]
    icntl[1] = -1  # output stream for error messages
    icntl[2] = -1  # output stream for diagnostic printing and statistics local to each MPI process
    icntl[3] = -1  # output stream for global information, collected on the host
    icntl[4] = -1  # level of printing for error, warning, and diagnostic messages
    mumps0 = MUMPS.Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)
    if MPI.Comm_rank(comm) == root
        MUMPS.associate_matrix!(mumps0, Aglo)
        MUMPS.associate_rhs!(mumps0, Bglo')
    end
    MUMPS.factorize!(mumps0)
    MUMPS.solve!(mumps0)
    MPI.Barrier(comm)
    if MPI.Comm_rank(comm) == root
        Schu_comp=Bglo*MUMPS.get_solution(mumps0)
    end
    finalize(mumps0)

    mumps1 = MUMPS.Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)
    if MPI.Comm_rank(comm) == root
        MUMPS.associate_matrix!(mumps1, Aglo)
        MUMPS.associate_rhs!(mumps1, Fglo)
    end
    MUMPS.factorize!(mumps1)
    MUMPS.solve!(mumps1)
    MPI.Barrier(comm)
    if MPI.Comm_rank(comm) == root
        RHS = Bglo*MUMPS.get_solution(mumps1)
    end
    finalize(mumps1)
    
    mumps2 = MUMPS.Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)
    if MPI.Comm_rank(comm) == root
        MUMPS.associate_matrix!(mumps2, Schu_comp)
        MUMPS.associate_rhs!(mumps2, RHS)
    end
    MUMPS.factorize!(mumps2)
    MUMPS.solve!(mumps2)
    MPI.Barrier(comm)
    if MPI.Comm_rank(comm) == root
        p = MUMPS.get_solution(mumps2)
    end
    finalize(mumps2)

    mumps3 = MUMPS.Mumps{Float64}(mumps_symmetric, icntl, default_cntl64)
    if MPI.Comm_rank(comm) == root
        MUMPS.associate_matrix!(mumps3, Aglo)
        RHS1=Fglo-Bglo'*p
        MUMPS.associate_rhs!(mumps3, RHS1)
    end
    MUMPS.factorize!(mumps3)
    MUMPS.solve!(mumps3)
    MPI.Barrier(comm)
    if MPI.Comm_rank(comm) == root
        UV = MUMPS.get_solution(mumps3)
    end
    finalize(mumps3)
    return UV, p
    MPI.Finalize()
end