for (fname, elty) in ((:sqrm_spfct_get_rp_c, :Float32   ),
                      (:dqrm_spfct_get_rp_c, :Float64   ),
                      (:cqrm_spfct_get_rp_c, :ComplexF32),
                      (:zqrm_spfct_get_rp_c, :ComplexF64))
    @eval begin
        function qrm_spfct_get_rp(spfct :: qrm_spfct{$elty})
            err = $fname(spfct, spfct.ptr_rp)
            # Fix it with the release 3.0.5
            # qrm_check(err)
            rp = unsafe_wrap(Array, spfct.ptr_rp[], spfct.fct.m)
            return rp
        end
    end
end

for (fname, elty) in ((:sqrm_spfct_get_cp_c, :Float32   ),
                      (:dqrm_spfct_get_cp_c, :Float64   ),
                      (:cqrm_spfct_get_cp_c, :ComplexF32),
                      (:zqrm_spfct_get_cp_c, :ComplexF64))
    @eval begin
        function qrm_spfct_get_cp(spfct :: qrm_spfct{$elty})
            err = $fname(spfct, spfct.ptr_cp)
            # Fix it with the release 3.0.5
            # qrm_check(err)
            cp = unsafe_wrap(Array, spfct.ptr_cp[], spfct.fct.n)
            return cp
        end
    end
end

for (fname, elty) in ((:sqrm_spfct_get_r_c, :Float32   ),
                      (:dqrm_spfct_get_r_c, :Float64   ),
                      (:cqrm_spfct_get_r_c, :ComplexF32),
                      (:zqrm_spfct_get_r_c, :ComplexF64))
    @eval begin
        function qrm_spfct_get_r(spfct :: qrm_spfct{$elty})
            spmat = qrm_spmat_init($elty)
            err = $fname(spfct, spmat)
            qrm_check(err)
            I = unsafe_wrap(Array, spmat.mat.irn, spmat.mat.nz)
            J = unsafe_wrap(Array, spmat.mat.jcn, spmat.mat.nz)
            V = unsafe_wrap(Array, spmat.mat.val, spmat.mat.nz)
            R = sparse(I, J, V, spmat.mat.m, spmat.mat.n)
            return R
        end
    end
end

for (fname, elty) in ((:sqrm_spmat_init_c, :Float32   ),
                      (:dqrm_spmat_init_c, :Float64   ),
                      (:cqrm_spmat_init_c, :ComplexF32),
                      (:zqrm_spmat_init_c, :ComplexF64))
    @eval begin
        function qrm_spmat_init!(spmat :: qrm_spmat{$elty})
            err = $fname(spmat)
            qrm_check(err)
            return nothing
        end

        function qrm_spmat_init(::Type{$elty})
            spmat = qrm_spmat{$elty}()
            qrm_spmat_init!(spmat)
            return spmat
        end

        function qrm_spmat_init!(spmat :: qrm_spmat{$elty}, m :: Integer, n :: Integer, irn :: AbstractVector{I}, jcn :: AbstractVector{I}, val :: AbstractVector{$elty}; sym :: Bool=false) where I <: Integer
            nz = length(irn)
            @assert nz == length(jcn)
            @assert nz == length(val)
            err = $fname(spmat)
            qrm_check(err)
            spmat.irn = irn
            spmat.jcn = jcn
            spmat.val = val
            spmat.mat.m = m
            spmat.mat.n = n
            spmat.mat.nz  = nz
            spmat.mat.sym = sym
            spmat.mat.irn = pointer(spmat.irn)
            spmat.mat.jcn = pointer(spmat.jcn)
            spmat.mat.val = pointer(spmat.val)
            return nothing
        end

        function qrm_spmat_init(m :: Integer, n :: Integer, irn :: AbstractVector{I}, jcn :: AbstractVector{I}, val :: AbstractVector{$elty}; sym :: Bool=false) where I <: Integer
            spmat = qrm_spmat{$elty}()
            qrm_spmat_init!(spmat, m, n, irn, jcn, val; sym = sym)
            return spmat
        end
        
        function qrm_spmat_init!(spmat :: qrm_spmat{$elty}, A :: SparseMatrixCSC{$elty,I}; sym :: Bool=false) where I <: Integer
            qrm_spmat_init!(spmat, size(A)..., findnz(A)...; sym = sym)
        end

        function qrm_spmat_init(A :: SparseMatrixCSC{$elty,I}; sym :: Bool=false) where I <: Integer
            spmat = qrm_spmat{$elty}()
            qrm_spmat_init!(spmat, A, sym=sym)
            return spmat
        end

        @inline qrm_spmat_init!(spmat :: qrm_spmat{$elty}, A :: Symmetric{$elty, SparseMatrixCSC{$elty,I}}) where I <: Integer = qrm_spmat_init!(spmat, A.data, sym=true)
        @inline qrm_spmat_init!(spmat :: qrm_spmat{$elty}, A :: Hermitian{$elty, SparseMatrixCSC{$elty,I}}) where I <: Integer = qrm_spmat_init!(spmat, A.data, sym=true)

        @inline qrm_spmat_init(A :: Symmetric{$elty, SparseMatrixCSC{$elty,I}}) where I <: Integer = qrm_spmat_init(A.data, sym=true)
        @inline qrm_spmat_init(A :: Hermitian{$elty, SparseMatrixCSC{$elty,I}}) where I <: Integer = qrm_spmat_init(A.data, sym=true)
    end
end

for (fname, elty) in ((:sqrm_spmat_destroy_c, :Float32   ),
                      (:dqrm_spmat_destroy_c, :Float64   ),
                      (:cqrm_spmat_destroy_c, :ComplexF32),
                      (:zqrm_spmat_destroy_c, :ComplexF64))
    @eval begin
        function qrm_spmat_destroy!(spmat :: qrm_spmat{$elty})
            err = $fname(spmat)
            # qrm_check(err)
            return nothing
        end
    end
end

for (fname, elty) in ((:sqrm_spfct_init_c, :Float32   ),
                      (:dqrm_spfct_init_c, :Float64   ),
                      (:cqrm_spfct_init_c, :ComplexF32),
                      (:zqrm_spfct_init_c, :ComplexF64))
    @eval begin
        function qrm_spfct_init!(spfct :: qrm_spfct{$elty}, spmat :: qrm_spmat{$elty})
            err = $fname(spfct, spmat)
            qrm_check(err)
            return nothing
        end

        function qrm_spfct_init(spmat :: qrm_spmat{$elty})
            spfct = qrm_spfct{$elty}()
            qrm_spfct_init!(spfct, spmat)
            return spfct
        end
    end
end

for (fname, elty) in ((:sqrm_spfct_destroy_c, :Float32   ),
                      (:dqrm_spfct_destroy_c, :Float64   ),
                      (:cqrm_spfct_destroy_c, :ComplexF32),
                      (:zqrm_spfct_destroy_c, :ComplexF64))
    @eval begin
        function qrm_spfct_destroy!(spfct :: qrm_spfct{$elty})
            err = $fname(spfct)
            qrm_check(err)
            return nothing
        end
    end
end

for (fname, elty) in ((:sqrm_analyse_c, :Float32   ),
                      (:dqrm_analyse_c, :Float64   ),
                      (:cqrm_analyse_c, :ComplexF32),
                      (:zqrm_analyse_c, :ComplexF64))
    @eval begin
        function qrm_analyse!(spmat :: qrm_spmat{$elty}, spfct :: qrm_spfct{$elty}; transp :: Char='n')
            err = $fname(spmat, spfct, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_analyse(spmat :: qrm_spmat{$elty}; transp :: Char='n')
            spfct = qrm_spfct_init(spmat)
            err = $fname(spmat, spfct, transp)
            qrm_check(err)
            return spfct
        end

        @inline qrm_analyse!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, spfct :: qrm_spfct{$elty}) = qrm_analyse!(spmat.parent, spfct, transp='t')
        @inline qrm_analyse!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}  , spfct :: qrm_spfct{$elty}) = qrm_analyse!(spmat.parent, spfct, transp='c')

        @inline qrm_analyse(spmat :: Transpose{$elty,qrm_spmat{$elty}}) = qrm_analyse(spmat.parent, transp='t')
        @inline qrm_analyse(spmat :: Adjoint{$elty,qrm_spmat{$elty}})   = qrm_analyse(spmat.parent, transp='c')
    end
end

for (fname, elty) in ((:sqrm_factorize_c, :Float32   ),
                      (:dqrm_factorize_c, :Float64   ),
                      (:cqrm_factorize_c, :ComplexF32),
                      (:zqrm_factorize_c, :ComplexF64))
    @eval begin
        function qrm_factorize!(spmat :: qrm_spmat{$elty}, spfct :: qrm_spfct{$elty}; transp :: Char='n')
            err = $fname(spmat, spfct, transp)
            qrm_check(err)
            return nothing
        end

        @inline qrm_factorize!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, spfct :: qrm_spfct{$elty}) = qrm_factorize!(spmat.parent, spfct, transp='t')
        @inline qrm_factorize!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}  , spfct :: qrm_spfct{$elty}) = qrm_factorize!(spmat.parent, spfct, transp='c')
    end
end

for (fname, elty) in ((:sqrm_solve_c, :Float32   ),
                      (:dqrm_solve_c, :Float64   ),
                      (:cqrm_solve_c, :ComplexF32),
                      (:zqrm_solve_c, :ComplexF64))
    @eval begin
        function qrm_solve!(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}, x :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                @assert length(x) == spfct.fct.n
            else
                @assert length(x) == spfct.fct.m
            end
            err = $fname(spfct, transp, b, x, nrhs)
            qrm_check(err)
            return nothing
        end

        function qrm_solve!(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                @assert size(x) == (spfct.fct.n, nrhs)
            else
                @assert size(x) == (spfct.fct.m, nrhs)
            end
            err = $fname(spfct, transp, b, x, nrhs)
            qrm_check(err)
            return nothing
        end

        function qrm_solve(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                x = zeros($elty, spfct.fct.n)
            else
                x = zeros($elty, spfct.fct.m)
            end
            err = $fname(spfct, transp, b, x, nrhs)
            qrm_check(err)
            return x
        end

        function qrm_solve(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                x = zeros($elty, spfct.fct.n, nrhs)
            else
                x = zeros($elty, spfct.fct.m, nrhs)
            end
            err = $fname(spfct, transp, b, x, nrhs)
            qrm_check(err)
            return x
        end

        @inline qrm_solve!(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_solve!(spfct.parent, b, x, transp='t')
        @inline qrm_solve!(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_solve!(spfct.parent, b, x, transp='t')

        @inline qrm_solve(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_solve(spfct.parent, b, transp='t')
        @inline qrm_solve(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_solve(spfct.parent, b, transp='t')

        @inline qrm_solve!(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_solve!(spfct.parent, b, x, transp='c')
        @inline qrm_solve!(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_solve!(spfct.parent, b, x, transp='c')

        @inline qrm_solve(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_solve(spfct.parent, b, transp='c')
        @inline qrm_solve(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_solve(spfct.parent, b, transp='c')
    end
end

for (fname, elty) in ((:sqrm_apply_c, :Float32   ),
                      (:dqrm_apply_c, :Float64   ),
                      (:cqrm_apply_c, :ComplexF32),
                      (:zqrm_apply_c, :ComplexF64))
    @eval begin
        function qrm_apply!(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            err = $fname(spfct, transp, b, nrhs)
            qrm_check(err)
            return nothing
        end

        function qrm_apply!(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            err = $fname(spfct, transp, b, nrhs)
            qrm_check(err)
            return nothing
        end

        function qrm_apply(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            z = copy(b)
            err = $fname(spfct, transp, z, nrhs)
            qrm_check(err)
            return z
        end

        function qrm_apply(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            z = copy(b)
            err = $fname(spfct, transp, z, nrhs)
            qrm_check(err)
            return z
        end

        @inline qrm_apply!(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_apply!(spfct.parent, b, transp='t')
        @inline qrm_apply!(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_apply!(spfct.parent, b, transp='t')

        @inline qrm_apply(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_apply(spfct.parent, b, transp='t')
        @inline qrm_apply(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_apply(spfct.parent, b, transp='t')

        @inline qrm_apply!(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_apply!(spfct.parent, b, transp='c')
        @inline qrm_apply!(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_apply!(spfct.parent, b, transp='c')

        @inline qrm_apply(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_apply(spfct.parent, b, transp='c')
        @inline qrm_apply(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_apply(spfct.parent, b, transp='c')
    end
end

for (fname, elty) in ((:sqrm_spmat_mv_c, :Float32   ),
                      (:dqrm_spmat_mv_c, :Float64   ),
                      (:cqrm_spmat_mv_c, :ComplexF32),
                      (:zqrm_spmat_mv_c, :ComplexF64))
    @eval begin
        function qrm_spmat_mv!(spmat :: qrm_spmat{$elty}, alpha :: $elty, x :: Vector{$elty}, beta :: $elty, y :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            err = $fname(spmat, transp, alpha, x, beta, y, nrhs)
            qrm_check(err)
            return nothing
        end

        function qrm_spmat_mv!(spmat :: qrm_spmat{$elty}, alpha :: $elty, x :: Matrix{$elty}, beta :: $elty, y :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(y, 2)
            err = $fname(spmat, transp, alpha, x, beta, y, nrhs)
            qrm_check(err)
            return nothing
        end

        @inline qrm_spmat_mv!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, alpha :: $elty, x :: Vector{$elty}, beta :: $elty, y :: Vector{$elty}) = qrm_spmat_mv!(spmat.parent, alpha, x, beta, y, transp='t')
        @inline qrm_spmat_mv!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, alpha :: $elty, x :: Matrix{$elty}, beta :: $elty, y :: Matrix{$elty}) = qrm_spmat_mv!(spmat.parent, alpha, x, beta, y, transp='t')

        @inline qrm_spmat_mv!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, alpha :: $elty, x :: Vector{$elty}, beta :: $elty, y :: Vector{$elty}) = qrm_spmat_mv!(spmat.parent, alpha, x, beta, y, transp='c')
        @inline qrm_spmat_mv!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, alpha :: $elty, x :: Matrix{$elty}, beta :: $elty, y :: Matrix{$elty}) = qrm_spmat_mv!(spmat.parent, alpha, x, beta, y, transp='c')

        @inline mul!(y :: Vector{$elty}, spmat :: qrm_spmat{$elty}, x :: Vector{$elty}, alpha :: $elty, beta :: $elty) = qrm_spmat_mv!(spmat, alpha, x, beta, y)
        @inline mul!(y :: Matrix{$elty}, spmat :: qrm_spmat{$elty}, x :: Matrix{$elty}, alpha :: $elty, beta :: $elty) = qrm_spmat_mv!(spmat, alpha, x, beta, y)

        @inline mul!(y :: Vector{$elty}, spmat :: Transpose{$elty,qrm_spmat{$elty}}, x :: Vector{$elty}, alpha :: $elty, beta :: $elty) = qrm_spmat_mv!(spmat, alpha, x, beta, y)
        @inline mul!(y :: Matrix{$elty}, spmat :: Transpose{$elty,qrm_spmat{$elty}}, x :: Matrix{$elty}, alpha :: $elty, beta :: $elty) = qrm_spmat_mv!(spmat, alpha, x, beta, y)

        @inline mul!(y :: Vector{$elty}, spmat :: Adjoint{$elty,qrm_spmat{$elty}}, x :: Vector{$elty}, alpha :: $elty, beta :: $elty) = qrm_spmat_mv!(spmat, alpha, x, beta, y)
        @inline mul!(y :: Matrix{$elty}, spmat :: Adjoint{$elty,qrm_spmat{$elty}}, x :: Matrix{$elty}, alpha :: $elty, beta :: $elty) = qrm_spmat_mv!(spmat, alpha, x, beta, y)
    end
end

for (fname, elty, subty) in ((:sqrm_spmat_nrm_c, :Float32   , :Float32),
                             (:dqrm_spmat_nrm_c, :Float64   , :Float64),
                             (:cqrm_spmat_nrm_c, :ComplexF32, :Float32),
                             (:zqrm_spmat_nrm_c, :ComplexF64, :Float64))
    @eval begin
        function qrm_spmat_nrm(spmat :: qrm_spmat{$elty}; ntype :: Char='f')
            nrm = Ref{$subty}(0)
            err = $fname(spmat, ntype, nrm)
            qrm_check(err)
            return nrm[]
        end
    end
end

for (fname, elty, subty) in ((:sqrm_vecnrm_c, :Float32   , :Float32),
                             (:dqrm_vecnrm_c, :Float64   , :Float64),
                             (:cqrm_vecnrm_c, :ComplexF32, :Float32),
                             (:zqrm_vecnrm_c, :ComplexF64, :Float64))
    @eval begin
        function qrm_vecnrm(x :: Vector{$elty}; ntype :: Char='2')
            n = length(x)
            nrhs = 1
            nrm = Ref{$subty}(0)
            err = $fname(x, n, nrhs, ntype, nrm)
            qrm_check(err)
            return nrm[]
        end

        function qrm_vecnrm(x :: Matrix{$elty}; ntype :: Char='2')
            n, nrhs = size(x)
            nrm = zeros($subty, nrhs)
            err = $fname(x, n, nrhs, ntype, nrm)
            qrm_check(err)
            return nrm
        end

        function qrm_vecnrm!(x :: Matrix{$elty}, nrm :: Vector{$subty}; ntype :: Char='2')
            n, nrhs = size(x)
            err = $fname(x, n, nrhs, ntype, nrm)
            qrm_check(err)
            return nothing
        end
    end
end

for (fname, elty) in ((:sqrm_spbackslash_c, :Float32   ),
                      (:dqrm_spbackslash_c, :Float64   ),
                      (:cqrm_spbackslash_c, :ComplexF32),
                      (:zqrm_spbackslash_c, :ComplexF64))
    @eval begin
        function qrm_spbackslash!(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}, x :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                @assert length(x) == spmat.mat.n
            else
                @assert length(x) == spmat.mat.m
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_spbackslash!(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                @assert size(x) == (spmat.mat.n, nrhs)
            else
                @assert size(x) == (spmat.mat.m, nrhs)
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end


        function qrm_spbackslash(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                x = zeros($elty, spmat.mat.n)
            else
                x = zeros($elty, spmat.mat.m)
            end
            bcopy = (spmat.mat.m ≥ spmat.mat.n) ? copy(b) : b
            err = $fname(spmat, bcopy, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        function qrm_spbackslash(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                x = zeros($elty, spmat.mat.n, nrhs)
            else
                x = zeros($elty, spmat.mat.m, nrhs)
            end
            bcopy = (spmat.mat.m ≥ spmat.mat.n) ? copy(b) : b
            err = $fname(spmat, bcopy, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        @inline qrm_spbackslash!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_spbackslash!(spmat.parent, b, x, transp='t')
        @inline qrm_spbackslash!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_spbackslash!(spmat.parent, b, x, transp='t')

        @inline qrm_spbackslash(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spmat.parent, b, transp='t')
        @inline qrm_spbackslash(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spmat.parent, b, transp='t')

        @inline qrm_spbackslash!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_spbackslash!(spmat.parent, b, x, transp='c')
        @inline qrm_spbackslash!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_spbackslash!(spmat.parent, b, x, transp='c')

        @inline qrm_spbackslash(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spmat.parent, b, transp='c')
        @inline qrm_spbackslash(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spmat.parent, b, transp='c')

        
        @inline (\)(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}) = qrm_spbackslash(spmat, b, transp='n')
        @inline (\)(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}) = qrm_spbackslash(spmat, b, transp='n')

        @inline (\)(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spmat.parent, b, transp='t')
        @inline (\)(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spmat.parent, b, transp='t')

        @inline (\)(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spmat.parent, b, transp='c')
        @inline (\)(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spmat.parent, b, transp='c')
                 
                 
    end
end

for (fname, elty) in ((:sqrm_spfct_backslash_c, :Float32   ),
                      (:dqrm_spfct_backslash_c, :Float64   ),
                      (:cqrm_spfct_backslash_c, :ComplexF32),
                      (:zqrm_spfct_backslash_c, :ComplexF64))
    @eval begin
        function qrm_spbackslash!(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}, x :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                @assert length(x) == spfct.fct.n
            else
                @assert length(x) == spfct.fct.m
            end
            err = $fname(spfct, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_spbackslash!(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                @assert size(x) == (spfct.fct.n, nrhs)
            else
                @assert size(x) == (spfct.fct.m, nrhs)
            end
            err = $fname(spfct, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end


        function qrm_spbackslash(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                x = zeros($elty, spfct.fct.n)
            else
                x = zeros($elty, spfct.fct.m)
            end
            bcopy = (spfct.fct.m ≥ spfct.fct.n) ? copy(b) : b
            err = $fname(spfct, bcopy, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        function qrm_spbackslash(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                x = zeros($elty, spfct.fct.n, nrhs)
            else
                x = zeros($elty, spfct.fct.m, nrhs)
            end
            bcopy = (spfct.fct.m ≥ spfct.fct.n) ? copy(b) : b
            err = $fname(spfct, bcopy, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        @inline qrm_spbackslash!(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_spbackslash!(spfct.parent, b, x, transp='t')
        @inline qrm_spbackslash!(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_spbackslash!(spfct.parent, b, x, transp='t')

        @inline qrm_spbackslash(spfct  :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spfct.parent, b, transp='t')
        @inline qrm_spbackslash(spfct  :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spfct.parent, b, transp='t')

        @inline qrm_spbackslash!(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_spbackslash!(spfct.parent, b, x, transp='c')
        @inline qrm_spbackslash!(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_spbackslash!(spfct.parent, b, x, transp='c')

        @inline qrm_spbackslash(spfct  :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spfct.parent, b, transp='c')
        @inline qrm_spbackslash(spfct  :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spfct.parent, b, transp='c')

                 
        @inline (\)(spfct :: qrm_spfct{$elty}, b :: Vector{$elty}) = qrm_spbackslash(spfct, b)
        @inline (\)(spfct :: qrm_spfct{$elty}, b :: Matrix{$elty}) = qrm_spbackslash(spfct, b)

        @inline (\)(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spfct.parent, b, transp='t')
        @inline (\)(spfct :: Transpose{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spfct.parent, b, transp='t')

        @inline (\)(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Vector{$elty}) = qrm_spbackslash(spfct.parent, b, transp='c')
        @inline (\)(spfct :: Adjoint{$elty,qrm_spfct{$elty}}, b :: Matrix{$elty}) = qrm_spbackslash(spfct.parent, b, transp='c')
    end
end

for (fname, elty) in ((:sqrm_spposv_c, :Float32   ),
                      (:dqrm_spposv_c, :Float64   ),
                      (:cqrm_spposv_c, :ComplexF32),
                      (:zqrm_spposv_c, :ComplexF64))
    @eval begin
        function qrm_spposv!(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}, x :: Vector{$elty})
            nrhs = 1
            err = $fname(spmat, b, x, nrhs)
            qrm_check(err)
            return nothing
        end
        
        function qrm_spposv!(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty})
            nrhs = size(b, 2)
            err = $fname(spmat, b, x, nrhs)
            qrm_check(err)
            return nothing
        end

        function qrm_spposv(spmat :: qrm_spmat{$elty}, b :: Vector{$elty})
            nrhs = 1
            x = zeros($elty, spmat.mat.n)
            bcopy = copy(b)
            err = $fname(spmat, bcopy, x, nrhs)
            qrm_check(err)
            return x
        end

        function qrm_spposv(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty})
            nrhs = size(b, 2)
            x = zeros($elty, spmat.mat.n, nrhs)
            bcopy = copy(b)
            err = $fname(spmat, bcopy, x, nrhs)
            qrm_check(err)
            return x
        end
    end
end

for (fname, elty) in ((:sqrm_least_squares_c, :Float32   ),
                      (:dqrm_least_squares_c, :Float64   ),
                      (:cqrm_least_squares_c, :ComplexF32),
                      (:zqrm_least_squares_c, :ComplexF64))
    @eval begin
        function qrm_least_squares!(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}, x :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                @assert length(x) == spmat.mat.n
            else
                @assert length(x) == spmat.mat.m
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_least_squares!(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                @assert size(x) == (spmat.mat.n, nrhs)
            else
                @assert size(x) == (spmat.mat.m, nrhs)
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_least_squares(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                x = zeros($elty, spmat.mat.n)
            else
                x = zeros($elty, spmat.mat.m)
            end
            bcopy = copy(b)
            err = $fname(spmat, bcopy, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        function qrm_least_squares(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                x = zeros($elty, spmat.mat.n, nrhs)
            else
                x = zeros($elty, spmat.mat.n, nrhs)
            end
            bcopy = copy(b)
            err = $fname(spmat, bcopy, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        @inline qrm_least_squares!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_least_squares!(spmat.parent, b, x, transp='t')
        @inline qrm_least_squares!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_least_squares!(spmat.parent, b, x, transp='t')

        @inline qrm_least_squares(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_least_squares(spmat.parent, b, transp='t')
        @inline qrm_least_squares(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_least_squares(spmat.parent, b, transp='t')

        @inline qrm_least_squares!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_least_squares!(spmat.parent, b, x, transp='c')
        @inline qrm_least_squares!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_least_squares!(spmat.parent, b, x, transp='c')

        @inline qrm_least_squares(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_least_squares(spmat.parent, b, transp='c')
        @inline qrm_least_squares(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_least_squares(spmat.parent, b, transp='c')

    end
end

for (fname, elty) in ((:sqrm_min_norm_c, :Float32   ),
                      (:dqrm_min_norm_c, :Float64   ),
                      (:cqrm_min_norm_c, :ComplexF32),
                      (:zqrm_min_norm_c, :ComplexF64))
    @eval begin
        function qrm_min_norm!(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}, x :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                @assert length(x) == spmat.mat.n
            else
                @assert length(x) == spmat.mat.m
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_min_norm!(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                @assert size(x) == (spmat.mat.n, nrhs)
            else
                @assert size(x) == (spmat.mat.m, nrhs)
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return nothing
        end

        function qrm_min_norm(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            if transp == 'n'
                x = zeros($elty, spmat.mat.n)
            else
                x = zeros($elty, spmat.mat.m)
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        function qrm_min_norm(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(b, 2)
            if transp == 'n'
                x = zeros($elty, spmat.mat.n, nrhs)
            else
                x = zeros($elty, spmat.mat.n, nrhs)
            end
            err = $fname(spmat, b, x, nrhs, transp)
            qrm_check(err)
            return x
        end

        @inline qrm_min_norm!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_min_norm!(spmat.parent, b, x, transp='t')
        @inline qrm_min_norm!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_min_norm!(spmat.parent, b, x, transp='t')

        @inline qrm_min_norm(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_min_norm(spmat.parent, b, transp='t')
        @inline qrm_min_norm(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_min_norm(spmat.parent, b, transp='t')

        @inline qrm_min_norm!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_min_norm!(spmat.parent, b, x, transp='c')
        @inline qrm_min_norm!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_min_norm!(spmat.parent, b, x, transp='c')

        @inline qrm_min_norm(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}) = qrm_min_norm(spmat.parent, b, transp='c')
        @inline qrm_min_norm(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}) = qrm_min_norm(spmat.parent, b, transp='c')

    end
end

for (fname, elty, subty) in ((:sqrm_residual_norm_c, :Float32   , :Float32),
                             (:dqrm_residual_norm_c, :Float64   , :Float64),
                             (:cqrm_residual_norm_c, :ComplexF32, :Float32),
                             (:zqrm_residual_norm_c, :ComplexF64, :Float64))
    @eval begin
        function qrm_residual_norm(spmat :: qrm_spmat{$elty}, b :: Vector{$elty}, x :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            nrm = Ref{$subty}(0)
            err = $fname(spmat, b, x, nrhs, nrm, transp)
            qrm_check(err)
            return nrm[]
        end

        function qrm_residual_norm(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(x, 2)
            nrm = zeros($subty, nrhs)
            err = $fname(spmat, b, x, nrhs, nrm, transp)
            qrm_check(err)
            return nrm
        end

        function qrm_residual_norm!(spmat :: qrm_spmat{$elty}, b :: Matrix{$elty}, x :: Matrix{$elty}, nrm :: Vector{$subty}; transp :: Char='n')
            nrhs = size(x, 2)
            err = $fname(spmat, b, x, nrhs, nrm, transp)
            qrm_check(err)
            return nothing
        end

        @inline qrm_residual_norm!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}, nrm :: Vector{$subty}) = qrm_residual_norm!(spmat.parent, b, x, transp='t')
        @inline qrm_residual_norm!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}  , b :: Matrix{$elty}, x :: Matrix{$elty}, nrm :: Vector{$subty}) = qrm_residual_norm!(spmat.parent, b, x, transp='c')

        @inline qrm_residual_norm(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_residual_norm(spmat.parent, b, x, transp='t')
        @inline qrm_residual_norm(spmat :: Transpose{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_residual_norm(spmat.parent, b, x, transp='t')

        @inline qrm_residual_norm(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Vector{$elty}, x :: Vector{$elty}) = qrm_residual_norm(spmat.parent, b, x, transp='c')
        @inline qrm_residual_norm(spmat :: Adjoint{$elty,qrm_spmat{$elty}}, b :: Matrix{$elty}, x :: Matrix{$elty}) = qrm_residual_norm(spmat.parent, b, x, transp='c')


    end
end

for (fname, elty, subty) in ((:sqrm_residual_orth_c, :Float32   , :Float32),
                             (:dqrm_residual_orth_c, :Float64   , :Float64),
                             (:cqrm_residual_orth_c, :ComplexF32, :Float32),
                             (:zqrm_residual_orth_c, :ComplexF64, :Float64))
    @eval begin
        function qrm_residual_orth(spmat :: qrm_spmat{$elty}, r :: Vector{$elty}; transp :: Char='n')
            nrhs = 1
            nrm = Ref{$subty}(0)
            err = $fname(spmat, r, nrhs, nrm, transp)
            qrm_check(err)
            return nrm[]
        end

        function qrm_residual_orth(spmat :: qrm_spmat{$elty}, r :: Matrix{$elty}; transp :: Char='n')
            nrhs = size(r, 2)
            nrm = zeros($subty, nrhs)
            err = $fname(spmat, r, nrhs, nrm, transp)
            qrm_check(err)
            return nrm
        end

        function qrm_residual_orth!(spmat :: qrm_spmat{$elty}, r :: Matrix{$elty}, nrm :: Vector{$subty}; transp :: Char='n')
            nrhs = size(r, 2)
            err = $fname(spmat, r, nrhs, nrm, transp)
            qrm_check(err)
            return nothing
        end

        @inline qrm_residual_orth!(spmat :: Transpose{$elty,qrm_spmat{$elty}}, r :: Matrix{$elty}, nrm :: Vector{$subty}) = qrm_residual_orth!(spmat.parent, r, transp='t')
        @inline qrm_residual_orth!(spmat :: Adjoint{$elty,qrm_spmat{$elty}}  , r :: Matrix{$elty}, nrm :: Vector{$subty}) = qrm_residual_orth!(spmat.parent, r, transp='c')

        @inline qrm_residual_orth(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, r :: Vector{$elty}) = qrm_residual_orth(spmat.parent, r, transp='t')
        @inline qrm_residual_orth(spmat  :: Transpose{$elty,qrm_spmat{$elty}}, r :: Matrix{$elty}) = qrm_residual_orth(spmat.parent, r, transp='t')

        @inline qrm_residual_orth(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}  , r :: Vector{$elty}) = qrm_residual_orth(spmat.parent, r, transp='c')
        @inline qrm_residual_orth(spmat  :: Adjoint{$elty,qrm_spmat{$elty}}  , r :: Matrix{$elty}) = qrm_residual_orth(spmat.parent, r, transp='c')
    end
end

function qrm_set(str :: String, val :: Number)
    if (str ∈ GICNTL) || (str ∈ PICNTL)
        err = qrm_glob_set_i4_c(str, val)
    elseif str ∈ RCNTL
        err = qrm_glob_set_r4_c(str, val)
    else
        err = Int32(23)
    end
    qrm_check(err)
    return nothing
end

for (finame, frname, elty) in ((:sqrm_spfct_set_i4_c, :sqrm_spfct_set_r4_c, :Float32   ),
                               (:dqrm_spfct_set_i4_c, :dqrm_spfct_set_r4_c, :Float64   ),
                               (:cqrm_spfct_set_i4_c, :cqrm_spfct_set_r4_c, :ComplexF32),
                               (:zqrm_spfct_set_i4_c, :zqrm_spfct_set_r4_c, :ComplexF64))
    @eval begin
        function qrm_set(spfct :: qrm_spfct{$elty}, str :: String, val :: Number)
            if str ∈ PICNTL
                err = $finame(spfct, str, val)
            elseif str ∈ RCNTL
                err = $frname(spfct, str, val)
            else
                err = Int32(23)
            end
            qrm_check(err)
            return nothing
        end
    end
end

function qrm_get(str :: String)
    if (str ∈ GICNTL) || (str ∈ PICNTL)
        val = Ref{Clonglong}(0)
        err = qrm_glob_get_i8_c(str, val)
    elseif str ∈ RCNTL
        val = Ref{Float32}(0)
        err = qrm_glob_get_r4_c(str, val)
    else
        err = Int32(23)
    end
    qrm_check(err)
    return val[]
end

for (finame, frname, elty) in ((:sqrm_spfct_get_i8_c, :sqrm_spfct_get_r4_c, :Float32   ),
                               (:dqrm_spfct_get_i8_c, :dqrm_spfct_get_r4_c, :Float64   ),
                               (:cqrm_spfct_get_i8_c, :cqrm_spfct_get_r4_c, :ComplexF32),
                               (:zqrm_spfct_get_i8_c, :zqrm_spfct_get_r4_c, :ComplexF64))
    @eval begin
        function qrm_get(spfct :: qrm_spfct{$elty}, str :: String)
            if (str ∈ PICNTL) || (str ∈ STATS)
                val = Ref{Clonglong}(0)
                err = $finame(spfct, str, val)
            elseif str ∈ RCNTL
                val = Ref{Float32}(0)
                err = $frname(spfct, str, val)
            else
                err = Int32(23)
            end
            qrm_check(err)
            return val[]
        end
    end
end

function qrm_init(ncpu :: Integer = 1, ngpu :: Integer = 0)
    err = qrm_init_c(ncpu, ngpu)
    qrm_check(err)
    return nothing
end

function qrm_finalize()
    qrm_finalize_c()
end

function qrm_update!(spmat :: qrm_spmat{T}, val :: AbstractVector{T}) where T
    spmat.val .= val
    return nothing
end

function qrm_update!(spmat :: qrm_spmat{T}, A :: SparseMatrixCSC{T,I}) where {T, I <: Integer}
    qrm_update!(spmat, A.nzval)
end
