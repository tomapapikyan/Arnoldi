
module cppBLASwrapper
    use :: iso_c_binding
    implicit none
    contains

    subroutine dgemv_wrapper(N, a, q_prev, q) bind (C, name = "dgemv_cpp")

        type(C_PTR), value :: a, q, q_prev
        real(8), pointer :: A_Matrix(:, :), vec(:), vec_prev(:)
        integer(C_INT), value :: N
        call C_F_POINTER(a, A_Matrix, [N, N])
        call C_F_POINTER(q_prev, vec_prev, [N])
        call C_F_POINTER(q, vec, [N])
        call dgemv('N', N, N, 1.d0, A_Matrix, N, vec_prev, 1, 0.d0, vec, 1)
    end subroutine

    real(8) function ddot_wrapper(N, v, w) bind (C, name = "ddot_cpp")

        real(8) :: ddot
        type(C_PTR), value :: v, w
        integer(C_INT), value :: N        
        real(8), pointer :: x(:), y(:)

        call C_F_POINTER(v, x, [N])
        call C_F_POINTER(w, y, [N])

        ddot_wrapper = ddot(n, x, 1, y, 1)

    end function

    subroutine daxpy_wrapper(N, a, x, y) bind (C, name = "daxpy_cpp")

        type(C_PTR), value :: x, y
        integer(C_INT), value :: N
        real(C_DOUBLE), value :: a
        real(8), pointer :: v(:), w(:)

        call C_F_POINTER(x, v, [N])
        call C_F_POINTER(y, w, [N])

        call daxpy(N, a, v, 1, w, 1)
    end subroutine

    real(8) function dnrm2_wrapper(N, x) bind (C, name = "dnrm2_cpp")
        real(8)::dnrm2
        type(C_PTR), value :: x
        integer(C_INT), value :: N
        real(8), pointer :: v(:)

        call C_F_POINTER(x, v, [N])
        dnrm2_wrapper = dnrm2(N, v, 1)
    end function

    subroutine dscal_wrapper(N, a, x) bind (C, name = "dscal_cpp")

        integer(C_INT), value :: N        
        type(C_PTR), value :: x
        real(c_DOUBLE), value :: a
        real(8), pointer :: v(:)

        call C_F_POINTER(x, v, [N])
        call dscal(N, a, v, 1)        
    end subroutine

    subroutine dhseqr_wrapper(N, H, wr, wi) bind (C, name = "dhseqr_cpp")
        type(C_PTR), value :: H, wr, wi
        integer(C_INT), value :: N
        real(8), pointer :: H_Matrix(:, :), wr_vec(:), wi_vec(:)
        integer :: info
        real(8), dimension(:), allocatable :: wrk

        allocate(wrk(N))
        call C_F_POINTER(H, H_Matrix, [N, N])
        call C_F_POINTER(wr, wr_vec, [N])
        call C_F_POINTER(wi, wi_vec, [N])

        call dhseqr('E', 'N', N, 1, N, H_Matrix, N, wr_vec, wi_vec, wi, 1, wrk, N, info)
    end subroutine
end module CppBLASwrapper
