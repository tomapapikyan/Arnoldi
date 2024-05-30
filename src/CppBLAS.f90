
module cppBLASwrapper
    use :: iso_c_binding
    implicit none
    contains

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
        deallocate(wrk)
    end subroutine
end module CppBLASwrapper
