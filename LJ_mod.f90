module LJ_mod   !Author: Miguel Benito de Lama (January 2023)
use general_mod
implicit none


contains


    subroutine LJ_potential_diff(coor, trial, atom, delta_v,r_cut,N,hL,L)
        !Calculates the Lennard Jones potential energy difference of two given structures
        implicit none
        integer(kind=precision),intent(in)::atom,N
        integer(kind=precision)::k
        real(kind=precision),intent(in)::coor(:,:), trial(:,:),r_cut,hL,L
        real(kind=precision)::r_ij_initial,r_ij_new
        real(kind=precision),intent(out)::delta_v
        integer(kind=precision):: i

        delta_v = 0
        do i = 1, N
            if (i .ne. atom) then
                r_ij_initial=0
                r_ij_new=0

                call distance_squared(coor(:,i),coor(:,atom),r_ij_initial,hL,L)
                call distance_squared(trial(:,i),trial(:,atom),r_ij_new,hL,L)

                if (r_ij_initial .lt. r_cut**(2.)) then
                    if(r_ij_new .lt. r_cut**(2.)) then
                        delta_v=delta_v+(((1/r_ij_new)**(6.)-(1/r_ij_new)**3.)-((1/r_ij_initial)**(6.)-(1/r_ij_initial)**3.))
                    else
                        delta_v=delta_v+(-((1/r_ij_initial)**(6.)-(1/r_ij_initial)**3.))
                    endif
                else
                    if(r_ij_new .lt. r_cut**(2.)) then
                        delta_v=delta_v+(((1/r_ij_new)**(6.)-(1/r_ij_new)**3.))
                    else
                        delta_v=delta_v
                    endif
                endif
            endif
        enddo
        delta_v=4*delta_v
    end subroutine LJ_potential_diff


    subroutine LJ_potential_diff_NL(coor, trial, neigh_mat_coor, atom, n_max, delta_v, r_cut, hL, L)
        !Calculates the Lennard Jones potential energy difference of two given structures using neighbour lists
        implicit none
        integer(kind=precision),intent(in)::atom, n_max
        integer,intent(in)::neigh_mat_coor(:,:)
        integer(kind=precision)::i
        real(kind=precision),intent(in)::coor(:,:), trial(:,:),r_cut,hL, L
        real(kind=precision)::r_ij_initial,r_ij_new
        real(kind=precision),intent(out)::delta_v


        delta_v=0
        do i=1,n_max
            if (neigh_mat_coor(i,atom).eq. 0) exit

            r_ij_initial=0
            r_ij_new=0

            call distance_squared(coor(:,neigh_mat_coor(i,atom)),coor(:,atom),r_ij_initial,hL,L)
            call distance_squared(trial(:,neigh_mat_coor(i,atom)),trial(:,atom),r_ij_new,hL,L)

            if (r_ij_initial .lt. r_cut**(2.)) then
                if(r_ij_new .lt. r_cut**(2.)) then
                    delta_v=delta_v+(((1/r_ij_new)**(6.)-(1/r_ij_new)**3.)-((1/r_ij_initial)**(6.)-(1/r_ij_initial)**3.))
                else
                    delta_v=delta_v+(-((1/r_ij_initial)**(6.)-(1/r_ij_initial)**3.))
                endif
            else
                if(r_ij_new .lt. r_cut**(2.)) then
                    delta_v=delta_v+(((1/r_ij_new)**(6.)-(1/r_ij_new)**3.))
                else
                    delta_v=delta_v
                endif
            endif
        enddo
        delta_v=4*delta_v
    end subroutine LJ_potential_diff_NL

end module LJ_mod