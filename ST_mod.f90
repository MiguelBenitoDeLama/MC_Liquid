module ST_mod   !Author: Miguel Benito de Lama (January 2023)
use general_mod
implicit none


contains


    subroutine ST_potential_diff(ini_coor, trial_coor, hL, L, r_cut, N, atom, delta_v)
        !Calculates the Stillinger potential energy difference of two given structures
        implicit none
        integer(kind=precision),intent(in):: N, atom
        integer(kind=precision):: i, j, k
        real(kind=precision),intent(in)::ini_coor(:,:), trial_coor(:,:), hL, L, r_cut
        real(kind=precision):: result, r_atom_i_2,r_atom_k_2, r_i_k_2, theta_jik, theta_ijk
        real(kind=precision),intent(out)::delta_v

        delta_v = 0

        do i=1,N
            if (i .ne. atom) then
                call distance_squared(ini_coor(:,atom),ini_coor(:,i), r_atom_i_2, hL, L) 
                if (r_atom_i_2 .lt. r_cut**2.) then
                    call f_2(r_atom_i_2,r_cut,result)
                    delta_v=delta_v-result
                    do k=i+1,N
                        if (k .ne. atom) then
                            call distance_squared(ini_coor(:,atom),ini_coor(:,k), r_atom_k_2, hL, L) 
                            if (r_atom_k_2 .lt. r_cut**2.) then
                                call angle_jik(ini_coor(:,atom),ini_coor(:,i), ini_coor(:,k), hL, L, theta_jik) ! angle theta_jik
                                call h_f_3(r_atom_i_2, r_atom_k_2, r_cut, theta_jik, result)
                                delta_v=delta_v-result
                            endif
                        endif
                    enddo

                    do k=1,N
                        if ((k .ne. atom) .and. (k .ne. i)) then
                            call distance_squared(ini_coor(:,i),ini_coor(:,k), r_i_k_2, hL, L) ! distance rjk
                            if (r_i_k_2 .lt. r_cut**2.) then
                                call angle_jik(ini_coor(:,i),ini_coor(:,atom), ini_coor(:,k), hL, L, theta_ijk) ! angle theta_ijk
                                call h_f_3(r_atom_i_2, r_i_k_2, r_cut, theta_ijk, result)
                                delta_v=delta_v-result
                            endif
                        endif
                    enddo
                endif

                call distance_squared(trial_coor(:,atom),trial_coor(:,i), r_atom_i_2, hL, L) 
                if (r_atom_i_2 .lt. r_cut**2.) then
                    call f_2(r_atom_i_2,r_cut,result)
                    delta_v=delta_v+result
                    do k=i+1,N
                        if (k .ne. atom) then
                            call distance_squared(trial_coor(:,atom),trial_coor(:,k), r_atom_k_2, hL, L) 
                            if (r_atom_k_2 .lt. r_cut**2.) then
                                call angle_jik(trial_coor(:,atom),trial_coor(:,i), trial_coor(:,k), hL, L, theta_jik) ! angle theta_jik
                                call h_f_3(r_atom_i_2, r_atom_k_2, r_cut, theta_jik, result)
                                delta_v=delta_v+result
                            endif
                        endif
                    enddo

                    do k=1,N
                        if ((k .ne. atom) .and. (k .ne. i)) then
                            call distance_squared(trial_coor(:,i),trial_coor(:,k), r_i_k_2, hL, L) ! distance rjk
                            if (r_i_k_2 .lt. r_cut**2.) then
                                call angle_jik(trial_coor(:,i),trial_coor(:,atom), trial_coor(:,k), hL, L, theta_ijk) ! angle theta_ijk
                                call h_f_3(r_atom_i_2, r_i_k_2, r_cut, theta_ijk, result)
                                delta_v=delta_v+result
                            endif
                        endif
                    enddo
                endif
            endif
        enddo
    end subroutine ST_potential_diff


    subroutine f_2(r_ij_2,r_cut,result)
        !Calculates the f_2 function used to obtain the V2 term of the potential
        implicit none
        real(kind=precision), intent(in)::r_cut, r_ij_2
        real(kind=precision), parameter:: A=7.049556277, B=0.6022245584
        real(kind=precision),intent(out):: result

        result=A* ((B * r_ij_2**(-2.)) - 1) * EXP(1 / (sqrt(r_ij_2) - r_cut))
    end subroutine f_2


    subroutine angle_jik(atom_i, atom_j, atom_k, hL,L, theta)
        !Calculates the angle jik between three given atoms, atom_i, atom_j and atom_k
        implicit none
        real(kind=precision), intent(in):: atom_i(3,1), atom_j(3,1), atom_k(3,1),hL,L
        real(kind=precision)::prod_v_ij_v_ik, mod_v_ij, mod_v_ik, vec_ij(3), vec_ik(3)
        real(kind=precision),intent(out):: theta
        integer(kind=precision)::i
        
        prod_v_ij_v_ik=0
        mod_v_ij=0
        mod_v_ik=0

        do i=1,3
            vec_ij(i)=atom_j(i,1)-atom_i(i,1)
            if (vec_ij(i) .gt. hL) then
                vec_ij(i)=vec_ij(i)-L
            elseif(vec_ij(i) .lt. -hL) then
                vec_ij(i)=vec_ij(i)+L
            endif
            vec_ik(i)=atom_k(i,1)-atom_i(i,1)
            if (vec_ik(i) .gt. hL) then
                vec_ik(i)=vec_ik(i)-L
            elseif(vec_ik(i) .lt. -hL) then
                vec_ik(i)=vec_ik(i)+L
            endif
            mod_v_ij=mod_v_ij+vec_ij(i)**(2.)
            mod_v_ik=mod_v_ik+vec_ik(i)**(2.)
            prod_v_ij_v_ik=prod_v_ij_v_ik+vec_ij(i)*vec_ik(i)
        enddo

        theta=ACOS((prod_v_ij_v_ik)/(sqrt(mod_v_ij)*sqrt(mod_v_ik)))
    end subroutine angle_jik

    
    subroutine h_f_3(r_ij_2, r_ik_2, r_cut, theta, result)
        !Calculates the h function used to obtain the V3 term of the potential
        implicit none
        real(kind=precision), intent(in):: r_ij_2, r_ik_2, r_cut, theta
        real(kind=precision), parameter:: lambda=21.0, gamma=1.20
        real(kind=precision):: exp_term, cos_term
        real(kind=precision),intent(out):: result

        exp_term = (gamma/(sqrt(r_ij_2)-r_cut))+(gamma/(sqrt(r_ik_2)-r_cut))
        cos_term = (COS(theta)+(1./3.))**(2.)
        result=lambda *EXP(exp_term) * cos_term

    end subroutine h_f_3


    subroutine ST_potential_diff_NL(ini_coor, trial_coor, neighbour_matrix, hL, L, r_cut, atom, n_max, delta_v)
        !Calculates the Stillinger potential energy difference of two given structures using neighbour lists
        implicit none
        integer(kind=precision),intent(in):: atom,n_max
        integer,intent(in)::neighbour_matrix(:,:)
        integer(kind=precision):: i, j, k
        real(kind=precision),intent(in)::ini_coor(:,:), trial_coor(:,:), hL, L, r_cut
        real(kind=precision)::f_2_result, h_f_3_result,r_atom_i_2,r_atom_k_2, r_i_k_2, theta_jik, theta_ijk
        real(kind=precision),intent(out)::delta_v

        delta_v = 0
        
        do i=1,n_max
            if (neighbour_matrix(i,atom).eq. 0) exit
            call distance_squared(ini_coor(:,atom),ini_coor(:,neighbour_matrix(i,atom)), r_atom_i_2, hL, L) 
            if (r_atom_i_2 .lt. r_cut**2.) then
                call f_2(r_atom_i_2,r_cut,f_2_result)
                delta_v=delta_v-f_2_result !V2 term
                do k=i+1,n_max
                    if (neighbour_matrix(k,atom).eq. 0) exit
                    call distance_squared(ini_coor(:,atom),ini_coor(:,neighbour_matrix(k,atom)), r_atom_k_2, hL, L) 
                    if (r_atom_k_2 .lt. r_cut**2.) then
                        call angle_jik(ini_coor(:,atom),ini_coor(:,neighbour_matrix(i,atom)), &
                        ini_coor(:,neighbour_matrix(k,atom)), hL, L, theta_jik) 
                        call h_f_3(r_atom_i_2, r_atom_k_2, r_cut, theta_jik, h_f_3_result)
                        delta_v=delta_v-h_f_3_result !V3 term
                    endif
                enddo

                do k=1,n_max
                    if (neighbour_matrix(k,neighbour_matrix(i,atom)).eq.0) exit
                    if (neighbour_matrix(k,neighbour_matrix(i,atom)) .ne. atom) then

                        call distance_squared(ini_coor(:,neighbour_matrix(i,atom)),&
                        ini_coor(:,neighbour_matrix(k,neighbour_matrix(i,atom))), r_i_k_2, hL, L)

                        if (r_i_k_2 .lt. r_cut**2.) then
                        
                            call angle_jik(ini_coor(:,neighbour_matrix(i,atom)),ini_coor(:,atom), &
                            ini_coor(:,neighbour_matrix(k,neighbour_matrix(i,atom))), hL, L, theta_ijk)

                            call h_f_3(r_atom_i_2, r_i_k_2, r_cut, theta_ijk, h_f_3_result)

                            delta_v=delta_v-h_f_3_result !V3 term
                        endif
                    endif
                enddo
            endif

            call distance_squared(trial_coor(:,atom),trial_coor(:,neighbour_matrix(i,atom)), r_atom_i_2, hL, L) 
            if (r_atom_i_2 .lt. r_cut**2.) then
                call f_2(r_atom_i_2,r_cut,f_2_result)
                delta_v=delta_v+f_2_result !V2 term
                do k=i+1,n_max
                    if (neighbour_matrix(k,atom).eq. 0) exit
                    call distance_squared(trial_coor(:,atom),trial_coor(:,neighbour_matrix(k,atom)), r_atom_k_2, hL, L) 
                    if (r_atom_k_2 .lt. r_cut**2.) then
                        call angle_jik(trial_coor(:,atom),trial_coor(:,neighbour_matrix(i,atom)),&
                            trial_coor(:,neighbour_matrix(k,atom)), hL, L, theta_jik) 
                        call h_f_3(r_atom_i_2, r_atom_k_2, r_cut, theta_jik, h_f_3_result)
                        delta_v=delta_v+h_f_3_result !V3 term
                    endif
                enddo

                do k=1,n_max
                    if (neighbour_matrix(k,neighbour_matrix(i,atom)).eq.0) exit
                    if (neighbour_matrix(k,neighbour_matrix(i,atom)) .ne. atom) then
                        call distance_squared(trial_coor(:,neighbour_matrix(i,atom)),&
                        trial_coor(:,neighbour_matrix(k,neighbour_matrix(i,atom))), r_i_k_2, hL, L)
                        if (r_i_k_2 .lt. r_cut**2.) then
                            call angle_jik(trial_coor(:,neighbour_matrix(i,atom)),trial_coor(:,atom), &
                            trial_coor(:,neighbour_matrix(k,neighbour_matrix(i,atom))), hL, L, theta_ijk)
                            call h_f_3(r_atom_i_2, r_i_k_2, r_cut, theta_ijk, h_f_3_result)
                            delta_v=delta_v+h_f_3_result !V3 term
                        endif
                    endif
                enddo
            endif
        enddo
    end subroutine ST_potential_diff_NL


end module ST_mod