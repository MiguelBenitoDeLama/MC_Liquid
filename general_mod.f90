module general_mod  !Author: Miguel Benito de Lama (January 2023)
implicit none
integer,parameter::precision=8 !SETS THE PRECISION OF THE PROGRAM


contains


    subroutine MC_movement(initial_coor,r_max,L,N,atom,trial_coor)
        !Changes at random the position of one atom, generating the trial coordinates
        implicit none
        real(kind=precision),intent(in)::initial_coor(:,:),r_max,L
        integer(kind=precision), intent(in)::N
        integer(kind=precision)::m
        real(kind=precision)::r,aux_coor
        integer(kind=precision), intent(out)::atom
        real(kind=precision),intent(out)::trial_coor(:,:)

        call random_number(r)
        atom = FLOOR( N * r ) + 1 !from 1 to N, takes one atom index, that is going to change its position    
        trial_coor = initial_coor
        do m = 1, 3 !Changes at random the position of the atom
            call random_number(r)
            aux_coor = trial_coor( m, atom ) + 2 * r_max * ( r - 0.5 ) 
            if (aux_coor .gt. L) then ! If the new position is greater than L 
                trial_coor(m,atom) = aux_coor - L
            elseif (aux_coor .lt. 0) then ! If the new position is lower than 0
                trial_coor(m,atom) = aux_coor + L
            else ! If the new position is inside the box
                trial_coor(m,atom) = trial_coor(m,atom) + 2 * r_max * ( r - 0.5) 
            endif
        enddo
    end subroutine MC_movement


    subroutine distance_squared(atom_i,atom_j,r_sq,hL,L)
        !Calculates the squared distance between two given atoms, taking
        !into account the boundary conditions
        implicit none
        real(kind=precision), intent(in):: atom_i(3,1), atom_j(3,1),hL,L
        integer(kind=precision)::k
        real(kind=precision)::dist(3)
        real(kind=precision),intent(out)::r_sq

        r_sq = 0
        dist = 0
        do k = 1, 3
            dist(k) = atom_i(k,1) - atom_j(k,1)
            if (dist(k) .gt. hL) then
                dist(k) = dist(k) - L
            else if(dist(k) .lt. -hL) then
                dist(k) = dist(k) + L
            endif
            r_sq = r_sq + dist( k ) ** ( 2. )
        enddo
    end subroutine distance_squared


    subroutine distances_matrix(coor,N,hL,L,dist_mat)
        !Calculates the atom-atom distance matrix for all the atoms
        implicit none
        real(kind=precision), intent(in):: coor(:,:),hL,L
        integer(kind=precision), intent(in):: N
        real(kind=precision)::r_ij
        real(kind=precision), intent(inout)::dist_mat(N,N)
        integer(kind=precision):: i, j

        do i = 1, N
            do j = i+1, N
                call distance_squared(coor(:,i),coor(:,j),r_ij,hL,L)
                dist_mat(i,j) = r_ij
                dist_mat(j,i) = dist_mat(i,j)
            enddo
        enddo
    end subroutine distances_matrix


    subroutine update_distances_matrix(coor, atom ,N,hL,L,dist_mat)
        !Updates the atom-atom distance matrix for a given atom
        implicit none
        real(kind=precision), intent(in):: coor(:,:),hL,L
        integer(kind=precision), intent(in):: N, atom
        real(kind=precision)::r_ij
        real(kind=precision), intent(inout)::dist_mat(N,N)
        integer(kind=precision):: i

        do i = 1, N
            call distance_squared(coor(:,i),coor(:,atom),r_ij,hL,L)
            dist_mat(i,atom) = r_ij
            dist_mat(atom,i) = dist_mat(i,atom)
        enddo
    end subroutine update_distances_matrix


    subroutine decision(delta_v, temp, atom, initial_coor, trial_coor,success_r)
        !Decides whether the trial structure is accepted or not depending on the potential energy difference 
        implicit none
        integer(kind=precision),intent(in)::atom
        real(kind=precision), intent(in)::delta_v, temp, trial_coor(:,:)
        real(kind=precision), intent(inout)::initial_coor(:,:),success_r
        real(kind=precision):: alpha, r

        if (delta_v .le. 0.) then  !If negative or 0 the new structure is accepted
            initial_coor(:,atom) = trial_coor(:,atom)
            success_r = success_r + 1
        else !If positive the structure is accepted or not depending on the Boltzmann factor
            alpha = exp( - delta_v / temp )
            call random_number(r)
            if (alpha .ge. r) then
                initial_coor(:,atom) = trial_coor(:,atom)
                success_r = success_r + 1
            else
                continue
            endif
        endif
    end subroutine decision


    subroutine radial_distribution(N,r_points,g_rad,rho,N_accum,N_hist)
        !Calculates the radial distribution function
        implicit none
        real(kind=precision), intent(in):: rho
        integer(kind=precision), intent(in)::r_points,N,N_hist(:),N_accum
        real(kind=precision)::div
        real(kind=precision), intent(out)::g_rad(2,r_points)
        integer(kind=precision):: i
	
	div=real(N_accum)*real(N)/2.*4./3.*4.D0 * DATAN(1.D0)*rho
        do i = 1, r_points
            g_rad(2,i) = real(N_hist(i))/(div*((g_rad(1,i) + g_rad(1,2))**3.-g_rad(1,i)**3.))
        enddo
    end subroutine radial_distribution


    subroutine N_histogram(dist_mat, N,r_points,delta_r,N_accum,N_hist)
        !Calculates the histogram used for the calculations of g(r)
        implicit none
        real(kind=precision), intent(in):: delta_r,dist_mat(:,:)
        integer(kind=precision), intent(in)::r_points,N
        integer(kind=precision), intent(inout)::N_hist(:),N_accum
        integer(kind=precision):: i, j, k
        
        if (N_accum .eq. 0) N_hist=0
        do i = 1, r_points
            do j = 1, N
                do k = j+1, N
                    if ((sqrt(dist_mat(j,k)) .gt. (delta_r*(i-1))) .and. (sqrt(dist_mat(j,k)) .le. (delta_r*i))) then
                        N_hist(i) = N_hist(i) + 1
                    endif
                enddo
            enddo
        enddo
        N_accum = N_accum + 1
    end subroutine N_histogram


    subroutine neighbour_list(dist_mat, neighbour_mat, r_n_max, atom, N)
        !Calculates the neighbour list of a given atom (updates the column of the matrix for that atom)
        implicit none
        real(kind=precision), intent(in):: dist_mat(:,:),r_n_max
        integer, intent(inout)::neighbour_mat(:,:)
        integer(kind=precision):: i, counter
        integer(kind=precision), intent(in):: N, atom

        neighbour_mat(:,atom) = 0

        counter = 0
        do i = 1, N
            if ((sqrt(dist_mat(i,atom)) .le. r_n_max) .and. (i .ne. atom)) then
                counter = counter + 1
                neighbour_mat(counter,atom) = i
            endif
        enddo
    end subroutine neighbour_list


    subroutine gnuplot_radial()
        !Plots the radial distribution function using gnuplot
        open(19,file='radial_plot.plt')                        
        write(19,*)"set nokey"                             
        write(19,*)"set grid"                             
        write(19,*)'set xlabel "r"'                       
        write(19,*)'set ylabel "g(r)"'                        
        write(19,*)"set title 'Radial Distribution function g(r)'"                    
        write(19,*)'plot "radial.dat" u 1:2 w l'
        close(19)
        call system('gnuplot -p radial_plot.plt')
    end subroutine gnuplot_radial


end module general_mod