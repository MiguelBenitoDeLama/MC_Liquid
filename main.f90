program main         !Author: Miguel Benito de Lama (January 2023)
use general_mod             !Contains general use subroutines for the program
use LJ_mod                  !Contains the subroutines related with the Lennard Jones potential
use ST_mod                  !Contains the subroutines related with the Stillinger potential
implicit none
real :: time_1, time_2
integer ( kind = precision ) :: i, j, N, iterations, k, r_points, it, equil_iter, prod_iter, N_accum, atom, n_max
real ( kind = precision ) :: r, delta_v, L, hL, success_r, T, r_max, r_cut, rho, max_radial, delta_r
real ( kind = precision ) :: r_min, r_dist, r_skin, r_n
real ( kind = precision ), allocatable :: initial_coor ( :, : ), trial_coor ( :, : ), g_rad ( :, : )
real ( kind = precision ), allocatable :: dist_mat ( :, : ), positions_vec ( :, : )
integer, allocatable :: neighbour_mat ( :, : )
integer ( kind = precision ), allocatable :: N_hist ( : )
character ( len = 2 ) :: mode, trash
character ( len = 30 ) :: structure_file
logical :: control, NL_mode, structure_given

!________#_#_#_#_#_#________ Parameters setup and initial structure ________#_#_#_#_#_#________

call cpu_time ( time_1 )    !Stores current time in time_1, in order to measure the cpu time of the calculation


!Reads the input parameters from the input file "input.in"
open( unit = 17, file = 'input.in', action = 'read') 

read (17, *)
read (17,'(t37,I8)') N                  !Number of atoms
read (17,'(t37,f8.7)') rho              !Density
read (17,'(t37,f8.5)') T                !Temperature
read (17,'(t37,I8)') equil_iter         !Equilibrium sweeps
read (17,'(t37,I8)') prod_iter          !Production sweeps
read (17,'(t37,f8.7)') r_max            !Maximum MC displacement
read (17,'(t37,f8.5)') r_min            !Minimim interatomic distance for a random structure
read (17,'(t37,f8.5)') r_cut            !Potential cutoff distance
read (17,'(t37,f8.5)') max_radial       !Maximum value of distance for which g(r) is calculated
read (17,'(t37,I8)') r_points           !Number of points used for the calculations of g(r)
read (17,'(t37,A)') mode                !Potential type: Lennard Jones (LJ) or Stillinger (ST)
read (17,'(t37,L)') NL_mode             !Use of neighbour lists to optimize the calculation, can be True of False
read (17,'(t37,f8.5)') r_skin           !r_{skin} value, used for the calculation of neighbour lists
read (17,'(t37,L)') structure_given     !True if an initial structure is given, False if it is randomly generated
if ( structure_given ) then 
    read (17,'(t37,A)') structure_file  !Name of the initial structure file
    open ( 18, file = structure_file, action = 'read')
    read (18,*) N                       !Number of atoms from the initial structure file
    read (18,*)
endif
close (17)

L = ( real( N ) / rho ) ** ( 1. /3. )              !Length of the box used in the simulation

open(unit = 17, file = 'MC_simulation.log', action = 'write')
!A .log file that stores the parameters of the simulation and some results

write (17,*) "MC simulation parameters"
write (17,*)
write (17,*) "Number of atoms:                ", N
write (17,*) "Density:                        ", rho 
write (17,*) "Box length:                     ", L
write (17,*) "Temperature:                    ", T
write (17,*) "Equilibrium sweeps:             ", equil_iter
write (17,*) "Production sweeps:              ", prod_iter
write (17,*) "Maximum MC displacement:        ", r_max
write (17,*) "Minimum atom-atom distance:     ", r_min
write (17,*) "Potential cutoff distance:      ", r_cut
write (17,*) "Max g(r) distance:              ", max_radial
write (17,*) "Number of g(r) points:          ", r_points

if ( mode .eq. "LJ" ) then !Writes which potential type is being used to the log file
    write (17,*) "Lenard Jones Potential"
elseif ( mode .eq. "ST" ) then
    write (17,*) "Stillinger Potential"
endif

open( unit = 16, file = 'initial_structure.xyz', action = 'write')
!An .xyz file which stores the initial structure used in the simulation

write(16, *) N 
write(16, *)

!If the calculation is performed using a given structure the program reads the structure from "structure_file"
allocate( initial_coor ( 3, N ) )
if ( structure_given ) then
    write (17,*) "Using a given initial structure of ",N,"atoms"
    do i = 1, N
        read (18,*) trash, initial_coor (:, i)
        do j=1,3
            initial_coor ( j, i ) = initial_coor ( j, i ) / 2.1
        enddo
    enddo
    close(18)
else !If there is no a given initial structure a random initial structure is generated

    write (17,*) "Using a random initial structure of ",N,"atoms"
    do i = 1, N
        control = .True.
        do while (control)
            do j = 1, 3
                call random_number( r )
                initial_coor ( j, i ) = r * L
            enddo
            control = .False.
            do k = 1, i-1
                if ( control ) exit
                call distance_squared( initial_coor ( :, i ), initial_coor ( :, k ), r_dist, hL, L )
                if ( sqrt(r_dist) .lt. r_min ) control = .True.        
            enddo      
        enddo
    enddo
endif

!The initial structure is written in "initial_structure.xyz"
do i = 1, N
    write(16,*) initial_coor ( :, i )
enddo
close(16)

!________#_#_#_#_#_#________END Parameters setup and initial structure ________#_#_#_#_#_#________


!________#_#_#_#_#_#________ Montecarlo simulation ________#_#_#_#_#_#________

!Allocation of needed variables
allocate(trial_coor(3,N),dist_mat(N,N),N_hist(r_points),g_rad(2,r_points))
dist_mat = 0                    !Matrix of distances between atoms
N_accum = 0                         !Number of times the histogram of g(r) has been updated
success_r = 0                       !Success rate

equil_iter = equil_iter*N           !Number of equilibrium iterations
prod_iter = prod_iter * N           !Number of production iterations
hL = L/2.                           !Half length of the box used in the simulation
iterations = equil_iter + prod_iter !Total number of MC iterations 
delta_r = max_radial/r_points       !Distance step for the radial distribution function calculation 

do i=1,r_points
    g_rad(1,i) = delta_r * ( i - 1 )!g_rad stores r_points of distance and g(r) for each of them
enddo

if ( NL_mode ) then     !WITH NEIGHBOUR LISTS
    write (17,*) "Using neighbour lists."
    r_n = r_cut + r_skin                !Maximum distance of the atoms included in a neighbour list
    !n_max Number of elements that contains each neighbour list
    n_max = int( rho * 4 * 4.D0 * DATAN( 1.D0 ) * r_n ** ( 3. ) * 2. / 3. ) 
    allocate(positions_vec(3,N))    
    positions_vec = 0   !This vector stores the position of a given atom when its neighbour list was updated for the last time
    allocate(neighbour_mat(n_max,N))
    neighbour_mat = 0   !A matrix that contains the neighbour list of an atom on each column
    do it = 1, iterations
        if ( mod(it,N).eq.0) write(*,*) "Iteration", it / N, "of", iterations / N
        
        !MC_movement moves at random the position of a random atom
        call MC_movement(initial_coor, r_max, L, N, atom, trial_coor)

        if ( it .eq. 1 ) then
            !In the first iteration the whole distances matrix is calculated
            call distances_matrix(initial_coor,N,hL,L,dist_mat)
        else
            !In the following iterations only the row and column which correspond to the atom
            !moved at random by MC_movement is updated as the rest remains the same
            call update_distances_matrix(initial_coor, atom, N,hL,L,dist_mat)
        endif

        !Neighbour lists creation
        if (neighbour_mat(1,atom) .eq. 0) then
            !If the neighbour list is empty the position of the atom is stored in positions_vec
            !and the neighbour list is created
            positions_vec(:,atom) = initial_coor(:,atom)
            call neighbour_list(dist_mat,neighbour_mat,r_n, atom, N)
        else
            !If the neighbour list is not empty the distance between the atom position and the position of
            !the atom when the neighbour matrix was updated last time is calculated, if it is larger than
            !r_skin/2 the position and the neighbour list is updated
            call distance_squared(trial_coor(:,atom),positions_vec(:,atom),r_dist,hL,L)
            if ( sqrt(r_dist) .gt. r_skin / 2. ) then
                positions_vec(:,atom) = initial_coor(:,atom)
                call neighbour_list(dist_mat,neighbour_mat,r_n, atom, N)
            endif
        endif

        !Choose potential
        if (mode .eq. "LJ") then
            !If Lenard Jones potential is used the following subroutine calculated the potential
            !energy difference between the trial and the previous structure
            call LJ_potential_diff_NL(initial_coor, trial_coor, neighbour_mat, atom, n_max, delta_v, r_cut,hL,L)
        
        elseif (mode .eq. "ST") then
            !If Stillinger potential is used the neighbour list of the neighbours of the atom
            !that has been moved in the trial coordinates must be updated if they are outdated
            ! or empty
            do i = 1, n_max
                if (neighbour_mat(1,neighbour_mat(i,atom)) .eq. 0) then
                    positions_vec(:,neighbour_mat(i,atom)) = initial_coor(:,neighbour_mat(i,atom))
                    call neighbour_list(dist_mat,neighbour_mat,r_n, int(neighbour_mat(i,atom),8), N)
                else
                    call distance_squared(trial_coor(:,neighbour_mat(i,atom)),positions_vec(:,neighbour_mat(i,atom)),r_dist,hL,L)
                    if ( sqrt(r_dist) .gt. r_skin / 2. )then
                        positions_vec(:,neighbour_mat(i,atom)) = initial_coor(:,neighbour_mat(i,atom))
                        call neighbour_list(dist_mat,neighbour_mat,r_n, int(neighbour_mat(i,atom),8), N)
                    endif
                endif
            enddo
            !The following subroutine calculates the potential energy difference between
            !the trial and the previous structure
            call ST_potential_diff_NL(initial_coor, trial_coor, neighbour_mat, hL, L, r_cut, atom, n_max, delta_v)
        endif  

        !The following subroutine decides if the trial structure is accepted or not
        !depending on the value of the potential energy difference
        call decision(delta_v, T, atom, initial_coor, trial_coor, success_r)

        !The histogram of g(r) is updated every sweep (1st condition) only 
        !if the iterations are in the production iterations
        if ( (mod(it,N).eq.0 ) .and. (it .gt. equil_iter)) &
        call N_histogram(dist_mat, N,r_points,delta_r,N_accum,N_hist)
    enddo


else        !WITHOUT NEIGHBOUR LISTS
    write (17,*) "Not using neighbour lists."
    do it = 1, iterations
        if ( mod(it,N).eq.0 ) write(*,*)"Iteration", it / N, "of", iterations / N
        
        !MC_movement moves at random the position of a random atom
        call MC_movement(initial_coor,r_max,L,N,atom,trial_coor)

        if (it .eq. 1) then
            !In the first iteration the whole distances matrix is calculated
            call distances_matrix(initial_coor,N,hL,L,dist_mat)
        else
            !In the following iterations only the row and column which correspond to the atom
            !moved at random by MC_movement is updated as the rest remains the same
            call update_distances_matrix(initial_coor, atom, N,hL,L,dist_mat)
        endif
        
        !Choose potential
        if (mode .eq. "LJ") then
            !If Lenard Jones potential is used the following subroutine calculated the potential
            !energy difference between the trial and the previous structure
            call LJ_potential_diff(initial_coor, trial_coor, atom, delta_v, r_cut,N,hL,L)
        elseif (mode .eq. "ST") then
            !The following subroutine calculates the potential energy difference between
            !the trial and the previous structure
            call ST_potential_diff(initial_coor, trial_coor, hL, L, r_cut, N, atom, delta_v)
        endif  

        !The following subroutine decides if the trial structure is accepted or not
        !depending on the value of the potential energy difference
        call decision(delta_v, T, atom, initial_coor, trial_coor, success_r)

        !The histogram of g(r) is updated every sweep (1st condition) only 
        !if the iterations are in the production iterations
        if ( (mod(it,N).eq.0 ) .and. ( it .gt. equil_iter )) &
        call N_histogram(dist_mat, N,r_points,delta_r,N_accum,N_hist)
    enddo
endif

!________#_#_#_#_#_#________ END Montecarlo simulation ________#_#_#_#_#_#________

!________#_#_#_#_#_#________ Results ________#_#_#_#_#_#________

write(17,*)"__________________Results__________________"

!The success ratio is calculated and written to the log file
write (17,*) "Success ratio:                  ",success_r/iterations 

!The radial distribution function is calculated
call radial_distribution(N,r_points,g_rad,rho,N_accum,N_hist)

!The radial distribution function is written in "radial.dat"
open(unit = 19, file = 'radial.dat', action = 'write')
do i=1,r_points
    write(19,*)g_rad(:,i)
enddo
close(19)

!Finally a gnuplot file is created in order to plot the radial distribution function
call gnuplot_radial()

!The cpu time of the calculation is obtained and written to the log file
call cpu_time(time_2)
write (17,*) "Total CPU time                  ", time_2-time_1
close(17)

!________#_#_#_#_#_#________ END Results ________#_#_#_#_#_#________

end program main
