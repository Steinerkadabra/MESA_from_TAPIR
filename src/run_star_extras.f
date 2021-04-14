      ! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      
      real(dp), save ::  Lacc, Ladd


      
      ! these routines are called by the standard run_star check_model
      contains


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_energy => energy_routine
         s% other_adjust_mdot => other_adjust_mdot

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         s% job% initial_h2 = s% x_ctrl(1)*1d-6
         s% job% initial_he3 = s% x_ctrl(2)*1d-6
         s% job% initial_he4 = 0.224 + 2*s% initial_z - s% x_ctrl(2)*1d-6
         s% initial_y = 0.224 + 2*s% initial_z 
         s% job% initial_h1 = 1 - s% job% initial_h2 - s% job% initial_he3 -s% job% initial_he4 - s % initial_z

         s% accretion_h2 = s% x_ctrl(1)*1d-6
         s% accretion_he3 = s% x_ctrl(2)*1d-6
         s% accretion_he4 = 0.224 + 2*s% initial_z - s% x_ctrl(2)*1d-6
         s% accretion_h1 = 1 - s% job% initial_h2 - s% job% initial_he3 -s% job% initial_he4 - s % initial_z


         s% xa_central_lower_limit(1) = s% job% initial_h1 - 0.001

      end subroutine extras_controls
      
      
       subroutine other_adjust_mdot(id, ierr)
            use star_def
            implicit none
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            integer :: step, i
            type(star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr) ! retrieve star id


            if (.not. s% doing_relax) then
            s% mstar_dot = (s% x_ctrl(3)* (s% x_ctrl(7) - s% star_age)+ &
                            s% x_ctrl(4)* (s% star_age - s% x_ctrl(6) ))/&
                            (s% x_ctrl(7) - s% x_ctrl(6)) *Msun/secyer
            end if



        end subroutine other_adjust_mdot

subroutine energy_routine(id, ierr)
         use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         real(dp) :: mass, mdot, mdot_msun, alpha, radius, mcumul
         real(dp) :: xposk, dif, xpos, thingie, thingie2, Mstar, mass_dif
         integer :: k, xposk_orig

            real(dp), parameter :: mdot_break = 6.2e-6
            real(dp), parameter :: alpha_high = 0.2
            real(dp), parameter :: alpha_low = 0.005
            real(dp), parameter :: delta_break = 5.95e-6
            real(dp) :: numerator, denominator


         ierr = 0
         call star_ptr(id, s, ierr)

            mass = s% mstar
            mdot = s% mstar_dot
            mdot_msun = s% mstar_dot / Msun* secyer
            alpha = 0.1



            mass = s% mstar
            mdot = s% mstar_dot

            if (mdot == 0.) return

            radius = s% r(1)

            mcumul = 0.
            k = 0
            Lacc = 0.5 * (standard_cgrav * (s% mstar_dot) * (mass)) / radius

            !Ladd = min(1.0, float(s% model_number)/1000.)*alpha/2*standard_cgrav * (s% star_mass *Msun) * (s% mstar_dot) / (s% r(1)) ![erg/s]
            Ladd = alpha/2*standard_cgrav * (s% star_mass *Msun) * (s% mstar_dot) / (s% r(1)) ![erg/s]
            Lacc =  safe_log10((1-alpha)/2*standard_cgrav * (s% star_mass *Msun) *&
     (s% mstar_dot) / (s% r(1))  /Lsun) ![erg/s]
      Mstar=(s% star_mass *Msun) ![g]

       !do while ((Ladd / max(mcumul,1.) > 1d3) .and. mcumul < 0.8 * mass)
         !   k = k + 1
         !   mcumul = mcumul + s% dm(k)
         !enddo
   
  xpos = s% x_ctrl(5)
  do k = 1, s% nz
            if (s% star_mdot > 0.0d0 .and. s% m(k)>=(1.0d0-xpos)*(s% star_mass)*msol ) then
                  s% extra_heat(k) = 2.0d0*Ladd/(Mstar*xpos**2)*(s% m(k)/Mstar-(1.0d0-xpos))

            end if
         end do

  
       end subroutine energy_routine



      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         integer :: nlines, io, i
         type (star_info), pointer :: s
         real(kind=8) :: offset
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

 
      end subroutine extras_startup




      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
   real(dp) :: mass_dif 
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0


         s% hydro_save_photo = .false.
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
   character(len = 13) :: name_base
   character(len = 5) :: mass_name
   character(len = 18) :: final_name
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if






         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: max_temp
         integer :: max_pos, k
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         names(1) = 'log_opacity'
         names(2) = 'grad_temperature'
         names(3) = 'grad_L'
         do k = 1, s% nz
            vals(k,1) = safe_log10(s% opacity(k))
            vals(k,2) = safe_log10(s% grad_temperature(k))
            vals(k,3) = safe_log10(s% gradL(k))
         end do
 
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         open (unit = 14, file = "MESA_for_TAPIR.txt")
         write(unit=14, fmt="(F28.15)") s% Teff, s% star_age, s% r(1), s% photosphere_L, Lacc
         close (unit = 14)



         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      
      end module run_star_extras
      
