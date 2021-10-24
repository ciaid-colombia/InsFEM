module Mod_nsm_Output
   use MPI
   use Mod_iofile
   use Mod_int2str
   use Mod_PointerSetter
   use Mod_nsm_BaseBouope
   implicit none
   private
   public SetPointersOutput

   type, extends(PointerSetter) :: SPOutput
contains
      procedure :: SpecificSet => SpecificSetOutput
   end type
   type(SPOutput) :: SetPointersOutput

   real(rp), allocatable, dimension(:,:) :: force, momen
   integer(ip) :: ierr
 
contains

   subroutine SpecificSetOutput(d)
      implicit none
      class(SPOutput) :: d
      if (a%kfl_outfm==1) then
         if(a%kfl_openforce) call ConcatenateProcedures(ProcHook_Initializations,OpenFiles)
         call ConcatenateProcedures(ProcHook_Initializations,AllocForces)
         call ConcatenateProcedures(ProcHook_Initializations,ToZero)
         call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,ForceGauss)
         call ConcatenateProcedures(ProcHook_PostLoop,ComputeForce)
         call ConcatenateProcedures(ProcHook_PostLoop,OutputForce)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocForces)
      end if
   end subroutine

   subroutine AllocForces
      implicit none

      call a%Mesh%GetNbody(nbody)
      call a%Memor%alloc(ndime,nbody,force,'force','nsi_forces')
      call a%Memor%alloc(3,nbody,momen,'momen','nsi_forces')
   end subroutine

   subroutine DeallocForces
      implicit none

      call a%Mesh%GetNbody(nbody)
      call a%Memor%dealloc(ndime,nbody,force,'force','nsi_forces')
      call a%Memor%dealloc(3,nbody,momen,'momen','nsi_forces')   
   end subroutine

   subroutine ToZero
      implicit none

      a%force=0.0_rp
      a%momen=0.0_rp
      force=0.0_rp
      momen=0.0_rp
   end subroutine

   subroutine OpenFiles
      implicit none
      character(150) :: fil_force

      call a%Mesh%GetNbody(nbody)
      if (a%MPIrank == a%MPIroot) then 
         do ibody=1,nbody
            fil_force = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//adjustl(trim(int2str(ibody)))//'.frc'
            call iofile(zero,a%lun_force(ibody),fil_force,adjustl(trim(a%exmod))//'FORCES & MOMENTS')
            if (ndime==2) write (a%lun_force(ibody),12)
            if (ndime==3) write (a%lun_force(ibody),13)
         end do
         a%kfl_openforce = .false.
      endif      

     !Formats
     12   format( &
         '$        Forces and Moments over bodies          ',/,          &
         '$    Time',15x,'Fx',15x,'Fy',15x,'Mz')
     13   format( &
         '$          Forces and Moments over bodies            ',/,      &
         '$ Time',15x,'Fx',15x,'Fy',15x,'Fz',15x,'Mx',15x,'My'15x,'Mz')
   end subroutine

   subroutine ForceGauss
      implicit none

      if (ibody>0) then
         !Forces
         do idime=1,ndime
            force(idime,ibody)=force(idime,ibody) + tract(idime)*dsurf*weightfactor
         end do

         !Moments
         if (ndime==2) then
            momen(3,ibody) = momen(3,ibody) &
               + tract(2)*(gpcod(1)-a%origm(1))*dsurf*weightfactor &
               - tract(1)*(gpcod(2)-a%origm(2))*dsurf*weightfactor 
         elseif (ndime==3) then
            momen(1,ibody) = momen(1,ibody) &
               + tract(3)*(gpcod(2)-a%origm(2))*dsurf*weightfactor &
               - tract(2)*(gpcod(3)-a%origm(3))*dsurf*weightfactor          
            momen(2,ibody) = momen(1,ibody) &
               + tract(1)*(gpcod(3)-a%origm(3))*dsurf*weightfactor &
               - tract(3)*(gpcod(1)-a%origm(1))*dsurf*weightfactor  
            momen(3,ibody) = momen(3,ibody) &
               + tract(2)*(gpcod(1)-a%origm(1))*dsurf*weightfactor &
               - tract(1)*(gpcod(2)-a%origm(2))*dsurf*weightfactor
         end if
      end if
   end subroutine

   subroutine ComputeForce
      implicit none
      if (nbody>0) then
        !call MPI_ALLREDUCE(force,a%force,ndime*nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
        !call MPI_ALLREDUCE(momen,a%momen,3*nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
        call MPI_REDUCE(force,a%force,ndime*nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
        call MPI_REDUCE(momen,a%momen,3*nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
      end if
   end subroutine

   subroutine OutputForce
      implicit none

      call a%Mesh%GetNbody(nbody)
      if (a%MPIrank == a%MPIroot) then
         if (ndime==2) then
            do ibody = 1,nbody
               write(a%lun_force(ibody),20) a%ctime,&
                     -a%force(1,ibody)/a%adimf(1),  &
                     -a%force(2,ibody)/a%adimf(2),  &
                     -a%momen(3,ibody)/a%adimm(3)
            end do
         else if(ndime==3) then
            do ibody = 1,nbody
               write(a%lun_force(ibody),30) a%ctime,&
                     -a%force(1,ibody)/a%adimf(1),  &
                     -a%force(2,ibody)/a%adimf(2),  &
                     -a%force(3,ibody)/a%adimf(3),  &
                     -a%momen(1,ibody)/a%adimm(1),  &
                     -a%momen(2,ibody)/a%adimm(2),  &
                     -a%momen(3,ibody)/a%adimm(3)
            end do
         end if
      end if

      if (a%kfl_flush == 1) then
         do ibody = 1,nbody
            call flush(a%lun_force(ibody))
         end do
      endif
   

      !Formats
      20   format(1x,e12.6,2x,3(2x,e15.8))
      30   format(1x,e12.6,2x,6(2x,e15.8))
   end subroutine

end module
