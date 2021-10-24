   subroutine nsc_GhostCommunicate(a)
      use typre
      use Mod_NSCompressible
      implicit none
      class(NSCompressibleProblem)  :: a

      integer(ip):: ndime
   
      call a%Mesh%GetNdime(ndime)

      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%densf(:,1))
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%momen(:,:,1))
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%energ(:,1))

   end subroutine
