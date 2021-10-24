subroutine lev_zigzag(a,e,mnode,StartElementLayer,auxPointList,PointDone,PointListaux,PointList,auxlevel)
!DESCRIPTION
!   This get the type of element: cut=1 and non cut=0
!-----------------------------------------------------------------------
   use typre
   use Mod_Mesh
   use Mod_Memor   
   use Mod_LevelSet
   use Mod_Element
   use Mod_CutMesh
   use def_parame
   implicit none
   class(LevelSetProblem) :: a   
   class(FiniteElement) :: e
   
   integer(ip), intent(in)    :: mnode  
   integer(ip)                :: nelem,ndime,npoin
   real(rp), intent(inout)    :: auxlevel(*)
   logical, intent(inout)     :: PointDone(*),PointListaux(*)
   integer(ip), intent(inout) :: PointList(*),auxPointList
   integer(ip)                :: StartElementLayer
   
   !auxiliar variables
   integer(ip) :: elemj,jelem,jpoin,kpoin,knode,nodej,nodek,jnode,idime
   integer(ip) :: ToDoList(mnode),DoneList(mnode)   
   integer(ip) :: nNodeDone,nNodeToDo,poinStatus
   real(rp)    :: distance,mindist
   integer(ip) :: auxsign
   
   call a%Mesh%GetNpoin(npoin)   
   call a%Mesh%GetNelem(nelem)  
   call a%Mesh%GetNdime(ndime)    
   

   
   if (StartElementLayer > -1) then
      elemj=1
      do elemj=1+StartElementLayer,nelem   
         jelem=a%ElementListByLayers(elemj)      
         call a%Mesh%ElementLoad(jelem,e)             
               
         nNodeDone=0
         nNodeToDo=0
                  
         do jnode=1,e%pnode
            jpoin=e%lnods(jnode)
            if(PointDone(jpoin) .eqv. .false.)then
               nNodeToDo=nNodeToDo+1
               ToDoList(nNodeToDo)=jnode
            else
               nNodeDone=nNodeDone+1
               DoneList(nNodeDone)=jnode
            end if
         end do   
            
         do nodej=1,nNodeToDo
            mindist=1e12  
            jnode=ToDoList(nodej)
            jpoin=e%lnods(jnode)
                  
            do nodek=1,nNodeDone
               knode = DoneList(nodek)
               kpoin = e%lnods(knode)
                     
               !mimimum distance
               distance = 0.0_rp
               do idime=1,ndime
                  distance = distance + (e%elcod(idime,jnode) - e%elcod(idime,knode))**2.0_rp
               enddo  
               !sign definition
               call a%CutMesh%GetPointType(jpoin,poinStatus) 
               if(poinStatus<0)then
                  auxsign=-1_ip
               elseif(poinStatus>0)then
                  auxsign=1_ip
               end if
               distance = sqrt(distance)*auxsign
               
               distance = distance + auxlevel(kpoin)
               if(abs(distance) < abs(mindist))then          
                  mindist =distance
               end if
            end do         

            auxlevel(jpoin) = mindist  
            PointDone(jpoin)= .true.
   
            if(PointListaux(jpoin) .eqv. .false.)then
               PointListaux(jpoin) = .true.
               !Auxiliar List
               auxPointList = auxPointList + 1
               PointList(auxPointList) = jpoin            
            end if

         end do             
      end do   
   endif
   
end subroutine