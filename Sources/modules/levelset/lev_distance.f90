subroutine lev_distance(a,e,mnode,StartElementLayer,auxPointList,PointDone,PointListaux,PointList,auxlevel)
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
   integer(ip) :: startElementLayer
   
   !auxiliar variables
   integer(ip) :: elemi,elemk,elemj
   integer(ip) :: ToDoList(mnode),DoneList(mnode)   
   integer(ip) :: nNodeDone,nNodeToDo,nelemList,nPointDone
   real(rp)    :: distance,hipo,hipo2,auxhipo
   integer(ip) :: jpoin,jnode,kpoin,knode,nodej,poini,nodek,jelem,nipoin,poinj
   integer(ip) :: hpoin,hnode,idime
   
   real(rp)    :: grad(3),diference(3),diference2(3),auxdiference(3)
   real(rp)    :: localx1(3),localx2(3),vecprojx1(3),vecprojx2(3), &
                  localx3(3),vecprojx3(3),vecprojx21(3)                               
   integer(ip) :: gnode,nodeg,fnode,nodef,jdime,idir,fpoin,gpoin,poinStatus   
   real(rp)    :: projx1,projx2,projx3,normlocalx2,normlocalx3,normvecprojx2,projx21,&
                  normvecprojx1,grad2projx2,perpendicularcheck,ponderation 
   integer(ip) :: auxsign
   integer(ip) :: kfl_Hanging
   logical :: ParentsAreDone
   
   procedure(), pointer :: CheckHanging => NULL() 
   
   !External Procedures
   procedure() :: NULLSUB  

   call a%Mesh%GetNpoin(npoin)   
   call a%Mesh%GetNelem(nelem)  
   call a%Mesh%GetNdime(ndime)    
   
   !Set pointer for hangingnodes
   CheckHanging => NULLSUB
   call a%Mesh%GetHanging(kfl_Hanging)
   if (kfl_Hanging == 1) then
      CheckHanging => CheckAndUpdateHangingNodes
   endif

   if (StartElementLayer > -1) then
      elemj=1
      if(ndime==2)then
         do elemj=1+StartElementLayer,nelem   
            jelem=a%ElementListByLayers(elemj)      
            call a%Mesh%ElementLoad(jelem,e)             
                  
            nNodeDone=0
            nNodeToDo=0
                     
            do jnode=1,e%pnode
               jpoin=e%lnods(jnode)
               
               !HangingNodes
               call CheckHanging
               
               if(PointDone(jpoin) .eqv. .false.)then
                  nNodeToDo=nNodeToDo+1
                  ToDoList(nNodeToDo)=jnode
               else
                  nNodeDone=nNodeDone+1
                  DoneList(nNodeDone)=jnode
               end if
            end do      
            
            if (nnodetodo > 0) then
            !Initialization
            diference = 0.0_rp
            grad = 0.0_rp   
            hipo = 0.0_rp
               
            !all between the done nodes
            fnode=DoneList(1)
            fpoin=e%lnods(fnode)
            gnode=DoneList(2)
            gpoin=e%lnods(gnode)
            do jdime=1,ndime                    
               diference(jdime) = diference(jdime) + (e%elcod(jdime,gnode) - e%elcod(jdime,fnode))                    
               hipo = hipo + diference(jdime)**2.0_rp                  
            end do              
            hipo=sqrt(hipo)      
               
            !diference between the reference done-node and the todo-node
            auxdiference=0
            hipo2=0
            !diference between the reference done-node and the todo-node
            do nodek =1,nNodeToDo
               knode=ToDoList(nodek)
               kpoin=e%lnods(knode) 
               do jdime=1,ndime              
                     
                  auxdiference(jdime) = auxdiference(jdime) + (e%elcod(jdime,knode) - e%elcod(jdime,fnode))  
                  hipo2 = hipo2 + auxdiference(jdime)**2.0_rp
                  
               end do         
            end do 
            hipo2=sqrt(hipo2)

            !known gradient direction
            localx1   = diference/hipo
            projx1    = dot_product(auxdiference,localx1)
            vecprojx1 = projx1*localx1
                  
            vecprojx2 = auxdiference - vecprojx1 
            normlocalx2 = 0.0_rp
            do idime=1,ndime
               normlocalx2 = normlocalx2 + vecprojx2(idime)*vecprojx2(idime)
            end do
            normlocalx2 = sqrt(normlocalx2)          
            
            localx2 = vecprojx2/normlocalx2
            projx2    = dot_product(auxdiference,localx2)            

            grad(1) = (auxlevel(gpoin)-auxlevel(fpoin))/hipo         
            if(abs(grad(1))>=1.0_rp)then
               grad(2) = 0.0_rp
            else
               grad(2) = sqrt(1.0_rp-grad(1)**2.0_rp)          
            endif
               
            do nodej=1,nNodeToDo 
               jnode=ToDoList(nodej)
               jpoin=e%lnods(jnode)
                        
               !distance
               distance = 0.0_rp 
               call a%CutMesh%GetPointType(jpoin,poinStatus) 
               if(poinStatus<0)then
                  auxsign=-1_ip
               elseif(poinStatus>0)then
                  auxsign=1_ip
               end if
               
               distance = auxlevel(fpoin) + grad(1)*projx1 + grad(2)*auxsign*projx2            
               
               auxlevel(jpoin) = distance  
               PointDone(jpoin)= .true.
         
               if(PointListaux(jpoin) .eqv. .false.)then
                  PointListaux(jpoin) = .true.
                  !Auxiliar List
                  auxPointList = auxPointList + 1
                  PointList(auxPointList) = jpoin            
               end if

            end do   
            endif
         end do

      elseif(ndime==3)then
         
         do elemj=1+StartElementLayer,nelem   
            jelem=a%ElementListByLayers(elemj)      
            call a%Mesh%ElementLoad(jelem,e)             
                  
            nNodeDone=0
            nNodeToDo=0
                     
            do jnode=1,e%pnode
               jpoin=e%lnods(jnode)
               
               !HangingNodes
               call CheckHanging
               
               if(PointDone(jpoin) .eqv. .false.)then
                  nNodeToDo=nNodeToDo+1
                  ToDoList(nNodeToDo)=jnode
               else
                  nNodeDone=nNodeDone+1
                  DoneList(nNodeDone)=jnode
               end if
            end do 
            
            
            !Initialization
            diference    = 0.0_rp
            diference2   = 0.0_rp
            auxdiference = 0.0_rp   
            hipo    = 0.0_rp
            hipo2   = 0.0_rp
            auxhipo = 0.0_rp
            grad = 0.0_rp 
               
            !done nodes
            fnode=DoneList(1)
            fpoin=e%lnods(fnode)
            gnode=DoneList(2)
            gpoin=e%lnods(gnode)
            hnode=DoneList(3)
            hpoin=e%lnods(hnode)         
               
            do jdime=1,ndime               
                        
               diference(jdime) = diference(jdime) + (e%elcod(jdime,gnode) - e%elcod(jdime,fnode)) 
               diference2(jdime) = diference2(jdime) + (e%elcod(jdime,hnode) - e%elcod(jdime,fnode))
               hipo = hipo + diference(jdime)**2.0_rp   
               hipo2 = hipo2 + diference2(jdime)**2.0_rp 
                  
            end do         
            hipo=sqrt(hipo) 
            hipo2=sqrt(hipo2)
               
            !diference between the reference done-node and the todo-node
            do nodek =1,nNodeToDo
               knode=ToDoList(nodek)
               kpoin=e%lnods(knode) 
               do jdime=1,ndime               
                        
                  auxdiference(jdime) = auxdiference(jdime) + (e%elcod(jdime,knode) - e%elcod(jdime,fnode))  
                  auxhipo = auxhipo + auxdiference(jdime)**2.0_rp
                     
               end do         
            end do 
            auxhipo=sqrt(auxhipo)

            !------------------------------------------------------------------------------------------------
            !projection and definition of the distance function            
            localx1   = diference/hipo
            projx21   = dot_product(diference2,localx1)
            vecprojx21 = projx21*localx1
            !diference2 is not necessary perpendicular to diference
            vecprojx2  = diference2 - vecprojx21
               
            normvecprojx2=0.0_rp
            do idime=1,ndime         
               normvecprojx2 = normvecprojx2 + vecprojx2(idime)*vecprojx2(idime)            
            end do
            normvecprojx2 = sqrt(normvecprojx2)        
            localx2 = vecprojx2/normvecprojx2
               
      !          perpendicularcheck=dot_product(localx1,localx2)
      !          if(perpendicularcheck/=0.0_rp) write(*,*) 1,jelem,perpendicularcheck         
               
            !now we can project the unknown value in x1 and x2
            projx1 = dot_product(auxdiference,localx1)
            projx2 = dot_product(auxdiference,localx2)
               
            vecprojx3 = auxdiference - projx1*localx1 - projx2*localx2         
               
            normlocalx3 = 0.0_rp
            do idime=1,ndime
               normlocalx3 = normlocalx3 + vecprojx3(idime)*vecprojx3(idime)
            end do
            normlocalx3 = sqrt(normlocalx3) 
               
            !normlocalx3 can be zero when nNodeToDo=0  
            if(normlocalx3==0.0_rp)then
               localx3 = 0.0_rp
            else         
               localx3 = vecprojx3/normlocalx3
            end if
            projx3= dot_product(auxdiference,localx3)

            !Grad(1) is always in the x1 direction by definition
            grad(1) = (auxlevel(gpoin)-auxlevel(fpoin))/hipo 
            !Grad(2) is not necessary in the x2 direction
            grad(2) = (auxlevel(hpoin)-(auxlevel(fpoin)+grad(1)*projx21))/normvecprojx2
               
            
      !          if(abs(grad(2)*projx2)<=1e-10)then
      !             grad(2)=0.0_rp
      !          end if
      !          
      !          if(abs((grad(1)*projx1))<=1e-10)then
      !             grad(1)=0.0_rp
      !          end if
               
   !          if((abs(grad(1))>=1.0_rp) .or. (abs(grad(2))>=1.0_rp))then
   !             if(abs(grad(1))>abs(grad(2)))then
   !                grad(2)=0.0_rp
   !                grad(1)=1.0_rp*(grad(1)/abs(grad(1)))
   !             else
   !                grad(1)=0.0_rp
   !                grad(2)=1.0_rp*(grad(2)/abs(grad(2)))
   !             end if
   !          elseif((abs(grad(1))>=1.0_rp) .and. (abs(grad(2))<1.0_rp))then
   !             grad(2)=0.0_rp
   !             grad(1)=1.0_rp*(grad(1)/abs(grad(1)))
   !          elseif((abs(grad(2))>=1.0_rp) .and. (abs(grad(1))<1.0_rp))then
   !             grad(1)=0.0_rp
   !             grad(2)=1.0_rp*(grad(2)/abs(grad(2)))
   !          end if         
               

            if((grad(1)**2.0_rp+grad(2)**2.0_rp)>1.0_rp)then
               grad(3) = 0.0_rp
   !             ponderation= sqrt(grad(1)**2.0_rp + grad(2)**2.0_rp)
   !             grad(1) = grad(1)/ponderation
   !             grad(2) = grad(2)/ponderation
               
            else
               grad(3) = sqrt(1.0_rp-(grad(1)**2.0_rp+grad(2)**2.0_rp))          
            endif
                  
               
            do nodej=1,nNodeToDo 
               jnode=ToDoList(nodej)
               jpoin=e%lnods(jnode)
                        
               !distance
               distance = 0.0_rp
               call a%CutMesh%GetPointType(jpoin,poinStatus) 
               if(poinStatus<0)then
                  auxsign=-1_ip
               elseif(poinStatus>0)then
                  auxsign=1_ip
               end if
               distance = auxlevel(fpoin) + grad(1)*projx1 + grad(2)*projx2 + grad(3)*projx3*auxsign
               auxlevel(jpoin) = distance  
               PointDone(jpoin)= .true.
         
               if(PointListaux(jpoin) .eqv. .false.)then
                  PointListaux(jpoin) = .true.
                  !Auxiliar List
                  auxPointList = auxPointList + 1
                  PointList(auxPointList) = jpoin            
               end if

            end do             
         end do
      end if
   endif
   
contains
   
   subroutine CheckAndUpdateHangingNodes
      implicit none
      !Check for Hanging Nodes
         if (PointDone(jpoin) .eqv. .false.) then
            call CheckLogicalParents(a%Mesh,jpoin,PointDone,ParentsAreDone)
            if (ParentsAreDone) then
               call interpolateHangingValueSingle(a%Mesh,jpoin,1_ip,auxlevel)
               PointDone(jpoin) = .true.
            endif
         endif
   end subroutine

   

end subroutine
