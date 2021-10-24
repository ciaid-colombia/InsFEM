subroutine ManufacturedMesh(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip) :: ManufacturedNprocX, ManufacturedProcPosition(3)=0
   
   logical :: kfl_isPowerOf,kfl_InternalRadiusAtZero
   integer(ip) :: UpGhost0,DownGhost0,LeftGhost0,RightGhost0
   integer(ip) :: UpElem0, DownElem0,LeftElem0,RightElem0
   integer(ip) :: UpLeftGhost,UpRightGhost,DownLeftGhost,DownRightGhost
   integer(ip) :: UpLeftElem0,UpRightElem0,DownLeftElem0,DownRightElem0
   real(rp) :: ProcCoordSpanX,ElemCoordSpanX
   real(rp) :: coordx,basecoordx,basecoordy,coordy
   integer(ip) :: baseGhost0,TotalNpoinX
   
   integer(ip) :: basepointcol, jcol,irow,ipoin,groupi
   integer(ip) :: BlockDown,BlockUp,BlockRight,BlockLeft,rowdown,rowup,colright,colleft
   integer(ip) :: ielem,ispos
   integer(ip), allocatable :: LocalToGlobal(:), ProcessorList(:)

   integer(ip) :: TopGhost0,BottomGhost0
   integer(ip) :: TopElem0,BottomElem0
   integer(ip) :: TopUpGhost,TopRightGhost,TopDownGhost,TopLeftGhost
   integer(ip) :: BottomUpGhost,BottomRightGhost,BottomDownGhost,BottomLeftGhost
   integer(ip) :: TopUpElem0,TopRightElem0,TopDownElem0,TopLeftElem0
   integer(ip) :: BottomUpElem0,BottomRightElem0,BottomDownElem0,BottomLeftElem0

   integer(ip) :: TopUpLeftGhost,TopUpRightGhost,TopDownLeftGhost,TopDownRightGhost
   integer(ip) :: BottomUpLeftGhost,BottomUpRightGhost,BottomDownLeftGhost,BottomDownRightGhost
   integer(ip) :: TopUpLeftElem0,TopUpRightElem0,TopDownLeftElem0,TopDownRightElem0
   integer(ip) :: BottomUpLeftElem0,BottomUpRightElem0,BottomDownLeftElem0,BottomDownRightElem0

   real(rp) :: basecoordz,coordz
   integer(ip) :: basepointflr,kflr,jpoin
   integer(ip) :: BlockTop,BlockBottom,floortop,floorbottom
   
   real(rp) :: AxialLength,ExternalRadius,InternalRadius
   integer(ip) :: PointsInLength,PointsInRadius,DivisionsInAngle
   real(rp) :: DeltaAngleProc,DeltaAngle
   real(rp) :: SpanInRadialCoord,SpanInLengthCoord
   integer(ip) :: LocalAngleDivisionNpoin,LocalAngleNelem,TimesAngleIsDivided
   integer(ip) :: InnerGhost0,OtherPcsGhost0
   integer(ip) :: InnerElem0,UpInnerElem0,DownInnerElem0,OtherInnerElem0

   real(rp) :: BaseProcessorAngle,BaseAngle,BaseGroupAngle,BaseDivisionRadius
   integer(ip) :: ZeroPointAngle,FirstInnerPoint
   integer(ip) :: bottomleftpoint,bottomrightpoint,topleftpoint,toprightpoint,rightpoint
   integer(ip) :: floorbottomleftpoint,floorrightpoint,floorbottomrightpoint,floortopleftpoint,floortoprightpoint
   integer(ip) :: ceilingbottomleftpoint,ceilingrightpoint,ceilingbottomrightpoint,ceilingtopleftpoint,ceilingtoprightpoint

   if (a%ndime == 2 .and. a%ManufacturedType == 'LINQU') then
      !call isPowerOf(a%MPIsize,4_ip,kfl_isPowerOf)
      if ( nint(sqrt(real(a%MPIsize)))**2 /= a%MPIsize) then
         kfl_isPowerOf = .false.
      else
         kfl_isPowerOf = .true.
      endif
      
      if (kfl_isPowerOf .eqv. .false.) then
         call runend('2d ManufacturedMesh, MPIsize the square of a natural number')
      endif
      
      ManufacturedNprocX = sqrt(real(a%MPIsize,rp))
      
      ManufacturedProcPosition(1) = a%MPIrank/ManufacturedNprocX
      ManufacturedProcPosition(2) = mod(a%MPIrank,ManufacturedNprocX)
      
      a%gnpoin = a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%MPIsize
      a%npoinLocal = a%ManufacturedNpoinX*a%ManufacturedNpoinX
      
      a%nelem = (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      a%npoinGhost = 0
      UpGhost0 = -1
      DownGhost0 = -1
      LeftGhost0 = -1
      RightGhost0 = -1

      if (ManufacturedProcPosition(2) /= 0) then
         UpGhost0 = a%npoinLocal + a%npoinGhost
         UpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + a%ManufacturedNpoinX-1
      endif
      if (ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownGhost0 = a%npoinLocal + a%npoinGhost
         DownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + a%ManufacturedNpoinX-1
      endif
      if (ManufacturedProcPosition(1) /= 0) then
         LeftGhost0 = a%npoinLocal + a%npoinGhost
         LeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + a%ManufacturedNpoinX-1
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         RightGhost0 = a%npoinLocal + a%npoinGhost
         RightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + a%ManufacturedNpoinX-1
      endif

      UpRightGhost = -1
      DownRightGhost = -1
      UpLeftGhost = -1
      DownLeftGhost = -1
      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
         a%npoinGhost = a%npoinGhost + 1_ip
         UpLeftGhost = a%npoinLocal + a%npoinGhost
         UpLeftElem0 = a%nelem
         a%nelem = a%nelem + 1
      endif
      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         a%npoinGhost = a%npoinGhost + 1_ip
         DownLeftGhost = a%npoinLocal + a%npoinGhost
         DownLeftElem0 = a%nelem
         a%nelem = a%nelem + 1
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
         a%npoinGhost = a%npoinGhost + 1_ip
         UpRightGhost = a%npoinLocal + a%npoinGhost
         UpRightElem0 = a%nelem
         a%nelem = a%nelem + 1
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         a%npoinGhost = a%npoinGhost + 1_ip
         DownRightGhost = a%npoinLocal + a%npoinGhost
         DownRightElem0 = a%nelem
         a%nelem = a%nelem + 1
      endif

      a%npoin = a%npoinLocal+a%npoinGhost
      a%gnelem = (a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*a%MPIsize
      
      a%poin0 = a%MPIrank*a%npoinLocal
      a%elem0 = a%MPIrank*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)
      
      !Fill local nodes and elements
      call a%Memor%alloc(1_ip,a%pnods,'pnods','ManufacturedMesh')
      call a%Memor%alloc(a%nelem*4_ip,a%lnods,'lnods','ManufacturedMesh')
      call a%Memor%alloc(2_ip,a%npoin,a%coord,'coord','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','ManufacturedMesh')
      a%geoblock=0
      
      call a%Memor%alloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !Build local elements, nodes and coord
      TotalNpoinX = ManufacturedNprocX*a%ManufacturedNpoinX
      ElemCoordSpanX = 1_rp/real(TotalNpoinX-1,rp)
      
      basecoordx = ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(1))
      basecoordy = 1- ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(2))
      
      !Local Coordinates
      do jcol = 1,a%ManufacturedNpoinX
         coordx = basecoordx + (jcol-1)*ElemCoordSpanX
         basepointcol = (jcol-1)*a%ManufacturedNpoinX
         do irow = 1,a%ManufacturedNpoinX
            ipoin = basepointcol+irow   
            coordy = basecoordy - (irow-1)*ElemCoordSpanX
            a%coord(1,ipoin) = coordx
            a%coord(2,ipoin) = coordy
         enddo
      enddo

      BlockDown = a%npoinLocal
      BlockUp   = -BlockDown
      BlockRight = a%npoinLocal*ManufacturedNprocX
      BlockLeft = -BlockRight
      rowdown = 1
      rowup = -rowdown
      colright = a%ManufacturedNpoinX
      colleft = -colright

      !Local Elements
      a%pnods(1) = 4
      ielem = 0
      do jcol = 1,a%ManufacturedNpoinX-1
         basepointcol = (jcol-1)*a%ManufacturedNpoinX
         do irow = 1,a%ManufacturedNpoinX-1
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + irow + colright + rowdown
            a%lnods(ispos+2) = basepointcol + irow + colright
            a%lnods(ispos+3) = basepointcol + irow
            a%lnods(ispos+4) = basepointcol + irow + rowdown
         enddo
      enddo

      !LocalToGlobal for local points
      do ipoin = 1,a%npoinLocal
         LocalToGlobal(ipoin) = a%poin0+ipoin
         ProcessorList(ipoin) = a%MPIrank
      enddo

      !Ghost points and elements
      if (UpGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockUp + (colright-1)
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(UpGhost0 + ipoin) = baseGhost0 + 1 + colright*(ipoin-1)
            ProcessorList(UpGhost0 + ipoin) = a%MPIrank-1
         enddo
         do jcol = 1,a%ManufacturedNpoinX-1
            basepointcol = (jcol-1)*a%ManufacturedNpoinX
            
            ielem = UpElem0 + jcol
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + 1 + colright
            a%lnods(ispos+2) = upGhost0 + jcol + 1
            a%lnods(ispos+3) = upGhost0 + jcol 
            a%lnods(ispos+4) = basepointcol + 1
            
         enddo
      endif
      if (DownGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockDown
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(DownGhost0 + ipoin) = baseGhost0 + 1 + colright*(ipoin-1)
            ProcessorList(DownGhost0 + ipoin) = a%MPIrank+1
         enddo
         
         do jcol = 1,a%ManufacturedNpoinX-1
            basepointcol = (jcol-1)*a%ManufacturedNpoinX + a%ManufacturedNpoinX - 1
            
            ielem = DownElem0 + jcol
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = DownGhost0 + jcol + 1
            a%lnods(ispos+2) = basepointcol + 1 + colright
            a%lnods(ispos+3) = basepointcol + 1
            a%lnods(ispos+4) = DownGhost0 + jcol
            
         enddo
      endif
      if (LeftGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockLeft + colright*(a%ManufacturedNpoinX-1)
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(LeftGhost0 + ipoin) = baseGhost0 + ipoin
            ProcessorList(LeftGhost0 + ipoin) = a%MPIrank - ManufacturedNprocX
         enddo
         basepointcol = 0
         do irow = 1,a%ManufacturedNpoinX-1
            ielem = LeftElem0 + irow 
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + irow + rowdown 
            a%lnods(ispos+2) = basepointcol + irow
            a%lnods(ispos+3) = LeftGhost0 + irow
            a%lnods(ispos+4) = LeftGhost0 + irow + 1
                        
         enddo
      endif
      if (RightGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockRight
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(RightGhost0 + ipoin) = baseGhost0 + ipoin
            ProcessorList(RightGhost0 + ipoin) = a%MPIrank + ManufacturedNprocX
         enddo
         basepointcol = colright*(a%ManufacturedNpoinX-1)
         do irow = 1,a%ManufacturedNpoinX-1
            ielem = RightElem0 + irow
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = RightGhost0 + irow + 1
            a%lnods(ispos+2) = RightGhost0 + irow
            a%lnods(ispos+3) = basepointcol + irow
            a%lnods(ispos+4) = basepointcol + irow + rowdown
         enddo
      endif

      if (UpLeftGhost /= -1) then
         LocalToGlobal(UpLeftGhost) = a%poin0 + 1 + BlockUp + BlockLeft + colright*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) 
         ProcessorList(UpLeftGhost) = a%MPIrank - 1 - ManufacturedNprocX
         
         ielem = UpLeftElem0 + 1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = 1
         a%lnods(ispos+2) = UpGhost0 + 1
         a%lnods(ispos+3) = UpLeftGhost
         a%lnods(ispos+4) = LeftGhost0 + 1
      endif   
      
      if (UpRightGhost /= -1) then
         LocalToGlobal(UpRightGhost) = a%poin0 + 1 + BlockUp + BlockRight + rowdown*(a%ManufacturedNpoinX-1)
         ProcessorList(UpRightGhost) = a%MPIrank - 1 + ManufacturedNprocX
         
         ielem = UpRightElem0 + 1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = RightGhost0 + 1
         a%lnods(ispos+2) = UpRightGhost
         a%lnods(ispos+3) = UpGhost0 + a%ManufacturedNpoinX
         a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 
      endif   
      
      if (DownLeftGhost /= -1) then
         LocalToGlobal(DownLeftGhost) = a%poin0 + 1 + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1)
         ProcessorList(DownLeftGhost) = a%MPIrank + 1 - ManufacturedNprocX
         
         ielem = DownLeftElem0 + 1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = DownGhost0 + 1
         a%lnods(ispos+2) = a%ManufacturedNpoinX
         a%lnods(ispos+3) = LeftGhost0 + a%ManufacturedNpoinX
         a%lnods(ispos+4) = DownLeftGhost
         
      endif
         
      if (DownRightGhost /= -1) then
         LocalToGlobal(DownRightGhost) = a%poin0 + 1 + BlockDown + BlockRight 
         ProcessorList(DownRightGhost) = a%MPIrank + 1 + ManufacturedNprocX
         
         ielem = DownRightElem0 + 1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = DownRightGhost
         a%lnods(ispos+2) = RightGhost0 + a%ManufacturedNpoinX
         a%lnods(ispos+3) = BlockDown
         a%lnods(ispos+4) = DownGhost0 + a%ManufacturedNpoinX
         
      endif   
         
         
      
      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      call a%Memor%dealloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%dealloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !External normal, we simply allocate it
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ManufacturedMesh')
      
   endif

   if (a%ndime == 2 .and. a%ManufacturedType == 'LINTR') then
      
      !call isPowerOf(a%MPIsize,4_ip,kfl_isPowerOf)
      if ( nint(sqrt(real(a%MPIsize)))**2 /= a%MPIsize) then
         kfl_isPowerOf = .false.
      else
         kfl_isPowerOf = .true.
      endif
      
      if (kfl_isPowerOf .eqv. .false.) then
         call runend('2d ManufacturedMesh, MPIsize the square of a natural number')
      endif
      
      ManufacturedNprocX = sqrt(real(a%MPIsize,rp))
      
      ManufacturedProcPosition(1) = a%MPIrank/ManufacturedNprocX
      ManufacturedProcPosition(2) = mod(a%MPIrank,ManufacturedNprocX)
      
      a%gnpoin = a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%MPIsize
      a%npoinLocal = a%ManufacturedNpoinX*a%ManufacturedNpoinX
      
      a%nelem = 2*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      a%npoinGhost = 0
      UpGhost0 = -1
      DownGhost0 = -1
      LeftGhost0 = -1
      RightGhost0 = -1
      if (ManufacturedProcPosition(2) /= 0) then
         UpGhost0 = a%npoinLocal + a%npoinGhost
         UpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 2*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownGhost0 = a%npoinLocal + a%npoinGhost
         DownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 2*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= 0) then
         LeftGhost0 = a%npoinLocal + a%npoinGhost
         LeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 2*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         RightGhost0 = a%npoinLocal + a%npoinGhost
         RightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 2*(a%ManufacturedNpoinX-1)
      endif
      
      UpRightGhost = -1
      DownRightGhost = -1
      UpLeftGhost = -1
      DownLeftGhost = -1
      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
         a%npoinGhost = a%npoinGhost + 1_ip
         UpLeftGhost = a%npoinLocal + a%npoinGhost
         UpLeftElem0 = a%nelem
         a%nelem = a%nelem + 1
      endif
      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         a%npoinGhost = a%npoinGhost + 1_ip
         DownLeftGhost = a%npoinLocal + a%npoinGhost
         DownLeftElem0 = a%nelem
         a%nelem = a%nelem + 2
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
         a%npoinGhost = a%npoinGhost + 1_ip
         UpRightGhost = a%npoinLocal + a%npoinGhost
         UpRightElem0 = a%nelem
         a%nelem = a%nelem + 2
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         a%npoinGhost = a%npoinGhost + 1_ip
         DownRightGhost = a%npoinLocal + a%npoinGhost
         DownRightElem0 = a%nelem
         a%nelem = a%nelem + 1
      endif
      
      a%npoin = a%npoinLocal+a%npoinGhost
      a%gnelem = 2*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*a%MPIsize
      
      a%poin0 = a%MPIrank*a%npoinLocal
      a%elem0 = a%MPIrank*2*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)
      
      !Fill local nodes and elements
      call a%Memor%alloc(1_ip,a%pnods,'pnods','ManufacturedMesh')
      call a%Memor%alloc(a%nelem*3_ip,a%lnods,'lnods','ManufacturedMesh')
      call a%Memor%alloc(2_ip,a%npoin,a%coord,'coord','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','ManufacturedMesh')
      a%geoblock=0
      
      call a%Memor%alloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !Build local elements, nodes and coord
      TotalNpoinX = ManufacturedNprocX*a%ManufacturedNpoinX
      ElemCoordSpanX = 1_rp/real(TotalNpoinX-1,rp)
      
      basecoordx = ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(1))
      basecoordy = 1- ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(2))
      
      !Local Coordinates
      do jcol = 1,a%ManufacturedNpoinX
         coordx = basecoordx + (jcol-1)*ElemCoordSpanX
         basepointcol = (jcol-1)*a%ManufacturedNpoinX
         do irow = 1,a%ManufacturedNpoinX
            ipoin = basepointcol+irow   
            coordy = basecoordy - (irow-1)*ElemCoordSpanX
            a%coord(1,ipoin) = coordx
            a%coord(2,ipoin) = coordy
         enddo
      enddo
      
      BlockDown = a%npoinLocal
      BlockUp   = -BlockDown
      BlockRight = a%npoinLocal*ManufacturedNprocX
      BlockLeft = -BlockRight
      rowdown = 1
      rowup = -rowdown
      colright = a%ManufacturedNpoinX
      colleft = -colright
      
      !Local Elements
      a%pnods(1) = 3
      ielem = 0
      do jcol = 1,a%ManufacturedNpoinX-1
         basepointcol = (jcol-1)*a%ManufacturedNpoinX
         do irow = 1,a%ManufacturedNpoinX-1
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + irow
            a%lnods(ispos+2) = basepointcol + irow + rowdown
            a%lnods(ispos+3) = basepointcol + irow + colright
            
            ielem = ielem +1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + irow + rowdown
            a%lnods(ispos+2) = basepointcol + irow + colright + rowdown
            a%lnods(ispos+3) = basepointcol + irow + colright
         enddo
      enddo
      
      !LocalToGlobal for local points
      do ipoin = 1,a%npoinLocal
         LocalToGlobal(ipoin) = a%poin0+ipoin
         ProcessorList(ipoin) = a%MPIrank
      enddo
      !Ghost points and elements
      if (UpGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockUp + (colright-1)
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(UpGhost0 + ipoin) = baseGhost0 + 1 + colright*(ipoin-1)
            ProcessorList(UpGhost0 + ipoin) = a%MPIrank-1
         enddo
         do jcol = 1,a%ManufacturedNpoinX-1
            basepointcol = (jcol-1)*a%ManufacturedNpoinX
            
            ielem = UpElem0 + (jcol-1)*2 + 1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + 1
            a%lnods(ispos+2) = upGhost0+jcol + 1
            a%lnods(ispos+3) = upGhost0+jcol 
            
            ielem = UpElem0 + (jcol-1)*2 + 2
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + 1
            a%lnods(ispos+2) = basepointcol + 1 + colright
            a%lnods(ispos+3) = upGhost0+jcol + 1
         enddo
      endif
      if (DownGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockDown
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(DownGhost0 + ipoin) = baseGhost0 + 1 + colright*(ipoin-1)
            ProcessorList(DownGhost0 + ipoin) = a%MPIrank+1
         enddo
         
         do jcol = 1,a%ManufacturedNpoinX-1
            basepointcol = (jcol-1)*a%ManufacturedNpoinX + a%ManufacturedNpoinX - 1
            
            ielem = DownElem0 + (jcol-1)*2 + 1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = DownGhost0 + jcol
            a%lnods(ispos+2) = basepointcol + 1 + colright
            a%lnods(ispos+3) = basepointcol + 1
            
            ielem = DownElem0 + (jcol-1)*2 + 2
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = DownGhost0 + jcol
            a%lnods(ispos+2) = DownGhost0 + jcol + 1
            a%lnods(ispos+3) = basepointcol + 1 + colright 
         enddo
      endif
      if (LeftGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockLeft + colright*(a%ManufacturedNpoinX-1)
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(LeftGhost0 + ipoin) = baseGhost0 + ipoin
            ProcessorList(LeftGhost0 + ipoin) = a%MPIrank - ManufacturedNprocX
         enddo
         basepointcol = 0
         do irow = 1,a%ManufacturedNpoinX-1
            ielem = LeftElem0 + (irow-1)*2 + 1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = LeftGhost0+irow +1
            a%lnods(ispos+2) = basepointcol + irow
            a%lnods(ispos+3) = LeftGhost0+irow
            
                        
            ielem = LeftElem0 + (irow-1)*2 + 2
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = LeftGhost0+irow +1
            a%lnods(ispos+2) = basepointcol + irow+1
            a%lnods(ispos+3) = basepointcol + irow
         enddo
      endif
      if (RightGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockRight
         do ipoin = 1,a%ManufacturedNpoinX
            LocalToGlobal(RightGhost0 + ipoin) = baseGhost0 + ipoin
            ProcessorList(RightGhost0 + ipoin) = a%MPIrank + ManufacturedNprocX
         enddo
         basepointcol = colright*(a%ManufacturedNpoinX-1)
         do irow = 1,a%ManufacturedNpoinX-1
            ielem = RightElem0 + (irow-1)*2 + 1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + irow+1
            a%lnods(ispos+2) = RightGhost0+irow
            a%lnods(ispos+3) = basepointcol + irow
                        
            ielem = RightElem0 + (irow-1)*2 + 2
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = basepointcol + irow+1
            a%lnods(ispos+2) = RightGhost0+irow + 1
            a%lnods(ispos+3) = RightGhost0+irow 
         enddo
      endif
      
      if (UpLeftGhost /= -1) then
         LocalToGlobal(UpLeftGhost) = a%poin0 + 1 + BlockUp + BlockLeft + colright*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) 
         ProcessorList(UpLeftGhost) = a%MPIrank - 1 - ManufacturedNprocX
         
!          ielem = UpLeftElem0+1
!          ispos = (ielem-1)*a%pnods(1)
!          a%lnods(ispos+1) = LeftGhost0 + 1
!          a%lnods(ispos+2) = UpGhost0+1
!          a%lnods(ispos+3) = UpLeftGhost
         
         ielem = UpLeftElem0+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = LeftGhost0 + 1
         a%lnods(ispos+2) = 1
         a%lnods(ispos+3) = UpGhost0+1
      endif   
      
      if (UpRightGhost /= -1) then
         LocalToGlobal(UpRightGhost) = a%poin0 +1 + BlockUp + BlockRight + rowdown*(a%ManufacturedNpoinX-1)
         ProcessorList(UpRightGhost) = a%MPIrank - 1 + ManufacturedNprocX
         
         ielem = UpRightElem0+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = 1 + colright*(a%ManufacturedNPoinX-1) 
         a%lnods(ispos+2) = UpRightGhost
         a%lnods(ispos+3) = UpGhost0 + a%ManufacturedNpoinX
         
         ielem = UpRightElem0+2
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = 1 + colright*(a%ManufacturedNPoinX-1)
         a%lnods(ispos+2) = RightGhost0+1
         a%lnods(ispos+3) = UpRightGhost
      endif   
      
      if (DownLeftGhost /= -1) then
         LocalToGlobal(DownLeftGhost) = a%poin0 +1 + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1)
         ProcessorList(DownLeftGhost) = a%MPIrank + 1 - ManufacturedNprocX
         
         ielem = DownLeftElem0+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = DownLeftGhost
         a%lnods(ispos+2) = a%ManufacturedNpoinX
         a%lnods(ispos+3) = LeftGhost0 + a%ManufacturedNpoinX
         
         ielem = DownLeftElem0+2
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = DownLeftGhost
         a%lnods(ispos+2) = DownGhost0+1
         a%lnods(ispos+3) = a%ManufacturedNpoinX
      endif
         
      if (DownRightGhost /= -1) then
         LocalToGlobal(DownRightGhost) = a%poin0 +1 + BlockDown + BlockRight 
         ProcessorList(DownRightGhost) = a%MPIrank + 1 + ManufacturedNprocX
         
         ielem = DownRightElem0+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = DownGhost0 + a%ManufacturedNpoinX
         a%lnods(ispos+2) = RightGhost0 + a%ManufacturedNpoinX
         a%lnods(ispos+3) = BlockDown
         
!          ielem = DownRightElem0+2
!          ispos = (ielem-1)*a%pnods(1)
!          a%lnods(ispos+1) = DownGhost0 + a%ManufacturedNpoinX
!          a%lnods(ispos+2) = DownRightGhost
!          a%lnods(ispos+3) = RightGhost0 + a%ManufacturedNpoinX
      endif   
         
         
     
      
            
      
      
      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      call a%Memor%dealloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%dealloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !External normal, we simply allocate it
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ManufacturedMesh')
      
   endif
   
   if (a%ndime == 2 .and. a%ManufacturedType == 'POLAR') then

      InternalRadius = a%ManufacturedInternalRadius 
      ExternalRadius = a%ManufacturedExternalRadius
      PointsInRadius = a%ManufacturedNpoinX 
      DivisionsInAngle  = a%ManufacturedDivisionsInAngle 
      
      DeltaAngleProc = 2_rp*dacos(-1.D0)/real(a%MPIsize,rp)
      DeltaAngle     = DeltaAngleProc/DivisionsInAngle

      kfl_InternalRadiusAtZero = .false.
      if (abs(InternalRadius) <= 1.0D-6) then
          kfl_InternalRadiusAtZero = .true.
          InternalRadius = ExternalRadius/(PointsInRadius-1)
          PointsInRadius = PointsInRadius - 1
      endif
          
      SpanInRadialCoord = (ExternalRadius-InternalRadius)/(PointsInRadius-1)

      LocalAngleDivisionNpoin = 0
      LocalAngleNelem = 0

      do irow = 1,PointsInRadius
         TimesAngleIsDivided = 2_ip**(irow-1)
         LocalAngleNelem = LocalAngleNelem + 3_ip*TimesAngleIsDivided 
         LocalAngleDivisionNpoin = LocalAngleDivisionNpoin + TimesAngleIsDivided
      enddo
      LocalAngleNelem = LocalAngleNelem - 3_ip*2_ip**(PointsInRadius-1) 
      
      a%npoinLocal = LocalAngleDivisionNpoin*DivisionsInAngle
      a%gnpoin = a%npoinLocal*a%MPIsize

      a%nelem = (LocalAngleNelem*DivisionsInAngle) - 2_ip*(PointsInRadius-1) 
      a%gnelem = LocalAngleNelem*DivisionsInAngle*a%MPIsize

      a%poin0 = a%MPIrank*a%npoinLocal
      a%elem0 = a%MPIrank*(LocalAngleNelem*DivisionsInAngle) 

      if (kfl_InternalRadiusAtZero) then 
           a%gnpoin = a%gnpoin + 1 
           a%gnelem = a%gnelem + DivisionsInAngle*a%MPIsize
           a%elem0 = a%MPIrank*(LocalAngleNelem+1)*DivisionsInAngle 
           if (a%MPIrank == a%MPIsize-1) then
               !Im last pc with inner point as local
               a%npoinLocal = a%npoinLocal + 1
               a%nelem = a%nelem + DivisionsInAngle - 1_ip 
           endif
      endif

      a%npoinGhost = 0

      UpGhost0 = -1
      DownGhost0 = -1
      InnerGhost0 = -1
      OtherPcsGhost0 = -1

      !First the upper and lower processors element connections (required for r=0)
      if (a%MPIsize == 1) then
         !Special case of serial run (Local elements at teta=0 self-connection) 
         a%nelem = a%nelem + 2_ip*(PointsInRadius-1)
      else
         !Parallel run
         !Upper processor conection
         UpGhost0 = a%npoinLocal + a%npoinGhost
         UpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + PointsInRadius
         a%nelem = a%nelem + 2_ip*(PointsInRadius-1)
         !Lower processor conection
         DownGhost0 = a%npoinLocal + a%npoinGhost
         DownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + PointsInRadius
         a%nelem = a%nelem + 2_ip*(PointsInRadius-1)
      endif

      !Then the elements at radius equal zero 
      if (kfl_InternalRadiusAtZero) then
         if (a%MPIsize == 1) then
            !Serial run (local element at r=0 self-connection)
            a%nelem = a%nelem + 1_ip
         else
            !Parallel run
            if (a%MPIrank /= a%MPIsize-1) then
               !Im not last pc: Inner processor conection
               !I use localpoints and InnerGhost0
               InnerGhost0 = a%npoinLocal + a%npoinGhost 
               a%npoinGhost = a%npoinGhost + 1
               InnerElem0 = a%nelem
               a%nelem = a%nelem + (DivisionsInAngle-1)
               !Im not last pc: Inner-Upper processor conection
               !I use UpGhost0 and InnerGhost0
               UpInnerElem0 = a%nelem
               a%nelem = a%nelem + 1
               !Im not last pc: Inner-Lower processor conection
               !I use DownGhost0 and InnerGhost0
               DownInnerElem0 = a%nelem
               a%nelem = a%nelem + 1
            else
               !Im last pc: All inner elements of other pcs except upper and lower
               OtherPcsGhost0 = a%npoinLocal + a%npoinGhost
               a%npoinGhost = a%npoinGhost + DivisionsInAngle*(a%MPIsize-1) - 2
               OtherInnerElem0 = a%nelem
               a%nelem = a%nelem + DivisionsInAngle*(a%MPIsize-1) - 1
               !Im last pc: Inner-Upper processor conection
               !I use UpGhost0 and InnerLocalpoint
               UpInnerElem0 = a%nelem
               a%nelem = a%nelem + 1
               !Im last pc: Inner-Lower processor conection
               !I use DownGhost0 and InnerLocalpoint
               DownInnerElem0 = a%nelem
               a%nelem = a%nelem + 1
            endif
         endif
      endif


      a%npoin = a%npoinLocal+a%npoinGhost

      
      !Fill local nodes and elements
      call a%Memor%alloc(1_ip,a%pnods,'pnods','ManufacturedMesh')
      call a%Memor%alloc(a%nelem*3_ip,a%lnods,'lnods','ManufacturedMesh')
      call a%Memor%alloc(2_ip,a%npoin,a%coord,'coord','ManufacturedMesh')
      
      call a%Memor%alloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !Build local elements, nodes and coord

      !Local Coordinates
      BaseProcessorAngle = (a%MPIrank-1)*DeltaAngleProc 
      ipoin = 0
      do jcol=1,DivisionsInAngle
        BaseAngle = BaseProcessorAngle + (jcol-1)*DeltaAngle
        !Radial nodes at BaseAngle (to locate the angular nodes)
        do irow=1,PointsInRadius
          BaseDivisionRadius =  InternalRadius + (irow-1)*SpanInRadialCoord
          coordx = BaseDivisionRadius*cos(BaseAngle)
          coordy = BaseDivisionRadius*sin(BaseAngle)
          ipoin = ipoin + 1
          a%coord(1,ipoin) = coordx
          a%coord(2,ipoin) = coordy
        enddo
        !Angular nodes at each radius
        do irow=1,PointsInRadius
           TimesAngleIsDivided = 2_ip**(irow-1)
           BaseDivisionRadius =  InternalRadius + (irow-1)*SpanInRadialCoord
           do groupi = 1,TimesAngleIsDivided-1
              BaseGroupAngle = BaseAngle + groupi*DeltaAngle/TimesAngleIsDivided
              coordx = BaseDivisionRadius*cos(BaseGroupAngle)
              coordy = BaseDivisionRadius*sin(BaseGroupAngle)
              ipoin = ipoin + 1
              a%coord(1,ipoin) = coordx
              a%coord(2,ipoin) = coordy
           enddo
        enddo
      enddo
   
      if (kfl_InternalRadiusAtZero) then
        if (a%MPIrank == a%MPIsize-1) then
          ipoin = ipoin + 1
          a%coord(1,ipoin) = 0.0_rp
          a%coord(2,ipoin) = 0.0_rp
        endif
      endif

      !LocalToGlobal for local points
      do ipoin = 1,a%npoinLocal
         LocalToGlobal(ipoin) = a%poin0+ipoin
         ProcessorList(ipoin) = a%MPIrank
      enddo

      !Local Elements
      a%pnods(1) = 3
      ielem = 0

      do jcol=1,DivisionsInAngle-1
         ZeroPointAngle = (jcol-1)*LocalAngleDivisionNpoin
         !First radial group
         irow = 1
         bottomleftpoint = ZeroPointAngle + irow 
         rightpoint = bottomleftpoint + PointsInRadius
         bottomrightpoint = bottomleftpoint + 1
         topleftpoint = ZeroPointAngle + LocalAngleDivisionNpoin + 1
         toprightpoint = topleftpoint + irow

         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = bottomleftpoint
         a%lnods(ispos+2) = bottomrightpoint
         a%lnods(ispos+3) = rightpoint
         
         ielem = ielem +1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = bottomleftpoint
         a%lnods(ispos+2) = rightpoint
         a%lnods(ispos+3) = topleftpoint

         ielem = ielem +1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = topleftpoint
         a%lnods(ispos+2) = rightpoint
         a%lnods(ispos+3) = toprightpoint

         !Second radial groups
         irow = 2
         FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
         !Second radial Lower group
         bottomleftpoint = ZeroPointAngle + irow 
         rightpoint = FirstInnerPoint + 2**(irow-1) - 1
         bottomrightpoint = bottomleftpoint + 1
         topleftpoint =  FirstInnerPoint
         toprightpoint = rightpoint + 1

         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = bottomleftpoint
         a%lnods(ispos+2) = bottomrightpoint
         a%lnods(ispos+3) = rightpoint
         
         ielem = ielem +1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = bottomleftpoint
         a%lnods(ispos+2) = rightpoint
         a%lnods(ispos+3) = topleftpoint

         ielem = ielem +1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = topleftpoint
         a%lnods(ispos+2) = rightpoint
         a%lnods(ispos+3) = toprightpoint

         !Second radial Upper group  
         bottomleftpoint = FirstInnerPoint
         rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
         bottomrightpoint = rightpoint -1
         topleftpoint =  ZeroPointAngle + irow + LocalAngleDivisionNpoin
         toprightpoint = topleftpoint + 1

         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = bottomleftpoint
         a%lnods(ispos+2) = bottomrightpoint
         a%lnods(ispos+3) = rightpoint
         
         ielem = ielem +1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = bottomleftpoint
         a%lnods(ispos+2) = rightpoint
         a%lnods(ispos+3) = topleftpoint

         ielem = ielem +1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = topleftpoint
         a%lnods(ispos+2) = rightpoint
         a%lnods(ispos+3) = toprightpoint

         !Other radial groups
         do irow=3,PointsInRadius-1
           TimesAngleIsDivided = 2_ip**(irow-1)
           FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
           !Lower group
           bottomleftpoint = ZeroPointAngle + irow
           rightpoint = FirstInnerPoint + 2**(irow-1) - 1
           bottomrightpoint = bottomleftpoint + 1
           topleftpoint =  FirstInnerPoint
           toprightpoint = rightpoint + 1

           ielem = ielem+1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = bottomrightpoint
           a%lnods(ispos+3) = rightpoint
           
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = topleftpoint
   
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = topleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = toprightpoint

           !Inner groups
           do groupi = 2,TimesAngleIsDivided-1
              bottomleftpoint = FirstInnerPoint + groupi - 2
              rightpoint = FirstInnerPoint + 2**(irow-1) - 1 + (groupi-1)*2
              bottomrightpoint = rightpoint - 1
              topleftpoint =  bottomleftpoint + 1
              toprightpoint = rightpoint + 1
 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = bottomleftpoint
              a%lnods(ispos+2) = bottomrightpoint
              a%lnods(ispos+3) = rightpoint
              
              ielem = ielem +1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = bottomleftpoint
              a%lnods(ispos+2) = rightpoint
              a%lnods(ispos+3) = topleftpoint
      
              ielem = ielem +1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = topleftpoint
              a%lnods(ispos+2) = rightpoint
              a%lnods(ispos+3) = toprightpoint

           enddo
           !Upper group  
           bottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
           rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
           bottomrightpoint = rightpoint - 1
           topleftpoint =  ZeroPointAngle + irow + LocalAngleDivisionNpoin
           toprightpoint = topleftpoint + 1

           ielem = ielem+1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = bottomrightpoint
           a%lnods(ispos+3) = rightpoint
           
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = topleftpoint
   
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = topleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = toprightpoint

         enddo
      enddo

      !This part cares for not connecting top (ghost) points in local elements
      ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
      !First radial group
      irow = 1
      bottomleftpoint = ZeroPointAngle + irow 
      rightpoint = bottomleftpoint + PointsInRadius
      bottomrightpoint = bottomleftpoint + 1

      ielem = ielem+1
      ispos = (ielem-1)*a%pnods(1)
      a%lnods(ispos+1) = bottomleftpoint
      a%lnods(ispos+2) = bottomrightpoint
      a%lnods(ispos+3) = rightpoint
      
      !Second radial groups
      irow = 2
      FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
      !Second radial Lower group
      bottomleftpoint = ZeroPointAngle + irow 
      rightpoint = FirstInnerPoint + 2**(irow-1) - 1
      bottomrightpoint = bottomleftpoint + 1
      topleftpoint =  FirstInnerPoint
      toprightpoint = rightpoint + 1

      ielem = ielem+1
      ispos = (ielem-1)*a%pnods(1)
      a%lnods(ispos+1) = bottomleftpoint
      a%lnods(ispos+2) = bottomrightpoint
      a%lnods(ispos+3) = rightpoint
      
      ielem = ielem +1
      ispos = (ielem-1)*a%pnods(1)
      a%lnods(ispos+1) = bottomleftpoint
      a%lnods(ispos+2) = rightpoint
      a%lnods(ispos+3) = topleftpoint

      ielem = ielem +1
      ispos = (ielem-1)*a%pnods(1)
      a%lnods(ispos+1) = topleftpoint
      a%lnods(ispos+2) = rightpoint
      a%lnods(ispos+3) = toprightpoint

      !Second radial Upper group  
      bottomleftpoint = FirstInnerPoint
      rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
      bottomrightpoint = rightpoint -1

      ielem = ielem+1
      ispos = (ielem-1)*a%pnods(1)
      a%lnods(ispos+1) = bottomleftpoint
      a%lnods(ispos+2) = bottomrightpoint
      a%lnods(ispos+3) = rightpoint
      
      !Other radial groups
      do irow=3,PointsInRadius-1
        TimesAngleIsDivided = 2_ip**(irow-1)
        FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
        !Lower group
        bottomleftpoint = ZeroPointAngle + irow
        rightpoint = FirstInnerPoint + 2**(irow-1) - 1
        bottomrightpoint = bottomleftpoint + 1
        topleftpoint =  FirstInnerPoint
        toprightpoint = rightpoint + 1

        ielem = ielem+1
        ispos = (ielem-1)*a%pnods(1)
        a%lnods(ispos+1) = bottomleftpoint
        a%lnods(ispos+2) = bottomrightpoint
        a%lnods(ispos+3) = rightpoint
        
        ielem = ielem +1
        ispos = (ielem-1)*a%pnods(1)
        a%lnods(ispos+1) = bottomleftpoint
        a%lnods(ispos+2) = rightpoint
        a%lnods(ispos+3) = topleftpoint
   
        ielem = ielem +1
        ispos = (ielem-1)*a%pnods(1)
        a%lnods(ispos+1) = topleftpoint
        a%lnods(ispos+2) = rightpoint
        a%lnods(ispos+3) = toprightpoint

        !Inner groups
        do groupi = 2,TimesAngleIsDivided-1
           bottomleftpoint = FirstInnerPoint + groupi - 2
           rightpoint = FirstInnerPoint + 2**(irow-1) - 1 + (groupi-1)*2
           bottomrightpoint = rightpoint - 1
           topleftpoint =  bottomleftpoint + 1
           toprightpoint = rightpoint + 1
 
           ielem = ielem+1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = bottomrightpoint
           a%lnods(ispos+3) = rightpoint
           
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = topleftpoint
      
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = topleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = toprightpoint

        enddo
        !Upper group  
        bottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
        rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
        bottomrightpoint = rightpoint - 1

        ielem = ielem+1
        ispos = (ielem-1)*a%pnods(1)
        a%lnods(ispos+1) = bottomleftpoint
        a%lnods(ispos+2) = bottomrightpoint
        a%lnods(ispos+3) = rightpoint
        
      enddo

      !This part creates the zero radius local elements belonging to the last pc
      if (kfl_InternalRadiusAtZero) then
        if (a%MPIrank == a%MPIsize-1) then 
          do jcol=1,DivisionsInAngle-1
             ZeroPointAngle = (jcol-1)*LocalAngleDivisionNpoin
             !First radial group
             irow = 1
             bottomleftpoint = a%npoinLocal 
             bottomrightpoint = ZeroPointAngle + 1
             toprightpoint = jcol*LocalAngleDivisionNpoin + 1

             ielem = ielem+1
             ispos = (ielem-1)*a%pnods(1)
             a%lnods(ispos+1) = bottomleftpoint
             a%lnods(ispos+2) = bottomrightpoint
             a%lnods(ispos+3) = toprightpoint
          enddo
        endif
      endif


      !Ghost points and elements 
        
      !First the upper and lower processors element connections (required for r=0)
      if (a%MPIsize == 1) then
           !Special case of serial run (Local elements self-connection at teta=0) 
           ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
           !First radial group
           irow = 1
           bottomleftpoint = ZeroPointAngle + irow 
           rightpoint = bottomleftpoint + PointsInRadius
           topleftpoint = irow
           toprightpoint = topleftpoint + 1
  
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = topleftpoint
  
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = topleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = toprightpoint

           !Second radial groups
           irow = 2
           FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
           !Second radial Upper group  
           bottomleftpoint = FirstInnerPoint
           rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
           topleftpoint =  irow
           toprightpoint = topleftpoint + 1
      
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = bottomleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = topleftpoint
      
           ielem = ielem +1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = topleftpoint
           a%lnods(ispos+2) = rightpoint
           a%lnods(ispos+3) = toprightpoint

           !Other radial groups
           do irow=3,PointsInRadius-1
             TimesAngleIsDivided = 2_ip**(irow-1)
             FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
             !Upper group  
             bottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
             rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
             bottomrightpoint = rightpoint - 1
             topleftpoint =  irow
             toprightpoint = topleftpoint + 1
      
             ielem = ielem +1
             ispos = (ielem-1)*a%pnods(1)
             a%lnods(ispos+1) = bottomleftpoint
             a%lnods(ispos+2) = rightpoint
             a%lnods(ispos+3) = topleftpoint
            
             ielem = ielem +1
             ispos = (ielem-1)*a%pnods(1)
             a%lnods(ispos+1) = topleftpoint
             a%lnods(ispos+2) = rightpoint
             a%lnods(ispos+3) = toprightpoint
            
           enddo

      else
            !Parallel run

            !Upper processor conection

            !GhostPoints Global numbering
            if (a%MPIrank < a%MPIsize-1) then
               baseGhost0  = a%poin0 + a%npoinLocal
               do ipoin = 1,PointsInRadius
                  LocalToGlobal(UpGhost0 + ipoin) = baseGhost0 + ipoin
                  ProcessorList(UpGhost0 + ipoin) = a%MPIrank + 1
               enddo
            else 
               !Special case of being the last processor 
               baseGhost0 = 0
               do ipoin = 1,PointsInRadius
                  LocalToGlobal(UpGhost0 + ipoin) = baseGhost0 + ipoin
                  ProcessorList(UpGhost0 + ipoin) = 1
               enddo
            endif

            ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin

            !First radial group
            irow = 1
            bottomleftpoint = ZeroPointAngle + irow 
            rightpoint = bottomleftpoint + PointsInRadius
            topleftpoint = UpGhost0 + irow
            toprightpoint = topleftpoint + 1
   
            ispos = UpElem0*a%pnods(1)
            a%lnods(ispos+1) = bottomleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = topleftpoint
   
            ispos = (UpElem0+1)*a%pnods(1)
            a%lnods(ispos+1) = topleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = toprightpoint

            !Second radial group
            irow = 2
            FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)

            !Second radial Upper group  
            bottomleftpoint = FirstInnerPoint
            rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
            topleftpoint =  UpGhost0 + irow 
            toprightpoint = topleftpoint + 1
   
            ispos = (UpElem0+2)*a%pnods(1)
            a%lnods(ispos+1) = bottomleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = topleftpoint
   
            ispos = (UpElem0+3)*a%pnods(1)
            a%lnods(ispos+1) = topleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = toprightpoint
   

            !Other radial groups
            do irow=3,PointsInRadius-1
              TimesAngleIsDivided = 2_ip**(irow-1)
              FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
              !Upper group  
              bottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
              rightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
              topleftpoint = UpGhost0 + irow;
              toprightpoint = topleftpoint + 1;
      
              ispos = (UpElem0 + (irow-1)*2)*a%pnods(1)
              a%lnods(ispos+1) = bottomleftpoint
              a%lnods(ispos+2) = rightpoint
              a%lnods(ispos+3) = topleftpoint
      
              ispos = (UpElem0 + (irow-1)*2 +1)*a%pnods(1)
              a%lnods(ispos+1) = topleftpoint
              a%lnods(ispos+2) = rightpoint
              a%lnods(ispos+3) = toprightpoint
              
            enddo

            !Lower processor conection

            !GhostPoints Global numbering
            if (a%MPIrank > 0) then
               baseGhost0  = a%poin0 - LocalAngleDivisionNpoin
               do ipoin = 1,PointsInRadius
                  ProcessorList(DownGhost0 + ipoin) = a%MPIrank - 1
               enddo
            else 
               !Special case of being the first processor 
               baseGhost0 = a%npoinLocal*a%MPIsize - LocalAngleDivisionNpoin
               do ipoin = 1,PointsInRadius
                  ProcessorList(DownGhost0 + ipoin) = a%MPIsize - 1 
               enddo
            endif

            !First radial group
            irow = 1
            LocalToGlobal(DownGhost0 + irow) = baseGhost0 + irow
            bottomleftpoint = DownGhost0 + irow 
            rightpoint = bottomleftpoint + 1
            topleftpoint = irow
            toprightpoint = topleftpoint + 1
        
            ispos = (DownElem0)*a%pnods(1)
            a%lnods(ispos+1) = bottomleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = topleftpoint
        
            ispos = (DownElem0+1)*a%pnods(1)
            a%lnods(ispos+1) = topleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = toprightpoint
        
            !Second radial groups
            irow = 2
            FirstInnerPoint = baseGhost0 + PointsInRadius + 2**(irow-2)
            LocalToGlobal(DownGhost0 + irow) = FirstInnerPoint
        
            !Second radial Upper group  
            bottomleftpoint = DownGhost0 + irow 
            rightpoint = bottomleftpoint + 1
            topleftpoint = irow 
            toprightpoint = topleftpoint + 1
        
            ispos = (DownElem0+2)*a%pnods(1)
            a%lnods(ispos+1) = bottomleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = topleftpoint
        
            ispos = (DownElem0+3)*a%pnods(1)
            a%lnods(ispos+1) = topleftpoint
            a%lnods(ispos+2) = rightpoint
            a%lnods(ispos+3) = toprightpoint
        
            !Other radial groups
            do irow=3,PointsInRadius-1
              TimesAngleIsDivided = 2_ip**(irow-1)
              FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
              LocalToGlobal(DownGhost0 + irow) = FirstInnerPoint + 2**(irow-1) - 2
              !Upper group  
              bottomleftpoint = DownGhost0 + irow 
              rightpoint = bottomleftpoint + 1
              topleftpoint = irow 
              toprightpoint = topleftpoint + 1
        
              ispos = (DownElem0+(irow-1)*2)*a%pnods(1)
              a%lnods(ispos+1) = bottomleftpoint
              a%lnods(ispos+2) = rightpoint
              a%lnods(ispos+3) = topleftpoint
           
              ispos = (DownElem0+(irow-1)*2 +1)*a%pnods(1)
              a%lnods(ispos+1) = topleftpoint
              a%lnods(ispos+2) = rightpoint
              a%lnods(ispos+3) = toprightpoint
              
            enddo

            TimesAngleIsDivided = 2_ip**(PointsInRadius-1)
            FirstInnerPoint = FirstInnerPoint + 2**(PointsInRadius-2) - 1
            LocalToGlobal(DownGhost0 + PointsInRadius) = FirstInnerPoint + 2**(PointsInRadius-1) - 2
      
      endif

      !Then the ghost points of elements at radius equal zero
      if (kfl_InternalRadiusAtZero) then
         if (a%MPIsize == 1) then
               !Serial run: Inside the zero radius local elements self-connection
               !Special case of serial run self-connection, anyway inside this logic 
               ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
               !First radial group
               irow = 1
               bottomleftpoint = a%npoinLocal 
               bottomrightpoint = ZeroPointAngle + 1
               toprightpoint = irow

               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = bottomleftpoint
               a%lnods(ispos+2) = bottomrightpoint
               a%lnods(ispos+3) = toprightpoint
          else
            !Parallel run !Therefore r_0 = 0 (inner, upper and lower) element connections between processors
            
            !First the inner elements conection between processors 
            if (a%MPIrank /= a%MPIsize-1) then
               !Not last pc: Inside the inner elements parallel pc connection to last pc 
               !GhostPoints Global numbering
               LocalToGlobal(InnerGhost0 + 1) = a%gnpoin 
               ProcessorList(InnerGhost0 + 1) = a%MPIsize-1
               
               do jcol=1,DivisionsInAngle-1
                  ZeroPointAngle = (jcol-1)*LocalAngleDivisionNpoin
                  !First radial group
                  irow = 1
                  bottomleftpoint = InnerGhost0 + 1 
                  bottomrightpoint = ZeroPointAngle + 1
                  toprightpoint = jcol*LocalAngleDivisionNpoin + 1
      
                  ielem = InnerElem0 + jcol
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = bottomleftpoint
                  a%lnods(ispos+2) = bottomrightpoint
                  a%lnods(ispos+3) = toprightpoint
               enddo

               !Not last pc: Inside the inner element parallel connection with upper pc and last pc 
               ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
               !First radial group
               irow = 1
               bottomleftpoint = InnerGhost0 + 1 
               bottomrightpoint = ZeroPointAngle + 1
               toprightpoint = UpGhost0 + 1
      
               ielem = UpInnerElem0 + 1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = bottomleftpoint
               a%lnods(ispos+2) = bottomrightpoint
               a%lnods(ispos+3) = toprightpoint

               !Not last pc: Inside the inner element parallel connection with lower pc and last pc 
               !First radial group
               irow = 1
               bottomleftpoint = InnerGhost0 + 1 
               bottomrightpoint = DownGhost0 + 1
               toprightpoint = 1
      
               ielem = DownInnerElem0 + 1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = bottomleftpoint
               a%lnods(ispos+2) = bottomrightpoint
               a%lnods(ispos+3) = toprightpoint
            else
               !Last pc: Inside the inner element connection to inner ghosts in other pcs
               !OtherPcsGhost0  
               !Last pc: Connection with all inner elements of other pcs
               do groupi = 1,a%MPIsize-1
                  baseGhost0 = LocalAngleDivisionNpoin*DivisionsInAngle*(groupi-1)
                  !All inner elements 
                  !care for the lower group that is connected with the lower to the inner other pc 
                  if (groupi == 1) then
                     bottomleftpoint = a%npoinLocal
                     bottomrightpoint = UpGhost0 + 1 
                     toprightpoint = OtherPcsGhost0 + 1

                     ielem = OtherInnerElem0 + 1
                  else
                     LocalToGlobal(OtherPcsGhost0 + DivisionsInAngle*(groupi-1)) = baseGhost0 + 1
                     ProcessorList(OtherPcsGhost0 + DivisionsInAngle*(groupi-1)) = groupi

                     bottomleftpoint = a%npoinLocal
                     bottomrightpoint = OtherPcsGhost0 + DivisionsInAngle*(groupi-1) 
                     toprightpoint = OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + 1

                     ielem = OtherInnerElem0 + DivisionsInAngle*(groupi-1) + 1
                  endif
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = bottomleftpoint
                  a%lnods(ispos+2) = bottomrightpoint
                  a%lnods(ispos+3) = toprightpoint
                  !All other inner elements 
                  do jcol=2,DivisionsInAngle-1
                     ZeroPointAngle = (jcol-1)*LocalAngleDivisionNpoin
                     if ((groupi == a%MPIsize-1) .and. (jcol == DivisionsInAngle-1)) then
                        LocalToGlobal(OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol - 1) = baseGhost0 + ZeroPointAngle + 1
                        ProcessorList(OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol - 1) = groupi
                        bottomleftpoint = a%npoinLocal
                        bottomrightpoint = OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol - 1 
                        toprightpoint = DownGhost0 + 1 
      
                     else
                        LocalToGlobal(OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol - 1) = baseGhost0 + ZeroPointAngle + 1
                        ProcessorList(OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol - 1) = groupi
                        bottomleftpoint = a%npoinLocal
                        bottomrightpoint = OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol - 1 
                        toprightpoint = OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + jcol
      
                     endif
                     ielem = OtherInnerElem0 + DivisionsInAngle*(groupi-1) + jcol
                     ispos = (ielem-1)*a%pnods(1)
                     a%lnods(ispos+1) = bottomleftpoint
                     a%lnods(ispos+2) = bottomrightpoint
                     a%lnods(ispos+3) = toprightpoint
                  enddo
                  !care for the top group that is connected with the upper to the inner other pc 
                  ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
                  if (groupi /= a%MPIsize-1) then
                        LocalToGlobal(OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + DivisionsInAngle - 1) = baseGhost0 + ZeroPointAngle + 1
                        ProcessorList(OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + DivisionsInAngle - 1) = groupi
                        bottomleftpoint = a%npoinLocal
                        bottomrightpoint = OtherPcsGhost0 + DivisionsInAngle*(groupi-1) + DivisionsInAngle - 1
                        toprightpoint = OtherPcsGhost0 + DivisionsInAngle*groupi

                        ielem = OtherInnerElem0 + DivisionsInAngle*(groupi-1) + DivisionsInAngle
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = bottomleftpoint
                        a%lnods(ispos+2) = bottomrightpoint
                        a%lnods(ispos+3) = toprightpoint
                  endif
               enddo
               !Im last pc: Inner-Upper processor conection
               !I use UpGhost0 and InnerLocalpoint
               ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
               !First radial group
               irow = 1
               bottomleftpoint = a%npoinLocal
               bottomrightpoint = ZeroPointAngle + 1  
               toprightpoint = UpGhost0 + 1
   
               ielem = UpInnerElem0 + 1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = bottomleftpoint
               a%lnods(ispos+2) = bottomrightpoint
               a%lnods(ispos+3) = toprightpoint

               !Im last pc: Inner-Lower processor conection
               !I use DownGhost0 and InnerLocalpoint
               !First radial group
               irow = 1
               bottomleftpoint = a%npoinLocal
               bottomrightpoint = DownGhost0 + 1  
               toprightpoint = 1
   
               ielem = DownInnerElem0 + 1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = bottomleftpoint
               a%lnods(ispos+2) = bottomrightpoint
               a%lnods(ispos+3) = toprightpoint
            endif
             
          endif
     endif
  

      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      call a%Memor%dealloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%dealloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !External normal, we simply allocate it
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ManufacturedMesh')

   endif
   
   if (a%ndime == 3 .and. a%ManufacturedType == 'LINHE') then
      
      if ( nint((real(a%MPIsize))**(1.0/3.0))**3 /= a%MPIsize) then
         kfl_isPowerOf = .false.
      else
         kfl_isPowerOf = .true.
      endif
      
      if (kfl_isPowerOf .eqv. .false.) then
         call runend('3d ManufacturedMesh, MPIsize the cube of a natural number')
      endif
      
      ManufacturedNprocX = (real(a%MPIsize,rp))**(1.0/3.0)
      
      ManufacturedProcPosition(3) = a%MPIrank/(ManufacturedNprocX*ManufacturedNprocX)
      ManufacturedProcPosition(2) = mod(a%MPIrank,ManufacturedNprocX)
      ManufacturedProcPosition(1) = (a%MPIrank-ManufacturedProcPosition(3)*ManufacturedNprocX*ManufacturedNprocX)/ManufacturedNprocX

      
      a%gnpoin = a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%MPIsize
      a%npoinLocal = a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%ManufacturedNpoinX
      
      a%nelem = (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      a%npoinGhost = 0

      BottomGhost0 = -1
      TopGhost0 = -1
      UpGhost0 = -1
      DownGhost0 = -1
      LeftGhost0 = -1
      RightGhost0 = -1

      if (ManufacturedProcPosition(3) /= 0) then
         BottomGhost0 = a%npoinLocal + a%npoinGhost
         BottomElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1) then
         TopGhost0 = a%npoinLocal + a%npoinGhost
         TopElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(2) /= 0) then
         UpGhost0 = a%npoinLocal + a%npoinGhost
         UpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownGhost0 = a%npoinLocal + a%npoinGhost
         DownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= 0) then
         LeftGhost0 = a%npoinLocal + a%npoinGhost
         LeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         RightGhost0 = a%npoinLocal + a%npoinGhost
         RightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      
      UpRightGhost = -1
      DownRightGhost = -1
      UpLeftGhost = -1
      DownLeftGhost = -1
      TopLeftGhost = -1
      TopRightGhost = -1
      TopUpGhost = -1
      TopDownGhost = -1
      BottomLeftGhost = -1
      BottomRightGhost = -1
      BottomUpGhost = -1
      BottomDownGhost = -1

      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
         UpLeftGhost = a%npoinLocal + a%npoinGhost
         UpLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownLeftGhost = a%npoinLocal + a%npoinGhost
         DownLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
         UpRightGhost = a%npoinLocal + a%npoinGhost
         UpRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownRightGhost = a%npoinLocal + a%npoinGhost
         DownRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(1) /= 0) then
         TopLeftGhost = a%npoinLocal + a%npoinGhost
         TopLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         TopRightGhost = a%npoinLocal + a%npoinGhost
         TopRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
         TopUpGhost = a%npoinLocal + a%npoinGhost
         TopUpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         TopDownGhost = a%npoinLocal + a%npoinGhost
         TopDownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(1) /= 0) then
         BottomLeftGhost = a%npoinLocal + a%npoinGhost
         BottomLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         BottomRightGhost = a%npoinLocal + a%npoinGhost
         BottomRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
         BottomUpGhost = a%npoinLocal + a%npoinGhost
         BottomUpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         BottomDownGhost = a%npoinLocal + a%npoinGhost
         BottomDownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + (a%ManufacturedNpoinX-1)
      endif

      TopUpLeftGhost = -1
      TopUpRightGhost = -1
      TopDownLeftGhost = -1
      TopDownRightGhost = -1
      BottomUpLeftGhost = -1
      BottomUpRightGhost = -1
      BottomDownLeftGhost = -1
      BottomDownRightGhost = -1

      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1) then
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      TopUpLeftGhost = a%npoinLocal + a%npoinGhost
		      TopUpLeftElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      TopDownLeftGhost = a%npoinLocal + a%npoinGhost
		      TopDownLeftElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      TopUpRightGhost = a%npoinLocal + a%npoinGhost
		      TopUpRightElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      TopDownRightGhost = a%npoinLocal + a%npoinGhost
		      TopDownRightElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
      endif
      if (ManufacturedProcPosition(3) /= 0) then
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      BottomUpLeftGhost = a%npoinLocal + a%npoinGhost
		      BottomUpLeftElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      BottomDownLeftGhost = a%npoinLocal + a%npoinGhost
		      BottomDownLeftElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      BottomUpRightGhost = a%npoinLocal + a%npoinGhost
		      BottomUpRightElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      BottomDownRightGhost = a%npoinLocal + a%npoinGhost
		      BottomDownRightElem0 = a%nelem
		      a%nelem = a%nelem + 1
	      endif
      endif
      
      a%npoin = a%npoinLocal+a%npoinGhost
      a%gnelem = (a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*a%MPIsize
      
      a%poin0 = a%MPIrank*a%npoinLocal
      a%elem0 = a%MPIrank*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)
      
      !Fill local nodes and elements
      call a%Memor%alloc(1_ip,a%pnods,'pnods','ManufacturedMesh')
      call a%Memor%alloc(a%nelem*8_ip,a%lnods,'lnods','ManufacturedMesh')
      call a%Memor%alloc(3_ip,a%npoin,a%coord,'coord','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','ManufacturedMesh')
      a%geoblock=0
      
      call a%Memor%alloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !Build local elements, nodes and coord
      TotalNpoinX = ManufacturedNprocX*a%ManufacturedNpoinX
      ElemCoordSpanX = 1_rp/real(TotalNpoinX-1,rp)
      
      basecoordx = ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(1))
      basecoordy = 1- ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(2))
      basecoordz = ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(3))
      
      !Local Coordinates
      do kflr = 1,a%ManufacturedNpoinX
	      coordz = basecoordz + (kflr-1)*ElemCoordSpanX
	      basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      do jcol = 1,a%ManufacturedNpoinX
		      coordx = basecoordx + (jcol-1)*ElemCoordSpanX
		      basepointcol = (jcol-1)*a%ManufacturedNpoinX
		      do irow = 1,a%ManufacturedNpoinX
			      ipoin = basepointflr+basepointcol+irow   
			      coordy = basecoordy - (irow-1)*ElemCoordSpanX
			      a%coord(1,ipoin) = coordx
			      a%coord(2,ipoin) = coordy
			      a%coord(3,ipoin) = coordz
		      enddo
	      enddo
      enddo
      
      BlockDown = a%npoinLocal
      BlockUp   = -BlockDown
      BlockRight = a%npoinLocal*ManufacturedNprocX
      BlockLeft = -BlockRight
      BlockTop = a%npoinLocal*ManufacturedNprocX*ManufacturedNprocX
      BlockBottom = -BlockTop

      rowdown = 1
      rowup = -rowdown
      colright = a%ManufacturedNpoinX
      colleft = -colright
      floortop = a%ManufacturedNpoinX*a%ManufacturedNpoinX
      floorbottom = -floortop
      
      !Local Elements
      a%pnods(1) = 8
      ielem = 0
      do kflr = 1,a%ManufacturedNpoinX-1
	      basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      do jcol = 1,a%ManufacturedNpoinX-1
		      basepointcol = (jcol-1)*a%ManufacturedNpoinX
		      do irow = 1,a%ManufacturedNpoinX-1
			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + floortop + basepointcol  + irow + rowdown
			      a%lnods(ispos+2) = basepointflr + floortop + basepointcol  + irow + colright + rowdown
			      a%lnods(ispos+3) = basepointflr + basepointcol + irow + colright + rowdown
			      a%lnods(ispos+4) = basepointflr + basepointcol + irow + rowdown
			      a%lnods(ispos+5) = basepointflr + floortop + basepointcol  + irow 
			      a%lnods(ispos+6) = basepointflr + floortop + basepointcol  + irow + colright
			      a%lnods(ispos+7) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+8) = basepointflr + basepointcol + irow

		      enddo
	      enddo
      enddo
      
      !LocalToGlobal for local points
      do ipoin = 1,a%npoinLocal
         LocalToGlobal(ipoin) = a%poin0+ipoin
         ProcessorList(ipoin) = a%MPIrank
      enddo
      !Ghost points and elements
      if (BottomGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockBottom + floortop*(a%ManufacturedNpoinX-1)
	 do ipoin = 1,a%ManufacturedNpoinX*a%ManufacturedNpoinX
		 LocalToGlobal(BottomGhost0 + ipoin) = baseGhost0 + ipoin
		 ProcessorList(BottomGhost0 + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX
	 enddo
	 do jcol = 1,a%ManufacturedNpoinX-1
		 basepointcol = (jcol-1)*a%ManufacturedNpoinX
		 do irow = 1,a%ManufacturedNpoinX-1

			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1) + irow
			 ispos = (ielem-1)*a%pnods(1)
			 a%lnods(ispos+1) = basepointcol + irow + rowdown
			 a%lnods(ispos+2) = basepointcol + colright + irow + rowdown 
			 a%lnods(ispos+3) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+4) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1  
			 a%lnods(ispos+5) = basepointcol + irow
			 a%lnods(ispos+6) = basepointcol + colright + irow
			 a%lnods(ispos+7) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
			 a%lnods(ispos+8) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow

		 enddo
	 enddo
      endif
      if (TopGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockTop 
	 do ipoin = 1,a%ManufacturedNpoinX*a%ManufacturedNpoinX
		 LocalToGlobal(TopGhost0 + ipoin) = baseGhost0 + ipoin
		 ProcessorList(TopGhost0 + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX
	 enddo
         basepointflr = floortop*(a%ManufacturedNpoinX-1)
	 do jcol = 1,a%ManufacturedNpoinX-1
		 basepointcol = (jcol-1)*a%ManufacturedNpoinX
		 do irow = 1,a%ManufacturedNpoinX-1

			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1) + irow
			 ispos = (ielem-1)*a%pnods(1)
			 a%lnods(ispos+1) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+2) = TopGhost0 + jcol*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+3) = basepointflr + basepointcol + colright + irow + rowdown 
			 a%lnods(ispos+4) = basepointflr + basepointcol + irow + rowdown
			 a%lnods(ispos+5) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow
			 a%lnods(ispos+6) = TopGhost0 + jcol*a%ManufacturedNpoinX + irow
			 a%lnods(ispos+7) = basepointflr + basepointcol + colright + irow
			 a%lnods(ispos+8) = basepointflr + basepointcol + irow

		 enddo
	 enddo
      endif
      if (UpGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + 1 + BlockUp + rowdown*(a%ManufacturedNpoinX-1) 
	 do jpoin = 1,a%ManufacturedNpoinX
		 do ipoin = 1,a%ManufacturedNpoinX
			 LocalToGlobal(UpGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + colright*(ipoin-1)
			 ProcessorList(UpGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank-1
		 enddo
	 enddo
	 do kflr = 1,a%ManufacturedNpoinX-1
		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
		 do jcol = 1,a%ManufacturedNpoinX-1
			 basepointcol = (jcol-1)*a%ManufacturedNpoinX

			 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1) + jcol
			 ispos = (ielem-1)*a%pnods(1)
			 a%lnods(ispos+1) = basepointflr + floortop + basepointcol + 1
			 a%lnods(ispos+2) = basepointflr + floortop + basepointcol + 1 + colright 
			 a%lnods(ispos+3) = basepointflr + basepointcol + 1 + colright 
			 a%lnods(ispos+4) = basepointflr + basepointcol + 1
			 a%lnods(ispos+5) = UpGhost0 + kflr*a%ManufacturedNpoinX + jcol
			 a%lnods(ispos+6) = UpGhost0 + kflr*a%ManufacturedNpoinX + jcol + 1
			 a%lnods(ispos+7) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
			 a%lnods(ispos+8) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol

		 enddo
	 enddo
      endif
      if (DownGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + 1 + BlockDown
	 do jpoin = 1,a%ManufacturedNpoinX
		 do ipoin = 1,a%ManufacturedNpoinX
			 LocalToGlobal(DownGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + colright*(ipoin-1)
			 ProcessorList(DownGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank+1
		 enddo
	 enddo
	 do kflr = 1,a%ManufacturedNpoinX-1
		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
		 do jcol = 1,a%ManufacturedNpoinX-1
			 basepointcol = (jcol-1)*a%ManufacturedNpoinX + a%ManufacturedNpoinX - 1

			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1) + jcol
			 ispos = (ielem-1)*a%pnods(1)
			 a%lnods(ispos+1) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
			 a%lnods(ispos+2) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol + 1
			 a%lnods(ispos+3) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
			 a%lnods(ispos+4) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol
			 a%lnods(ispos+5) = basepointflr + floortop + basepointcol + 1 
			 a%lnods(ispos+6) = basepointflr + floortop + basepointcol + 1 + colright 
			 a%lnods(ispos+7) = basepointflr + basepointcol + 1 + colright 
			 a%lnods(ispos+8) = basepointflr + basepointcol + 1 
		 enddo
	 enddo
      endif
      if (LeftGhost0 /= -1) then
         !GhostPoints Global numbering
	 baseGhost0 = a%poin0 + BlockLeft + colright*(a%ManufacturedNpoinX-1)
	 do jpoin = 1,a%ManufacturedNpoinX
		 do ipoin = 1,a%ManufacturedNpoinX
			 LocalToGlobal(LeftGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + ipoin
			 ProcessorList(LeftGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank - ManufacturedNprocX
		 enddo
	 enddo
         basepointcol = 0
	 do kflr = 1,a%ManufacturedNpoinX-1
		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
		 do irow = 1,a%ManufacturedNpoinX-1
			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1) + irow
			 ispos = (ielem-1)*a%pnods(1)
			 a%lnods(ispos+1) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+2) = basepointflr + floortop + irow + rowdown 
			 a%lnods(ispos+3) = basepointflr + irow + rowdown 
			 a%lnods(ispos+4) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+5) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow
			 a%lnods(ispos+6) = basepointflr + floortop + irow
			 a%lnods(ispos+7) = basepointflr + irow
			 a%lnods(ispos+8) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
		 enddo
	 enddo
      endif
      if (RightGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockRight
	 do jpoin = 1,a%ManufacturedNpoinX
		 do ipoin = 1,a%ManufacturedNpoinX
			 LocalToGlobal(RightGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + ipoin
			 ProcessorList(RightGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank + ManufacturedNprocX
		 enddo
	 enddo
         basepointcol = colright*(a%ManufacturedNpoinX-1)
	 do kflr = 1,a%ManufacturedNpoinX-1
		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
		 do irow = 1,a%ManufacturedNpoinX-1
			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1) + irow
			 ispos = (ielem-1)*a%pnods(1)
			 a%lnods(ispos+1) = basepointflr + floortop + basepointcol + irow + rowdown  
			 a%lnods(ispos+2) = RightGhost0 + kflr*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+3) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow + 1
			 a%lnods(ispos+4) = basepointflr + basepointcol + irow + rowdown  
			 a%lnods(ispos+5) = basepointflr + floortop + basepointcol + irow 
			 a%lnods(ispos+6) = RightGhost0 + kflr*a%ManufacturedNpoinX + irow
			 a%lnods(ispos+7) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
			 a%lnods(ispos+8) = basepointflr + basepointcol + irow 

		 enddo
	 enddo
      endif



      
      if (UpLeftGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + 1 + BlockUp + BlockLeft + colright*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(UpLeftGhost + ipoin) = baseGhost0 + (ipoin-1)*floortop 
		      ProcessorList(UpLeftGhost + ipoin) = a%MPIrank - 1 - ManufacturedNprocX
	      enddo

	      do kflr = 1,a%ManufacturedNpoinX-1
		      ielem = UpLeftElem0 + kflr
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+2) = kflr*floortop + 1
		      a%lnods(ispos+3) = (kflr-1)*floortop + 1
		      a%lnods(ispos+4) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+5) = UpLeftGhost + kflr + 1
		      a%lnods(ispos+6) = UpGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+7) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+8) = UpLeftGhost + kflr 
	      enddo
      endif   
      
      if (DownLeftGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + 1 + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1)  
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(DownLeftGhost + ipoin) = baseGhost0 + (ipoin-1)*floortop 
		      ProcessorList(DownLeftGhost + ipoin) = a%MPIrank + 1 - ManufacturedNprocX
	      enddo

	      do kflr = 1,a%ManufacturedNpoinX-1
		      ielem = DownLeftElem0 + kflr
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownLeftGhost + kflr + 1 
		      a%lnods(ispos+2) = DownGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+4) = DownLeftGhost + kflr 
		      a%lnods(ispos+5) = LeftGhost0 + (kflr+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+6) = kflr*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+7) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+8) = LeftGhost0 + kflr*a%ManufacturedNpoinX
	      enddo
         
      endif
      
      if (UpRightGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + 1 + BlockUp + BlockRight + rowdown*(a%ManufacturedNpoinX-1) 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(UpRightGhost + ipoin) = baseGhost0 + (ipoin-1)*floortop
		      ProcessorList(UpRightGhost + ipoin) = a%MPIrank - 1 + ManufacturedNprocX
	      enddo

	      do kflr = 1,a%ManufacturedNpoinX-1
		      ielem = UpRightElem0 + kflr
		      ispos = (ielem-1)*a%pnods(1)
                      a%lnods(ispos+1) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
		      a%lnods(ispos+2) = RightGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
                      a%lnods(ispos+4) = (kflr-1)*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
                      a%lnods(ispos+5) = UpGhost0 + (kflr+1)*a%ManufacturedNpoinX
                      a%lnods(ispos+6) = UpRightGhost + kflr + 1
                      a%lnods(ispos+7) = UpRightGhost + kflr
                      a%lnods(ispos+8) = UpGhost0 + kflr*a%ManufacturedNpoinX
	      enddo
         
      endif   
         
      if (DownRightGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + 1 + BlockDown + BlockRight 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(DownRightGhost + ipoin) = baseGhost0 + (ipoin-1)*floortop
		      ProcessorList(DownRightGhost + ipoin) = a%MPIrank + 1 + ManufacturedNprocX
	      enddo
         
	      do kflr = 1,a%ManufacturedNpoinX-1
		      ielem = DownRightElem0 + kflr
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownGhost0 + (kflr+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+2) = DownRightGhost + kflr + 1
		      a%lnods(ispos+3) = DownRightGhost + kflr
		      a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+5) = (kflr+1)*a%ManufacturedNPoinX*a%ManufacturedNPoinX
		      a%lnods(ispos+6) = RightGhost0 + (kflr+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+7) = RightGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+8) = kflr*a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      enddo

      endif   
         
      if (TopLeftGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + BlockTop + BlockLeft + colright*(a%ManufacturedNpoinX-1)
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(TopLeftGhost + ipoin) = baseGhost0 + ipoin
		      ProcessorList(TopLeftGhost + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX - ManufacturedNprocX
	      enddo
         
	      do irow = 1,a%ManufacturedNpoinX-1
		      ielem = TopLeftElem0 + irow
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = TopLeftGhost + irow + 1
		      a%lnods(ispos+2) = TopGhost0 + irow + 1 
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+4) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1
		      a%lnods(ispos+5) = TopLeftGhost + irow
		      a%lnods(ispos+6) = TopGhost0 + irow
		      a%lnods(ispos+7) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+8) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
	      enddo

      endif   
         
      if (TopRightGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + BlockTop + BlockRight 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(TopRightGhost + ipoin) = baseGhost0 + ipoin
		      ProcessorList(TopRightGhost + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX + ManufacturedNprocX
	      enddo
         
	      do irow = 1,a%ManufacturedNpoinX-1
		      ielem = TopRightElem0 + irow
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = TopGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 
		      a%lnods(ispos+2) = TopRightGhost + irow + 1
		      a%lnods(ispos+3) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1
		      a%lnods(ispos+4) = floortop*(a%ManufacturedNpoinX-1) + colright*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+5) = TopGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+6) = TopRightGhost + irow
		      a%lnods(ispos+7) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+8) = floortop*(a%ManufacturedNpoinX-1) + colright*(a%ManufacturedNpoinX-1) + irow
	      enddo

      endif   
     
      if (TopUpGhost /= -1) then
	      baseGhost0 = a%poin0 + 1 + BlockTop + BlockUp + rowdown*(a%ManufacturedNpoinX-1) 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(TopUpGhost + ipoin) = baseGhost0 + colright*(ipoin-1)
		      ProcessorList(TopUpGhost + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX -1
	      enddo

	      do jcol = 1,a%ManufacturedNpoinX-1
		      ielem = TopUpElem0 + jcol
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+2) = TopGhost0 + jcol*a%ManufacturedNpoinX + 1 
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + colright*jcol + 1 
		      a%lnods(ispos+4) = floortop*(a%ManufacturedNpoinX-1) + colright*(jcol-1) + 1
		      a%lnods(ispos+5) = TopUpGhost + jcol
		      a%lnods(ispos+6) = TopUpGhost + jcol + 1
		      a%lnods(ispos+7) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1
		      a%lnods(ispos+8) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol 

	      enddo
      endif   
            
      if (TopDownGhost /= -1) then
	      baseGhost0 = a%poin0 + 1 + BlockTop + BlockDown  
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(TopDownGhost + ipoin) = baseGhost0 + colright*(ipoin-1)
		      ProcessorList(TopDownGhost + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX +1
	      enddo

	      do jcol = 1,a%ManufacturedNpoinX-1
		      ielem = TopDownElem0 + jcol
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = TopDownGhost + jcol
		      a%lnods(ispos+2) = TopDownGhost + jcol + 1
		      a%lnods(ispos+3) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1
		      a%lnods(ispos+4) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol 
		      a%lnods(ispos+5) = TopGhost0 + jcol*a%ManufacturedNpoinX
		      a%lnods(ispos+6) = TopGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+7) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+8) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*(jcol-1) + 1

	      enddo
      endif   

      if (BottomLeftGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + BlockBottom + BlockLeft + (a%ManufacturedNpoinX-1)*floortop + colright*(a%ManufacturedNpoinX-1)
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(BottomLeftGhost + ipoin) = baseGhost0 + ipoin
		      ProcessorList(BottomLeftGhost + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX - ManufacturedNprocX
	      enddo
         
	      do irow = 1,a%ManufacturedNpoinX-1
		      ielem = BottomLeftElem0 + irow
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + irow + 1
		      a%lnods(ispos+2) = irow + rowdown
		      a%lnods(ispos+3) = BottomGhost0 + irow + 1 
		      a%lnods(ispos+4) = BottomLeftGhost + irow + 1
		      a%lnods(ispos+5) = LeftGhost0 + irow
		      a%lnods(ispos+6) = irow
		      a%lnods(ispos+7) = BottomGhost0 + irow
		      a%lnods(ispos+8) = BottomLeftGhost + irow
	      enddo

      endif   

      if (BottomRightGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + BlockBottom + BlockRight + (a%ManufacturedNpoinX-1)*floortop
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(BottomRightGhost + ipoin) = baseGhost0 + ipoin
		      ProcessorList(BottomRightGhost + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX + ManufacturedNprocX
	      enddo
         
	      do irow = 1,a%ManufacturedNpoinX-1
		      ielem = BottomRightElem0 + irow
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+2) = RightGhost0 + irow + 1
		      a%lnods(ispos+3) = BottomRightGhost + irow + 1
		      a%lnods(ispos+4) = BottomGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 
		      a%lnods(ispos+5) = colright*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+6) = RightGhost0 + irow
		      a%lnods(ispos+7) = BottomRightGhost + irow
		      a%lnods(ispos+8) = BottomGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
	      enddo

      endif   
      
      if (BottomUpGhost /= -1) then
	      baseGhost0 = a%poin0 + 1 + BlockBottom + BlockUp + (a%ManufacturedNpoinX-1)*floortop + rowdown*(a%ManufacturedNpoinX-1) 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(BottomUpGhost + ipoin) = baseGhost0 + colright*(ipoin-1)
		      ProcessorList(BottomUpGhost + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX -1
	      enddo

	      do jcol = 1,a%ManufacturedNpoinX-1
		      ielem = BottomUpElem0 + jcol
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = colright*(jcol-1) + 1
		      a%lnods(ispos+2) = colright*jcol + 1 
		      a%lnods(ispos+3) = BottomGhost0 + jcol*a%ManufacturedNpoinX + 1 
		      a%lnods(ispos+4) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+5) = UpGhost0 + jcol 
		      a%lnods(ispos+6) = UpGhost0 + jcol + 1
		      a%lnods(ispos+7) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+8) = BottomUpGhost + jcol

	      enddo
      endif   
            
      if (BottomDownGhost /= -1) then
	      baseGhost0 = a%poin0 + 1 + BlockBottom + BlockDown + (a%ManufacturedNpoinX-1)*floortop
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(BottomDownGhost + ipoin) = baseGhost0 + colright*(ipoin-1)
		      ProcessorList(BottomDownGhost + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX +1
	      enddo

	      do jcol = 1,a%ManufacturedNpoinX-1
		      ielem = BottomDownElem0 + jcol
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownGhost0 + jcol
		      a%lnods(ispos+2) = DownGhost0 + jcol + 1
		      a%lnods(ispos+3) = BottomDownGhost + jcol + 1
		      a%lnods(ispos+4) = BottomDownGhost + jcol
		      a%lnods(ispos+5) = a%ManufacturedNpoinX + colright*(jcol-1) 
		      a%lnods(ispos+6) = a%ManufacturedNpoinX + colright*jcol
		      a%lnods(ispos+7) = BottomGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+8) = BottomGhost0 + jcol*a%ManufacturedNpoinX

	      enddo
      endif   

      if (TopUpLeftGhost /= -1) then
	      LocalToGlobal(TopUpLeftGhost) = a%poin0 + BlockTop + BlockUp + BlockLeft + a%ManufacturedNpoinX*a%ManufacturedNpoinX 
	      ProcessorList(TopUpLeftGhost) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX - 1 - ManufacturedNprocX

	      ielem = TopUpLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = TopLeftGhost + 1
	      a%lnods(ispos+2) = TopGhost0 + 1
	      a%lnods(ispos+3) = (a%ManufacturedNpoinX-1)*floortop + 1
	      a%lnods(ispos+4) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+5) = TopUpLeftGhost
	      a%lnods(ispos+6) = TopUpGhost + 1
	      a%lnods(ispos+7) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+8) = UpLeftGhost + a%ManufacturedNpoinX
      endif   
      
      if (TopDownLeftGhost /= -1) then
	      LocalToGlobal(TopDownLeftGhost) = a%poin0 + 1 + BlockTop + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1) 
	      ProcessorList(TopDownLeftGhost) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX + 1 - ManufacturedNprocX

	      ielem = TopDownLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = TopDownLeftGhost  
	      a%lnods(ispos+2) = TopDownGhost + 1
	      a%lnods(ispos+3) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+4) = DownLeftGhost + a%ManufacturedNpoinX 
	      a%lnods(ispos+5) = TopLeftGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+6) = TopGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+7) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+8) = LeftGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         
      endif

      if (TopUpRightGhost /= -1) then
	      LocalToGlobal(TopUpRightGhost) = a%poin0 + BlockTop + BlockUp + BlockRight + a%ManufacturedNpoinX 
	      ProcessorList(TopUpRightGhost) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX  - 1 + ManufacturedNprocX

	      ielem = TopUpRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = TopGhost0 + 1 + a%ManufacturedNPoinX*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+2) = TopRightGhost + 1
	      a%lnods(ispos+3) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+4) = (a%ManufacturedNpoinX-1)*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+5) = TopUpGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+6) = TopUpRightGhost 
	      a%lnods(ispos+7) = UpRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+8) = UpGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX

      endif   
         
      if (TopDownRightGhost /= -1) then
	      LocalToGlobal(TopDownRightGhost) =  a%poin0 + 1 + BlockTop + BlockDown + BlockRight 
	      ProcessorList(TopDownRightGhost) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX + 1 + ManufacturedNprocX

	      ielem = TopDownRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = TopDownGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = TopDownRightGhost 
	      a%lnods(ispos+3) = DownRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = DownGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+5) = TopGhost0 + a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+6) = TopRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+7) = RightGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+8) = a%ManufacturedNPoinX*a%ManufacturedNPoinX*a%ManufacturedNPoinX

      endif   

      if (BottomUpLeftGhost /= -1) then
	      LocalToGlobal(BottomUpLeftGhost) = a%poin0 + BlockBottom + BlockUp + BlockLeft + a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%ManufacturedNpoinX 
	      ProcessorList(BottomUpLeftGhost) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX - 1 - ManufacturedNprocX

	      ielem = BottomUpLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = LeftGhost0 + 1
	      a%lnods(ispos+2) = 1
	      a%lnods(ispos+3) = BottomGhost0 + 1
	      a%lnods(ispos+4) = BottomLeftGhost + 1
	      a%lnods(ispos+5) = UpLeftGhost + 1 
	      a%lnods(ispos+6) = UpGhost0 + 1
	      a%lnods(ispos+7) = BottomUpGhost + 1
	      a%lnods(ispos+8) = BottomUpLeftGhost 
      endif   
      
      if (BottomDownLeftGhost /= -1) then
	      LocalToGlobal(BottomDownLeftGhost) = a%poin0 + 1 + BlockBottom + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1) + (a%ManufacturedNpoinX-1)*floortop 
	      ProcessorList(BottomDownLeftGhost) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX + 1 - ManufacturedNprocX

	      ielem = BottomDownLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = DownLeftGhost + 1 
	      a%lnods(ispos+2) = DownGhost0 + 1
	      a%lnods(ispos+3) = BottomDownGhost + 1
	      a%lnods(ispos+4) = BottomDownLeftGhost 
	      a%lnods(ispos+5) = LeftGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+6) = a%ManufacturedNpoinX
	      a%lnods(ispos+7) = BottomGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+8) = BottomLeftGhost + a%ManufacturedNpoinX

      endif
      
      if (BottomUpRightGhost /= -1) then
	      LocalToGlobal(BottomUpRightGhost) = a%poin0 + 1 + BlockBottom + BlockUp + BlockRight + rowdown*(a%ManufacturedNpoinX-1) + (a%ManufacturedNpoinX-1)*floortop
	      ProcessorList(BottomUpRightGhost) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX - 1 + ManufacturedNprocX

	      ielem = BottomUpRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = 1 + colright*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+2) = RightGhost0 + 1
	      a%lnods(ispos+3) = BottomRightGhost + 1
	      a%lnods(ispos+4) = BottomGhost0 + 1 + a%ManufacturedNPoinX*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+5) = UpGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+6) = UpRightGhost + 1
	      a%lnods(ispos+7) = BottomUpRightGhost
	      a%lnods(ispos+8) = BottomUpGhost + a%ManufacturedNpoinX
         
      endif   
         
      if (BottomDownRightGhost /= -1) then
	      LocalToGlobal(BottomDownRightGhost) = a%poin0 + 1 + BlockBottom + BlockDown + BlockRight + (a%ManufacturedNpoinX-1)*floortop
	      ProcessorList(BottomDownRightGhost) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX + 1 + ManufacturedNprocX

	      ielem = BottomDownRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = DownGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = DownRightGhost + 1
	      a%lnods(ispos+3) = BottomDownRightGhost
	      a%lnods(ispos+4) = BottomDownGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+5) = a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+6) = RightGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+7) = BottomRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+8) = BottomGhost0 + a%ManufacturedNPoinX*a%ManufacturedNPoinX

      endif   

      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      call a%Memor%dealloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%dealloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !External normal, we simply allocate it
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ManufacturedMesh')
      
   endif
   
   if (a%ndime == 3 .and. a%ManufacturedType == 'LINTE') then
      
      if ( nint((real(a%MPIsize))**(1.0/3.0))**3 /= a%MPIsize) then
         kfl_isPowerOf = .false.
      else
         kfl_isPowerOf = .true.
      endif
      
      if (kfl_isPowerOf .eqv. .false.) then
         call runend('3d ManufacturedMesh, MPIsize the cube of a natural number')
      endif
      
      ManufacturedNprocX = (real(a%MPIsize,rp))**(1.0/3.0)
      
      ManufacturedProcPosition(3) = a%MPIrank/(ManufacturedNprocX*ManufacturedNprocX)
      ManufacturedProcPosition(2) = mod(a%MPIrank,ManufacturedNprocX)
      ManufacturedProcPosition(1) = (a%MPIrank-ManufacturedProcPosition(3)*ManufacturedNprocX*ManufacturedNprocX)/ManufacturedNprocX

      
      a%gnpoin = a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%MPIsize
      a%npoinLocal = a%ManufacturedNpoinX*a%ManufacturedNpoinX*a%ManufacturedNpoinX
      
      a%nelem = 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      a%npoinGhost = 0

      BottomGhost0 = -1
      TopGhost0 = -1
      UpGhost0 = -1
      DownGhost0 = -1
      LeftGhost0 = -1
      RightGhost0 = -1

      if (ManufacturedProcPosition(3) /= 0) then
         BottomGhost0 = a%npoinLocal + a%npoinGhost
         BottomElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1) then
         TopGhost0 = a%npoinLocal + a%npoinGhost
         TopElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(2) /= 0) then
         UpGhost0 = a%npoinLocal + a%npoinGhost
         UpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownGhost0 = a%npoinLocal + a%npoinGhost
         DownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= 0) then
         LeftGhost0 = a%npoinLocal + a%npoinGhost
         LeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         RightGhost0 = a%npoinLocal + a%npoinGhost
         RightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX*a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)*(a%ManufacturedNpoinX-1)
      endif

      UpRightGhost = -1
      DownRightGhost = -1
      UpLeftGhost = -1
      DownLeftGhost = -1
      TopLeftGhost = -1
      TopRightGhost = -1
      TopUpGhost = -1
      TopDownGhost = -1
      BottomLeftGhost = -1
      BottomRightGhost = -1
      BottomUpGhost = -1
      BottomDownGhost = -1

      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
         UpLeftGhost = a%npoinLocal + a%npoinGhost
         UpLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost
         a%nelem = a%nelem + 3*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownLeftGhost = a%npoinLocal + a%npoinGhost
         DownLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
         UpRightGhost = a%npoinLocal + a%npoinGhost
         UpRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         DownRightGhost = a%npoinLocal + a%npoinGhost
         DownRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost
         a%nelem = a%nelem + 3*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(1) /= 0) then
         TopLeftGhost = a%npoinLocal + a%npoinGhost
         TopLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         TopRightGhost = a%npoinLocal + a%npoinGhost
         TopRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost
         a%nelem = a%nelem + 3*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
         TopUpGhost = a%npoinLocal + a%npoinGhost
         TopUpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost
         a%nelem = a%nelem + 3*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         TopDownGhost = a%npoinLocal + a%npoinGhost
         TopDownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(1) /= 0) then
         BottomLeftGhost = a%npoinLocal + a%npoinGhost
         BottomLeftElem0 = a%nelem
         a%npoinGhost = a%npoinGhost 
         a%nelem = a%nelem + 3*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(1) /= ManufacturedNprocX-1) then
         BottomRightGhost = a%npoinLocal + a%npoinGhost
         BottomRightElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
         BottomUpGhost = a%npoinLocal + a%npoinGhost
         BottomUpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + a%ManufacturedNpoinX
         a%nelem = a%nelem + 6*(a%ManufacturedNpoinX-1)
      endif
      if (ManufacturedProcPosition(3) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
         BottomDownGhost = a%npoinLocal + a%npoinGhost
         BottomDownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost 
         a%nelem = a%nelem + 3*(a%ManufacturedNpoinX-1)
      endif

      TopUpLeftGhost = -1
      TopUpRightGhost = -1
      TopDownLeftGhost = -1
      TopDownRightGhost = -1
      BottomUpLeftGhost = -1
      BottomUpRightGhost = -1
      BottomDownLeftGhost = -1
      BottomDownRightGhost = -1

      if (ManufacturedProcPosition(3) /= ManufacturedNprocX-1) then
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost
		      TopUpLeftGhost = a%npoinLocal + a%npoinGhost
		      TopUpLeftElem0 = a%nelem
		      a%nelem = a%nelem + 2
	      endif
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      TopDownLeftGhost = a%npoinLocal + a%npoinGhost
		      TopDownLeftElem0 = a%nelem
		      a%nelem = a%nelem + 6
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost
		      TopUpRightGhost = a%npoinLocal + a%npoinGhost
		      TopUpRightElem0 = a%nelem
		      a%nelem = a%nelem + 2
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost 
		      TopDownRightGhost = a%npoinLocal + a%npoinGhost
		      TopDownRightElem0 = a%nelem
		      a%nelem = a%nelem + 2
	      endif
      endif
      if (ManufacturedProcPosition(3) /= 0) then
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost
		      BottomUpLeftGhost = a%npoinLocal + a%npoinGhost
		      BottomUpLeftElem0 = a%nelem
		      a%nelem = a%nelem + 2
	      endif
	      if (ManufacturedProcPosition(1) /= 0 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost
		      BottomDownLeftGhost = a%npoinLocal + a%npoinGhost
		      BottomDownLeftElem0 = a%nelem
		      a%nelem = a%nelem + 2
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= 0) then
		      a%npoinGhost = a%npoinGhost + 1_ip
		      BottomUpRightGhost = a%npoinLocal + a%npoinGhost
		      BottomUpRightElem0 = a%nelem
		      a%nelem = a%nelem + 6
	      endif
	      if (ManufacturedProcPosition(1) /= ManufacturedNprocX-1 .and. ManufacturedProcPosition(2) /= ManufacturedNprocX-1) then
		      a%npoinGhost = a%npoinGhost 
		      BottomDownRightGhost = a%npoinLocal + a%npoinGhost
		      BottomDownRightElem0 = a%nelem
		      a%nelem = a%nelem + 2
	      endif
      endif


      a%npoin = a%npoinLocal+a%npoinGhost
      a%gnelem = 6*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*a%MPIsize
      
      a%poin0 = a%MPIrank*a%npoinLocal
      a%elem0 = a%MPIrank*6*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)*(a%ManufacturedNpoinX+1)
      
      !Fill local nodes and elements
      call a%Memor%alloc(1_ip,a%pnods,'pnods','ManufacturedMesh')
      call a%Memor%alloc(a%nelem*4_ip,a%lnods,'lnods','ManufacturedMesh')
      call a%Memor%alloc(3_ip,a%npoin,a%coord,'coord','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','ManufacturedMesh')
      a%geoblock=0
      
      call a%Memor%alloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !Build local elements, nodes and coord
      TotalNpoinX = ManufacturedNprocX*a%ManufacturedNpoinX
      ElemCoordSpanX = 1_rp/real(TotalNpoinX-1,rp)
      
      basecoordx = ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(1))
      basecoordy = 1- ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(2))
      basecoordz = ElemCoordSpanX*(a%ManufacturedNpoinX*ManufacturedProcPosition(3))
      
      !Local Coordinates
      do kflr = 1,a%ManufacturedNpoinX
	      coordz = basecoordz + (kflr-1)*ElemCoordSpanX
	      basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      do jcol = 1,a%ManufacturedNpoinX
		      coordx = basecoordx + (jcol-1)*ElemCoordSpanX
		      basepointcol = (jcol-1)*a%ManufacturedNpoinX
		      do irow = 1,a%ManufacturedNpoinX
			      ipoin = basepointflr+basepointcol+irow   
			      coordy = basecoordy - (irow-1)*ElemCoordSpanX
			      a%coord(1,ipoin) = coordx
			      a%coord(2,ipoin) = coordy
			      a%coord(3,ipoin) = coordz
		      enddo
	      enddo
      enddo

      BlockDown = a%npoinLocal
      BlockUp   = -BlockDown
      BlockRight = a%npoinLocal*ManufacturedNprocX
      BlockLeft = -BlockRight
      BlockTop = a%npoinLocal*ManufacturedNprocX*ManufacturedNprocX
      BlockBottom = -BlockTop

      rowdown = 1
      rowup = -rowdown
      colright = a%ManufacturedNpoinX
      colleft = -colright
      floortop = a%ManufacturedNpoinX*a%ManufacturedNpoinX
      floorbottom = -floortop
      
      !Local Elements
      a%pnods(1) = 4
      ielem = 0
      do kflr = 1,a%ManufacturedNpoinX-1
	      basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      do jcol = 1,a%ManufacturedNpoinX-1
		      basepointcol = (jcol-1)*a%ManufacturedNpoinX
		      do irow = 1,a%ManufacturedNpoinX-1
			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + basepointcol + irow
			      a%lnods(ispos+2) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+3) = basepointflr + floortop + basepointcol  + irow 
			      a%lnods(ispos+4) = basepointflr + floortop + basepointcol  + irow + rowdown

			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+2) = basepointflr + floortop + basepointcol  + irow + colright
			      a%lnods(ispos+3) = basepointflr + floortop + basepointcol  + irow 
			      a%lnods(ispos+4) = basepointflr + floortop + basepointcol  + irow + rowdown

			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + basepointcol + irow 
			      a%lnods(ispos+2) = basepointflr + basepointcol + irow + rowdown
			      a%lnods(ispos+3) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+4) = basepointflr + floortop + basepointcol  + irow + rowdown

			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + floortop + basepointcol  + irow + colright + rowdown
			      a%lnods(ispos+2) = basepointflr + floortop + basepointcol  + irow + colright
			      a%lnods(ispos+3) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+4) = basepointflr + floortop + basepointcol  + irow + rowdown

			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + basepointcol + irow + rowdown
			      a%lnods(ispos+2) = basepointflr + basepointcol + irow + colright + rowdown
			      a%lnods(ispos+3) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+4) = basepointflr + floortop + basepointcol  + irow + rowdown

			      ielem = ielem+1
			      ispos = (ielem-1)*a%pnods(1)
			      a%lnods(ispos+1) = basepointflr + basepointcol + irow + colright 
			      a%lnods(ispos+2) = basepointflr + basepointcol + irow + colright + rowdown
			      a%lnods(ispos+3) = basepointflr + floortop + basepointcol  + irow + colright + rowdown
			      a%lnods(ispos+4) = basepointflr + floortop + basepointcol  + irow + rowdown
		      enddo
	      enddo
      enddo
      
      !LocalToGlobal for local points
      do ipoin = 1,a%npoinLocal
         LocalToGlobal(ipoin) = a%poin0+ipoin
         ProcessorList(ipoin) = a%MPIrank
      enddo


      !Ghost points and elements
      if (BottomGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockBottom + floortop*(a%ManufacturedNpoinX-1)
	       do ipoin = 1,a%ManufacturedNpoinX*a%ManufacturedNpoinX
		        LocalToGlobal(BottomGhost0 + ipoin) = baseGhost0 + ipoin
		        ProcessorList(BottomGhost0 + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX
	       enddo
      	 do jcol = 1,a%ManufacturedNpoinX-1
      		 basepointcol = (jcol-1)*a%ManufacturedNpoinX
      		 do irow = 1,a%ManufacturedNpoinX-1
      
      			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 1
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+2) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
      			 a%lnods(ispos+3) = basepointcol + irow
      			 a%lnods(ispos+4) = basepointcol + irow + rowdown

      			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 2
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
      			 a%lnods(ispos+2) = basepointcol + colright + irow
      			 a%lnods(ispos+3) = basepointcol + irow
      			 a%lnods(ispos+4) = basepointcol + irow + rowdown
      
      			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 3
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+2) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1  
      			 a%lnods(ispos+3) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
      			 a%lnods(ispos+4) = basepointcol + irow + rowdown
      
      			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 4
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointcol + colright + irow + rowdown 
      			 a%lnods(ispos+2) = basepointcol + colright + irow
      			 a%lnods(ispos+3) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
      			 a%lnods(ispos+4) = basepointcol + irow + rowdown

      			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 5
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1  
      			 a%lnods(ispos+2) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+3) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
      			 a%lnods(ispos+4) = basepointcol + irow + rowdown

      			 ielem = BottomElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 6
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow  
      			 a%lnods(ispos+2) = BottomGhost0 + jcol*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+3) = basepointcol + colright + irow + rowdown 
      			 a%lnods(ispos+4) = basepointcol + irow + rowdown

      		 enddo
      	 enddo
      endif
      if (TopGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockTop 
	       do ipoin = 1,a%ManufacturedNpoinX*a%ManufacturedNpoinX
		       LocalToGlobal(TopGhost0 + ipoin) = baseGhost0 + ipoin
		       ProcessorList(TopGhost0 + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX
	       enddo
         basepointflr = floortop*(a%ManufacturedNpoinX-1)
      	 do jcol = 1,a%ManufacturedNpoinX-1
      		 basepointcol = (jcol-1)*a%ManufacturedNpoinX
      		 do irow = 1,a%ManufacturedNpoinX-1
      
      			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 1
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + irow
      			 a%lnods(ispos+2) = basepointflr + basepointcol + colright + irow
      			 a%lnods(ispos+3) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1

      			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 2 
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + colright + irow
      			 a%lnods(ispos+2) = TopGhost0 + jcol*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+3) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1

      			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 3
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + irow
      			 a%lnods(ispos+2) = basepointflr + basepointcol + irow + rowdown
      			 a%lnods(ispos+3) = basepointflr + basepointcol + colright + irow
      			 a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1

      			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 4
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = TopGhost0 + jcol*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+2) = TopGhost0 + jcol*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+3) = basepointflr + basepointcol + colright + irow
      			 a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1

      			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 5
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + irow + rowdown
      			 a%lnods(ispos+2) = basepointflr + basepointcol + colright + irow + rowdown 
      			 a%lnods(ispos+3) = basepointflr + basepointcol + colright + irow
      			 a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1

      			 ielem = TopElem0 + (jcol-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 6
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + colright + irow
      			 a%lnods(ispos+2) = basepointflr + basepointcol + colright + irow + rowdown 
      			 a%lnods(ispos+3) = TopGhost0 + jcol*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + irow + 1
      
      		 enddo
      	 enddo
      endif
      if (UpGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + 1 + BlockUp + rowdown*(a%ManufacturedNpoinX-1) 
	       do jpoin = 1,a%ManufacturedNpoinX
		        do ipoin = 1,a%ManufacturedNpoinX
			         LocalToGlobal(UpGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + colright*(ipoin-1)
			         ProcessorList(UpGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank-1
		        enddo
         enddo
         do kflr = 1,a%ManufacturedNpoinX-1
           basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
           do jcol = 1,a%ManufacturedNpoinX-1
          	 basepointcol = (jcol-1)*a%ManufacturedNpoinX
         
          	 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 1
          	 ispos = (ielem-1)*a%pnods(1)
          	 a%lnods(ispos+1) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol
          	 a%lnods(ispos+2) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+3) = UpGhost0 + kflr*a%ManufacturedNpoinX + jcol
          	 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + 1
         
          	 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 2
          	 ispos = (ielem-1)*a%pnods(1)
          	 a%lnods(ispos+1) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+2) = UpGhost0 + kflr*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+3) = UpGhost0 + kflr*a%ManufacturedNpoinX + jcol
          	 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + 1

          	 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 3
          	 ispos = (ielem-1)*a%pnods(1)
          	 a%lnods(ispos+1) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol
          	 a%lnods(ispos+2) = basepointflr + basepointcol + 1
          	 a%lnods(ispos+3) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + 1

          	 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 4
          	 ispos = (ielem-1)*a%pnods(1)
          	 a%lnods(ispos+1) = basepointflr + floortop + basepointcol + 1 + colright 
          	 a%lnods(ispos+2) = UpGhost0 + kflr*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+3) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + 1

          	 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 5
          	 ispos = (ielem-1)*a%pnods(1)
          	 a%lnods(ispos+1) = basepointflr + basepointcol + 1
          	 a%lnods(ispos+2) = basepointflr + basepointcol + 1 + colright 
          	 a%lnods(ispos+3) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + 1

          	 ielem = UpElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 6
          	 ispos = (ielem-1)*a%pnods(1)
          	 a%lnods(ispos+1) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
          	 a%lnods(ispos+2) = basepointflr + basepointcol + 1 + colright 
          	 a%lnods(ispos+3) = basepointflr + floortop + basepointcol + 1 + colright 
          	 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + 1
           enddo
         enddo
      endif
      if (DownGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + 1 + BlockDown
      	 do jpoin = 1,a%ManufacturedNpoinX
      		 do ipoin = 1,a%ManufacturedNpoinX
      			 LocalToGlobal(DownGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + colright*(ipoin-1)
      			 ProcessorList(DownGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank+1
      		 enddo
      	 enddo
      	 do kflr = 1,a%ManufacturedNpoinX-1
      		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
      		 do jcol = 1,a%ManufacturedNpoinX-1
      			 basepointcol = (jcol-1)*a%ManufacturedNpoinX + a%ManufacturedNpoinX - 1
      
      			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 1 
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + 1 
      			 a%lnods(ispos+2) = basepointflr + basepointcol + 1 + colright 
      			 a%lnods(ispos+3) = basepointflr + floortop + basepointcol + 1 
      			 a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
     
      			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 2
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + 1 + colright 
      			 a%lnods(ispos+2) = basepointflr + floortop + basepointcol + 1 + colright 
      			 a%lnods(ispos+3) = basepointflr + floortop + basepointcol + 1 
      			 a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
      
      			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 3
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + 1 
      			 a%lnods(ispos+2) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol
      			 a%lnods(ispos+3) = basepointflr + basepointcol + 1 + colright 
      			 a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
      
      			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 4
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol + 1
      			 a%lnods(ispos+2) = basepointflr + floortop + basepointcol + 1 + colright 
      			 a%lnods(ispos+3) = basepointflr + basepointcol + 1 + colright 
      			 a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
      
      			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 5
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol
      			 a%lnods(ispos+2) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
      			 a%lnods(ispos+3) = basepointflr + basepointcol + 1 + colright 
      			 a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
      
      			 ielem = DownElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (jcol-1)*6 + 6
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + 1 + colright 
      			 a%lnods(ispos+2) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + jcol + 1
      			 a%lnods(ispos+3) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol + 1
      			 a%lnods(ispos+4) = DownGhost0 + kflr*a%ManufacturedNpoinX + jcol
      		 enddo
      	 enddo
      endif
      if (LeftGhost0 /= -1) then
         !GhostPoints Global numbering
      	 baseGhost0 = a%poin0 + BlockLeft + colright*(a%ManufacturedNpoinX-1)
      	 do jpoin = 1,a%ManufacturedNpoinX
      		 do ipoin = 1,a%ManufacturedNpoinX
      			 LocalToGlobal(LeftGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + ipoin
      			 ProcessorList(LeftGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank - ManufacturedNprocX
      		 enddo
      	 enddo
               basepointcol = 0
      	 do kflr = 1,a%ManufacturedNpoinX-1
      		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
      		 do irow = 1,a%ManufacturedNpoinX-1

      			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 1 
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+2) = basepointflr + irow
      			 a%lnods(ispos+3) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1

      			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 2
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + irow
      			 a%lnods(ispos+2) = basepointflr + floortop + irow
      			 a%lnods(ispos+3) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1

      			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 3
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+2) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+3) = basepointflr + irow
      			 a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1

      			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 4
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + floortop + irow + rowdown 
      			 a%lnods(ispos+2) = basepointflr + floortop + irow
      			 a%lnods(ispos+3) = basepointflr + irow
      			 a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1

      			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 5
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+2) = basepointflr + irow + rowdown 
      			 a%lnods(ispos+3) = basepointflr + irow
      			 a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1

      			 ielem = LeftElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 6
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + irow
      			 a%lnods(ispos+2) = basepointflr + irow + rowdown 
      			 a%lnods(ispos+3) = basepointflr + floortop + irow + rowdown 
      			 a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + irow + 1
      		 enddo
      	 enddo
      endif
      if (RightGhost0 /= -1) then
         !GhostPoints Global numbering
         baseGhost0 = a%poin0 + BlockRight
      	 do jpoin = 1,a%ManufacturedNpoinX
      		 do ipoin = 1,a%ManufacturedNpoinX
      			 LocalToGlobal(RightGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = baseGhost0 + floortop*(jpoin-1) + ipoin
      			 ProcessorList(RightGhost0 + (jpoin-1)*a%ManufacturedNpoinX + ipoin) = a%MPIrank + ManufacturedNprocX
      		 enddo
      	 enddo
               basepointcol = colright*(a%ManufacturedNpoinX-1)
      	 do kflr = 1,a%ManufacturedNpoinX-1
      		 basepointflr = (kflr-1)*a%ManufacturedNpoinX*a%ManufacturedNpoinX
      		 do irow = 1,a%ManufacturedNpoinX-1

      			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 1 
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + irow 
      			 a%lnods(ispos+2) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+3) = basepointflr + floortop + basepointcol + irow 
      			 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + irow + rowdown  

      			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 2
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+2) = RightGhost0 + kflr*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+3) = basepointflr + floortop + basepointcol + irow 
      			 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + irow + rowdown  

      			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 3
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + irow 
      			 a%lnods(ispos+2) = basepointflr + basepointcol + irow + rowdown 
      			 a%lnods(ispos+3) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + irow + rowdown  

      			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 4
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = RightGhost0 + kflr*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+2) = RightGhost0 + kflr*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+3) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + irow + rowdown  

      			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 5
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = basepointflr + basepointcol + irow + rowdown  
      			 a%lnods(ispos+2) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+3) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + irow + rowdown  

      			 ielem = RightElem0 + (kflr-1)*(a%ManufacturedNpoinX-1)*6 + (irow-1)*6 + 6
      			 ispos = (ielem-1)*a%pnods(1)
      			 a%lnods(ispos+1) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow
      			 a%lnods(ispos+2) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+3) = RightGhost0 + kflr*a%ManufacturedNpoinX + irow + 1
      			 a%lnods(ispos+4) = basepointflr + floortop + basepointcol + irow + rowdown    
      		 enddo
      	 enddo
      endif
      if (UpLeftGhost /= -1) then
	      do kflr = 1,a%ManufacturedNpoinX-1

		      ielem = UpLeftElem0 + (kflr-1)*3 + 1 
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = kflr*floortop + 1
		      a%lnods(ispos+2) = UpGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + 1

		      ielem = UpLeftElem0 + (kflr-1)*3 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+2) = (kflr-1)*floortop + 1
		      a%lnods(ispos+3) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + 1

		      ielem = UpLeftElem0 + (kflr-1)*3 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = UpGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+2) = (kflr-1)*floortop + 1
		      a%lnods(ispos+3) = kflr*floortop + 1
		      a%lnods(ispos+4) = LeftGhost0 + kflr*a%ManufacturedNpoinX + 1
	      enddo
      endif   
      if (DownLeftGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + 1 + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1)  
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(DownLeftGhost + ipoin) = baseGhost0 + (ipoin-1)*floortop 
		      ProcessorList(DownLeftGhost + ipoin) = a%MPIrank + 1 - ManufacturedNprocX
	      enddo

	      do kflr = 1,a%ManufacturedNpoinX-1

		      ielem = DownLeftElem0 + (kflr-1)*6 + 1
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+2) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+3) = LeftGhost0 + (kflr+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownLeftGhost + kflr + 1 

		      ielem = DownLeftElem0 + (kflr-1)*6 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+2) = kflr*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+3) = LeftGhost0 + (kflr+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownLeftGhost + kflr + 1 

		      ielem = DownLeftElem0 + (kflr-1)*6 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+2) = DownLeftGhost + kflr 
		      a%lnods(ispos+3) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownLeftGhost + kflr + 1 

		      ielem = DownLeftElem0 + (kflr-1)*6 + 4
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+2) = kflr*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+3) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownLeftGhost + kflr + 1 

		      ielem = DownLeftElem0 + (kflr-1)*6 + 5
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownLeftGhost + kflr 
		      a%lnods(ispos+2) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownLeftGhost + kflr + 1 

		      ielem = DownLeftElem0 + (kflr-1)*6 + 6
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = (kflr-1)*floortop + a%ManufacturedNpoinX
		      a%lnods(ispos+2) = DownGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = DownGhost0 + kflr*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+4) = DownLeftGhost + kflr + 1 
	      enddo
         
      endif

      if (UpRightGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + 1 + BlockUp + BlockRight + rowdown*(a%ManufacturedNpoinX-1) 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(UpRightGhost + ipoin) = baseGhost0 + (ipoin-1)*floortop
		      ProcessorList(UpRightGhost + ipoin) = a%MPIrank - 1 + ManufacturedNprocX
	      enddo

	      do kflr = 1,a%ManufacturedNpoinX-1

		      ielem = UpRightElem0 + (kflr-1)*6 + 1 
		      ispos = (ielem-1)*a%pnods(1)
          a%lnods(ispos+1) = UpGhost0 + kflr*a%ManufacturedNpoinX
          a%lnods(ispos+2) = UpRightGhost + kflr
          a%lnods(ispos+3) = UpGhost0 + (kflr+1)*a%ManufacturedNpoinX
          a%lnods(ispos+4) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 

		      ielem = UpRightElem0 + (kflr-1)*6 + 2 
		      ispos = (ielem-1)*a%pnods(1)
          a%lnods(ispos+1) = UpRightGhost + kflr
          a%lnods(ispos+2) = UpRightGhost + kflr + 1
          a%lnods(ispos+3) = UpGhost0 + (kflr+1)*a%ManufacturedNpoinX
          a%lnods(ispos+4) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 

		      ielem = UpRightElem0 + (kflr-1)*6 + 3
		      ispos = (ielem-1)*a%pnods(1)
          a%lnods(ispos+1) = UpGhost0 + kflr*a%ManufacturedNpoinX
          a%lnods(ispos+2) = (kflr-1)*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
          a%lnods(ispos+3) = UpRightGhost + kflr
          a%lnods(ispos+4) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 

		      ielem = UpRightElem0 + (kflr-1)*6 + 4
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = RightGhost0 + kflr*a%ManufacturedNpoinX + 1
          a%lnods(ispos+2) = UpRightGhost + kflr + 1
          a%lnods(ispos+3) = UpRightGhost + kflr
          a%lnods(ispos+4) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 

		      ielem = UpRightElem0 + (kflr-1)*6 + 5
		      ispos = (ielem-1)*a%pnods(1)
          a%lnods(ispos+1) = (kflr-1)*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
		      a%lnods(ispos+2) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
          a%lnods(ispos+3) = UpRightGhost + kflr
          a%lnods(ispos+4) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 

		      ielem = UpRightElem0 + (kflr-1)*6 + 6
		      ispos = (ielem-1)*a%pnods(1)
          a%lnods(ispos+1) = UpRightGhost + kflr
		      a%lnods(ispos+2) = RightGhost0 + (kflr-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = RightGhost0 + kflr*a%ManufacturedNpoinX + 1
          a%lnods(ispos+4) = kflr*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
	      enddo
         
      endif   
         
      if (DownRightGhost /= -1) then
         
	      do kflr = 1,a%ManufacturedNpoinX-1

		      ielem = DownRightElem0 + (kflr-1)*3 + 1
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = kflr*a%ManufacturedNPoinX*a%ManufacturedNPoinX
		      a%lnods(ispos+2) = RightGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+3) = (kflr+1)*a%ManufacturedNPoinX*a%ManufacturedNPoinX
		      a%lnods(ispos+4) = DownGhost0 + (kflr+1)*a%ManufacturedNpoinX

		      ielem = DownRightElem0 + (kflr-1)*3 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = RightGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+2) = RightGhost0 + (kflr+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+3) = (kflr+1)*a%ManufacturedNPoinX*a%ManufacturedNPoinX
		      a%lnods(ispos+4) = DownGhost0 + (kflr+1)*a%ManufacturedNpoinX

		      ielem = DownRightElem0 + (kflr-1)*3 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = kflr*a%ManufacturedNPoinX*a%ManufacturedNPoinX
		      a%lnods(ispos+2) = DownGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+3) = RightGhost0 + kflr*a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownGhost0 + (kflr+1)*a%ManufacturedNpoinX
	      enddo

      endif   
         
      if (TopLeftGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + BlockTop + BlockLeft + colright*(a%ManufacturedNpoinX-1)
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(TopLeftGhost + ipoin) = baseGhost0 + ipoin
		      ProcessorList(TopLeftGhost + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX - ManufacturedNprocX
	      enddo
         
	      do irow = 1,a%ManufacturedNpoinX-1

		      ielem = TopLeftElem0 + (irow-1)*6 + 1
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+3) = TopLeftGhost + irow
		      a%lnods(ispos+4) = TopLeftGhost + irow + 1

		      ielem = TopLeftElem0 + (irow-1)*6 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+2) = TopGhost0 + irow
		      a%lnods(ispos+3) = TopLeftGhost + irow
		      a%lnods(ispos+4) = TopLeftGhost + irow + 1

		      ielem = TopLeftElem0 + (irow-1)*6 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+2) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+4) = TopLeftGhost + irow + 1

		      ielem = TopLeftElem0 + (irow-1)*6 + 4
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = TopGhost0 + irow + 1 
		      a%lnods(ispos+2) = TopGhost0 + irow
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+4) = TopLeftGhost + irow + 1

		      ielem = TopLeftElem0 + (irow-1)*6 + 5
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+4) = TopLeftGhost + irow + 1

		      ielem = TopLeftElem0 + (irow-1)*6 + 6
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+3) = TopGhost0 + irow + 1 
		      a%lnods(ispos+4) = TopLeftGhost + irow + 1
	      enddo

      endif   
         
      if (TopRightGhost /= -1) then
         
	      do irow = 1,a%ManufacturedNpoinX-1

		      ielem = TopRightElem0 + (irow-1)*3 + 1
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + colright*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+2) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+3) = TopGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+4) = TopGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 

		      ielem = TopRightElem0 + (irow-1)*3 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + colright*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + colright*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+3) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+4) = TopGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 

		      ielem = TopRightElem0 + (irow-1)*3 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + colright*(a%ManufacturedNpoinX-1) + irow + rowdown
		      a%lnods(ispos+2) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1
		      a%lnods(ispos+3) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+4) = TopGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 
	      enddo

      endif   
     
      if (TopUpGhost /= -1) then

	      do jcol = 1,a%ManufacturedNpoinX-1

		      ielem = TopUpElem0 + (jcol-1)*3 + 1
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol 
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + colright*(jcol-1) + 1
		      a%lnods(ispos+3) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1
		      a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1

		      ielem = TopUpElem0 + (jcol-1)*3 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + colright*(jcol-1) + 1
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + colright*jcol + 1 
		      a%lnods(ispos+3) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1 
		      a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1

		      ielem = TopUpElem0 + (jcol-1)*3 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + colright*jcol + 1 
		      a%lnods(ispos+3) = TopGhost0 + jcol*a%ManufacturedNpoinX + 1 
		      a%lnods(ispos+4) = TopGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1

	      enddo
      endif   
            
      if (TopDownGhost /= -1) then
	      baseGhost0 = a%poin0 + 1 + BlockTop + BlockDown  
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(TopDownGhost + ipoin) = baseGhost0 + colright*(ipoin-1)
		      ProcessorList(TopDownGhost + ipoin) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX +1
	      enddo

	      do jcol = 1,a%ManufacturedNpoinX-1

		      ielem = TopDownElem0 + (jcol-1)*6 + 1 
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*(jcol-1) + 1
		      a%lnods(ispos+2) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+3) = TopGhost0 + jcol*a%ManufacturedNpoinX
		      a%lnods(ispos+4) = TopDownGhost + jcol

		      ielem = TopDownElem0 + (jcol-1)*6 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+2) = TopGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+3) = TopGhost0 + jcol*a%ManufacturedNpoinX
		      a%lnods(ispos+4) = TopDownGhost + jcol

		      ielem = TopDownElem0 + (jcol-1)*6 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*(jcol-1) + 1
		      a%lnods(ispos+2) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol 
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+4) = TopDownGhost + jcol

		      ielem = TopDownElem0 + (jcol-1)*6 + 4
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = TopDownGhost + jcol + 1
		      a%lnods(ispos+2) = TopGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+4) = TopDownGhost + jcol

		      ielem = TopDownElem0 + (jcol-1)*6 + 5
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol 
		      a%lnods(ispos+2) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1
		      a%lnods(ispos+3) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+4) = TopDownGhost + jcol

		      ielem = TopDownElem0 + (jcol-1)*6 + 6
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = floortop*(a%ManufacturedNpoinX-1) + rowdown*(a%ManufacturedNpoinX-1) + colright*jcol + 1
		      a%lnods(ispos+2) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + jcol + 1
		      a%lnods(ispos+3) = TopDownGhost + jcol + 1
		      a%lnods(ispos+4) = TopDownGhost + jcol

	      enddo
      endif   

      if (BottomLeftGhost /= -1) then
         
	      do irow = 1,a%ManufacturedNpoinX-1

		      ielem = BottomLeftElem0 + (irow-1)*3 + 1 
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + irow
		      a%lnods(ispos+2) = irow
		      a%lnods(ispos+3) = LeftGhost0 + irow
		      a%lnods(ispos+4) = LeftGhost0 + irow + 1

		      ielem = BottomLeftElem0 + (irow-1)*3 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = irow + rowdown
		      a%lnods(ispos+2) = irow
		      a%lnods(ispos+3) = BottomGhost0 + irow
		      a%lnods(ispos+4) = LeftGhost0 + irow + 1

		      ielem = BottomLeftElem0 + (irow-1)*3 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + irow
		      a%lnods(ispos+2) = BottomGhost0 + irow + 1 
		      a%lnods(ispos+3) = irow + rowdown
		      a%lnods(ispos+4) = LeftGhost0 + irow + 1
	      enddo

      endif   

      if (BottomRightGhost /= -1) then
	      !GhostPoints Global numbering
	      baseGhost0 = a%poin0 + BlockBottom + BlockRight + (a%ManufacturedNpoinX-1)*floortop
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(BottomRightGhost + ipoin) = baseGhost0 + ipoin
		      ProcessorList(BottomRightGhost + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX + ManufacturedNprocX
	      enddo
         
	      do irow = 1,a%ManufacturedNpoinX-1

		      ielem = BottomRightElem0 + (irow-1)*6 + 1 
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+2) = BottomRightGhost + irow
		      a%lnods(ispos+3) = colright*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+4) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown

		      ielem = BottomRightElem0 + (irow-1)*6 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomRightGhost + irow
		      a%lnods(ispos+2) = RightGhost0 + irow
		      a%lnods(ispos+3) = colright*(a%ManufacturedNpoinX-1) + irow
		      a%lnods(ispos+4) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown

		      ielem = BottomRightElem0 + (irow-1)*6 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow
		      a%lnods(ispos+2) = BottomGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 
		      a%lnods(ispos+3) = BottomRightGhost + irow
		      a%lnods(ispos+4) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown

		      ielem = BottomRightElem0 + (irow-1)*6 + 4
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = RightGhost0 + irow + 1
		      a%lnods(ispos+2) = RightGhost0 + irow
		      a%lnods(ispos+3) = BottomRightGhost + irow
		      a%lnods(ispos+4) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown

		      ielem = BottomRightElem0 + (irow-1)*6 + 5
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + irow + 1 
		      a%lnods(ispos+2) = BottomRightGhost + irow + 1
		      a%lnods(ispos+3) = BottomRightGhost + irow
		      a%lnods(ispos+4) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown

		      ielem = BottomRightElem0 + (irow-1)*6 + 6
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomRightGhost + irow
		      a%lnods(ispos+2) = BottomRightGhost + irow + 1
		      a%lnods(ispos+3) = RightGhost0 + irow + 1
		      a%lnods(ispos+4) = colright*(a%ManufacturedNpoinX-1) + irow + rowdown
	      enddo

      endif   
      
      if (BottomUpGhost /= -1) then
	      baseGhost0 = a%poin0 + 1 + BlockBottom + BlockUp + (a%ManufacturedNpoinX-1)*floortop + rowdown*(a%ManufacturedNpoinX-1) 
	      do ipoin = 1,a%ManufacturedNpoinX
		      LocalToGlobal(BottomUpGhost + ipoin) = baseGhost0 + colright*(ipoin-1)
		      ProcessorList(BottomUpGhost + ipoin) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX -1
	      enddo

	      do jcol = 1,a%ManufacturedNpoinX-1

		      ielem = BottomUpElem0 + (jcol-1)*6 + 1 
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomUpGhost + jcol
		      a%lnods(ispos+2) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+3) = UpGhost0 + jcol 
		      a%lnods(ispos+4) = colright*(jcol-1) + 1

		      ielem = BottomUpElem0 + (jcol-1)*6 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+2) = UpGhost0 + jcol + 1
		      a%lnods(ispos+3) = UpGhost0 + jcol 
		      a%lnods(ispos+4) = colright*(jcol-1) + 1

		      ielem = BottomUpElem0 + (jcol-1)*6 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomUpGhost + jcol
		      a%lnods(ispos+2) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+3) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+4) = colright*(jcol-1) + 1

		      ielem = BottomUpElem0 + (jcol-1)*6 + 4
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = colright*jcol + 1 
		      a%lnods(ispos+2) = UpGhost0 + jcol + 1
		      a%lnods(ispos+3) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+4) = colright*(jcol-1) + 1

		      ielem = BottomUpElem0 + (jcol-1)*6 + 5
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + (jcol-1)*a%ManufacturedNpoinX + 1
		      a%lnods(ispos+2) = BottomGhost0 + jcol*a%ManufacturedNpoinX + 1 
		      a%lnods(ispos+3) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+4) = colright*(jcol-1) + 1

		      ielem = BottomUpElem0 + (jcol-1)*6 + 6
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomUpGhost + jcol + 1
		      a%lnods(ispos+2) = BottomGhost0 + jcol*a%ManufacturedNpoinX + 1 
		      a%lnods(ispos+3) = colright*jcol + 1 
		      a%lnods(ispos+4) = colright*(jcol-1) + 1

	      enddo
      endif   
            
      if (BottomDownGhost /= -1) then

	      do jcol = 1,a%ManufacturedNpoinX-1

		      ielem = BottomDownElem0 + (jcol-1)*3 + 1 
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + jcol*a%ManufacturedNpoinX
		      a%lnods(ispos+2) = BottomGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+3) = a%ManufacturedNpoinX + colright*(jcol-1) 
		      a%lnods(ispos+4) = DownGhost0 + jcol

		      ielem = BottomDownElem0 + (jcol-1)*3 + 2
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = BottomGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+2) = a%ManufacturedNpoinX + colright*jcol
		      a%lnods(ispos+3) = a%ManufacturedNpoinX + colright*(jcol-1) 
		      a%lnods(ispos+4) = DownGhost0 + jcol 

		      ielem = BottomDownElem0 + (jcol-1)*3 + 3
		      ispos = (ielem-1)*a%pnods(1)
		      a%lnods(ispos+1) = DownGhost0 + jcol + 1
		      a%lnods(ispos+2) = a%ManufacturedNpoinX + colright*jcol
		      a%lnods(ispos+3) = BottomGhost0 + (jcol+1)*a%ManufacturedNpoinX
		      a%lnods(ispos+4) = DownGhost0 + jcol 

	      enddo
      endif   

      if (TopUpLeftGhost /= -1) then

	      ielem = TopUpLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = LeftGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+2) = (a%ManufacturedNpoinX-1)*floortop + 1
	      a%lnods(ispos+3) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+4) = TopLeftGhost + 1

	      ielem = TopUpLeftElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = UpGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+2) = (a%ManufacturedNpoinX-1)*floortop + 1
	      a%lnods(ispos+3) = TopGhost0 + 1
	      a%lnods(ispos+4) = TopLeftGhost + 1
      endif   
      
      if (TopDownLeftGhost /= -1) then
	      LocalToGlobal(TopDownLeftGhost) = a%poin0 + 1 + BlockTop + BlockDown + BlockLeft + colright*(a%ManufacturedNpoinX-1) 
	      ProcessorList(TopDownLeftGhost) = a%MPIrank + ManufacturedNprocX*ManufacturedNprocX + 1 - ManufacturedNprocX

	      ielem = TopDownLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = LeftGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+2) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+3) = TopLeftGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopDownLeftGhost  

	      ielem = TopDownLeftElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = TopGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+3) = TopLeftGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopDownLeftGhost  

	      ielem = TopDownLeftElem0 + 3
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = LeftGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+2) = DownLeftGhost + a%ManufacturedNpoinX 
	      a%lnods(ispos+3) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopDownLeftGhost  

	      ielem = TopDownLeftElem0 + 4
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = TopDownGhost + 1
	      a%lnods(ispos+2) = TopGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+3) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopDownLeftGhost  

	      ielem = TopDownLeftElem0 + 5
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = DownLeftGhost + a%ManufacturedNpoinX 
	      a%lnods(ispos+2) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+3) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopDownLeftGhost  

	      ielem = TopDownLeftElem0 + 6
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = (a%ManufacturedNpoinX-1)*floortop + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = DownGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+3) = TopDownGhost + 1
	      a%lnods(ispos+4) = TopDownLeftGhost  
         
      endif

      if (TopUpRightGhost /= -1) then

	      ielem = TopUpRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = UpGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+2) = (a%ManufacturedNpoinX-1)*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+3) = UpRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopGhost0 + 1 + a%ManufacturedNPoinX*(a%ManufacturedNPoinX-1) 

	      ielem = TopUpRightElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = (a%ManufacturedNpoinX-1)*floortop + 1 + colright*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+2) = RightGhost0 + (a%ManufacturedNpoinX-1)*a%ManufacturedNpoinX + 1
	      a%lnods(ispos+3) = UpRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopGhost0 + 1 + a%ManufacturedNPoinX*(a%ManufacturedNPoinX-1) 

      endif   
         
      if (TopDownRightGhost /= -1) then

	      ielem = TopDownRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = a%ManufacturedNPoinX*a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+2) = RightGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+3) = TopGhost0 + a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+4) = TopDownGhost + a%ManufacturedNpoinX

	      ielem = TopDownRightElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = a%ManufacturedNPoinX*a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+2) = DownGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+3) = RightGhost0 + a%ManufacturedNpoinX*a%ManufacturedNpoinX
	      a%lnods(ispos+4) = TopDownGhost + a%ManufacturedNpoinX

      endif   

      if (BottomUpLeftGhost /= -1) then

	      ielem = BottomUpLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = 1
	      a%lnods(ispos+2) = UpGhost0 + 1
	      a%lnods(ispos+3) = BottomUpGhost + 1
	      a%lnods(ispos+4) = LeftGhost0 + 1

	      ielem = BottomUpLeftElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomUpGhost + 1
	      a%lnods(ispos+2) = BottomGhost0 + 1
	      a%lnods(ispos+3) = 1
	      a%lnods(ispos+4) = LeftGhost0 + 1
      endif   
      
      if (BottomDownLeftGhost /= -1) then

	      ielem = BottomDownLeftElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = a%ManufacturedNpoinX
	      a%lnods(ispos+3) = LeftGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = DownLeftGhost + 1 

	      ielem = BottomDownLeftElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = DownGhost0 + 1
	      a%lnods(ispos+2) = a%ManufacturedNpoinX
	      a%lnods(ispos+3) = BottomGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = DownLeftGhost + 1 

      endif
      
      if (BottomUpRightGhost /= -1) then
	      LocalToGlobal(BottomUpRightGhost) = a%poin0 + 1 + BlockBottom + BlockUp + BlockRight + rowdown*(a%ManufacturedNpoinX-1) + (a%ManufacturedNpoinX-1)*floortop
	      ProcessorList(BottomUpRightGhost) = a%MPIrank - ManufacturedNprocX*ManufacturedNprocX - 1 + ManufacturedNprocX

	      ielem = BottomUpRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomUpGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = BottomUpRightGhost
	      a%lnods(ispos+3) = UpGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 

	      ielem = BottomUpRightElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomUpRightGhost
	      a%lnods(ispos+2) = UpRightGhost + 1
	      a%lnods(ispos+3) = UpGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 

	      ielem = BottomUpRightElem0 + 3
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomUpGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = BottomGhost0 + 1 + a%ManufacturedNPoinX*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+3) = BottomUpRightGhost
	      a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 

	      ielem = BottomUpRightElem0 + 4
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = RightGhost0 + 1
	      a%lnods(ispos+2) = UpRightGhost + 1
	      a%lnods(ispos+3) = BottomUpRightGhost
	      a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 

	      ielem = BottomUpRightElem0 + 5
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomGhost0 + 1 + a%ManufacturedNPoinX*(a%ManufacturedNPoinX-1) 
	      a%lnods(ispos+2) = BottomRightGhost + 1
	      a%lnods(ispos+3) = BottomUpRightGhost
	      a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 

	      ielem = BottomUpRightElem0 + 6
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomUpRightGhost
	      a%lnods(ispos+2) = BottomRightGhost + 1
	      a%lnods(ispos+3) = RightGhost0 + 1
	      a%lnods(ispos+4) = 1 + colright*(a%ManufacturedNPoinX-1) 
         
      endif   
         
      if (BottomDownRightGhost /= -1) then

	      ielem = BottomDownRightElem0 + 1
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomGhost0 + a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+2) = BottomRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+3) = a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+4) = DownGhost0 + a%ManufacturedNpoinX

	      ielem = BottomDownRightElem0 + 2
	      ispos = (ielem-1)*a%pnods(1)
	      a%lnods(ispos+1) = BottomRightGhost + a%ManufacturedNpoinX
	      a%lnods(ispos+2) = RightGhost0 + a%ManufacturedNpoinX
	      a%lnods(ispos+3) = a%ManufacturedNPoinX*a%ManufacturedNPoinX
	      a%lnods(ispos+4) = DownGhost0 + a%ManufacturedNpoinX

      endif   

      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      call a%Memor%dealloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%dealloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !External normal, we simply allocate it
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ManufacturedMesh')

   endif

   if (a%ndime == 3 .and. a%ManufacturedType == 'CYLIN') then

      InternalRadius = a%ManufacturedInternalRadius 
      ExternalRadius = a%ManufacturedExternalRadius
      AxialLength =  a%ManufacturedAxialLength 
      PointsInRadius = a%ManufacturedNpoinX 
      DivisionsInAngle  = a%ManufacturedDivisionsInAngle 
      PointsInLength  = a%ManufacturedNpoinZ 
      
      DeltaAngleProc = 2_rp*dacos(-1.D0)/real(a%MPIsize,rp)
      DeltaAngle     = DeltaAngleProc/DivisionsInAngle

      kfl_InternalRadiusAtZero = .false.
      if (abs(InternalRadius) <= 1.0D-6) then
          kfl_InternalRadiusAtZero = .true.
          InternalRadius = ExternalRadius/(PointsInRadius-1)
          PointsInRadius = PointsInRadius - 1
      endif
          
      SpanInRadialCoord = (ExternalRadius-InternalRadius)/(PointsInRadius-1)
      SpanInLengthCoord = AxialLength/(PointsInLength-1)

      LocalAngleDivisionNpoin = 0
      LocalAngleNelem = 0

      do irow = 1,PointsInRadius
         TimesAngleIsDivided = 2_ip**(irow-1)
         LocalAngleNelem = LocalAngleNelem + 9_ip*TimesAngleIsDivided*(PointsInLength-1) 
         LocalAngleDivisionNpoin = LocalAngleDivisionNpoin + TimesAngleIsDivided
      enddo
      LocalAngleNelem = LocalAngleNelem - 9_ip*2_ip**(PointsInRadius-1)*(PointsInLength-1)   
      
      a%npoinLocal = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength
      a%gnpoin = a%npoinLocal*a%MPIsize

      a%nelem = (LocalAngleNelem*DivisionsInAngle) - 6_ip*(PointsInRadius-1)*(PointsInLength-1)    
      a%gnelem = LocalAngleNelem*DivisionsInAngle*a%MPIsize

      a%poin0 = a%MPIrank*a%npoinLocal
      a%elem0 = a%MPIrank*(LocalAngleNelem*DivisionsInAngle) 

      if (kfl_InternalRadiusAtZero) then 
           a%gnpoin = a%gnpoin + PointsInLength 
           a%gnelem = a%gnelem + 3_ip*DivisionsInAngle*a%MPIsize*(PointsInLength-1)  
           a%elem0 = a%MPIrank*(LocalAngleNelem+3*(PointsInLength-1))*DivisionsInAngle 
           if (a%MPIrank == a%MPIsize-1) then
               !Im last pc with inner point as local
               a%npoinLocal = a%npoinLocal + PointsInLength
               a%nelem = a%nelem + 3*(PointsInLength-1)*(DivisionsInAngle - 1_ip) 
           endif
      endif

      a%npoinGhost = 0
  
      UpGhost0 = -1
      DownGhost0 = -1
      InnerGhost0 = -1
      OtherPcsGhost0 = -1

      !First the upper and lower processors element connections (required for r=0)
      if (a%MPIsize == 1) then
         !Special case of serial run (Local elements at teta=0 self-connection) 
         a%nelem = a%nelem + 6_ip*(PointsInRadius-1)*(PointsInLength-1)
      else
         !Parallel run
         !Upper processor conection
         UpGhost0 = a%npoinLocal + a%npoinGhost
         UpElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + PointsInRadius*PointsInLength
         a%nelem = a%nelem + 6_ip*(PointsInRadius-1)*(PointsInLength-1)
         !Lower processor conection
         DownGhost0 = a%npoinLocal + a%npoinGhost
         DownElem0 = a%nelem
         a%npoinGhost = a%npoinGhost + PointsInRadius*PointsInLength
         a%nelem = a%nelem + 6_ip*(PointsInRadius-1)*(PointsInLength-1)
      endif

      !Then the elements at radius equal zero 
      if (kfl_InternalRadiusAtZero) then
         if (a%MPIsize == 1) then
            !Serial run (local element at r=0 self-connection)
            a%nelem = a%nelem + 3_ip*(PointsInLength-1)
         else
            !Parallel run
            if (a%MPIrank /= a%MPIsize-1) then
               !Im not last pc: Inner processor conection
               !I use localpoints and InnerGhost0
               InnerGhost0 = a%npoinLocal + a%npoinGhost
               a%npoinGhost = a%npoinGhost + PointsInLength
               InnerElem0 = a%nelem
               a%nelem = a%nelem + 3_ip*(DivisionsInAngle-1)*(PointsInLength-1)
               !Im not last pc: Inner-Upper processor conection
               !I use UpGhost0 and InnerGhost0
               UpInnerElem0 = a%nelem
               a%nelem = a%nelem + 3_ip*(PointsInLength-1)  
               !Im not last pc: Inner-Lower processor conection
               !I use DownGhost0 and InnerGhost0
               DownInnerElem0 = a%nelem
               a%nelem = a%nelem + 3_ip*(PointsInLength-1) 
            else
               !Im last pc: All inner elements of other pcs except upper and lower
               OtherPcsGhost0 = a%npoinLocal + a%npoinGhost
               a%npoinGhost = a%npoinGhost + (DivisionsInAngle*(a%MPIsize-1) - 2)*PointsInLength
               OtherInnerElem0 = a%nelem
               a%nelem = a%nelem + 3_ip*(DivisionsInAngle*(a%MPIsize-1) - 1)*(PointsInLength-1)
               !Im last pc: Inner-Upper processor conection
               !I use UpGhost0 and InnerLocalpoint
               UpInnerElem0 = a%nelem
               a%nelem = a%nelem + 3_ip*(PointsInLength-1)
               !Im last pc: Inner-Lower processor conection
               !I use DownGhost0 and InnerLocalpoint
               DownInnerElem0 = a%nelem
               a%nelem = a%nelem + 3_ip*(PointsInLength-1)
            endif
         endif
      endif


      a%npoin = a%npoinLocal+a%npoinGhost

      !Fill local nodes and elements
      call a%Memor%alloc(1_ip,a%pnods,'pnods','ManufacturedMesh')
      call a%Memor%alloc(a%nelem*4_ip,a%lnods,'lnods','ManufacturedMesh')
      call a%Memor%alloc(3_ip,a%npoin,a%coord,'coord','ManufacturedMesh')
      
      call a%Memor%alloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !Build local elements, nodes and coord

      !Local Coordinates
      BaseProcessorAngle = (a%MPIrank-1)*DeltaAngleProc 
      ipoin = 0
      do kflr = 1,PointsInLength
         coordz = (kflr-1)*SpanInLengthCoord
         do jcol=1,DivisionsInAngle
           BaseAngle = BaseProcessorAngle + (jcol-1)*DeltaAngle
           !Radial nodes at BaseAngle (to locate the angular nodes)
           do irow=1,PointsInRadius
             BaseDivisionRadius =  InternalRadius + (irow-1)*SpanInRadialCoord
             coordx = BaseDivisionRadius*cos(BaseAngle)
             coordy = BaseDivisionRadius*sin(BaseAngle)
             ipoin = ipoin + 1
             a%coord(1,ipoin) = coordx
             a%coord(2,ipoin) = coordy
             a%coord(3,ipoin) = coordz
           enddo
           !Angular nodes at each radius
           do irow=1,PointsInRadius
              TimesAngleIsDivided = 2_ip**(irow-1)
              BaseDivisionRadius =  InternalRadius + (irow-1)*SpanInRadialCoord
              do groupi = 1,TimesAngleIsDivided-1
                 BaseGroupAngle = BaseAngle + groupi*DeltaAngle/TimesAngleIsDivided
                 coordx = BaseDivisionRadius*cos(BaseGroupAngle)
                 coordy = BaseDivisionRadius*sin(BaseGroupAngle)
                 ipoin = ipoin + 1
                 a%coord(1,ipoin) = coordx
                 a%coord(2,ipoin) = coordy
                 a%coord(3,ipoin) = coordz
              enddo
           enddo
         enddo
      enddo
   
      if (kfl_InternalRadiusAtZero) then
        if (a%MPIrank == a%MPIsize-1) then
           do kflr = 1,PointsInLength
              coordz = (kflr-1)*SpanInLengthCoord
              ipoin = ipoin + 1
              a%coord(1,ipoin) = 0.0_rp
              a%coord(2,ipoin) = 0.0_rp
              a%coord(3,ipoin) = coordz
           enddo
        endif
      endif

      !LocalToGlobal for local points
      do ipoin = 1,a%npoinLocal
         LocalToGlobal(ipoin) = a%poin0+ipoin
         ProcessorList(ipoin) = a%MPIrank
      enddo

      BlockUp   = LocalAngleDivisionNpoin*DivisionsInAngle
      BlockDown = -BlockUp

      !Local Elements
      a%pnods(1) = 4
      ielem = 0

      do kflr = 1,PointsInLength-1
         basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
         do jcol=1,DivisionsInAngle-1
            ZeroPointAngle = basepointflr + (jcol-1)*LocalAngleDivisionNpoin
            !First radial group
            irow = 1
            floorbottomleftpoint = ZeroPointAngle + irow 
            floorrightpoint = floorbottomleftpoint + PointsInRadius
            floorbottomrightpoint = floorbottomleftpoint + 1
            floortopleftpoint = ZeroPointAngle + LocalAngleDivisionNpoin + 1
            floortoprightpoint = floortopleftpoint + irow
            ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
            ceilingrightpoint = floorrightpoint + BlockUp 
            ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
            ceilingtopleftpoint = floortopleftpoint + BlockUp
            ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
   
            !Lower tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomrightpoint
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorrightpoint 
            a%lnods(ispos+4) = floorbottomleftpoint  

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floorrightpoint 
            a%lnods(ispos+2) = floorbottomrightpoint   
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint

            !Middle tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint  
            a%lnods(ispos+2) = ceilingtopleftpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint 
            a%lnods(ispos+4) = floortopleftpoint 

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint 
            a%lnods(ispos+2) = floorrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint
            a%lnods(ispos+4) = floortopleftpoint   

            !Upper tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingtopleftpoint 
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = floortoprightpoint 
            a%lnods(ispos+4) = floortopleftpoint  

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floortoprightpoint 
            a%lnods(ispos+2) = floorrightpoint   
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            !Second radial groups
            irow = 2
            FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
            !Second radial Lower group
            floorbottomleftpoint = ZeroPointAngle + irow 
            floorrightpoint = FirstInnerPoint + 2**(irow-1) - 1
            floorbottomrightpoint = floorbottomleftpoint + 1
            floortopleftpoint =  FirstInnerPoint
            floortoprightpoint = floorrightpoint + 1
            ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
            ceilingrightpoint = floorrightpoint + BlockUp 
            ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
            ceilingtopleftpoint = floortopleftpoint + BlockUp
            ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
            !Lower tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomrightpoint
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorrightpoint 
            a%lnods(ispos+4) = floorbottomleftpoint  

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floorrightpoint 
            a%lnods(ispos+2) = floorbottomrightpoint   
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint

            !Middle tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint  
            a%lnods(ispos+2) = ceilingtopleftpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint 
            a%lnods(ispos+4) = floortopleftpoint 

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint 
            a%lnods(ispos+2) = floorrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint
            a%lnods(ispos+4) = floortopleftpoint   

            !Upper tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingtopleftpoint 
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = floortoprightpoint 
            a%lnods(ispos+4) = floortopleftpoint  

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floortoprightpoint 
            a%lnods(ispos+2) = floorrightpoint   
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint
   
            !Second radial Upper group  
            floorbottomleftpoint = FirstInnerPoint
            floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
            floorbottomrightpoint = floorrightpoint -1
            floortopleftpoint =  ZeroPointAngle + irow + LocalAngleDivisionNpoin
            floortoprightpoint = floortopleftpoint + 1
            ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
            ceilingrightpoint = floorrightpoint + BlockUp 
            ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
            ceilingtopleftpoint = floortopleftpoint + BlockUp
            ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
            !Lower tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomrightpoint
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorrightpoint 
            a%lnods(ispos+4) = floorbottomleftpoint  

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floorrightpoint 
            a%lnods(ispos+2) = floorbottomrightpoint   
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint

            !Middle tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint  
            a%lnods(ispos+2) = ceilingtopleftpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint 
            a%lnods(ispos+4) = floortopleftpoint 

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint 
            a%lnods(ispos+2) = floorrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint
            a%lnods(ispos+4) = floortopleftpoint   

            !Upper tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingtopleftpoint 
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = floortoprightpoint 
            a%lnods(ispos+4) = floortopleftpoint  

            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floortoprightpoint 
            a%lnods(ispos+2) = floorrightpoint   
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint
 

   
            !Other radial groups
            do irow=3,PointsInRadius-1
              TimesAngleIsDivided = 2_ip**(irow-1)
              FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
              !Lower group
              floorbottomleftpoint = ZeroPointAngle + irow
              floorrightpoint = FirstInnerPoint + 2**(irow-1) - 1
              floorbottomrightpoint = floorbottomleftpoint + 1
              floortopleftpoint =  FirstInnerPoint
              floortoprightpoint = floorrightpoint + 1
              ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
              ceilingrightpoint = floorrightpoint + BlockUp 
              ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
              ceilingtopleftpoint = floortopleftpoint + BlockUp
              ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
              !Lower tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint 
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = ceilingbottomrightpoint
              a%lnods(ispos+4) = floorbottomleftpoint
  
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomrightpoint
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = floorrightpoint 
              a%lnods(ispos+4) = floorbottomleftpoint  
  
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = floorrightpoint 
              a%lnods(ispos+2) = floorbottomrightpoint   
              a%lnods(ispos+3) = ceilingbottomrightpoint
              a%lnods(ispos+4) = floorbottomleftpoint
  
              !Middle tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint  
              a%lnods(ispos+2) = ceilingtopleftpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
  
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint 
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint 
              a%lnods(ispos+4) = floortopleftpoint 
  
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint 
              a%lnods(ispos+2) = floorrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint
              a%lnods(ispos+4) = floortopleftpoint   
  
              !Upper tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingtopleftpoint 
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
  
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = floortoprightpoint 
              a%lnods(ispos+4) = floortopleftpoint  
  
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = floortoprightpoint 
              a%lnods(ispos+2) = floorrightpoint   
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
   
              !Inner groups
              do groupi = 2,TimesAngleIsDivided-1
                 floorbottomleftpoint = FirstInnerPoint + groupi - 2
                 floorrightpoint = FirstInnerPoint + 2**(irow-1) - 1 + (groupi-1)*2
                 floorbottomrightpoint = floorrightpoint - 1
                 floortopleftpoint =  floorbottomleftpoint + 1
                 floortoprightpoint = floorrightpoint + 1
                 ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
                 ceilingrightpoint = floorrightpoint + BlockUp 
                 ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
                 ceilingtopleftpoint = floortopleftpoint + BlockUp
                 ceilingtoprightpoint = floortoprightpoint + BlockUp 
      
              
                !Lower tetras 
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = ceilingbottomrightpoint
                a%lnods(ispos+4) = floorbottomleftpoint
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomrightpoint
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorrightpoint 
                a%lnods(ispos+4) = floorbottomleftpoint  
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = floorrightpoint 
                a%lnods(ispos+2) = floorbottomrightpoint   
                a%lnods(ispos+3) = ceilingbottomrightpoint
                a%lnods(ispos+4) = floorbottomleftpoint
    
                !Middle tetras 
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
    
                !Upper tetras 
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint

              enddo

              !Upper group  
              floorbottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
              floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
              floorbottomrightpoint = floorrightpoint - 1
              floortopleftpoint =  ZeroPointAngle + irow + LocalAngleDivisionNpoin
              floortoprightpoint = floortopleftpoint + 1
              ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
              ceilingrightpoint = floorrightpoint + BlockUp 
              ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
              ceilingtopleftpoint = floortopleftpoint + BlockUp
              ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
              !Lower tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint 
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = ceilingbottomrightpoint
              a%lnods(ispos+4) = floorbottomleftpoint
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomrightpoint
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = floorrightpoint 
              a%lnods(ispos+4) = floorbottomleftpoint  
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = floorrightpoint 
              a%lnods(ispos+2) = floorbottomrightpoint   
              a%lnods(ispos+3) = ceilingbottomrightpoint
              a%lnods(ispos+4) = floorbottomleftpoint
    
              !Middle tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint  
              a%lnods(ispos+2) = ceilingtopleftpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint 
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint 
              a%lnods(ispos+4) = floortopleftpoint 
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint 
              a%lnods(ispos+2) = floorrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint
              a%lnods(ispos+4) = floortopleftpoint   
    
              !Upper tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingtopleftpoint 
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = floortoprightpoint 
              a%lnods(ispos+4) = floortopleftpoint  
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = floortoprightpoint 
              a%lnods(ispos+2) = floorrightpoint   
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
   
            enddo
         enddo
   
         !This part cares for not connecting top (ghost) points in local elements
         ZeroPointAngle = basepointflr + (DivisionsInAngle-1)*LocalAngleDivisionNpoin
         !First radial group
         irow = 1
         floorbottomleftpoint = ZeroPointAngle + irow 
         floorrightpoint = floorbottomleftpoint + PointsInRadius
         floorbottomrightpoint = floorbottomleftpoint + 1
         ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
         ceilingrightpoint = floorrightpoint + BlockUp 
         ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
   
         !Lower tetras 
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomleftpoint 
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = ceilingbottomrightpoint
         a%lnods(ispos+4) = floorbottomleftpoint
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomrightpoint
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = floorrightpoint 
         a%lnods(ispos+4) = floorbottomleftpoint  
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = floorrightpoint 
         a%lnods(ispos+2) = floorbottomrightpoint   
         a%lnods(ispos+3) = ceilingbottomrightpoint
         a%lnods(ispos+4) = floorbottomleftpoint
         
         !Second radial groups
         irow = 2
         FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
         !Second radial Lower group
         floorbottomleftpoint = ZeroPointAngle + irow 
         floorrightpoint = FirstInnerPoint + 2**(irow-1) - 1
         floorbottomrightpoint = floorbottomleftpoint + 1
         floortopleftpoint =  FirstInnerPoint
         floortoprightpoint = floorrightpoint + 1
         ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
         ceilingrightpoint = floorrightpoint + BlockUp 
         ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
         ceilingtopleftpoint = floortopleftpoint + BlockUp
         ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
         !Lower tetras 
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomleftpoint 
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = ceilingbottomrightpoint
         a%lnods(ispos+4) = floorbottomleftpoint
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomrightpoint
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = floorrightpoint 
         a%lnods(ispos+4) = floorbottomleftpoint  
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = floorrightpoint 
         a%lnods(ispos+2) = floorbottomrightpoint   
         a%lnods(ispos+3) = ceilingbottomrightpoint
         a%lnods(ispos+4) = floorbottomleftpoint
      
         !Middle tetras 
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomleftpoint  
         a%lnods(ispos+2) = ceilingtopleftpoint
         a%lnods(ispos+3) = ceilingrightpoint
         a%lnods(ispos+4) = floortopleftpoint
    
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomleftpoint 
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = floorbottomleftpoint 
         a%lnods(ispos+4) = floortopleftpoint 
    
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingrightpoint 
         a%lnods(ispos+2) = floorrightpoint
         a%lnods(ispos+3) = floorbottomleftpoint
         a%lnods(ispos+4) = floortopleftpoint   
      
         !Upper tetras 
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingtopleftpoint 
         a%lnods(ispos+2) = ceilingtoprightpoint
         a%lnods(ispos+3) = ceilingrightpoint
         a%lnods(ispos+4) = floortopleftpoint
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingrightpoint
         a%lnods(ispos+2) = ceilingtoprightpoint
         a%lnods(ispos+3) = floortoprightpoint 
         a%lnods(ispos+4) = floortopleftpoint  
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = floortoprightpoint 
         a%lnods(ispos+2) = floorrightpoint   
         a%lnods(ispos+3) = ceilingrightpoint
         a%lnods(ispos+4) = floortopleftpoint
   
         !Second radial Upper group  
         floorbottomleftpoint = FirstInnerPoint
         floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
         floorbottomrightpoint = floorrightpoint -1
         ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
         ceilingrightpoint = floorrightpoint + BlockUp 
         ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
   
         !Lower tetras 
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomleftpoint 
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = ceilingbottomrightpoint
         a%lnods(ispos+4) = floorbottomleftpoint
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = ceilingbottomrightpoint
         a%lnods(ispos+2) = ceilingrightpoint
         a%lnods(ispos+3) = floorrightpoint 
         a%lnods(ispos+4) = floorbottomleftpoint  
      
         ielem = ielem+1
         ispos = (ielem-1)*a%pnods(1)
         a%lnods(ispos+1) = floorrightpoint 
         a%lnods(ispos+2) = floorbottomrightpoint   
         a%lnods(ispos+3) = ceilingbottomrightpoint
         a%lnods(ispos+4) = floorbottomleftpoint
         
         !Other radial groups
         do irow=3,PointsInRadius-1
           TimesAngleIsDivided = 2_ip**(irow-1)
           FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
           !Lower group
           floorbottomleftpoint = ZeroPointAngle + irow
           floorrightpoint = FirstInnerPoint + 2**(irow-1) - 1
           floorbottomrightpoint = floorbottomleftpoint + 1
           floortopleftpoint =  FirstInnerPoint
           floortoprightpoint = floorrightpoint + 1
           ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
           ceilingrightpoint = floorrightpoint + BlockUp 
           ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
           ceilingtopleftpoint = floortopleftpoint + BlockUp
           ceilingtoprightpoint = floortoprightpoint + BlockUp 
   
            !Lower tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint
         
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomrightpoint
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorrightpoint 
            a%lnods(ispos+4) = floorbottomleftpoint  
         
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floorrightpoint 
            a%lnods(ispos+2) = floorbottomrightpoint   
            a%lnods(ispos+3) = ceilingbottomrightpoint
            a%lnods(ispos+4) = floorbottomleftpoint
         
            !Middle tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint  
            a%lnods(ispos+2) = ceilingtopleftpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint
       
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingbottomleftpoint 
            a%lnods(ispos+2) = ceilingrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint 
            a%lnods(ispos+4) = floortopleftpoint 
       
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint 
            a%lnods(ispos+2) = floorrightpoint
            a%lnods(ispos+3) = floorbottomleftpoint
            a%lnods(ispos+4) = floortopleftpoint   
         
            !Upper tetras 
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingtopleftpoint 
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint
         
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = ceilingrightpoint
            a%lnods(ispos+2) = ceilingtoprightpoint
            a%lnods(ispos+3) = floortoprightpoint 
            a%lnods(ispos+4) = floortopleftpoint  
         
            ielem = ielem+1
            ispos = (ielem-1)*a%pnods(1)
            a%lnods(ispos+1) = floortoprightpoint 
            a%lnods(ispos+2) = floorrightpoint   
            a%lnods(ispos+3) = ceilingrightpoint
            a%lnods(ispos+4) = floortopleftpoint
   
           !Inner groups
           do groupi = 2,TimesAngleIsDivided-1
              floorbottomleftpoint = FirstInnerPoint + groupi - 2
              floorrightpoint = FirstInnerPoint + 2**(irow-1) - 1 + (groupi-1)*2
              floorbottomrightpoint = floorrightpoint - 1
              floortopleftpoint =  floorbottomleftpoint + 1
              floortoprightpoint = floorrightpoint + 1
              ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
              ceilingrightpoint = floorrightpoint + BlockUp 
              ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
              ceilingtopleftpoint = floortopleftpoint + BlockUp
              ceilingtoprightpoint = floortoprightpoint + BlockUp 
    
            
               !Lower tetras 
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingbottomleftpoint 
               a%lnods(ispos+2) = ceilingrightpoint
               a%lnods(ispos+3) = ceilingbottomrightpoint
               a%lnods(ispos+4) = floorbottomleftpoint
            
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingbottomrightpoint
               a%lnods(ispos+2) = ceilingrightpoint
               a%lnods(ispos+3) = floorrightpoint 
               a%lnods(ispos+4) = floorbottomleftpoint  
            
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = floorrightpoint 
               a%lnods(ispos+2) = floorbottomrightpoint   
               a%lnods(ispos+3) = ceilingbottomrightpoint
               a%lnods(ispos+4) = floorbottomleftpoint
            
               !Middle tetras 
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingbottomleftpoint  
               a%lnods(ispos+2) = ceilingtopleftpoint
               a%lnods(ispos+3) = ceilingrightpoint
               a%lnods(ispos+4) = floortopleftpoint
          
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingbottomleftpoint 
               a%lnods(ispos+2) = ceilingrightpoint
               a%lnods(ispos+3) = floorbottomleftpoint 
               a%lnods(ispos+4) = floortopleftpoint 
          
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingrightpoint 
               a%lnods(ispos+2) = floorrightpoint
               a%lnods(ispos+3) = floorbottomleftpoint
               a%lnods(ispos+4) = floortopleftpoint   
            
               !Upper tetras 
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingtopleftpoint 
               a%lnods(ispos+2) = ceilingtoprightpoint
               a%lnods(ispos+3) = ceilingrightpoint
               a%lnods(ispos+4) = floortopleftpoint
            
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = ceilingrightpoint
               a%lnods(ispos+2) = ceilingtoprightpoint
               a%lnods(ispos+3) = floortoprightpoint 
               a%lnods(ispos+4) = floortopleftpoint  
            
               ielem = ielem+1
               ispos = (ielem-1)*a%pnods(1)
               a%lnods(ispos+1) = floortoprightpoint 
               a%lnods(ispos+2) = floorrightpoint   
               a%lnods(ispos+3) = ceilingrightpoint
               a%lnods(ispos+4) = floortopleftpoint
   
           enddo
           !Upper group  
           floorbottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
           floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
           floorbottomrightpoint = floorrightpoint - 1
           ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
           ceilingrightpoint = floorrightpoint + BlockUp 
           ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
   
           !Lower tetras 
           ielem = ielem+1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = ceilingbottomleftpoint 
           a%lnods(ispos+2) = ceilingrightpoint
           a%lnods(ispos+3) = ceilingbottomrightpoint
           a%lnods(ispos+4) = floorbottomleftpoint
           
           ielem = ielem+1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = ceilingbottomrightpoint
           a%lnods(ispos+2) = ceilingrightpoint
           a%lnods(ispos+3) = floorrightpoint 
           a%lnods(ispos+4) = floorbottomleftpoint  
           
           ielem = ielem+1
           ispos = (ielem-1)*a%pnods(1)
           a%lnods(ispos+1) = floorrightpoint 
           a%lnods(ispos+2) = floorbottomrightpoint   
           a%lnods(ispos+3) = ceilingbottomrightpoint
           a%lnods(ispos+4) = floorbottomleftpoint
           
         enddo
      enddo

      !This part creates the zero radius local elements belonging to the last pc
      if (kfl_InternalRadiusAtZero) then
        if (a%MPIrank == a%MPIsize-1) then 
          do kflr = 1,PointsInLength-1
             basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
             do jcol=1,DivisionsInAngle-1
                ZeroPointAngle =  basepointflr + (jcol-1)*LocalAngleDivisionNpoin
                floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                floorbottomrightpoint = ZeroPointAngle + 1
                floortoprightpoint = basepointflr + jcol*LocalAngleDivisionNpoin + 1
                ceilingbottomleftpoint = floorbottomleftpoint + 1
                ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
                ceilingtoprightpoint = floortoprightpoint + BlockUp
   
                !Tetras 

                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingbottomrightpoint
                a%lnods(ispos+4) = floorbottomleftpoint
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floorbottomleftpoint  
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorbottomrightpoint   
                a%lnods(ispos+3) = ceilingbottomrightpoint
                a%lnods(ispos+4) = floorbottomleftpoint


             enddo
           enddo
        endif
      endif

      !Ghost points and elements 
        
      !First the upper and lower processors element connections (required for r=0)
      if (a%MPIsize == 1) then
           !Special case of serial run (Local elements self-connection at teta=0) 
           do kflr = 1,PointsInLength-1
              basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
              ZeroPointAngle = basepointflr + (DivisionsInAngle-1)*LocalAngleDivisionNpoin
              !First radial group
              irow = 1
              floorbottomleftpoint = ZeroPointAngle + irow 
              floorrightpoint = floorbottomleftpoint + PointsInRadius
              floortopleftpoint = basepointflr + irow
              floortoprightpoint = floortopleftpoint + 1
              ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
              ceilingrightpoint = floorrightpoint + BlockUp 
              ceilingtopleftpoint = floortopleftpoint + BlockUp
              ceilingtoprightpoint = floortoprightpoint + BlockUp 
        
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint  
              a%lnods(ispos+2) = ceilingtopleftpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint 
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint 
              a%lnods(ispos+4) = floortopleftpoint 
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint 
              a%lnods(ispos+2) = floorrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint
              a%lnods(ispos+4) = floortopleftpoint   
    
              !Upper tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingtopleftpoint 
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = floortoprightpoint 
              a%lnods(ispos+4) = floortopleftpoint  
    
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = floortoprightpoint 
              a%lnods(ispos+2) = floorrightpoint   
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
      
              !Second radial groups
              irow = 2
              FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
              !Second radial Upper group  
              floorbottomleftpoint = FirstInnerPoint
              floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
              floortopleftpoint =  basepointflr + irow
              floortoprightpoint = floortopleftpoint + 1
              ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
              ceilingrightpoint = floorrightpoint + BlockUp 
              ceilingtopleftpoint = floortopleftpoint + BlockUp
              ceilingtoprightpoint = floortoprightpoint + BlockUp 
            
              !Middle tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint  
              a%lnods(ispos+2) = ceilingtopleftpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
          
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingbottomleftpoint 
              a%lnods(ispos+2) = ceilingrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint 
              a%lnods(ispos+4) = floortopleftpoint 
          
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint 
              a%lnods(ispos+2) = floorrightpoint
              a%lnods(ispos+3) = floorbottomleftpoint
              a%lnods(ispos+4) = floortopleftpoint   
  
              !Upper tetras 
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingtopleftpoint 
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
   
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = ceilingrightpoint
              a%lnods(ispos+2) = ceilingtoprightpoint
              a%lnods(ispos+3) = floortoprightpoint 
              a%lnods(ispos+4) = floortopleftpoint  
   
              ielem = ielem+1
              ispos = (ielem-1)*a%pnods(1)
              a%lnods(ispos+1) = floortoprightpoint 
              a%lnods(ispos+2) = floorrightpoint   
              a%lnods(ispos+3) = ceilingrightpoint
              a%lnods(ispos+4) = floortopleftpoint
      
              !Other radial groups
              do irow=3,PointsInRadius-1
                TimesAngleIsDivided = 2_ip**(irow-1)
                FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
                !Upper group  
                floorbottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
                floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
                floorbottomrightpoint = rightpoint - 1
                floortopleftpoint = basepointflr + irow
                floortoprightpoint = floortopleftpoint + 1
                ceilingbottomleftpoint = floorbottomleftpoint + BlockUp 
                ceilingrightpoint = floorrightpoint + BlockUp 
                ceilingtopleftpoint = floortopleftpoint + BlockUp
                ceilingtoprightpoint = floortoprightpoint + BlockUp 

                 !Middle tetras 
                 ielem = ielem+1
                 ispos = (ielem-1)*a%pnods(1)
                 a%lnods(ispos+1) = ceilingbottomleftpoint  
                 a%lnods(ispos+2) = ceilingtopleftpoint
                 a%lnods(ispos+3) = ceilingrightpoint
                 a%lnods(ispos+4) = floortopleftpoint
             
                 ielem = ielem+1
                 ispos = (ielem-1)*a%pnods(1)
                 a%lnods(ispos+1) = ceilingbottomleftpoint 
                 a%lnods(ispos+2) = ceilingrightpoint
                 a%lnods(ispos+3) = floorbottomleftpoint 
                 a%lnods(ispos+4) = floortopleftpoint 
             
                 ielem = ielem+1
                 ispos = (ielem-1)*a%pnods(1)
                 a%lnods(ispos+1) = ceilingrightpoint 
                 a%lnods(ispos+2) = floorrightpoint
                 a%lnods(ispos+3) = floorbottomleftpoint
                 a%lnods(ispos+4) = floortopleftpoint   
      
                 !Upper tetras 
                 ielem = ielem+1
                 ispos = (ielem-1)*a%pnods(1)
                 a%lnods(ispos+1) = ceilingtopleftpoint 
                 a%lnods(ispos+2) = ceilingtoprightpoint
                 a%lnods(ispos+3) = ceilingrightpoint
                 a%lnods(ispos+4) = floortopleftpoint
      
                 ielem = ielem+1
                 ispos = (ielem-1)*a%pnods(1)
                 a%lnods(ispos+1) = ceilingrightpoint
                 a%lnods(ispos+2) = ceilingtoprightpoint
                 a%lnods(ispos+3) = floortoprightpoint 
                 a%lnods(ispos+4) = floortopleftpoint  
      
                 ielem = ielem+1
                 ispos = (ielem-1)*a%pnods(1)
                 a%lnods(ispos+1) = floortoprightpoint 
                 a%lnods(ispos+2) = floorrightpoint   
                 a%lnods(ispos+3) = ceilingrightpoint
                 a%lnods(ispos+4) = floortopleftpoint
               
              enddo
           enddo
      else
            !Parallel run

            !Upper processor conection

            !GhostPoints Global numbering
            if (a%MPIrank < a%MPIsize-1) then
               baseGhost0  = a%poin0 + a%npoinLocal
               do kflr = 1,PointsInLength
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                  do ipoin = 1,PointsInRadius
                     LocalToGlobal(UpGhost0 + (kflr-1)*PointsInRadius + ipoin) = baseGhost0 + basepointflr + ipoin
                     ProcessorList(UpGhost0 + (kflr-1)*PointsInRadius + ipoin) = a%MPIrank + 1
                  enddo
               enddo
            else 
               !Special case of being the last processor 
               baseGhost0 = 0
               do kflr = 1,PointsInLength
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                  do ipoin = 1,PointsInRadius
                     LocalToGlobal(UpGhost0 + (kflr-1)*PointsInRadius + ipoin) = baseGhost0 + basepointflr + ipoin
                     ProcessorList(UpGhost0 + (kflr-1)*PointsInRadius + ipoin) = 1
                  enddo
               enddo
            endif

            do kflr = 1,PointsInLength-1
               basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
               ZeroPointAngle = basepointflr + (DivisionsInAngle-1)*LocalAngleDivisionNpoin

               !First radial group
               irow = 1
               floorbottomleftpoint = ZeroPointAngle + irow 
               floorrightpoint = floorbottomleftpoint + PointsInRadius
               floortopleftpoint = UpGhost0 + (kflr-1)*PointsInRadius + irow
               floortoprightpoint = floortopleftpoint + 1
               ceilingbottomleftpoint = floorbottomleftpoint + BlockUp
               ceilingrightpoint = floorrightpoint + BlockUp
               ceilingtopleftpoint = floortopleftpoint + PointsInRadius
               ceilingtoprightpoint = floortoprightpoint + PointsInRadius

               !Middle tetras 
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 2)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
      
               !Upper tetras 
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 3)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 4)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 5)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint

               !Second radial group
               irow = 2
               FirstInnerPoint = ZeroPointAngle + PointsInRadius + 2**(irow-2)
   
               !Second radial Upper group  
               floorbottomleftpoint = FirstInnerPoint
               floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
               floortopleftpoint =  UpGhost0 + (kflr-1)*PointsInRadius + irow 
               floortoprightpoint = floortopleftpoint + 1
               ceilingbottomleftpoint = floorbottomleftpoint + BlockUp
               ceilingrightpoint = floorrightpoint + BlockUp
               ceilingtopleftpoint = floortopleftpoint + PointsInRadius
               ceilingtoprightpoint = floortoprightpoint + PointsInRadius

               !Middle tetras 
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 2)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
      
               !Upper tetras 
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 3)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 4)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 5)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint

               !Other radial groups
               do irow=3,PointsInRadius-1
                 TimesAngleIsDivided = 2_ip**(irow-1)
                 FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
                 !Upper group  
                 floorbottomleftpoint = FirstInnerPoint + 2**(irow-1) - 2
                 floorrightpoint = FirstInnerPoint + 2**(irow-1) + 2**(irow) - 3
                 floortopleftpoint = UpGhost0 + (kflr-1)*PointsInRadius + irow;
                 floortoprightpoint = floortopleftpoint + 1;
                 ceilingbottomleftpoint = floorbottomleftpoint + BlockUp
                 ceilingrightpoint = floorrightpoint + BlockUp
                 ceilingtopleftpoint = floortopleftpoint + PointsInRadius
                 ceilingtoprightpoint = floortoprightpoint + PointsInRadius
         
               !Middle tetras 
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 2)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
      
               !Upper tetras 
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 3)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 4)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
      
               ispos = (UpElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 5)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
                 
               enddo
            enddo

            !Lower processor conection
   
            !GhostPoints Global numbering
            if (a%MPIrank > 0) then
               baseGhost0  = a%poin0 - LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + LocalAngleDivisionNpoin*(DivisionsInAngle-1)
               do kflr = 1,PointsInLength
                  do ipoin = 1,PointsInRadius
                     ProcessorList(DownGhost0 + (kflr-1)*PointsInRadius + ipoin) = a%MPIrank - 1
                  enddo
               enddo
            else 
               !Special case of being the first processor 
               baseGhost0 = a%npoinLocal*(a%MPIsize-1) + LocalAngleDivisionNpoin*(DivisionsInAngle-1)
               do kflr = 1,PointsInLength
                  do ipoin = 1,PointsInRadius
                     ProcessorList(DownGhost0 + (kflr-1)*PointsInRadius + ipoin) = a%MPIsize - 1
                  enddo
               enddo
            endif

            do kflr = 1,PointsInLength-1
               basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
   
               !First radial group
               irow = 1
               LocalToGlobal(DownGhost0 + (kflr-1)*PointsInRadius + irow) = baseGhost0 + (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle + irow
               floorbottomleftpoint = DownGhost0 + (kflr-1)*PointsInRadius + irow 
               floorrightpoint = floorbottomleftpoint + 1
               floortopleftpoint = basepointflr + irow
               floortoprightpoint = floortopleftpoint + 1
               ceilingbottomleftpoint = floorbottomleftpoint + PointsInRadius 
               ceilingrightpoint = ceilingbottomleftpoint + 1
               ceilingtopleftpoint = floortopleftpoint + BlockUp
               ceilingtoprightpoint = ceilingtopleftpoint + 1

               !Middle tetras 
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 2)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
 
      
               !Upper tetras 
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 3)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 4)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 5)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint

               !Second radial groups
               irow = 2
               FirstInnerPoint = baseGhost0 + (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle + PointsInRadius + 2**(irow-2)
               LocalToGlobal(DownGhost0 + (kflr-1)*PointsInRadius + irow) = FirstInnerPoint
           
               !Second radial Upper group  
               floorbottomleftpoint = DownGhost0 + (kflr-1)*PointsInRadius + irow 
               floorrightpoint = floorbottomleftpoint + 1
               floortopleftpoint = basepointflr + irow
               floortoprightpoint = floortopleftpoint + 1
               ceilingbottomleftpoint = floorbottomleftpoint + PointsInRadius 
               ceilingrightpoint = ceilingbottomleftpoint + 1
               ceilingtopleftpoint = floortopleftpoint + BlockUp
               ceilingtoprightpoint = ceilingtopleftpoint + 1

               !Middle tetras 
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 2)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
 
      
               !Upper tetras 
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 3)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 4)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 5)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint

               !Other radial groups
               do irow=3,PointsInRadius-1
                 TimesAngleIsDivided = 2_ip**(irow-1)
                 FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
                 LocalToGlobal(DownGhost0 + (kflr-1)*PointsInRadius + irow ) = FirstInnerPoint + 2**(irow-1) - 2
                 !Upper group  
                 floorbottomleftpoint = DownGhost0 + (kflr-1)*PointsInRadius + irow 
                 floorrightpoint = floorbottomleftpoint + 1
                 floortopleftpoint = basepointflr + irow
                 floortoprightpoint = floortopleftpoint + 1
                 ceilingbottomleftpoint = floorbottomleftpoint + PointsInRadius 
                 ceilingrightpoint = ceilingbottomleftpoint + 1
                 ceilingtopleftpoint = floortopleftpoint + BlockUp
                 ceilingtoprightpoint = ceilingtopleftpoint + 1
           
               !Middle tetras 
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint  
                a%lnods(ispos+2) = ceilingtopleftpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint 
                a%lnods(ispos+4) = floortopleftpoint 
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 2)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint 
                a%lnods(ispos+2) = floorrightpoint
                a%lnods(ispos+3) = floorbottomleftpoint
                a%lnods(ispos+4) = floortopleftpoint   
 
      
               !Upper tetras 
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 3)*a%pnods(1)
                a%lnods(ispos+1) = ceilingtopleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 4)*a%pnods(1)
                a%lnods(ispos+1) = ceilingrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floortopleftpoint  
      
               ispos = (DownElem0 + (kflr-1)*(PointsInRadius-1)*6 + (irow-1)*6 + 5)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorrightpoint   
                a%lnods(ispos+3) = ceilingrightpoint
                a%lnods(ispos+4) = floortopleftpoint
                 
               enddo
               !Care for defining last radius ghosts
               TimesAngleIsDivided = 2_ip**(PointsInRadius-1)
               FirstInnerPoint = FirstInnerPoint + 2**(PointsInRadius-2) - 1
               LocalToGlobal(DownGhost0 + kflr*PointsInRadius) = FirstInnerPoint + 2**(PointsInRadius-1) - 2

            enddo
            
            !Care for defining last floor ghosts
            !First radial group
            irow = 1
            LocalToGlobal(DownGhost0 + (PointsInLength-1)*PointsInRadius + irow) = baseGhost0 + (PointsInLength-1)*LocalAngleDivisionNpoin*DivisionsInAngle + irow

            !Second radial groups
            irow = 2
            FirstInnerPoint = baseGhost0 + (PointsInLength-1)*LocalAngleDivisionNpoin*DivisionsInAngle + PointsInRadius + 2**(irow-2)
            LocalToGlobal(DownGhost0 + (PointsInLength-1)*PointsInRadius + irow) = FirstInnerPoint
           
            !Other radial groups
            do irow=3,PointsInRadius
              TimesAngleIsDivided = 2_ip**(irow-1)
              FirstInnerPoint = FirstInnerPoint + 2**(irow-2) - 1
              LocalToGlobal(DownGhost0 + (PointsInLength-1)*PointsInRadius + irow ) = FirstInnerPoint + 2**(irow-1) - 2
            enddo

      endif

      !Then the ghost points of elements at radius equal zero
      if (kfl_InternalRadiusAtZero) then
         if (a%MPIsize == 1) then
            !Serial run: Inside the zero radius local elements self-connection
            !Special case of serial run self-connection, anyway inside this logic 
            do kflr = 1,PointsInLength-1
               basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
               ZeroPointAngle = basepointflr + (DivisionsInAngle-1)*LocalAngleDivisionNpoin
               !First radial group
               irow = 1
               floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
               floorbottomrightpoint = ZeroPointAngle + 1
               floortoprightpoint = basepointflr + 1
               ceilingbottomleftpoint = floorbottomleftpoint + 1
               ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
               ceilingtoprightpoint = floortoprightpoint + BlockUp

               !Tetras 
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomleftpoint 
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = ceilingbottomrightpoint
                a%lnods(ispos+4) = floorbottomleftpoint
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = ceilingbottomrightpoint
                a%lnods(ispos+2) = ceilingtoprightpoint
                a%lnods(ispos+3) = floortoprightpoint 
                a%lnods(ispos+4) = floorbottomleftpoint  
    
                ielem = ielem+1
                ispos = (ielem-1)*a%pnods(1)
                a%lnods(ispos+1) = floortoprightpoint 
                a%lnods(ispos+2) = floorbottomrightpoint   
                a%lnods(ispos+3) = ceilingbottomrightpoint
                a%lnods(ispos+4) = floorbottomleftpoint

            enddo
         else
            !Parallel run !Therefore r_0 = 0 (inner, upper and lower) element connections between processors
            !First the inner elements conection between processors 
            if (a%MPIrank /= a%MPIsize-1) then
               !Not last pc: Inside the inner elements parallel pc connection to last pc 
               !GhostPoints Global numbering
               do kflr = 1,PointsInLength-1
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle

                  LocalToGlobal(InnerGhost0 + kflr) = a%npoinLocal*a%MPIsize + kflr
                  ProcessorList(InnerGhost0 + kflr) = a%MPIsize-1
                  
                  do jcol=1,DivisionsInAngle-1
                     ZeroPointAngle = basepointflr + (jcol-1)*LocalAngleDivisionNpoin
                     !First radial group
                     irow = 1
                     floorbottomleftpoint = InnerGhost0 + kflr 
                     floorbottomrightpoint = ZeroPointAngle + 1
                     floortoprightpoint = basepointflr + jcol*LocalAngleDivisionNpoin + 1
                     ceilingbottomleftpoint = floorbottomleftpoint + 1 
                     ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
                     ceilingtoprightpoint = floortoprightpoint + BlockUp
         
                     !Tetras 
                     ielem = InnerElem0 + (kflr-1)*3_ip*(DivisionsInAngle-1) + (jcol-1)*3_ip + 1
                     ispos = (ielem-1)*a%pnods(1)
                     a%lnods(ispos+1) = ceilingbottomleftpoint 
                     a%lnods(ispos+2) = ceilingtoprightpoint
                     a%lnods(ispos+3) = ceilingbottomrightpoint
                     a%lnods(ispos+4) = floorbottomleftpoint
                     
                     ielem = InnerElem0 + (kflr-1)*3_ip*(DivisionsInAngle-1) + (jcol-1)*3_ip + 2
                     ispos = (ielem-1)*a%pnods(1)
                     a%lnods(ispos+1) = ceilingbottomrightpoint
                     a%lnods(ispos+2) = ceilingtoprightpoint
                     a%lnods(ispos+3) = floortoprightpoint 
                     a%lnods(ispos+4) = floorbottomleftpoint  
                     
                     ielem = InnerElem0 + (kflr-1)*3_ip*(DivisionsInAngle-1) + (jcol-1)*3_ip + 3
                     ispos = (ielem-1)*a%pnods(1)
                     a%lnods(ispos+1) = floortoprightpoint 
                     a%lnods(ispos+2) = floorbottomrightpoint   
                     a%lnods(ispos+3) = ceilingbottomrightpoint
                     a%lnods(ispos+4) = floorbottomleftpoint
                  enddo

               enddo
               LocalToGlobal(InnerGhost0 + PointsInLength) = a%npoinLocal*a%MPIsize + PointsInLength
               ProcessorList(InnerGhost0 + PointsInLength) = a%MPIsize-1

               !Not last pc: Inside the inner element parallel connection with upper pc and last pc 
               do kflr = 1,PointsInLength-1
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                  ZeroPointAngle = basepointflr + (DivisionsInAngle-1)*LocalAngleDivisionNpoin
                  !First radial group
                  irow = 1
                  floorbottomleftpoint = InnerGhost0 + kflr 
                  floorbottomrightpoint = ZeroPointAngle + 1
                  floortoprightpoint = UpGhost0 + (kflr-1)*PointsInRadius + 1 
                  ceilingbottomleftpoint = floorbottomleftpoint + 1 
                  ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
                  ceilingtoprightpoint = floortoprightpoint + PointsInRadius
         
                  ielem = UpInnerElem0 + (kflr-1)*3_ip + 1
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomleftpoint 
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint

                  ielem = UpInnerElem0 + (kflr-1)*3_ip + 2
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomrightpoint
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = floortoprightpoint 
                  a%lnods(ispos+4) = floorbottomleftpoint  

                  ielem = UpInnerElem0 + (kflr-1)*3_ip + 3
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = floortoprightpoint 
                  a%lnods(ispos+2) = floorbottomrightpoint   
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint
               enddo

               !Not last pc: Inside the inner element parallel connection with lower pc and last pc 
               do kflr = 1,PointsInLength-1
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                  !First radial group
                  irow = 1
                  floorbottomleftpoint = InnerGhost0 + kflr 
                  floorbottomrightpoint = DownGhost0 + (kflr-1)*PointsInRadius + 1
                  floortoprightpoint = basepointflr + 1
                  ceilingbottomleftpoint = floorbottomleftpoint + 1 
                  ceilingbottomrightpoint = floorbottomrightpoint + PointsInRadius
                  ceilingtoprightpoint = floortoprightpoint + BlockUp
         
                  ielem = DownInnerElem0 + (kflr-1)*3_ip + 1
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomleftpoint 
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint

                  ielem = DownInnerElem0 + (kflr-1)*3_ip + 2
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomrightpoint
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = floortoprightpoint 
                  a%lnods(ispos+4) = floorbottomleftpoint  

                  ielem = DownInnerElem0 + (kflr-1)*3_ip + 3
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = floortoprightpoint 
                  a%lnods(ispos+2) = floorbottomrightpoint   
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint
               enddo
            else
               !Last pc: Inside the inner element connection to inner ghosts in other pcs
               !OtherPcsGhost0  
               !Last pc: Connection with all inner elements of other pcs
               do groupi = 1,a%MPIsize-1
                  baseGhost0 = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength*(groupi-1)
                  !All inner elements 
                  !care for the lower group that is connected with the lower to the inner other pc 
                  if (groupi == 1) then
                     do kflr = 1,PointsInLength-1
                        floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                        floorbottomrightpoint = UpGhost0 + (kflr-1)*PointsInRadius + 1   
                        floortoprightpoint = OtherPcsGhost0 + kflr
                        ceilingbottomleftpoint = floorbottomleftpoint + 1 
                        ceilingbottomrightpoint = floorbottomrightpoint + PointsInRadius
                        ceilingtoprightpoint = floortoprightpoint + 1
   
                        ielem = OtherInnerElem0 + (kflr-1)*3 + 1
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = ceilingbottomleftpoint 
                        a%lnods(ispos+2) = ceilingtoprightpoint
                        a%lnods(ispos+3) = ceilingbottomrightpoint
                        a%lnods(ispos+4) = floorbottomleftpoint
      
                        ielem = OtherInnerElem0 + (kflr-1)*3 + 2
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = ceilingbottomrightpoint
                        a%lnods(ispos+2) = ceilingtoprightpoint
                        a%lnods(ispos+3) = floortoprightpoint 
                        a%lnods(ispos+4) = floorbottomleftpoint  
      
                        ielem = OtherInnerElem0 + (kflr-1)*3 + 3
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = floortoprightpoint 
                        a%lnods(ispos+2) = floorbottomrightpoint   
                        a%lnods(ispos+3) = ceilingbottomrightpoint
                        a%lnods(ispos+4) = floorbottomleftpoint
                     enddo
                  else
                     do kflr = 1,PointsInLength-1
                        basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle

                        LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength - PointsInLength + kflr) = baseGhost0 + basepointflr + 1
                        ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength - PointsInLength + kflr) = groupi
      
                        floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                        floorbottomrightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength - PointsInLength + kflr 
                        floortoprightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + kflr
                        ceilingbottomleftpoint = floorbottomleftpoint + 1 
                        ceilingbottomrightpoint = floorbottomrightpoint + 1
                        ceilingtoprightpoint = floortoprightpoint + 1
      
                        ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (kflr-1)*3 + 1
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = ceilingbottomleftpoint 
                        a%lnods(ispos+2) = ceilingtoprightpoint
                        a%lnods(ispos+3) = ceilingbottomrightpoint
                        a%lnods(ispos+4) = floorbottomleftpoint
      
                        ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (kflr-1)*3 + 2
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = ceilingbottomrightpoint
                        a%lnods(ispos+2) = ceilingtoprightpoint
                        a%lnods(ispos+3) = floortoprightpoint 
                        a%lnods(ispos+4) = floorbottomleftpoint  
      
                        ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (kflr-1)*3 + 3
                        ispos = (ielem-1)*a%pnods(1)
                        a%lnods(ispos+1) = floortoprightpoint 
                        a%lnods(ispos+2) = floorbottomrightpoint   
                        a%lnods(ispos+3) = ceilingbottomrightpoint
                        a%lnods(ispos+4) = floorbottomleftpoint
                     enddo
                     basepointflr = (PointsInLength-1)*LocalAngleDivisionNpoin*DivisionsInAngle

                     LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength) = baseGhost0 + basepointflr + 1
                     ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength) = groupi
                  endif
                  !All other inner elements 
                  do jcol=2,DivisionsInAngle-1
                     ZeroPointAngle = (jcol-1)*LocalAngleDivisionNpoin
                     if ((groupi == a%MPIsize-1) .and. (jcol == DivisionsInAngle-1)) then
                        do kflr = 1,PointsInLength-1
                           basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                           LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr) = baseGhost0 + basepointflr + ZeroPointAngle + 1
                           ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr) = groupi
                           floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                           floorbottomrightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr 
                           floortoprightpoint = DownGhost0 + (kflr-1)*PointsInRadius + 1
                           ceilingbottomleftpoint = floorbottomleftpoint + 1 
                           ceilingbottomrightpoint = floorbottomrightpoint + 1
                           ceilingtoprightpoint = floortoprightpoint + PointsInRadius
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 1
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = ceilingbottomleftpoint 
                           a%lnods(ispos+2) = ceilingtoprightpoint
                           a%lnods(ispos+3) = ceilingbottomrightpoint
                           a%lnods(ispos+4) = floorbottomleftpoint
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 2
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = ceilingbottomrightpoint
                           a%lnods(ispos+2) = ceilingtoprightpoint
                           a%lnods(ispos+3) = floortoprightpoint 
                           a%lnods(ispos+4) = floorbottomleftpoint  
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 3
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = floortoprightpoint 
                           a%lnods(ispos+2) = floorbottomrightpoint   
                           a%lnods(ispos+3) = ceilingbottomrightpoint
                           a%lnods(ispos+4) = floorbottomleftpoint
                        enddo
                        basepointflr = (PointsInLength-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                        LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength) = baseGhost0 + basepointflr + ZeroPointAngle + 1
                        ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength) = groupi

                     else
                        do kflr = 1,PointsInLength-1
                           basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                           LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr) = baseGhost0 + basepointflr + ZeroPointAngle + 1
                           ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr) = groupi
                           floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                           floorbottomrightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr 
                           floortoprightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength + kflr
                           ceilingbottomleftpoint = floorbottomleftpoint + 1 
                           ceilingbottomrightpoint = floorbottomrightpoint + 1
                           ceilingtoprightpoint = floortoprightpoint + 1

                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 1
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = ceilingbottomleftpoint 
                           a%lnods(ispos+2) = ceilingtoprightpoint
                           a%lnods(ispos+3) = ceilingbottomrightpoint
                           a%lnods(ispos+4) = floorbottomleftpoint
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 2
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = ceilingbottomrightpoint
                           a%lnods(ispos+2) = ceilingtoprightpoint
                           a%lnods(ispos+3) = floortoprightpoint 
                           a%lnods(ispos+4) = floorbottomleftpoint  
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 3
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = floortoprightpoint 
                           a%lnods(ispos+2) = floorbottomrightpoint   
                           a%lnods(ispos+3) = ceilingbottomrightpoint
                           a%lnods(ispos+4) = floorbottomleftpoint
                        enddo
                        basepointflr = (PointsInLength-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                        LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength) = baseGhost0 + basepointflr + ZeroPointAngle + 1
                        ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength) = groupi
                     endif
                  enddo
                  !care for the top group that is connected with the upper to the inner other pc 
                  ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
                  if (groupi /= a%MPIsize-1) then
                        jcol = DivisionsInAngle
                        do kflr = 1,PointsInLength-1
                           basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                           LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr) = baseGhost0 + basepointflr + ZeroPointAngle + 1
                           ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr) = groupi
                           floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                           floorbottomrightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 2)*PointsInLength + kflr 
                           floortoprightpoint = OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength + kflr
                           ceilingbottomleftpoint = floorbottomleftpoint + 1 
                           ceilingbottomrightpoint = floorbottomrightpoint + 1
                           ceilingtoprightpoint = floortoprightpoint + 1

                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 1
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = ceilingbottomleftpoint 
                           a%lnods(ispos+2) = ceilingtoprightpoint
                           a%lnods(ispos+3) = ceilingbottomrightpoint
                           a%lnods(ispos+4) = floorbottomleftpoint
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 2
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = ceilingbottomrightpoint
                           a%lnods(ispos+2) = ceilingtoprightpoint
                           a%lnods(ispos+3) = floortoprightpoint 
                           a%lnods(ispos+4) = floorbottomleftpoint  
         
                           ielem = OtherInnerElem0 + (groupi-1)*DivisionsInAngle*(PointsInLength-1)*3 + (jcol - 1)*(PointsInLength-1)*3 + (kflr-1)*3 + 3
                           ispos = (ielem-1)*a%pnods(1)
                           a%lnods(ispos+1) = floortoprightpoint 
                           a%lnods(ispos+2) = floorbottomrightpoint   
                           a%lnods(ispos+3) = ceilingbottomrightpoint
                           a%lnods(ispos+4) = floorbottomleftpoint
                        enddo
                        basepointflr = (PointsInLength-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                        LocalToGlobal(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength) = baseGhost0 + basepointflr + ZeroPointAngle + 1
                        ProcessorList(OtherPcsGhost0 + (groupi-1)*DivisionsInAngle*PointsInLength + (jcol - 1)*PointsInLength) = groupi
                  endif
               enddo
               !Im last pc: Inner-Upper processor conection
               !I use UpGhost0 and InnerLocalpoint
               ZeroPointAngle = (DivisionsInAngle-1)*LocalAngleDivisionNpoin
               do kflr = 1,PointsInLength-1
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                  floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                  floorbottomrightpoint = basepointflr + ZeroPointAngle + 1 
                  floortoprightpoint = UpGhost0 + (kflr-1)*PointsInRadius + 1   
                  ceilingbottomleftpoint = floorbottomleftpoint + 1 
                  ceilingbottomrightpoint = floorbottomrightpoint + BlockUp
                  ceilingtoprightpoint = floortoprightpoint + PointsInRadius

                  ielem = UpInnerElem0 + (kflr-1)*3 + 1
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomleftpoint 
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint
         
                  ielem = UpInnerElem0 + (kflr-1)*3 + 2
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomrightpoint
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = floortoprightpoint 
                  a%lnods(ispos+4) = floorbottomleftpoint  
         
                  ielem = UpInnerElem0 + (kflr-1)*3 + 3
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = floortoprightpoint 
                  a%lnods(ispos+2) = floorbottomrightpoint   
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint
               enddo
               !Im last pc: Inner-Lower processor conection
               !I use DownGhost0 and InnerLocalpoint
               do kflr = 1,PointsInLength-1
                  basepointflr = (kflr-1)*LocalAngleDivisionNpoin*DivisionsInAngle
                  floorbottomleftpoint = LocalAngleDivisionNpoin*DivisionsInAngle*PointsInLength + kflr
                  floorbottomrightpoint = DownGhost0 + (kflr-1)*PointsInRadius + 1    
                  floortoprightpoint = basepointflr + 1 
                  ceilingbottomleftpoint = floorbottomleftpoint + 1 
                  ceilingbottomrightpoint = floorbottomrightpoint + PointsInRadius
                  ceilingtoprightpoint = floortoprightpoint + BlockUp

                  ielem = DownInnerElem0 + (kflr-1)*3 + 1
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomleftpoint 
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint
         
                  ielem = DownInnerElem0 + (kflr-1)*3 + 2
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = ceilingbottomrightpoint
                  a%lnods(ispos+2) = ceilingtoprightpoint
                  a%lnods(ispos+3) = floortoprightpoint 
                  a%lnods(ispos+4) = floorbottomleftpoint  
         
                  ielem = DownInnerElem0 + (kflr-1)*3 + 3
                  ispos = (ielem-1)*a%pnods(1)
                  a%lnods(ispos+1) = floortoprightpoint 
                  a%lnods(ispos+2) = floorbottomrightpoint   
                  a%lnods(ispos+3) = ceilingbottomrightpoint
                  a%lnods(ispos+4) = floorbottomleftpoint
               enddo
            endif
         endif
      endif

      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      call a%Memor%dealloc(a%npoin,LocalToGlobal,'LocalToGlobal','ManufacturedMesh')
      call a%Memor%dealloc(a%npoin,ProcessorList,'ProcessorList','ManufacturedMesh')
      
      !External normal, we simply allocate it
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','ManufacturedMesh')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ManufacturedMesh')

   endif

end subroutine

subroutine ManufacturedMeshb(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   !List of bodies, we simply allocate it
   call a%Memor%alloc(a%nboun,a%lbody,'lbody','lbody') 




end subroutine
