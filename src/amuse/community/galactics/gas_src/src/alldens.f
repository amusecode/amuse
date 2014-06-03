! density except halo df component
      function alldensnohalo(r,z) result(totaldens)      
       real totaldens
      common /flags/ idiskflag, ibulgeflag, ihaloflag
      totaldens=0
      if(idiskflag.EQ.2.OR.idiskflag.EQ.3)totaldens=totaldens+
     &  thickdiskdens(r,z)
      if(idiskflag.EQ.1.OR.idiskflag.EQ.3)totaldens=totaldens+
     &  tdskdens(r,z)   
      if(ibulgeflag.EQ.1) then
       print*,'dont know how to get bulge density!'        
       stop
      endif 

      if(ihaloflag.EQ.2) totaldens=totaldens+halofixeddens(r,z)     
      if(ibulgeflag.EQ.2) totaldens=totaldens+bulgefixeddens(r,z)
      if(ibulgeflag.EQ.3) totaldens=totaldens+bulgedfdens(r,z)
      end function

! density except bulge df component
      function alldensnobulge(r,z) result(totaldens)      
       real totaldens
      common /flags/ idiskflag, ibulgeflag, ihaloflag
      totaldens=0
      if(idiskflag.EQ.2.OR.idiskflag.EQ.3)totaldens=totaldens+
     &  thickdiskdens(r,z)
      if(idiskflag.EQ.1.OR.idiskflag.EQ.3)totaldens=totaldens+
     &  tdskdens(r,z)   
      if(ihaloflag.EQ.1) then
       print*,'dont know how to get halo density!'        
       stop
      endif 
      if(ihaloflag.EQ.2) totaldens=totaldens+halofixeddens(r,z)
      if(ihaloflag.EQ.3) totaldens=totaldens+halodfdens(r,z)
      if(ibulgeflag.EQ.2) totaldens=totaldens+bulgefixeddens(r,z)
      end function


