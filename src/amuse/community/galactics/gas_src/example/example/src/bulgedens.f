
      function bulgedens(r,z)
      common /flags/ idiskflag, ibulgeflag, ihaloflag

      bulgedens=0
      psi = pot(r,z)
      psi0= pot(0.,0.)
      if(ibulgeflag.EQ.1) 
     &       bulgedens = bulgedens+ bulgedenspsi(psi,psi0)
      if(ibulgeflag.EQ.3) 
     &       bulgedens = bulgedens+ bulgedfdenspsi(psi,psi0)  
      return
      end
