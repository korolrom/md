subroutine BackIntegration(dimensionality,nsurfaces,natoms,nbeads,istep)
    use potential
    IMPLICIT NONE
    
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads,istep
    
    real(8)::hopprob(nsurfaces),velocityveti(dimensionality,natoms,nbeads)
    real(8)::positionnmsave(dimensionality,natoms,nbeads),velocitynmsave(dimensionality,natoms,nbeads)
    real(8)::prob,rtemp,sign,evotime
    integer::idim,ibead,jbead,iatom,isur,jsur,flag_hopattempt,regionnew
    complex(8)::comtemp
    real(8)::randomnumber,gaussrannum
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derisurface(idim,iatom,ibead)/mass(iatom)*timestep*0.5d0
    end do
    end do
    end do
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1    
    positionnmsave(idim,iatom,:)=0.0d0
    velocitynmsave(idim,iatom,:)=0.0d0
        do ibead=1,nbeads,1
        do jbead=1,nbeads,1
            positionnmsave(idim,iatom,ibead)=positionnmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*position(idim,iatom,jbead)
            velocitynmsave(idim,iatom,ibead)=velocitynmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*velocity(idim,iatom,jbead)
        end do
        end do
    end do
    end do
    
    evotime=timestep
    call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
    call PotentialRead(dimensionality,nsurfaces,natoms,nbeads,1)
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derisurface(idim,iatom,ibead)/mass(iatom)*timestep*0.5d0
    end do
    end do
    end do
    ! do isur=1,nsurfaces,1
        ! elecoeffobs(isur)=elecoeffevo(isur)*CDexp((0.0d0,-1.0d0)*phase(isur))
    ! end do
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1    
    do ibead=1,nbeads,1
        positionrecord(idim,iatom,ibead,istep-1)=position(idim,iatom,ibead)
        velocityrecord(idim,iatom,ibead,istep-1)=velocity(idim,iatom,ibead)
    end do
    end do
    end do
    
    if (nsurfaces.gt.1) then
        ! flag_hopattempt=0
        
        ! hopprob=0.0d0
        ! do isur=1,nsurfaces,1
            ! do idim=1,dimensionality,1
            ! do ibead=1,nbeads,1
                ! hopprob(isur)=hopprob(isur)+velocity(idim,ibead)*derivativecoupling(isur,surface,idim,ibead)
            ! end do
            ! end do
            ! hopprob(isur)=hopprob(isur)/dble(nbeads)
            ! ! comtemp=elecoeffobs(surface)*Dconjg(elecoeffobs(isur))
            ! ! hopprob(isur)=hopprob(isur)*Dreal(comtemp)*timestep
            ! ! hopprob(isur)=-2.0d0*hopprob(isur)/(Cdabs(elecoeffobs(surface))**2)
            ! ! if (hopprob(isur).lt.0.0d0) then
                ! ! hopprob(isur)=0.0d0
            ! ! end if
            ! hopprob(isur)=Dexp(-Dabs(hopprob(isur))*timestep)
            ! hopprob(isur)=0.5d0*(1.0d0-hopprob(isur))
        ! end do
        ! do isur=nsurfaces,2,-1
            ! do jsur=isur-1,1,-1
                ! hopprob(isur)=hopprob(isur)+hopprob(jsur)
            ! end do
        ! end do
        
        ! prob=randomnumber(randomseed)
        ! do isur=1,nsurfaces,1
            ! if (prob.lt.(hopprob(isur))) then

                ! ! write(*,*) velocity, derivativecoupling, potsurface            
                ! flag_hopattempt=1
                ! if (isur.gt.1) then
                ! weight=weight/(hopprob(isur)-hopprob(isur-1))
                ! ! write(*,*) istep, hopprob(isur)-hopprob(isur-1)
                ! else
                ! weight=weight/hopprob(1)
                ! ! write(*,*) istep, hopprob(1)
                ! end if
                
                ! rtemp=0.0d0
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! rtemp=rtemp+derivativecoupling(isur,surface,idim,ibead)**2
                ! end do
                ! end do
                ! rtemp=1.0d0/Dsqrt(rtemp)
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! derivativecoupling(isur,surface,idim,ibead)=derivativecoupling(isur,surface,idim,ibead)*rtemp
                ! end do
                ! end do
                
                ! rtemp=0.0d0
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! rtemp=rtemp+velocity(idim,ibead)*derivativecoupling(isur,surface,idim,ibead)
                ! end do
                ! end do
                ! if (rtemp.gt.0.0d0) then
                    ! sign=1.0d0
                ! else
                    ! sign=-1.0d0
                ! end if
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! velocityveti(idim,ibead)=velocity(idim,ibead)-derivativecoupling(isur,surface,idim,ibead)*rtemp
                ! end do
                ! end do
                
                ! hopposition(nsurfaceswitch+nfrustratedhops+1)=istep-1
                
                ! rtemp=rtemp**2+(potsurface(surface)-potsurface(isur))*dble(nbeads)*2.0d0/mass
                ! if (rtemp.gt.0.0d0) then
                    ! rtemp=Dsqrt(rtemp)
                    ! do idim=1,dimensionality,1
                    ! do ibead=1,nbeads,1
                        ! velocity(idim,ibead)=velocityveti(idim,ibead)+derivativecoupling(isur,surface,idim,ibead)*rtemp*sign
                    ! end do
                    ! end do
                    ! nsurfaceswitch=nsurfaceswitch+1

                    ! hopstate(nsurfaceswitch+nfrustratedhops)=surface
                    ! hopsuccess(nsurfaceswitch+nfrustratedhops)=1
                    ! hopvelocity(nsurfaceswitch+nfrustratedhops,:)=velocity(1,:)
                    
                    ! surface=isur
                    ! call PotentialRead (dimensionality,nsurfaces,nbeads,1)
                ! else    
                    ! nfrustratedhops=nfrustratedhops+1
                    ! hopstate(nsurfaceswitch+nfrustratedhops)=isur
                    ! hopsuccess(nsurfaceswitch+nfrustratedhops)=0
                    ! hopvelocity(nsurfaceswitch+nfrustratedhops,:)=velocity(1,:)
                    
                ! end if
                ! exit
            ! end if
        ! end do
    
        ! if (flag_hopattempt.eq.0) then
            ! weight=weight/(1.0d0-hopprob(nsurfaces))
        ! end if
        
    end if
    
    call centroid(dimensionality,nsurfaces,natoms,nbeads)
    call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,2)

    if (xi.lt.0.0d0) then
        regionnew=0
    else
        regionnew=1
    end if
    if(region-regionnew.gt.0) then
        crossing=crossing+1
    end if
    region=regionnew
    
end subroutine

subroutine ConstainedIntegration(dimensionality,nsurfaces,natoms,nbeads,nescaleratio,istep)
    use potential
    IMPLICIT NONE
    
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads,nescaleratio,istep
    real(8)::hopprob(nsurfaces),eigenvalue(nsurfaces),velocityveti(dimensionality,nbeads)
    real(8)::positionnmsave(dimensionality,natoms,nbeads),velocitynmsave(dimensionality,natoms,nbeads)
    complex(8)::transition(nsurfaces,nsurfaces),evolution(nsurfaces,nsurfaces)
    real(8)::eletimestep,prob,rtemp,sign,evotime
    integer::idim,ibead,jbead,iatom,isur,jsur,ielestep,iflag_hophappen
    complex(8)::comtemp
    
    complex(8)::work(50)
    real(8)::rwork(100)
    integer::iwork(40),info
        
    real(8)::randomnumber,gaussrannum
    eletimestep=timestep/dble(nescaleratio)
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        position(idim,iatom,ibead)=positionrecord(idim,iatom,ibead,istep)
    end do
    end do
    end do
    call PotentialRead(dimensionality,nsurfaces,natoms,nbeads,1)
    
    iflag_hophappen=0
    if (hopmarker.gt.0) then
    if (istep.eq.hopposition(hopmarker)) then

        iflag_hophappen=1
    
        ! velocity(1,:)=-hopvelocity(hopmarker,:)
        ! hopprob=0.0d0
        ! isur=hopstate(hopmarker)
        
        ! do idim=1,dimensionality,1
        ! do ibead=1,nbeads,1
            ! hopprob(isur)=hopprob(isur)+velocity(idim,ibead)*derivativecoupling(isur,surface,idim,ibead)
        ! end do
        ! end do
        ! hopprob(isur)=hopprob(isur)/dble(nbeads)
        ! comtemp=elecoeffobs(surface)*Dconjg(elecoeffobs(isur))
        ! hopprob(isur)=hopprob(isur)*Dreal(comtemp)*timestep
        ! hopprob(isur)=-2.0d0*hopprob(isur)/(Cdabs(elecoeffobs(surface))**2)
        ! if (hopprob(isur).lt.0.0d0) then
            ! hopprob(isur)=0.0d0
        ! end if
        
        ! weight=weight*hopprob(isur)
        
        ! if (hopsuccess(hopmarker).eq.1) then

            ! surface=isur
            ! hopmarker=hopmarker-1
            ! call PotentialRead(dimensionality,nsurfaces,nbeads,1)
        
        ! else if (hopsuccess(hopmarker).eq.0) then
           
            ! hopmarker=hopmarker-1
            
        ! end if
    
        ! do ibead=1,nbeads,1
            ! velocity(1,ibead)=-velocityrecord(ibead,istep)
        ! end do
    
    end if
    end if
    
    if (iflag_hophappen.eq.0) then
    
        ! do ibead=1,nbeads,1
            ! velocity(1,ibead)=-velocityrecord(ibead,istep)
        ! end do
        
        ! if (nsurfaces.gt.1) then
            ! hopprob=0.0d0
            ! do isur=1,nsurfaces,1
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! hopprob(isur)=hopprob(isur)+velocity(idim,ibead)*derivativecoupling(isur,surface,idim,ibead)
                ! end do
                ! end do
                ! hopprob(isur)=hopprob(isur)/dble(nbeads)
                ! comtemp=elecoeffobs(surface)*Dconjg(elecoeffobs(isur))
                ! hopprob(isur)=hopprob(isur)*Dreal(comtemp)*timestep
                ! hopprob(isur)=-2.0d0*hopprob(isur)/(Cdabs(elecoeffobs(surface))**2)
                ! if (hopprob(isur).lt.0.0d0) then
                    ! hopprob(isur)=0.0d0
                ! end if
            ! end do
            ! do jsur=nsurfaces-1,1,-1
                    ! hopprob(nsurfaces)=hopprob(nsurfaces)+hopprob(jsur)
            ! end do

            ! weight=weight*(1.0d0-hopprob(nsurfaces))
        ! end if
                
    end if
        
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derisurface(idim,iatom,ibead)/mass(iatom)*timestep*0.5d0
    end do
    end do
    end do
    do idim=1,dimensionality,1
    do iatom=1,natoms,1    
    positionnmsave(idim,iatom,:)=0.0d0
    velocitynmsave(idim,iatom,:)=0.0d0
        do ibead=1,nbeads,1
        do jbead=1,nbeads,1
            positionnmsave(idim,iatom,ibead)=positionnmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*position(idim,iatom,jbead)
            velocitynmsave(idim,iatom,ibead)=velocitynmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*velocity(idim,iatom,jbead)
        end do
        end do
    end do
    end do
    
    evotime=eletimestep/2.0d0
    call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
    call PotentialRead(dimensionality,nsurfaces,natoms,nbeads,1)
    do isur=1,nsurfaces,1
        phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
    end do

    if (nsurfaces.gt.1) then    
        ! do ielestep=1,nescaleratio,1
    
            ! transition=(0.0d0,0.0d0)
            ! do isur=1,nsurfaces,1
            ! do jsur=isur+1,nsurfaces,1
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! transition(isur,jsur)=transition(isur,jsur)+velocity(idim,ibead)*derivativecoupling(isur,jsur,idim,ibead)
                ! end do
                ! end do
                ! transition(isur,jsur)=transition(isur,jsur)/dble(nbeads)
                ! transition(isur,jsur)=transition(isur,jsur)*CDexp((0.0d0,1.0d0)*(phase(isur)-phase(jsur)))
                ! transition(isur,jsur)=-transition(isur,jsur)*(0.0d0,1.0d0)
                ! transition(jsur,isur)=Dconjg(transition(isur,jsur))
            ! end do
            ! end do
            
            ! call zheevd('V','U',nsurfaces,transition,nsurfaces,eigenvalue,work,50,rwork,100,iwork,40,info)
            ! evolution=(0.0d0,0.0d0)
            ! do isur=1,nsurfaces,1
                ! evolution(isur,isur)=CDexp(-(0.0d0,1.0d0)*eigenvalue(isur)*eletimestep)
            ! end do
            ! evolution=matmul(transition,evolution)
            ! evolution=matmul(evolution,Dconjg(transpose(transition)))
            
            ! elecoeffobs=(0.0d0,0.0d0)
            ! do isur=1,nsurfaces,1
            ! do jsur=1,nsurfaces,1
                ! elecoeffobs(isur)=elecoeffobs(isur)+evolution(isur,jsur)*elecoeffevo(jsur)
            ! end do
            ! end do
            ! elecoeffevo=elecoeffobs
        
            ! if (ielestep.ne.nescaleratio) then
            ! do isur=1,nsurfaces,1
                ! phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
            ! end do
            
            ! evotime=evotime+eletimestep
            ! call freerpevolution(dimensionality,nbeads,evotime,positionnmsave,velocitynmsave)
            ! call PotentialRead (dimensionality,nsurfaces,nbeads,1)
            ! do isur=1,nsurfaces,1
                ! phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
            ! end do
            ! end if
        
        ! end do            
    end if
    
    do isur=1,nsurfaces,1
        phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
    end do
    ! evotime=evotime+eletimestep/2.0d0
    evotime=timestep
    
    call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
    call PotentialRead(dimensionality,nsurfaces,natoms,nbeads,1)

    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derisurface(idim,iatom,ibead)/mass(iatom)*timestep*0.5d0
    end do
    end do
    end do
    do isur=1,nsurfaces,1
        elecoeffobs(isur)=elecoeffevo(isur)*CDexp((0.0d0,-1.0d0)*phase(isur))
    end do
    
end subroutine


subroutine Integration(dimensionality,nsurfaces,natoms,nbeads,nescaleratio)
    USE POTENTIAL
    IMPLICIT NONE
    
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads,nescaleratio
    real(8)::hopprob(nsurfaces),eigenvalue(nsurfaces),velocityveti(dimensionality,natoms,nbeads)
    ! real(8)::positionsave(dimensionality,natoms,nbeads),velocitysave(dimensionality,natoms,nbeads)
    real(8)::positionnmsave(dimensionality,natoms,nbeads),velocitynmsave(dimensionality,natoms,nbeads)
    complex(8)::transition(nsurfaces,nsurfaces),evolution(nsurfaces,nsurfaces)
   
    real(8)::eletimestep,prob,rtemp,rtemp2,sign,evotime,xisave,lagrangian,enetot,norm
    integer::idim,iatom,ibead,jbead,isur,jsur,ielestep,iterrattle,regionnew,rattleflag
    real(8)::dxisave(dimensionality,natoms)
    complex(8)::comtemp
    
    complex(8)::work(50)
    real(8)::rwork(100)
    integer::iwork(40),info
        
    real(8)::randomnumber,gaussrannum
    
    if (flag_constrain.eq.1) then
        xisave=xi
        dxisave=dxi_centroid
    end if
    
    eletimestep=timestep/dble(nescaleratio)
    
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)
    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "1",enetot,xi(1),xi(2),dxidt
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derisurface(idim,iatom,ibead)/mass(iatom)*timestep*0.5d0
    end do
    end do
    end do
    
    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "S7",enetot,xi(1),xi(2),dxidt
    
    if (flag_constrain.eq.1) then
        
        ! positionsave=position
        ! velocitysave=velocity
        
        do idim=1,dimensionality,1
        do iatom=1,natoms,1    
        positionnmsave(idim,iatom,:)=0.0d0
        velocitynmsave(idim,iatom,:)=0.0d0
            do ibead=1,nbeads,1
            do jbead=1,nbeads,1
                positionnmsave(idim,iatom,ibead)=positionnmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*position(idim,iatom,jbead)
                velocitynmsave(idim,iatom,ibead)=velocitynmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*velocity(idim,iatom,jbead)
            end do
            end do
        end do
        end do
        evotime=timestep
        call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
        call centroid(dimensionality,nsurfaces,natoms,nbeads)
        call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
        
        !FIRST RATTLE CORRECTION
        rattleflag=0
        iterrattle=1
        lagrangian=0.0d0
        do 
            
            rtemp=0.0d0
            do idim=1,dimensionality,1
            do iatom=1,natoms,1
                rtemp=rtemp+dxi_centroid(idim,iatom)*dxisave(idim,iatom)/mass(iatom)
            end do
            end do
            rtemp=rtemp*timestep**2
            rtemp=xi/rtemp
            
            lagrangian=lagrangian-rtemp
            
            rtemp=lagrangian*timestep**2
            do idim=1,dimensionality,1
            do iatom=1,natoms,1
                position_centroid(idim,iatom)=position_centroid(idim,iatom)+rtemp/mass(iatom)*dxisave(idim,iatom)
            end do
            end do
            
            call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
            
            if ((Dabs(xi).lt.1.0d-8).or.(Dabs(rtemp).lt.1.0d-10)) then
            ! if ((Dabs(xi).lt.1.0d-8)) then
                exit
            end if
            
            iterrattle=iterrattle+1
            
            if (iterrattle.eq.100) then
                ! if (rattleflag.eq.0) then
                    
                    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
                    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
                    
                    ! xiindex=3-xiindex
                    ! lagrangian=0.0d0
                    ! iterrattle=1
                    ! rattleflag=1
                    
                    ! write (*,*) "RATTLE SUGGESTS CHANGING REACTION CHANNEL"
                ! else
                    write (*,*) "WARNING IN RATTLE ITERATION, NOT CONVERGED"
                    exit
                ! end if

            end if

        end do

        ! write(*, *) "1, Lag=", lagrangian
        
        evotime=0.0d0
        call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
        
        do idim=1,dimensionality,1
        do iatom=1,natoms,1
        do ibead=1,nbeads,1
            ! velocity(idim,iatom,ibead)=velocitysave(idim,iatom,ibead)+lagrangian*timestep/mass(iatom)*dxisave(xiindex,idim,iatom)
            velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)+lagrangian*timestep/mass(iatom)*dxisave(idim,iatom)
        end do
        end do
        end do
        ! position=positionsave
    
    end if
    
    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "2",enetot,xi(1),xi(2),dxidt
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1    
    positionnmsave(idim,iatom,:)=0.0d0
    velocitynmsave(idim,iatom,:)=0.0d0
        do ibead=1,nbeads,1
        do jbead=1,nbeads,1
            positionnmsave(idim,iatom,ibead)=positionnmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*position(idim,iatom,jbead)
            velocitynmsave(idim,iatom,ibead)=velocitynmsave(idim,iatom,ibead)+normalmodectn(ibead,jbead)*velocity(idim,iatom,jbead)
        end do
        end do
    end do
    end do
    
    evotime=eletimestep/2.0d0
    call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
    call PotentialRead(dimensionality,nsurfaces,natoms,nbeads,1)
    if (flag_umbrella.eq.1) then
        call centroid(dimensionality,nsurfaces,natoms,nbeads)
        call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,3)
        call biaspotential(dimensionality,nsurfaces,natoms,nbeads)
    end if
    
    do isur=1,nsurfaces,1
        phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
    end do
    
    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "3",enetot,xi(1),xi(2),dxidt
    
    if (nsurfaces.gt.1) then
        ! do ielestep=1,nescaleratio,1
            ! transition=(0.0d0,0.0d0)
            ! do isur=1,nsurfaces,1
            ! do jsur=isur+1,nsurfaces,1
                ! do idim=1,dimensionality,1
                ! do iatom=1,natoms,1
                ! do ibead=1,nbeads,1
                    ! transition(isur,jsur)=transition(isur,jsur)+velocity(idim,iatom,ibead)*derivativecoupling(isur,jsur,idim,ibead)
                ! end do
                ! end do
                ! end do
                ! transition(isur,jsur)=transition(isur,jsur)/dble(nbeads)
                ! transition(isur,jsur)=transition(isur,jsur)*CDexp((0.0d0,1.0d0)*(phase(isur)-phase(jsur)))
                ! transition(isur,jsur)=-transition(isur,jsur)*(0.0d0,1.0d0)
                ! transition(jsur,isur)=Dconjg(transition(isur,jsur))
            ! end do
            ! end do
            
            ! call zheevd('V','U',nsurfaces,transition,nsurfaces,eigenvalue,work,50,rwork,100,iwork,40,info)
            ! evolution=(0.0d0,0.0d0)
            ! do isur=1,nsurfaces,1
                ! evolution(isur,isur)=CDexp(-(0.0d0,1.0d0)*eigenvalue(isur)*eletimestep)
            ! end do
            ! evolution=matmul(transition,evolution)
            ! evolution=matmul(evolution,Dconjg(transpose(transition)))
            
            ! elecoeffobs=(0.0d0,0.0d0)
            ! do isur=1,nsurfaces,1
            ! do jsur=1,nsurfaces,1
                ! elecoeffobs(isur)=elecoeffobs(isur)+evolution(isur,jsur)*elecoeffevo(jsur)
            ! end do
            ! end do
            ! elecoeffevo=elecoeffobs
        
            ! if (ielestep.ne.nescaleratio) then
            ! do isur=1,nsurfaces,1
                ! phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
            ! end do
            
            ! evotime=evotime+eletimestep
            ! call freerpevolution(dimensionality,nbeads,evotime,positionnmsave,velocitynmsave)
            ! call PotentialRead (dimensionality,nsurfaces,natoms,nbeads,1)
            ! do isur=1,nsurfaces,1
                ! phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
            ! end do
            ! end if
        ! end do            
    end if

    do isur=1,nsurfaces,1
        phase(isur)=phase(isur)+potsurface(isur)*eletimestep/2.0d0
    end do
    ! evotime=evotime+eletimestep/2.0d0
    evotime=timestep
    call freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
    call PotentialRead(dimensionality,nsurfaces,natoms,nbeads,1)
    if (flag_umbrella.eq.1) then
        call centroid(dimensionality,nsurfaces,natoms,nbeads)
        call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,3)
        call biaspotential(dimensionality,nsurfaces,natoms,nbeads)
    end if
    
    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)    
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "4",enetot,xi(1),xi(2),dxidt
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derisurface(idim,iatom,ibead)/mass(iatom)*timestep*0.5d0
    end do
    end do
    end do
    
    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "S8",enetot,xi(1),xi(2),dxidt
    
    if (flag_constrain.eq.1) then
        
        call centroid(dimensionality,nsurfaces,natoms,nbeads)
        call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
        
        !SECOND RATTLE
        rtemp=0.0d0
        do idim=1,dimensionality,1
        do iatom=1,natoms,1
        do ibead=1,nbeads,1    
            rtemp=rtemp+velocity(idim,iatom,ibead)*dxi_centroid(idim,iatom)
        end do
        end do
        end do
        rtemp=rtemp/dble(nbeads)
        rtemp2=0.0d0
        do idim=1,dimensionality,1
        do iatom=1,natoms,1  
            rtemp2=rtemp2+dxi_centroid(idim,iatom)*dxi_centroid(idim,iatom)/mass(iatom)
        end do
        end do    
        lagrangian=-rtemp/rtemp2
        
        ! write(*, *) "2, Lag=", lagrangian
        
        do idim=1,dimensionality,1
        do iatom=1,natoms,1
        do ibead=1,nbeads,1
            velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)+lagrangian*dxi_centroid(idim,iatom)/mass(iatom)
        end do
        end do
        end do
    
        call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    
    end if

    ! call centroid(dimensionality,nsurfaces,natoms,nbeads)
    ! call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,1)
    ! call energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,1)    
    ! write(*,"(A10,1pe20.6,1pe20.6,1pe20.6,1pe20.6)",advance='no') "5",enetot,xi(1),xi(2),dxidt
    ! write(*,"(A1)") " "
    
    if (nsurfaces.gt.1) then
        ! hopprob=0.0d0
        ! do isur=1,nsurfaces,1
            ! do idim=1,dimensionality,1
            ! do iatom=1,natoms,1
            ! do ibead=1,nbeads,1
                ! hopprob(isur)=hopprob(isur)+velocity(idim,iatom,ibead)*derivativecoupling(isur,surface,idim,ibead)
            ! end do
            ! end do
            ! end do
            ! hopprob(isur)=hopprob(isur)/dble(nbeads)
            ! comtemp=elecoeffobs(surface)*Dconjg(elecoeffobs(isur))
            ! hopprob(isur)=hopprob(isur)*Dreal(comtemp)*timestep
            ! hopprob(isur)=-2.0d0*hopprob(isur)/(Cdabs(elecoeffobs(surface))**2)
            ! if (hopprob(isur).lt.0.0d0) then
                ! hopprob(isur)=0.0d0
            ! end if
        ! end do
        ! do isur=nsurfaces,2,-1
            ! do jsur=isur-1,1,-1
                ! hopprob(isur)=hopprob(isur)+hopprob(jsur)
            ! end do
        ! end do
        
        ! prob=randomnumber(randomseed)
        ! do isur=1,nsurfaces,1
            ! if (prob.lt.(hopprob(isur))) then

                ! rtemp=0.0d0
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! rtemp=rtemp+derivativecoupling(isur,surface,idim,ibead)**2
                ! end do
                ! end do
                ! rtemp=1.0d0/Dsqrt(rtemp)
                ! do idim=1,dimensionality,1
                ! do ibead=1,nbeads,1
                    ! derivativecoupling(isur,surface,idim,ibead)=derivativecoupling(isur,surface,idim,ibead)*rtemp
                ! end do
                ! end do
                
                ! rtemp=0.0d0
                ! do idim=1,dimensionality,1
                ! do iatom=1,natoms,1
                ! do ibead=1,nbeads,1
                    ! rtemp=rtemp+velocity(idim,iatom,ibead)*derivativecoupling(isur,surface,idim,ibead)
                ! end do
                ! end do
                ! end do
                ! if (rtemp.gt.0.0d0) then
                    ! sign=1.0d0
                ! else
                    ! sign=-1.0d0
                ! end if
                ! do idim=1,dimensionality,1
                ! do iatom=1,natoms,1
                ! do ibead=1,nbeads,1
                    ! velocityveti(idim,iatom,ibead)=velocity(idim,iatom,ibead)-derivativecoupling(isur,surface,idim,ibead)*rtemp
                ! end do
                ! end do
                ! end do
                
                ! !TO MYSELF: RTEMP MIGHT NEED TO BE A MATRIX OF [NATOMS] SHAPE
                ! do iatom=1,natoms,1
                    ! rtemp=rtemp**2+(potsurface(surface)-potsurface(isur))*dble(nbeads)*2.0d0/mass(iatom)
                ! end do
                
                ! if (rtemp.gt.0.0d0) then
                    ! rtemp=Dsqrt(rtemp)
                    ! do idim=1,dimensionality,1
                    ! do iatom=1,natoms,1
                    ! do ibead=1,nbeads,1
                        ! velocity(idim,iatom,ibead)=velocityveti(idim,iatom,ibead)+derivativecoupling(isur,surface,idim,ibead)*rtemp*sign
                    ! end do
                    ! end do
                    ! end do
                    ! surface=isur
                    ! call PotentialRead(dimensionality,nsurfaces,nbeads,1)
                ! end if
                ! exit
            ! end if
        ! end do
    end if
    
    do isur=1,nsurfaces,1
        elecoeffobs(isur)=elecoeffevo(isur)*CDexp((0.0d0,-1.0d0)*phase(isur))
    end do
    
    call dividing_surface(dimensionality,nsurfaces,natoms,nbeads,2)
    
    if (xi.lt.0.0d0) then
        regionnew=0
    else
        regionnew=1
    end if
    if(regionnew-region.gt.0) then
        crossing=crossing+1
    end if
    region=regionnew
    
end subroutine

!USE INITIAL NORMAL MODE POSITION AND VELOCITY AS INPUTS,
!EVOLVE FREE RING-POLYMER FOR A PERIOD OF evotime; 
subroutine freerpevolution(dimensionality,natoms,nbeads,evotime,positionnmsave,velocitynmsave)
    use potential
    IMPLICIT NONE
    integer,intent(in)::dimensionality,natoms,nbeads
    real(8),intent(in)::positionnmsave(dimensionality,natoms,nbeads),velocitynmsave(dimensionality,natoms,nbeads),evotime
    integer::idim,ibead,jbead,iatom
    real(8)::costemp,sintemp,rtemp
    real(8)::positionnmevo(dimensionality,natoms,nbeads),velocitynmevo(dimensionality,natoms,nbeads)
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    velocitynmevo(idim,iatom,1)=velocitynmsave(idim,iatom,1)
    positionnmevo(idim,iatom,1)=positionnmsave(idim,iatom,1)+velocitynmsave(idim,iatom,1)*evotime
    do ibead=2,nbeads,1
        costemp=Dcos(normalmodefre(ibead)/betan*evotime)
        sintemp=Dsin(normalmodefre(ibead)/betan*evotime)
        velocitynmevo(idim,iatom,ibead)=costemp*velocitynmsave(idim,iatom,ibead)
        velocitynmevo(idim,iatom,ibead)=velocitynmevo(idim,iatom,ibead)-normalmodefre(ibead)/betan*sintemp*positionnmsave(idim,iatom,ibead)
        positionnmevo(idim,iatom,ibead)=costemp*positionnmsave(idim,iatom,ibead)
        positionnmevo(idim,iatom,ibead)=positionnmevo(idim,iatom,ibead)+betan/normalmodefre(ibead)*sintemp*velocitynmsave(idim,iatom,ibead)
    end do
    end do
    end do
    
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    position(idim,iatom,:)=0.0d0
    velocity(idim,iatom,:)=0.0d0
        do ibead=1,nbeads,1
        do jbead=1,nbeads,1
            position(idim,iatom,ibead)=position(idim,iatom,ibead)+normalmodentc(ibead,jbead)*positionnmevo(idim,iatom,jbead)
            velocity(idim,iatom,ibead)=velocity(idim,iatom,ibead)+normalmodentc(ibead,jbead)*velocitynmevo(idim,iatom,jbead)
        end do
        end do
    end do    
    end do    
   
end subroutine    

subroutine biaspotential(dimensionality,nsurfaces,natoms,nbeads)
    use potential
    IMPLICIT NONE
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads
    integer::idim,iatom,jdim,jatom,ibead
    real(8)::rtemp,rtemp2
    
    rtemp=phi-phiwindow
    potsurface(1)=potsurface(1)+0.5d0*umbrellak*rtemp**2
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
        derisurface(idim,iatom,:)=derisurface(idim,iatom,:)+umbrellak*rtemp*dphi_centroid(idim,iatom)
    end do
    end do

    rtemp=0.0d0
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
        rtemp=rtemp+dphi_centroid(idim,iatom)**2/mass(iatom)
    end do
    end do
    potsurface(1)=potsurface(1)-Dlog(rtemp)/inversetemp*0.5d0

    do idim=1,dimensionality,1
    do iatom=1,natoms,1   
        rtemp2=0.0d0
        do jdim=1,dimensionality,1
        do jatom=1,natoms,1
            rtemp2=rtemp2+d2phi_centroid(idim,iatom,jdim,jatom)*dphi_centroid(jdim,jatom)/mass(jatom)
        end do
        end do
        derisurface(idim,iatom,:)=derisurface(idim,iatom,:)-rtemp2/rtemp/inversetemp
    end do
    end do
    
end subroutine

subroutine centroid(dimensionality,nsurfaces,natoms,nbeads)
    USE POTENTIAL
    IMPLICIT NONE
    
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads
    integer::idim,iatom,ibead,ibond
    
    position_centroid=0.0d0
    do idim=1,dimensionality,1
    do iatom=1,natoms,1
    do ibead=1,nbeads,1
        position_centroid(idim,iatom)=position_centroid(idim,iatom)+position(idim,iatom,ibead)
    end do
    end do
    end do
    position_centroid=position_centroid/dble(nbeads)

end subroutine

!#######################################################################            
!#####  THREE MODES EXIST, WITH INTEGER "flag_enmode"
!#####  1 - RETURN ENERGY ONLY
!#####  2 - RETURN ELECTRONIC NORM ONLY
!#####  3 - RETURN BOTH
!#######################################################################            
        
subroutine energy_norm(dimensionality,nsurfaces,natoms,nbeads,enetot,norm,flag_enmode)
    USE POTENTIAL
    IMPLICIT NONE
    
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads,flag_enmode
    real(8),intent(out)::enetot,norm    
    
    real(8)::rtemp
    integer::idim,iatom,ibead,isur
    
    enetot=0.0d0
    norm=0.0d0
    
    if (flag_enmode.eq.2) then
        
        do isur=1,nsurfaces,1
            norm=norm+CDabs(elecoeffobs(isur))**2
        end do
    
    else
    
        do iatom=1,natoms,1
        rtemp=0.0d0
        do idim=1,dimensionality,1
        do ibead=1,nbeads-1,1
            rtemp=rtemp+(position(idim,iatom,ibead+1)-position(idim,iatom,ibead))**2
        end do
        rtemp=rtemp+(position(idim,iatom,1)-position(idim,iatom,nbeads))**2
        end do
        enetot=enetot+rtemp*mass(iatom)
        end do
        enetot=0.5d0*enetot/(betan**2)
        
        do iatom=1,natoms,1       
        rtemp=0.0d0
        do idim=1,dimensionality,1
        do ibead=1,nbeads,1
            rtemp=rtemp+velocity(idim,iatom,ibead)**2
        end do
        end do
        enetot=enetot+0.5d0*mass(iatom)*rtemp
        end do
        enetot=enetot+potsurface(surface)*dble(Nbeads)
        
        if (flag_enmode.eq.3) then 
            
            do isur=1,nsurfaces,1
                norm=norm+CDabs(elecoeffobs(isur))**2
            end do
            
        end if
        
    end if
    
    ! write(*,*) "PLEASE PROVIDE CORRECT FLAG_ENMODE"
    ! STOP
    
end subroutine


FUNCTION randomnumber(idum)
    INTEGER idum,IM1,IM2,IMM1,IA1,IA2
    INTEGER IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL(8) randomnumber,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1)
    PARAMETER (IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211)
    PARAMETER (IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB)
    PARAMETER (EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
    idum=max(-idum,1)
    idum2=idum
    do j=NTAB+8,1,-1
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    if (j.le.NTAB) iv(j)=idum
    enddo
    iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if(idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1) iy=iy+IMM1
    randomnumber=min(AM*iy,RNMX)
    RETURN
END FUNCTION

FUNCTION gaussrannum(sigma,idum)
    INTEGER idum
    REAL(8),intent(in)::sigma
    INTEGER iset
    REAL(8) fac,gset,rsq,v1,v2,randomnumber,gaussrannum
    SAVE iset,gset
    DATA iset/0/
    if (idum.lt.0) iset=0
    rsq=99.0d0
    if (iset.eq.0) then
    do while(rsq.ge.1..or.rsq.eq.0.)
    v1=2.0d0*randomnumber(idum)-1.0d0
    v2=2.0d0*randomnumber(idum)-1.0d0
    rsq=v1**2+v2**2
    end do
    fac=sqrt(-2.*log(rsq)/rsq)
    gset=v1*fac
    gaussrannum=v2*fac
    iset=1
    else
    gaussrannum=gset
    iset=0
    endif
    gaussrannum=DBLE(gaussrannum)*SIGMA
    return
END FUNCTION gaussrannum
