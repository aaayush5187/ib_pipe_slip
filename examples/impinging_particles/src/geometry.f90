!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
   real(WP), public :: bedbottom,bedheight,bedwidth
   

contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read,param_exists
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='impinging_particles')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
         integer :: i,j,k,nx
         real(WP) :: dx,Lx
         ! Form a box and add a wall on top at location of bed injection
         cfg%VF=1.0_WP
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  if (i.lt.cfg%imin) cfg%VF(i,j,k)=0.0_WP
                  if (i.gt.cfg%imax) cfg%VF(i,j,k)=0.0_WP
                  if (j.lt.cfg%jmin) cfg%VF(i,j,k)=0.0_WP
                  if (k.lt.cfg%kmin) cfg%VF(i,j,k)=0.0_WP
                  if (k.gt.cfg%kmax) cfg%VF(i,j,k)=0.0_WP
               end do
            end do
         end do
         ! Add side-walls to bed
         call param_read('Bed bottom',bedbottom)
         call param_read('Bed height',bedheight)
         call param_read('Bed width' ,bedwidth)
         call param_read('Lx',Lx); call param_read('nx',nx); dx=Lx/real(nx,WP)
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  if (abs(cfg%xm(i)).gt.0.5_WP*bedwidth          .and.&
                  &   abs(cfg%xm(i)).lt.0.5_WP*bedwidth+2.0_WP*dx.and.&
                  &   cfg%ym(j).gt.bedbottom                     .and.&
                  &   cfg%ym(j).lt.bedbottom+bedheight) cfg%VF(i,j,k)=0.0_WP
               end do
            end do
         end do
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
