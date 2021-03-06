SUBROUTINE setup_grid
    ! =============================================================
    ! Set up the grid
    ! =============================================================
    ! Subroutine for defining the grid of the GCM. Run once
    ! before the loop starts.
    ! -------------------------------------------------------------
    ! The following arrays has to be populated:
    !
    !  dxdy - Area of horizontal cell-walls.
    !  dz   - Height of k-cells in 3 dim (if z_timevar not used)
    !  dzt  - Height of k-cells in 4 dim (if z_timevar used)
    !  kmt  - Number of k-cells from surface to seafloor.
    !
    ! The following might be needed to calculate
    ! dxdy, uflux, and vflux
    !
    !  dzu - Height of each u-gridcell.
    !  dzv - Height of each v-gridcell.
    !  dxv -
    !  dyu -
    ! -------------------------------------------------------------

    USE mod_precdef
    USE mod_param

    USE mod_grid

    IMPLICIT none

    INTEGER               :: jj, ii, kk
    REAL(PP)              :: dlon, dlat

    ! Quarter degree resolution
    dlon = 360./imt; dlat = 180./(jmt-1)

    ! dx and dy in u and v points & grid area
    DO jj = 1, jmt
       DO ii = 1, imt

          dx         = dlon * radian * radius * COS( (-90.+dlat*(jj-0.5))*radian )
          dy         = dlat * radian * radius

          dxv(ii,jj) = dlon * radian * radius * COS( (-90.+dlat*jj)*radian )
          dyu(ii,jj) = dy

          dxdy(ii,jj) = dx * dy
       END DO
    END DO

    ! Read A(k) and B(k) to determine hybrid coordinate levels.
    ! p(k) = A(k) + B(k) * p_surface
    ALLOCATE ( aa(0:KM), bb(0:KM) )

    OPEN (12,FILE=trim(topoDataDir)//trim(zgridFile) )

    DO kk = 0, KM
      READ (12,'(10x,F12.6,4x,F10.8)') aa(kk), bb(kk)
    END DO

    CLOSE (12)


END SUBROUTINE setup_grid
