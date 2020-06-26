MODULE mod_getfile
  !!------------------------------------------------------------------------------
  !!
  !!       MODULE: mod_getfile
  !!
  !!          Reads data from netCDF files
  !!
  !!          Subroutines included:
  !!------------------------------------------------------------------------------

#ifndef no_netcdf

  USE mod_grid
  USE mod_log
  USE netcdf

  IMPLICIT NONE

  INTEGER                                    :: ierr, varid,ncid

  CONTAINS
  
  FUNCTION filledFileName(filePattern, currYear, currMon, currDay, currDmax, currHour, currMin, currSec)
  ! --------------------------------------------------
  !
  ! Purpose:
  ! Take filepattern and fill in year, month, day etc.
  !
  ! --------------------------------------------------
  
      CHARACTER (len=*)                       :: filePattern
      CHARACTER (len=LEN(filePattern))        :: filledFileName
      CHARACTER (len=8)                       :: timestamp_yyyymmdd
      
      INTEGER                                 :: currYear, currMon, currDay 
      INTEGER, OPTIONAL                       :: currDmax, currHour, currMin, currSec
      INTEGER                                 :: ichar
      
      filledFileName = filePattern      
      
      ! make a time step in format YYYYMMDD 
      ! sometimes this is used
      WRITE(timestamp_yyyymmdd(1:4),'(i4)') currYear
      WRITE(timestamp_yyyymmdd(5:6),'(i2)') currMon
      WRITE(timestamp_yyyymmdd(7:8),'(i2)') currDay 
      
      ! for every YYYY, replace with currYear
      ichar = INDEX(filledFileName,'YYYY')
      DO WHILE (ichar /= 0)
          WRITE(filledFileName(ichar:ichar+3),'(i4)') currYear
          ichar = INDEX(filledFileName,'YYYY')
      END DO
      
      ! for every MM, replace with currMon
      ichar = INDEX(filledFileName,'MM')
      DO WHILE (ichar /= 0)
          WRITE(filledFileName(ichar:ichar+1),'(i2.2)') currMon
          ichar = INDEX(filledFileName,'MM')
      END DO
      
      ! for every DD, replace with currDay
      ichar = INDEX(filledFileName,'DD')
      DO WHILE (ichar /= 0)
          write(filledFileName(ichar:ichar+1),'(i2.2)') currDay
          ichar = INDEX(filledFileName,'DD')
      END DO
      
      ! for every TSTSTSTS, replace with timestamp YYYYMMDD
      ichar = INDEX(filledFileName,'TSTSTSTS')
      DO WHILE (ichar /= 0)
          filledFileName = trim(filledFileName(:ichar-1))//trim(timestamp_yyyymmdd)//trim(filledFileName(ichar+8:))
          ichar = INDEX(filledFileName,'TSTSTSTS')
      END DO
      
      ! replace RUNID with RunID from namelist
      ichar = INDEX(filledFileName,'RUNID')
      DO WHILE (ichar /= 0)
          filledFileName = trim(filledFileName(:ichar-1))//trim(RunID)//trim(filledFileName(ichar+5:))
          ichar = INDEX(filledFileName,'RUNID')
      END DO
      
  END FUNCTION
  
  FUNCTION nferr(ierr,index)
      
      INTEGER                                 :: ierr,index
      CHARACTER (len=100)                     :: nferr
      
      IF(ierr .NE. 0) THEN
          PRINT*,NF90_STRERROR(ierr)
          STOP 
      ENDIF
      
  END FUNCTION
  
  
  FUNCTION get2DfieldNC(fieldFile ,varName, start2D, count2D, cextend)
  ! --------------------------------------------------
  !
  ! Purpose:
  ! Get 2D field data
  !
  ! --------------------------------------------------

      REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2DfieldNC
      REAL, ALLOCATABLE,   DIMENSION(:,:)     :: field
      REAL                                    :: scale_factor, add_offset

      INTEGER, DIMENSION(4)                   :: start2D  ,count2D

      CHARACTER (len=6), OPTIONAL             :: cextend

      CHARACTER (len=*)                       :: fieldFile ,varName

      ALLOCATE(field(count2D(1), count2D(2)))

      IF ( PRESENT(cextend)) THEN
            ALLOCATE(get2DfieldNC(1:imt, 1:jmt+1))
      ELSE
            ALLOCATE(get2DfieldNC(1:imt, 1:jmt))
      END IF
      
      IF (log_level >= 3) THEN
          PRINT*," get2DfieldNC: ",TRIM(fieldFile)
          PRINT*," get2DfieldNC: start2D, count2D", start2D, count2D
          PRINT*," get2DfieldNC: size(field,1), size(field,2) ", size(field,1), size(field,2)
      ENDIF
      
      IF (log_level >= 5) PRINT*,' get2DfieldNC: open ',TRIM(fieldFile)
      CALL err( NF90_OPEN(TRIM(fieldFile) ,NF90_NOWRITE ,ncid) )
      
      IF (log_level >= 5) PRINT*,' get2DfieldNC: inquire ',TRIM(varName)
      CALL err( NF90_INQ_VARID(ncid ,varName ,varid) )

      IF ( (l_subdom) .AND. (imindom > imaxdom) ) THEN
          
          IF (log_level >= 5) PRINT*,"get2DfieldNC: subdomain crosses end of grid"
          start2D(1) = imindom; count2D(1) = imthalf1
          ierr=NF90_GET_VAR(ncid ,varid ,field(1:imthalf1,:), start2D, count2D )
          IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
          
          start2D(1) = 1; count2D(1) = imthalf2
          ierr=NF90_GET_VAR(ncid ,varid ,field(imthalf1+1:imt,:), start2D, count2D )
          IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
 
      ELSE
          
          IF (log_level >= 5) PRINT*,"get2DfieldNC: ncid, varid, start2D, count2D",ncid, varid, start2D, count2D
          ierr=NF90_GET_VAR(ncid ,varid ,field, start2D, count2D )
          IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
          
      END IF
      
      IF (log_level >= 5) PRINT*,' get2DfieldNC: get_att scale_factor'
      ierr = NF90_GET_ATT(ncid, varid,"scale_factor", scale_factor) 
      IF (ierr .NE. 0) THEN
          IF (log_level >= 3) PRINT*,"get2DfieldNC: scale factor not found. set to 1.0" 
          scale_factor = 1.0
      END IF
      
      IF (log_level >= 5) PRINT*,' get2DfieldNC: get_att add_offset '
      ierr = NF90_GET_ATT(ncid, varid,"add_offset", add_offset)      
      IF (ierr .NE. 0) THEN
          IF (log_level >= 3) PRINT*,"get2DfieldNC: add_offset not found. set to 0.0"
          add_offset = 0.0
      END IF
      
      IF (log_level >= 5) PRINT*,' get2DfieldNC: close '
      CALL err( NF90_CLOSE(ncid) )
      
      IF (log_level >= 5) THEN
         PRINT*,"get2DfieldNC: put read data in array ",size(field,1),size(field,2),size(get2DfieldNC,1),size(get2DfieldNC,2)
      ENDIF
      get2DfieldNC(:,:) = field(:,:)*scale_factor + add_offset

   END FUNCTION get2DfieldNC

   FUNCTION get3DfieldNC(fieldFile ,varName, start3D, count3D, stcase, cextend)
   ! --------------------------------------------------
   !
   ! Purpose:
   ! Get 3D field data
   !
   ! --------------------------------------------------

       REAL, ALLOCATABLE,   DIMENSION(:,:,:)       :: get3DfieldNC
       REAL, ALLOCATABLE,   DIMENSION(:,:,:)       :: field
       REAL, ALLOCATABLE,   DIMENSION(:,:,:,:)     :: field4
       REAL                                        :: scale_factor, add_offset

       INTEGER, DIMENSION(4)                     :: start3D,count3D, ss, cc
       INTEGER                                   :: ii, kk

       CHARACTER (len=6), OPTIONAL               :: cextend
       CHARACTER (len=*)                         :: fieldFile ,varName, stcase

       IF (stcase == 'st')  THEN
            ALLOCATE(field(count3D(1), count3D(2),count3D(3)))
       ELSE IF (stcase == 'ts')  THEN
            ALLOCATE(field4(1,count3D(1), count3D(2),count3D(3)))
       ELSE IF (stcase == 'st_r') THEN
            ALLOCATE(field(count3D(3), count3D(2),count3D(1)))
       ELSE IF (stcase == 'ts_r') THEN
            ALLOCATE(field4(1,count3D(3), count3D(2),count3D(1)))
       END IF

       IF (stcase == 'st') THEN
          ss = start3D
          cc = count3D
       ELSE IF (stcase == 'ts') THEN
          ss(1) = start3D(4); ss(2:4) = start3D(1:3)
          cc(1) = count3D(4); cc(2:4) = count3D(1:3)
       ELSE IF (stcase == 'st_r') THEN
          ss(1) = start3D(3); ss(2) = start3D(2); ss(3) = start3D(1); ss(4) = start3D(4)
          cc(1) = count3D(3); cc(2) = count3D(2); cc(3) = count3D(1); cc(4) = count3D(4)
       ELSE IF (stcase == 'ts_r') THEN
          ss(1) = start3D(4); ss(2) = start3D(3); ss(3) = start3D(2); ss(4) = start3D(1)
          cc(1) = count3D(4); cc(2) = count3D(3); cc(3) = count3D(2); cc(4) = count3D(1)
       END IF

       IF ( PRESENT(cextend)) THEN
             ALLOCATE(get3DfieldNC(1:imt, 1:jmt+1, 1:km))
       ELSE
             ALLOCATE(get3DfieldNC(1:imt, 1:jmt, 1:km))
       END IF
       
       IF (log_level >= 5) PRINT*," get3DfieldNC open ",TRIM(fieldFile)
       CALL err( NF90_OPEN(TRIM(fieldFile), NF90_NOWRITE, ncid) )
       !IF(ierr .NE. 0) THEN
       !      PRINT*,"get3DfieldNC: fieldFile: TRIM(fieldFile)"
       !      IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
       !END IF
       
       IF (log_level >= 5) PRINT*," get3DfieldNC inquire varid ",varName
       CALL err( NF90_INQ_VARID(ncid, varName, varid) )
       !IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
       !IF(ierr .NE. 0) STOP 2
       
       IF ( (l_subdom) .AND. (imindom > imaxdom) ) THEN

           IF (stcase == 'st') THEN
              ss(1) = imindom; cc(1) = imthalf1
              ierr=NF90_GET_VAR(ncid, varid, field(1:imthalf1,:,:), ss, cc)
              IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
              ss(1) = 1; cc(1) = imthalf2
              ierr=NF90_GET_VAR(ncid ,varid ,field(imthalf1+1:imt,:,:), ss, cc )
              IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
              
           ELSE IF (stcase == 'st_r') THEN
              ss(3) = imindom; cc(3) = imthalf1
              ierr=NF90_GET_VAR(ncid, varid, field(:,:,1:imtdom), ss, cc)
              IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
              ss(3) = 1; cc(3) = imthalf2
              ierr=NF90_GET_VAR(ncid ,varid ,field(:,:,imtdom+1:imt), ss, cc )
              IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
              
           ELSE IF (stcase == 'ts') THEN
              ss(2) = imindom; cc(2) = imthalf1
              ierr=NF90_GET_VAR(ncid, varid, field4(:,1:imthalf1,:,:), ss, cc)
              IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
              ss(2) = 1; cc(2) = imthalf2
              ierr=NF90_GET_VAR(ncid ,varid ,field4(:,imthalf1+1:imt,:,:), ss, cc )
              IF (ierr .NE. 0) PRINT*,nferr(ierr,2)

           ELSE IF (stcase == 'ts_r') THEN
             ss(4) = imindom; cc(4) = imthalf1
             ierr=NF90_GET_VAR(ncid, varid, field4(:,:,:,1:imthalf1), ss, cc)
             IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
             ss(4) = 1; cc(4) = imthalf2
             ierr=NF90_GET_VAR(ncid ,varid ,field4(:,:,:,imthalf1+1:imt), ss, cc )
             IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
             
           END IF
           IF(ierr .NE. 0) STOP 3
           IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
           
       ELSE
           
           IF (stcase == 'st' .OR. stcase == 'st_r') ierr=NF90_GET_VAR(ncid, varid, field, ss, cc)
           IF (stcase == 'ts' .OR. stcase == 'ts_r') ierr=NF90_GET_VAR(ncid, varid, field4, ss, cc)
           IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
           
       END IF

       ierr = NF90_GET_ATT(ncid, varid,"scale_factor", scale_factor)
       IF (ierr .NE. 0) THEN 
           IF (log_level >= 3) PRINT*," get3DfieldNC: scale_factor not found. set to 1.0 "
           scale_factor = 1.0
       END IF
       ierr = NF90_GET_ATT(ncid, varid,"add_offset", add_offset)
       IF (ierr .NE. 0) THEN 
           IF (log_level >= 3) PRINT*," get3DfieldNC: add_offset not found. set to 0.0"
           add_offset = 0.0
       END IF
       
       CALL err( NF90_CLOSE(ncid) )
       !IF (ierr .NE. 0) PRINT*,nferr(ierr,2)
       !IF(ierr .NE. 0) STOP 4

       IF (stcase == 'st') THEN
          get3DfieldNC(:,:,:) = field(:,:,:)*scale_factor + add_offset
       ELSE IF (stcase == 'ts') THEN
          get3DfieldNC(:,:,:) = field4(1,:,:,:)*scale_factor + add_offset
       ELSE IF (stcase == 'st_r' .OR. stcase == 'ts_r') THEN
          DO ii = 1, imt
             DO kk = 1, km
               IF (stcase == 'st_r') get3DfieldNC(ii,:,kk) = field(kk,:,ii)*scale_factor + add_offset
               IF (stcase == 'ts_r') get3DfieldNC(ii,:,kk) = field4(1,kk,:,ii)*scale_factor + add_offset
             END DO
          END DO
       END IF

    end function get3DfieldNC
    
    
    SUBROUTINE err(ierr)
    
        INTEGER :: ierr
        
        IF( ierr /= 0) THEN
            PRINT*, NF90_STRERROR(ierr)
            STOP
        END IF
    
    END SUBROUTINE   

#endif
END MODULE mod_getfile
