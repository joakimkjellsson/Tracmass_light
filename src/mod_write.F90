MODULE mod_write
    !!------------------------------------------------------------------------------
    !!
    !!       MODULE: mod_write
    !!
    !!          Opening/ Closing outfiles
    !!          & Writing trajectories to files
    !!
    !!          Subroutines included:
    !!               - open_outfiles
    !!               - close_outfiles
    !!               - write_data
    !!
    !!------------------------------------------------------------------------------

    USE mod_precdef
    USE mod_traj
    USE mod_trajdef
    USE mod_time
    USE mod_loopvars
    USE mod_calendar

    IMPLICIT NONE

    CHARACTER(LEN=200)    :: fullWritePref
    CHARACTER(LEN=200)    :: outDataDir, outDataFile

    INTEGER               :: timeformat

    CONTAINS

    SUBROUTINE open_outfiles
    ! --------------------------------------------------
    !
    ! Purpose:
    ! Open outfiles ini, run, out, err
    !
    ! --------------------------------------------------

        fullWritePref =  TRIM(outDataDir)//TRIM(outDataFile)

        OPEN(UNIT=50, FILE = TRIM(fullWritePref)//'_ini.csv', STATUS='replace')
        OPEN(UNIT=51, FILE = TRIM(fullWritePref)//'_run.csv', STATUS='replace')
        OPEN(UNIT=52, FILE = TRIM(fullWritePref)//'_out.csv', STATUS='replace')

    END SUBROUTINE open_outfiles

    SUBROUTINE close_outfiles
    ! --------------------------------------------------
    !
    ! Purpose:
    ! Close outfiles ini, run, out, err
    !
    ! --------------------------------------------------

        CLOSE (50)
        CLOSE (51)
        CLOSE (52)

    END SUBROUTINE close_outfiles

    SUBROUTINE write_data(sel)
    ! --------------------------------------------------
    !
    ! Purpose:
    ! Write data to outfile
    !
    ! --------------------------------------------------

        CHARACTER(LEN=*)     :: sel

        SELECT CASE (TRIM(sel))


        ! INI file
        CASE ('ini')

            SELECT CASE(timeformat)

                CASE(0)
                ! Include time - tt in seconds
                WRITE(50,"(I8,4(',',F13.5),1(',',F20.5))")  ntrac, x1, y1, z1, subvol, tt
                RETURN

                CASE(1)
                ! Include time - Fraction ts
                WRITE(50,"(I8,5(',',F13.5))")  ntrac, x1, y1, z1, subvol, ts
                RETURN

                CASE(2)
                ! Include time - YYYY MM DD HH MM SS
                CALL tt_calendar
                WRITE(50,"(I8,4(',',F13.5),(',',I5),3(',',I3))")  ntrac, x1, y1, z1, subvol, &
                      dateYear, dateMon, dateDay, dateHour
                RETURN

            END SELECT

        ! RUN file
        CASE ('run')
            IF(  (write_frec == 1 .AND. trajectories(ntrac)%niter == niter-1   ) .OR. &
                 (write_frec == 2 .AND. scrivi                    ) .OR. &
                 (write_frec == 3                                 ) .OR. &
                 (write_frec == 5 .AND. .not.scrivi                  ) ) THEN

                SELECT CASE(timeformat)

                    CASE(0)
                    ! Include time - tt in seconds
                    WRITE(51,"(I8,4(',',F13.5),1(',',F20.5))")  ntrac, x1, y1, z1, subvol, tt
                    RETURN

                    CASE(1)
                    ! Include time - Fraction ts
                    WRITE(51,"(I8,5(',',F13.5))")  ntrac, x1, y1, z1, subvol, ts
                    RETURN

                    CASE(2)
                    ! Include time - YYYY MM DD HH MM SS
                    CALL tt_calendar
                    WRITE(51,"(I8,4(',',F13.5),(',',I5),3(',',I3))")  ntrac, x1, y1, z1, subvol,  &
                        dateYear, dateMon, dateDay, dateHour
                    RETURN

                END SELECT
            END IF

        ! OUT file
        CASE ('out')
            SELECT CASE(timeformat)

                CASE(0)
                ! Include time - tt in seconds
                WRITE(52,"(I8,4(',',F13.5),1(',',F20.5))")  ntrac, x1, y1, z1, subvol, tt
                RETURN

                CASE(1)
                ! Include time - Fraction ts
                WRITE(52,"(I8,5(',',F13.5))")  ntrac, x1, y1, z1, subvol, ts
                RETURN

                CASE(2)
                ! Include time - YYYY MM DD HH MM SS
                CALL tt_calendar
                WRITE(52,"(I8,4(',',F13.5),(',',I5),3(',',I3))")  ntrac, x1, y1, z1, subvol, dateYear, dateMon, dateDay, dateHour
                RETURN

            END SELECT
        END SELECT

    END SUBROUTINE write_data

END MODULE mod_write