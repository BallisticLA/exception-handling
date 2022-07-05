*> \brief \b CHECKCALL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CHECKCALL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/checkcall.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/checkcall.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/checkcall.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========    
*
*       SUBROUTINE CHECKCALL( FLAG_REPORT_INTERNAL, INFO_INTERNAL, 
*      $                      INFO_CALLARRAY, INFO_ARRAY,
*      $                      CALL_ID, LOC )
*
*       .. Scalar Arguments ..
*       INTEGER            CALL_ID, LOC
*       ..
*       .. Array Arguments ..
*       INTEGER            FLAG_REPORT_INTERNAL(2), INFO_INTERNAL(2)
*       INTEGER            INFO_CALLARRAY( * ), INFO_ARRAY( * )
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHECKCALL checks whether an internal LAPACK call signaled an
*> exception and updates reporting variables INFO_INTERNAL and
*> INFO_ARRAY accordingly.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FLAG_REPORT_INTERNAL
*> \verbatim
*>          FLAG_REPORT_INTERNAL is INTEGER array, dimension(2)
*>          FLAG_REPORT_INTERNAL(1) defines what kinds of exceptions 
*>          to report, using INFO and possibly also INFO_ARRAY for 
*>          more details.
*>          FLAG_REPORT_INTERNAL(1)
*>           <= -1 turns off all error checking
*>            =  0 standard error checks only (eg LDA < 0)
*>            =  1 also check inputs and outputs for Infs and NaN
*>           >=  2 also check input and output arguments of internal 
*>               LAPACK routines (not performed unless 
*>               FLAG_REPORT_INTERNAL(2) >= 1)
*>           FLAG_REPORT_INTERNAL(2) determines how errors are 
*>           reported:
*>           FLAG_REPORT_INTERNAL(2)
*>           <= 0 only returns INFO
*>            = 1 also return INFO_ARRAY with more details
*>            = 2 also calls REPORT_EXCEPTIONS if there are errors
*>           >= 3 also calls REPORT_EXCEPTIONS in internal LAPACK
*>                routines
*> \endverbatim
*>
*> \param[out] INFO_INTERNAL
*> \verbatim
*>          INFO_INTERNAL is INTEGER array, dimension(2)
*>          INFO_INTERNAL( 1 ) is not accessed.
*>          INFO_INTERNAL( 2 ) is used to track potential changes
*>          to INFO from errors detected by internal subroutine calls.
*> \endverbatim
*>
*> \param[in] INFO_CALLARRAY
*> \verbatim
*>          INFO_CALLARRAY is INTEGER array, dimension (*)
*>          When HOW >= 1, INFO_CALLARRAY is the INFO_ARRAY returned
*>          by the routine being checked.
*> \endverbatim
*>
*> \param[in,out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY is INTEGER ARRAY
*>          INFO_ARRAY can report error checks for each
*>          floating point input and output, and each 
*>          internal subroutine call.
*> \endverbatim
*>
*> \param[in] CALL_ID 
*> \verbatim
*>          CALL_ID is INTEGER
*>          CALL_ID is a unique identifier for each internal
*>          LAPACK call that may be checked for exceptions.
*> \endverbatim
*>
*> \param[in] LOC
*> \verbatim
*>          LOC is INTEGER
*>          LOC points to the entry of INFO_ARRAY used to report
*>          on the internal call identified by CALL_ID
*> \endverbatim
*>
*>
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEcomputational
*
*  =====================================================================
      SUBROUTINE CHECKCALL( FLAG_REPORT_INTERNAL, INFO_INTERNAL,
     $                      INFO_CALLARRAY, INFO_ARRAY,
     $                      CALL_ID, LOC )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            CALL_ID, LOC
*     ..
*     .. Array Arguments ..
      INTEGER            FLAG_REPORT_INTERNAL(2), INFO_INTERNAL(2)
      INTEGER            INFO_CALLARRAY( * ), INFO_ARRAY( * )
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            WHAT, HOW
      INTEGER            TMP, TMPINOUT, TMPCALLS, I, NUMARGS, NUMCALLS                
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      WHAT = FLAG_REPORT_INTERNAL( 1 )
      HOW = FLAG_REPORT_INTERNAL( 2 )
      IF (WHAT .GE. 2 .AND. HOW .GE. 1) THEN
*     Check inputs, outputs and internal calls of LAPACK call
*        Update INFO_ARRAY( LOC )
*        Determine whether called routine had any Inf or NaN
*        inputs or outputs
         TMP = INFO_ARRAY( LOC )
*        Determine whether called routine had any Infs or NaNs
*        reported by its own internal calls
         TMPCALLS = 0
*        NUMARGS = number of arguments reported in INFO_CALLARRAY
         NUMARGS = INFO_CALLARRAY(5)
*        NUMCALLS = number of internal subroutine calls reported in
*        INFO_ALLARRAY
         NUMCALLS = INFO_CALLARRAY(6)
         IF (NUMCALLS .GT. 0) THEN
            DO 10 I = 7+NUMARGS, 6+NUMARGS+NUMCALLS
               TMPCALLS = MAX( TMPCALLS, INFO_CALLARRAY( I ) )
10          CONTINUE
         ENDIF
         IF (TMPCALLS .GE. 1) TMP = MAX( TMP, 1 )
*        Determine whether called routine had any inputs or outputs
*        containing Infs or NaNs
         TMPINOUT = 0
         IF (NUMARGS .GT. 0) THEN
            DO 20 I = 7, 6+NUMARGS
               TMPINOUT = MAX( TMPINOUT, INFO_CALLARRAY( I ) )
20          CONTINUE
         ENDIF
         IF (TMPINOUT .GT. 0) TMP = MAX( TMP, TMPINOUT+1 )
         INFO_ARRAY( LOC ) = TMP
*        Update INFO_INTERNAL
         IF (INFO_INTERNAL( 2 ) .EQ. 0 .AND. TMP .GT. 0) 
     $      INFO_INTERNAL( 2 ) =  CALL_ID 
      ENDIF
*     
      RETURN
*
*     End of CHECKCALL
*
      END
