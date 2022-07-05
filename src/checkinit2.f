*> \brief \b CHECKINIT2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CHECKINIT2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/checkinit2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/checkinit2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/checkinit2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHECKINIT2( FLAG_REPORT_INTERNAL, INFO, 
*      $                       INFO_INTERNAL, INFO_ARRAY, NUMARGS, 
*      $                       NUMCALLS )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, NUMARGS, NUMCALLS
*       ..
*       .. Array Arguments ..
*       INTEGER            FLAG_REPORT_INTERNAL( 2 ) 
*       INTEGER            INFO_INTERNAL( 2 )
*       INTEGER            INFO_ARRAY( * )
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHECKINIT2 initializes variables used to report errors.
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
*>               LAPACK routines
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
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          INFO is used to report errors, and contains the value
*>          of INFO from standard input argument checking.
*> \endverbatim
*>
*> \param[out] INFO_INTERNAL
*> \verbatim
*>          INFO_INTERNAL is INTEGER array, dimension(2)
*>          INFO_INTERNAL( 1 ) is used to track potential changes
*>          to INFO from errors detected by checking input and output
*>          arguments.
*>          INFO_INTERNAL( 2 ) is used to track potential changes
*>          to INFO from errors detected by internal subroutine calls.
*>          Both entries are initialized to 0.
*> \endverbatim
*>
*> \param[in,out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY is INTEGER array, dimension( * )
*>          INFO_ARRAY is used to report errors, and is initialized.
*> \endverbatim
*>
*> \param[in] NUMARGS
*> \verbatim
*>          NUMARGS is INTEGER
*>          NUMARGS is the number of input and output arguments
*>          reported in INFO_ARRAY, in locations 
*>          INFO_ARRAY(7:6+NUMARGS). 
*>          NUMARGS is used to initialize INFO_ARRAY.
*> \endverbatim
*>
*> \param[in] NUMCALLS
*> \verbatim
*>          NUMCALLS is INTEGER
*>          NUMCALLS is the number of internal LAPACK calls
*>          reported in INFO_ARRAY, in locations
*>          INFO_ARRAY(7+NUMARGS:6+NUMARGS+NUMCALLS).
*>          NUMCALLS is used to initialize INFO_ARRAY.
*> \endverbatim
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
      SUBROUTINE CHECKINIT2( FLAG_REPORT_INTERNAL, INFO, INFO_INTERNAL, 
     $                       INFO_ARRAY, NUMARGS, NUMCALLS )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NUMARGS, NUMCALLS
*     ..
*     .. Array Arguments ..
      INTEGER            FLAG_REPORT_INTERNAL( 2 ), INFO_INTERNAL( 2 )
      INTEGER            INFO_ARRAY( * )
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            WHAT, HOW
*
*     ..
*     .. Executable Statements ..
*
      INFO_INTERNAL( 1 ) = 0
      INFO_INTERNAL( 2 ) = 0
      WHAT = FLAG_REPORT_INTERNAL( 1 )
      HOW = FLAG_REPORT_INTERNAL( 2 )
      IF (HOW .GE. 1) THEN
*        Initialize INFO_ARRAY
         INFO_ARRAY(1) = INFO
         INFO_ARRAY(2) = WHAT
         INFO_ARRAY(3) = HOW
         INFO_ARRAY(4) = INFO
         INFO_ARRAY(5) = NUMARGS
         INFO_ARRAY(6) = NUMCALLS
         IF (NUMARGS .GT. 0) THEN
            DO 10 I = 7, 6+NUMARGS
*              If argument already marked as checked, do not
*              reinitialize.
               IF (INFO_ARRAY(I).NE.0 .AND. INFO_ARRAY(I).NE.1)
     $            INFO_ARRAY(I) = -1    
10          CONTINUE
         ENDIF
         IF (NUMCALLS .GT. 0) THEN
            DO 20 I = 7+NUMARGS, 6+NUMARGS+NUMCALLS
               INFO_ARRAY(I) = -1
20          CONTINUE
         ENDIF
      END IF
*     
      RETURN
*
*     End of CHECKINIT2
*
      END