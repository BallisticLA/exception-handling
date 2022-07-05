*> \brief \b UPDATE_INFO
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download UPDATE_INFO + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/update_info.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/update_info.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/update_info.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE UPDATE_INFO( INFO, INFO_ARRAY, INFO_INTERNAL )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO 
*       ..
*       .. Array Arguments ..
*       INTEGER            INFO_ARRAY( * ), INFO_INTERNAL( 2 )
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> UPDATE_INFO updates INFO and INFO_ARRAY before an LAPACK routine 
*> returns.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          On input, INFO contains the standard value that LAPACK 
*>          would return.
*>          On output, INFO also reports any exceptional inputs or 
*>          outputs, or exceptions in internal subroutine calls, 
*>          as reported in INFO_INTERNAL
*> \endverbatim
*>
*> \param[in,out] INFO_ARRAY
*> \verbatim
*>          INFO_ARRAY(1) is set to the input value of INFO.
*>          INFO_ARRAY(4) is set to the updated value of INFO.
*>          
*> \endverbatim
*>
*> \param[in] INFO_INTERNAL
*> \verbatim
*>          INFO_INTERNAL is INTEGER array, dimension(2)
*>          INFO_INTERNAL( 1 ) is used to track potential changes
*>          to INFO from errors detected by checking input and output
*>          arguments.
*>          INFO_INTERNAL( 2 ) is used to track potential changes
*>          to INFO from errors detected by internal subroutine calls.
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
      SUBROUTINE UPDATE_INFO( INFO, INFO_ARRAY, INFO_INTERNAL )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO 
*     ..
*     .. Array Arguments ..
      INTEGER            INFO_ARRAY( * ), INFO_INTERNAL( 2 )
*
*  =====================================================================
*
*     ..
*     .. Executable Statements ..
*
      INFO_ARRAY(1) = INFO
      INFO_ARRAY(4) = INFO
      IF (INFO_INTERNAL(1) .GT. 0 .AND. INFO_ARRAY(4) .EQ. 0) THEN
*         Update INFO and INFO_ARRAY(4) to indicate exceptional
*         input or output
          INFO_ARRAY(4) = INFO_INTERNAL(1)
          INFO = INFO_INTERNAL(1)
      ENDIF
      IF (INFO_INTERNAL(2) .GT. 0 .AND. INFO_ARRAY(3) .GE. 1 .AND.
     $    INFO_ARRAY(4) .EQ. 0) THEN
*         Update INFO and INFO_ARRAY(4) to indicate exceptions
*         in internal LAPACK call
          INFO_ARRAY(4) = INFO_INTERNAL(2)
          INFO = INFO_INTERNAL(2)
      ENDIF
*     
      RETURN
*
*     End of UPDATE_INFO
*
      END