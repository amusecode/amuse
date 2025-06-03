      FUNCTION initialize_code()
      IMPLICIT NONE
      include "stopcond.inc"
      INTEGER :: initialize_code
      INTEGER :: set_support_for_condition
      INTEGER :: return
      initialize_code = 0
      return = set_support_for_condition(COLLISION_DETECTION)
      return = set_support_for_condition(PAIR_DETECTION)
      RETURN
      END FUNCTION
