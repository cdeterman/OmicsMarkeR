

#' @import assertive
assert_is_in_closed_range <- 
    function (x, lower = -Inf, upper = Inf) 
    {
        msg <- sprintf("%s is not in range %s to %s.", 
                       get_name_in_parent(x), lower, upper)
        assertive.base::assert_engine(is_in_closed_range, x, msg = msg, lower = lower, 
                      upper = upper)
    }

assert_is_positive <-
    function (x) 
    {
        msg <- sprintf("%s is non-positive.", get_name_in_parent(x))
        assertive.base::assert_engine(predicate = is_positive, x, msg=msg)
    }
