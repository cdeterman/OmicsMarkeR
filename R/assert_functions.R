
# Note, this function is directly duplicated from the
# assertive package.  After contacting the maintainer, 
# the function is currently not exported.  Until such a 
# time that it is I will keep this function here.  All credit
# goes to the original author of the assertive package.  I only
# wish to extend the use of this core function for other asserts.
# I will remove once the function is exported from the assertive package.
assert_engine <-
    function (x, predicate, msg, what = c("all", "any"), ...) 
    {
        handler <- match.fun(match.arg(getOption("assertive.severity"), 
                                       c("stop", "warning", "message")))
        what <- match.fun(match.arg(what))
        ok <- if (missing(x)) 
            predicate()
        else predicate(x, ...)
        if (!what(ok)) {
            if (missing(msg)) {
                if (is_scalar(ok)) {
                    msg <- cause(ok)
                }
                else {
                    stop("Bug in assertive; error message is missing")
                }
            }
            handler(msg, call. = FALSE)
        }
    }

assert_is_in_closed_range <- 
    function (x, lower = -Inf, upper = Inf) 
    {
        msg <- sprintf("%s is not in range %s to %s.", 
                       get_name_in_parent(x), lower, upper)
        assert_engine(x, is_in_closed_range, msg, lower = lower, 
                      upper = upper)
    }

assert_is_positive <-
    function (x) 
    {
        msg <- sprintf("%s is non-positive.", get_name_in_parent(x))
        assert_engine(x, is_positive, msg)
    }