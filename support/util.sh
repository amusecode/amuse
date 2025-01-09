# Trim excess whitespace off a set of words
#
# The words must not contain ?, * or [.
#
# Args:
#   A set of words, one per argument
#
# Example:
#   words="  a b  c "
#   words=$(trim ${words})
#   # words is now "a b c"
#
trim() {
     printf '%s' "$*"
}


# Determine whether a set of words is a subset of another set of words
# 
# Args:
#   needles: list of words (as a single argument, space separated) to find
#   haystack: list of words (as a single argument, space separated) to find in
#
# Example:
#
#   if is_subset "astropy amuse" "amuse astropy gnuplot" ; then
#       echo "AMUSE and Astropy found!"
#   fi
#
is_subset() {
    is_subset_all_found=1
    for needle in $1 ; do
        is_subset_found=0
        for hay in $2 ; do
            if [ "a${needle}" = "a${hay}" ] ; then
                is_subset_found=1
            fi
        done
        is_subset_all_found=$(( is_subset_all_found & is_subset_found ))
    done
    return $(( ! is_subset_all_found ))
}


# Filter out a set of words from another set of words
#
# Args:
#   remove: list of words (as a single argument, space separated) to remove
#   items: list of words (as a single argument, space separated) to remove them from
#
# Prints any words that are in "set" but not in "remove", space separated and on a
# single line, with an extra space at the end if non-empty.
#
filter_out() {
    filter_out_result=''
    for item in $2 ; do
        filter_out_remove=0
        for remove in $1 ; do
            if [ "a$item" = "a$remove" ] ; then
                filter_out_remove=1
            fi
        done
        if [ "$filter_out_remove" = "0"  ] ; then
            filter_out_result="${filter_out_result} ${item}"
        fi
    done
    filter_out_result=$(trim "$filter_out_result")
    printf '%s' "${filter_out_result}"
}


# Print the main directory of the code for a given package
#
# Codes may have multiple packages, e.g. with and without CUDA support, so different
# packages may have the same main directory.
#
code_directory() {
    package="$1"

    code=$(printf '%s' "${package}" | sed -e 's/^amuse-\([^-]*\).*/\1/')
    printf "src/amuse/community/${code}"
}

