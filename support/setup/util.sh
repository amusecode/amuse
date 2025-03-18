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


# Find any extra packages matching a base package
#
# Args:
#   base_name: Package base name of the form amuse-<code>
#   items: List of package names to search
#
# Given a base package name, e.g. amuse-ph4, this returns a list of all extra packages
# in the given list of the form amuse-ph4-*.
#
extra_packages() {
    extra_packages_result=''
    for item in $2 ; do
        if [ "a${1}" = "a${item%-*}" ] ; then
            extra_packages_result="${extra_packages_result} ${item}"
        fi
    done
    extra_packages_result=$(trim "${extra_packages_result}")
    printf '%s' "${extra_packages_result}"
}


# Check if any of the words in items starts with issue_
#
# Args:
#   items: list of words (as a single argument, space separated) to search
#
has_issue() {
    for item in $1 ; do
        if [ "a${item#issue_}" != "a${item}" ] ; then
            return 0
        fi
    done

    return 1
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


# Normalise a package name
#
# We allow the user to omit the amuse- prefix when referring to a package by name.
# This function adds it back on if it's been omitted.
#
# Args:
#   package: The given package name
#
# Returns:
#   The name prefixed with amuse- if it wasn't already
#
normalise_package_name() {
    package="$1"

    if [ "a${package#amuse-}" = "a${package}" ] ; then
        package="amuse-${package}"
    fi
    printf '%s' "${package}"
}


# Determine installed package name
#
# Our package names are directories, with underscores and dashes in them. The pip and
# conda packages we install have all dashes however, because that's the standard. This
# converts from our directory names to installed package names, so that we can match
# them correctly.
#
# Args:
#   package: The given package name
#
# Returns:
#   The installed name.
#
installed_package_name() {
    printf '%s' "$1" | tr '_' '-'
}


# Forward a command to a package's build system
#
# Args:
#   cmd: The command to forward, e.g. package, test, clean, distclean
#   package: The name of the package to forward to
#   brief: If set to "brief", print only a brief error rather than the full one
#
# Returns:
#   The exit code from the build system
#
forward_to_package() {
    cmd="$1"
    package="$2"
    brief="$3"

    announce_activity ${cmd} ${package}

    code_dir=$(code_directory "${package}")

    ec_file="$(exit_code_file ${cmd} ${package})"
    log_file="$(log_file ${cmd} ${package})"

    if [ "a${cmd}" = "aclean" ] || [ "a${cmd}" = "adistclean" ] ; then
        maketarget="${cmd}"
    else
        maketarget="${cmd}-${package}"
    fi

    (${GMAKE} -C "${code_dir}" "${maketarget}" ; echo "$?" >${ec_file}) 2>&1 | tee ${log_file}

    handle_result $(cat "$ec_file") "${cmd}" "${package}" "${log_file}" "${brief}"
    return $(cat "${ec_file}")
}

