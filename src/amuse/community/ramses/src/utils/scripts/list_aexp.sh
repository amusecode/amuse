#! /bin/bash
find $1 -name 'info*txt' -exec grep -H aexp '{}' + > dummy
sort dummy
rm dummy
