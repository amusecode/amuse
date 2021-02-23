#!/usr/bin/env bash

while read line; do
	pkgdir=`echo "print('amuse-${line}'.lower())" | python`;
	echo ${pkgdir}
	python setup_template.py ${line} > ${pkgdir}/setup.py;
done < community_package_names
