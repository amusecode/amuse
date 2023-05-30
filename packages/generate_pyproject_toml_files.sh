#!/usr/bin/env bash

while read -r line; do
	pkgdir=$(echo "print('amuse-${line}'.lower())" | python);
	echo ${pkgdir}
	cp pyproject.toml_template ${pkgdir}/pyproject.toml;
done < community_package_names
