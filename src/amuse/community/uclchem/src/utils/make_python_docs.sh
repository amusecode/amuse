#!/usr/bin/bash
#User must supply the path to their uclchem.github.io repo clone as the first argument
#don't forget to go add the id and title to the notebook markdowns

#Create markdown file for UCLCHEM docs from the docstrings in the python module
api_file=Tutorials/start-pythonapi.md
echo "---" > $api_file
echo "id: pythonapi" >> $api_file
echo "title: Python Reference" >> $api_file
echo "---" >> $api_file
pydoc-markdown -p uclchem --render-toc >> $api_file
sed -i 's|# Table of Contents|# Python API|g' $api_file
#get rid of version entries in the markdown
sed -i 's|^.*\_version.*||g' $api_file

#create default parameters file
python utils/generate_param_docs.py src/fortran_src/defaultparameters.f90 Tutorials/start-parameters.md

#create markdown file from notebooks then change image links to point to img folder
#of uclchem.github.io. Then move images and notebooks to the right place.
jupyter nbconvert --to markdown Tutorials/*.ipynb --output-dir Tutorials
#point images to paths in website
sed -i 's|\[png\](.*\/|[png](.\/assets\/|g' Tutorials/*.md
#remove styles from tables because they throw mdx errors
find Tutorials/*.md -exec python utils/remove_table_styles.py  {} +
mv Tutorials/*files/*.png $1/website/docs/assets
cp Tutorials/assets/* $1/website/docs/assets

mv Tutorials/*.md $1/website/docs
rm -r Tutorials/*files/
