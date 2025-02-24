#!/bin/bash

folder_path="."
find "$folder_path" -name '*out.maegz' | while read file; do
  base_name=$(basename "$file" .maegz)
  gzip -dc "$file" > temp_file.txt
  awk '
    BEGIN { inside_block=0; first_block_removed=0 }
    /f_m_ct {/ {
        if (first_block_removed == 0) {
            inside_block=1;
            first_block_removed=1;
            next;
        }
    }
    inside_block && /}/ { inside_block=0; next }
    !inside_block { print }
  ' temp_file.txt > temp_cleaned.txt
  new_file="${base_name}_ext.maegz"
  gzip -c temp_cleaned.txt > "$folder_path/$new_file"
  rm temp_file.txt temp_cleaned.txt
done
echo "RS results are successfully extracted."

python3 <<PYTHON
import re
import os
import gzip

directory = '.'

def process_file(filename):
    new_header = filename.split('-out_ext')[0]
    with gzip.open(os.path.join(directory, filename), 'rt', encoding='utf-8') as file:
        content = file.read()
    updated_content = re.sub(r'  "4K63_apo_Mutant_[^:]+:[^"]*"', f'  "{new_header}"', content)
    new_filename = filename.replace('out_ext.maegz', 'out_ext_name.maegz')
    with gzip.open(os.path.join(directory, new_filename), 'wt', encoding='utf-8') as file:
        file.write(updated_content)

for filename in os.listdir(directory):
    if filename.endswith('out_ext.maegz'):
        process_file(filename)
PYTHON
echo "RS result names are successfully changed."
