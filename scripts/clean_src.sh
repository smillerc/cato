#!/bin/bash

cd ../src
find . -name "*.f90"|while read fname; do
  echo "$fname"
  fprettify ${fname} --indent 2 --whitespace 3 --whitespace-intrinsics False --strict-indent --enable-replacements --c-relations
done

cd ../tests
find . -name "*.pf"|while read fname; do
  echo "$fname"
  fprettify ${fname} --indent 2 --whitespace 3 --whitespace-intrinsics False --strict-indent --enable-replacements --c-relations
done

cd ..
