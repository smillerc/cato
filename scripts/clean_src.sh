#!/bin/bash

fprettify *.f90 --indent 2 --whitespace 3 --whitespace-intrinsics False --strict-indent --enable-replacements --c-relations
