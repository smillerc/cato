---
project: CATO
summary: A 2D structured Euler solver
author: Sam Miller
src_dir: ./src
         ./build/generated
output_dir: ./doc
media_dir: ./media
exclude_dir: ./src/tests
project_github: https://github.com/smillerc/cato
github: https://github.com/smillerc
docmark: <
display: public
         protected
         private
source: true
graph: true
sort: alpha
coloured_edges: true
extra_filetypes: .inc !
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

--------------------

[TOC]

Brief description
-----------------

A flexible modern fortran code that solves the Euler fluid equations using the following flux solvers [AUSM+-up, (M-)AUSMPW+, SLAU(2), FVLEG]
