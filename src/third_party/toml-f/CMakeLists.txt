# This file is part of toml-f.
#
# Copyright (C) 2019-2020 Sebastian Ehlert
#
# Licensed under either of Apache License, Version 2.0 or MIT license at your option; you may not
# use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing permissions and limitations under
# the License.
add_library(
  toml-f
  src/tomlf/all.f90
  src/tomlf/constants.f90
  src/tomlf/de.f90
  src/tomlf/structure.f90
  src/tomlf/utils.f90
  src/tomlf/build.f90
  src/tomlf/datetime.f90
  src/tomlf/error.f90
  src/tomlf/ser.f90
  src/tomlf/type.f90
  src/tomlf/version.f90
  src/tomlf/build/array.f90
  src/tomlf/build/keyval.f90
  src/tomlf/build/table.f90
  src/tomlf/de/character.f90
  src/tomlf/de/tokenizer.f90
  src/tomlf/structure/base.f90
  src/tomlf/structure/vector.f90
  src/tomlf/type/array.f90
  src/tomlf/type/keyval.f90
  src/tomlf/type/table.f90
  src/tomlf/type/value.f90
  src/tomlf/utils/convert.f90
  src/tomlf/utils/verify.f90)
