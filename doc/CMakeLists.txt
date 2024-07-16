#
#  Copyright (C) 2024 by the authors of the RAYLEIGH code.
#
#  This file is part of RAYLEIGH.
#
#  RAYLEIGH is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  RAYLEIGH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RAYLEIGH; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.


CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

####################################################
# Generate the manual
####################################################

# Find the documentation sources
file(GLOB_RECURSE SOURCES
  source/*)

# Build HTML target
message(STATUS "===== Configuring documentation ===============")
message(STATUS "")

ADD_CUSTOM_COMMAND(
  OUTPUT
    ${PROJECT_BINARY_DIR}/doc/build/html/index.html
  COMMAND
    sphinx-build -M html ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/doc/build
  DEPENDS
    ${SOURCES}
  COMMENT
    "Generating the html documentation."
)

ADD_CUSTOM_TARGET(doc_html
  DEPENDS 
    ${PROJECT_BINARY_DIR}/doc/build/html/index.html
)

# Build pdf target
ADD_CUSTOM_COMMAND(
  OUTPUT
    ${PROJECT_BINARY_DIR}/doc/build/latex/rayleigh.pdf
  COMMAND
    sphinx-build -M latexpdf ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/doc/build
  DEPENDS
    ${SOURCES}
  COMMENT
    "Generating the latexpdf documentation."
)

ADD_CUSTOM_TARGET(doc_pdf
  DEPENDS
    ${PROJECT_BINARY_DIR}/doc/build/latex/rayleigh.pdf
)

################################################################
# The final target: Build it all
################################################################
message(STATUS "Setting up build information")
ADD_CUSTOM_TARGET(doc
  DEPENDS
    doc_html
    doc_pdf
)
