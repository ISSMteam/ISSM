#To generate list of directories
#
#find . -type d -exec ls -d {} \;

bin_SCRIPTS =
if WRAPPERS
if MATLAB
if !DEVELOPMENT
# TODO: There should be a better way of doing this so that manual updating is
#		not required when a new subdirectory or module is introduced (search
#		all subdirectories for file type, checking against an exclude list)
bin_SCRIPTS += \
	${ISSM_DIR}/src/m/archive/*.m \
	${ISSM_DIR}/src/m/array/*.m \
	${ISSM_DIR}/src/m/boundaryconditions/*.m \
	${ISSM_DIR}/src/m/classes/*.m \
	${ISSM_DIR}/src/m/classes/clusters/*.m \
	${ISSM_DIR}/src/m/classes/qmu/*.m \
	${ISSM_DIR}/src/m/classes/qmu/dakota_method/*.m \
	${ISSM_DIR}/src/m/consistency/*.m \
	${ISSM_DIR}/src/m/coordsystems/*.m \
	${ISSM_DIR}/src/m/dev/issmversion.m \
	${ISSM_DIR}/src/m/exp/*.m \
	${ISSM_DIR}/src/m/exp/operation/*.m \
	${ISSM_DIR}/src/m/extrusion/*.m \
	${ISSM_DIR}/src/m/geometry/*.m \
	${ISSM_DIR}/src/m/interp/*.m \
	${ISSM_DIR}/src/m/inversions/*.m \
	${ISSM_DIR}/src/m/io/*.m \
	${ISSM_DIR}/src/m/materials/*.m \
	${ISSM_DIR}/src/m/mech/*.m \
	${ISSM_DIR}/src/m/mesh/*.m \
	${ISSM_DIR}/src/m/mesh/planet/*.m \
	${ISSM_DIR}/src/m/mesh/rifts/*.m \
	${ISSM_DIR}/src/m/miscellaneous/*.m \
	${ISSM_DIR}/src/m/modules/*.m \
	${ISSM_DIR}/src/m/modeldata/*.m \
	${ISSM_DIR}/src/m/netcdf/*.m \
	${ISSM_DIR}/src/m/os/*.m \
	${ISSM_DIR}/src/m/parameterization/*.m \
	${ISSM_DIR}/src/m/partition/*.m \
	${ISSM_DIR}/src/m/plot/*.m \
	${ISSM_DIR}/src/m/plot/colormaps/*.m \
	${ISSM_DIR}/src/m/print/*.m \
	${ISSM_DIR}/src/m/qmu/*.m \
	${ISSM_DIR}/src/m/qmu/examples/*.m \
	${ISSM_DIR}/src/m/qmu/plot/*.m \
	${ISSM_DIR}/src/m/qmu/setupdesign/*.m \
	${ISSM_DIR}/src/m/shp/*.m \
	${ISSM_DIR}/src/m/solve/*.m \
	${ISSM_DIR}/src/m/solvers/*.m \
	${ISSM_DIR}/src/m/string/*.m
endif
endif

if PYTHON
if !DEVELOPMENT
# TODO: There should be a better way of doing this so that manual updating is
#		not required when a new subdirectory or module is introduced (search
#		all subdirectories for file type, checking against an exclude list)
bin_SCRIPTS += \
	${ISSM_DIR}/src/m/archive/*.py \
	${ISSM_DIR}/src/m/array/*.py \
	${ISSM_DIR}/src/m/boundaryconditions/*.py \
	${ISSM_DIR}/src/m/classes/*.py \
	${ISSM_DIR}/src/m/classes/clusters/*.py \
	${ISSM_DIR}/src/m/classes/qmu/*.py \
	${ISSM_DIR}/src/m/classes/qmu/dakota_method/*.py \
	${ISSM_DIR}/src/m/consistency/*.py \
	${ISSM_DIR}/src/m/coordsystems/*.py \
	${ISSM_DIR}/src/m/dev/issmversion.py \
	${ISSM_DIR}/src/m/exp/*.py \
	${ISSM_DIR}/src/m/extrusion/*.py \
	${ISSM_DIR}/src/m/geometry/*.py \
	${ISSM_DIR}/src/m/interp/*.py \
	${ISSM_DIR}/src/m/inversions/*.py \
	${ISSM_DIR}/src/m/io/*.py \
	${ISSM_DIR}/src/m/materials/*.py \
	${ISSM_DIR}/src/m/mech/*.py \
	${ISSM_DIR}/src/m/mesh/*.py \
	${ISSM_DIR}/src/m/mesh/planet/*.py \
	${ISSM_DIR}/src/m/mesh/rifts/*.py \
	${ISSM_DIR}/src/m/miscellaneous/*.py \
	${ISSM_DIR}/src/m/modules/*.py \
	${ISSM_DIR}/src/m/netcdf/*.py \
	${ISSM_DIR}/src/m/os/*.py \
	${ISSM_DIR}/src/m/parameterization/*.py \
	${ISSM_DIR}/src/m/partition/*.py \
	${ISSM_DIR}/src/m/plot/*.py \
	${ISSM_DIR}/src/m/plot/colormaps/*.py \
	${ISSM_DIR}/src/m/qmu/*.py \
	${ISSM_DIR}/src/m/qmu/setupdesign/*.py \
	${ISSM_DIR}/src/m/shp/*.py \
	${ISSM_DIR}/src/m/solve/*.py \
	${ISSM_DIR}/src/m/solvers/*.py
endif
endif

if JAVASCRIPT
if !DEVELOPMENT
# TODO: There should be a better way of doing this so that manual updating is
#		not required when a new subdirectory or module is introduced (search
#		all subdirectories for file type, checking against an exclude list)
js_scripts = \
	${ISSM_DIR}/src/m/array/*.js \
	${ISSM_DIR}/src/m/boundaryconditions/*.js \
	${ISSM_DIR}/src/m/classes/*.js \
	${ISSM_DIR}/src/m/classes/clusters/*.js \
	${ISSM_DIR}/src/m/consistency/*.js \
	${ISSM_DIR}/src/m/exp/*.js \
	${ISSM_DIR}/src/m/extrusion/*.js \
	${ISSM_DIR}/src/m/geometry/*.js \
	${ISSM_DIR}/src/m/inversions/*.js \
	${ISSM_DIR}/src/m/io/*.js \
	${ISSM_DIR}/src/m/materials/*.js \
	${ISSM_DIR}/src/m/mesh/*.js \
	${ISSM_DIR}/src/m/miscellaneous/*.js \
	${ISSM_DIR}/src/m/parameterization/*.js \
	${ISSM_DIR}/src/m/plot/*.js \
	${ISSM_DIR}/src/m/print/*.js \
	${ISSM_DIR}/src/m/shp/*.js \
	${ISSM_DIR}/src/m/solve/*.js \
	${ISSM_DIR}/src/m/solvers/*.js
bin_SCRIPTS += issm-bin.js

issm-bin.js: ${js_scripts}
	rm -f issm-bin.js
	cat ${js_scripts} > issm-bin.js

endif
endif

endif
