#!/usr/bin/env bash


read -d '' USAGE <<EOF
$(basename ${0}) [options] INPUT_CSV OUTPUT_NETCDF

Options:

  --run_res=RUN_RES, optional
    Spatial resolution of the grids at which we run the footprint model. A flux
    footprint is modeled and temporarily saved at this resolution per half-hour
    record (row) in the input CSV file. Default: 0.1m.
    
  --out_res=OUT_RES, optional
    Spatial resolution of the output NetCDF file where all the footprints of the
    half-hourly records in the input CSV file are saved. Default: 1m.

  --tmp_dir=TMP_DIR, optional 
    A temporary directory to save temporary output. If given, temporary data
    will be saved here but removed after the script finishes. However, this
    directory will NOT be removed. If not given, by default, the parent
    directory of the output NetCDF file.

  --njobs=NJOBS, optional
    Number of CPUs to run the footprint model over multiple half-hourly records
    in parallel. If 0, all available CPUs will be used. Default: 1.

  --grid_srs=SRS, optional
    Spatial reference system (SRS) of the output grid in GDAL-supported
    format. It must be a $(tput bold)projected SRS with map unit in meters$(tput sgr0). See
    https://gdal.org/programs/gdalsrsinfo.html?highlight=srs_def#cmdoption-gdalsrsinfo-arg-srs_def
    for details on GDAL-supported formats of SRS. The easiest format is to use
    an EPSG code, for example, "epsg:32633" is the SRS of WGS 84 / UTM zone
    33N. Go to https://epsg.io/ to find the EPSG code your desired SRS. Leave it
    empty for non-georeferenced simple planar cartesian coordinate system. 
    $(tput smul)Note for GFZ TEAM group's Zarnekow EC tower, we use epsg:32633$(tput sgr0).

  --grid_bounds="XMIN XMAX YMIN YMAX", optional
    Domain of the output grid on which the footprint to be estimated, given in
    the order of "xmin xmax ymin ymax", in meters in the coordinate system given
    by the option $(tput bold)grid_srs$(tput sgr0). Four values $(tput bold)must
    be double-quoted and separated by space$(tput sgr0). If not given, default
    bounds "-500 500 -500 500".
    $(tput smul)Note for GFZ TEAM group's Zarnekow EC tower, we use "360725.0 361725.0 5971285.0 5972285.0"$(tput sgr0).
     
    
  --receptor_xy="X Y", optional
    Location "x y" of receptor/measurement in the coordinate sysstem of the
    output grid. Two values $(tput bold)must be double-quoted and separated by
    space$(tput sgr0). If not given, default location (0, 0).
    $(tput smul)Note for GFZ TEAM group's Zarnekow EC tower, we use "361224.023952403 5971784.23062857"$(tput sgr0)

Arguments:

  INPUT_CSV, 
    Path to the CSV file where each row is a half-hourly record and at minimum
    contains the following columns with the $(tput bold)EXACT$(tput sgr0) names: Timestamp, qc_Tau,
    z_ensemble, z0_ensemble, wind_speed, u_, L, st_dev_v_, wind_dir.
    where, 
      Timestamp: timestamp of the half-hourly record.  
      qc_Tau: QC flag by turbulance condition of this half-hourly record. If
        qc_Tau >= 2, no footprint is modeled.
      z_ensemble: Aerodynamic/effective height of receptor/measurement (i.e.,
        height above zero-displacement plane), meter
      z0_ensemble: Roughness length, meter
      wind_speed: Mean wind speed, m*s^-1
      u_: Friction velocity, m*s^-1
      L: Monin-Obukhov length, meter
      st_dev_v_: Standard deviation of cross-wind speed, m*s^-1
      wind_dir: wind direction, degree, measured with regard to true north. 

  OUTPUT_NETCDF
    Path to the output NetCDF file that contains all the footprint rasters of
    all the half-hourly records in the input CSV file. 

e.g.

$(basename ${0}) --njobs=4 ./tests/data/test_03a_gapfill_z0d_knn_output.csv ./tests/data/test_04_gen_footprints_netcdf_output_by_user_run.nc

EOF

exe_dir=$(readlink -f ${0} | xargs dirname)

proc_res=0.1
out_res=1
tmp_dir=
njobs=1
grid_srs=
grid_bounds="-500 500 -500 500"
receptor_xy="0 0"

OPTS=`getopt -o h --long run_res:,out_res:,tmp_dir:,njobs:,grid_srs:,grid_bounds:,receptor_xy: -n "${0}" -- "$@"`
if [[ $? != 0 ]]; then echo "Failed parsing options" >&2 ; echo "${USAGE}" ; exit 1 ; fi
eval set -- "${OPTS}"
while true;
do
    case "${1}" in
        -h | --hep ) 
            echo "${USAGE}" ; exit 0 ;;
        --run_res )
            case "${2}" in
                "") shift 2 ;;
                *) proc_res=${2} ; shift 2 ;;
            esac ;;
        --out_res )
            case "${2}" in
                "") shift 2 ;;
                *) out_res=${2} ; shift 2 ;;
            esac ;;
        --tmp_dir )
            case "${2}" in
                "") shift 2 ;;
                *) tmp_dir=${2} ; shift 2 ;;
            esac ;;
        --njobs )
            case "${2}" in
                "") shift 2 ;;
                *) njobs=${2} ; shift 2 ;;
            esac ;;
        --grid_srs )
            case "${2}" in
                "") shift 2 ;;
                *) grid_srs=${2} ; shift 2 ;;
            esac ;;
        --grid_bounds )
            case "${2}" in
                "") shift 2 ;;
                *) grid_bounds=${2} ; shift 2 ;;
            esac ;;
        --receptor_xy )
            case "${2}" in
                "") shift 2 ;;
                *) receptor_xy=${2} ; shift 2 ;;
            esac ;;
        -- ) shift ; break ;;
        * ) break ;;
    esac
done

MINPARAMS=2
if [[ ${#} < ${MINPARAMS} ]]; then
    echo "Missing positional arguments"
    echo "${USAGE}"
    exit 1
fi

if [[ -z ${1} ]] || [[ -z ${2} ]] ; then
    echo "Missing parameters!"
    echo
    echo "${USAGE}"
    exit 1
fi

grid_bounds_arr=($(echo ${grid_bounds}))
if [[ ${#grid_bounds_arr[@]} -ne 4 ]]; then
    echo "--grid_bounds must be a string of four values separated by space and double-quoted."
    echo
    echo "${USAGE}"
    exit 1
fi
receptor_xy_arr=($(echo ${receptor_xy}))
if [[ ${#receptor_xy_arr[@]} -ne 2 ]]; then
    echo "--receptor_xy must be a string of four values separated by space and double-quoted."
    echo
    echo "${USAGE}"
    exit 1
fi

in_csv=${1}
out_nc=${2}

if [[ -z ${tmp_dir} ]]; then
    tmp_dir="$(dirname ${out_nc})"
fi

PYCMD="python ${exe_dir}/fluxfm/fluxfm/cli/estimate_footprint_km.py"

echo "Input CSV of half-hourly records = ${in_csv}"
echo "Output NetCDF of footprint rasters = ${out_nc}"
echo "Run the footprint model at resolution (meter) = ${proc_res}"
echo "Save footprint rasters in NetCDF at resolution (meter) = ${out_res}"
echo "Run and save footprint rasters in the spatial reference system = $([[ -z ${grid_srs} ]] && echo 'simple planar cartesian' || echo ${grid_srs})"
echo "The bounds of the footprint raster grid = ${grid_bounds_arr[@]}"
echo "The location of the receptor in the footprint raster grid = ${receptor_xy_arr[@]}"
echo "Temporary directory to save intermediate temporary data = ${tmp_dir}"
echo "Number of half-hourly records to run in parallel = $([[ ${njobs} -eq 0 ]] && echo 'all available CPUs' || echo ${njobs})"

# GDAL can not export scale_factor at the moment due to a bug
# scale_factor HAS to be 1.
scale_factor=1

# For export array variables for using `parallel`
import_array () {
  local func=$1; shift;
  export $func='() {
    '"$(for arr in $@; do
          declare -p $arr|sed '1s/declare -./&g/'
        done)"'
  }'
}

out_prefix=$(basename ${out_nc})
out_prefix="${tmp_dir}/${out_prefix/%.nc}"

# Below using awk script to extract specific columns 
# a string as awk script
read -r -d '' AWK_CMD <<- EOF
BEGIN {
    FS = ",";
    ncols=split(colstr, my_cols, ",");
}
{
    if (NR==1) {
        for (i=1; i<=NF; i++) {
            ix[\$i] = i;
        }
        for (i=1; i<ncols; i++) {
            printf "%s,", my_cols[i];
        }
        printf "%s\n", my_cols[i]
    }
    if (NR>1) {
        for (i=1; i<ncols; i++) {
            printf "%s,", \$ix[my_cols[i]]
        }
        printf "%s\n", \$ix[my_cols[i]]
    }
}
EOF

tfile_csv=$(mktemp -p ${tmp_dir} --suffix .csv)
awk -v colstr="Timestamp,qc_Tau,z_ensemble,z0_ensemble,wind_speed,u_,L,st_dev_v_,wind_dir" \
    "${AWK_CMD}" < "${in_csv}" | tail -n +2 > ${tfile_csv}

fp_file_list=$(mktemp -p ${tmp_dir} --suffix ".txt")
out_fp_file_list=$(mktemp -p ${tmp_dir} --suffix ".txt")

export PYCMD
export in_csv out_nc proc_res out_res tmp_dir out_prefix fp_file_list scale_factor
export exe_dir grid_srs
import_array my_importer grid_bounds_arr receptor_xy_arr
export out_fp_file_list

doit_ffm_func () {
    linestr="${1}"
    IFS=',' read -ra linearr <<< ${linestr}
    echo "<--- ${linearr[0]}"

    # If qc_Tau >=2, we do not model the footprint.
    [[ ${linearr[1]} -ge 2 ]] && echo "Skip due to QC" && return
    
    echo ${linestr} | grep -q "NaN" && echo "Skip due to NaN values" && return
   
    fp_label=$(date -u -d "${linearr[0]}" "+%Y%m%d%H%M%S")
    fp_file="${out_prefix}_${fp_label}.tif"

    read -r -d '' INI_STR <<- EOF
[meta_variables]
footprint_model = kormann and meixner
footprint_label = ${fp_label}

[input_variables]
; Aerodynamic/effective height of receptor/measurement (i.e., height above
; zero-displacement plane), meter
receptor_height = ${linearr[2]}
; Roughness length, meter
roughness_length = ${linearr[3]} 
; Mean wind speed, m*s^-1
alongwind_speed = ${linearr[4]}
; Friction velocity, m*s^-1
friction_velocity = ${linearr[5]} 
; Monin-Obukhov length, meter
obukhov_length = ${linearr[6]}
; Standard deviation of cross-wind speed, m*s^-1
crosswind_speed_sd = ${linearr[7]}

; Spatial reference system of the output grid in GDAL-supported format. Leave
; it empty for non-georeferenced simple coordinate system. 
grid_spatial_reference = ${grid_srs}

; Domain of the output grid on which the footprint to be estimated, 
; given in (xmin, xmax, ymin, ymax) in four lines, meter
; grid_domain = 360805.0
;               361645.0
;               5971365.0
;               5972205.0
grid_domain = ${grid_bounds_arr[@]}

; Resolution of the output grid, meter.
grid_resolution = ${proc_res}
; Location (x, y) of receptor/measurement in the coordinate system of the grid.
receptor_location = ${receptor_xy_arr[@]}

[optional_variables]
; Mean wind direction with regard to the north designated by
; "north_for_wind_direction", degree.  Leave it empty for default value 0
; degree, that is, the given north aligns with mean wind direction.
wind_direction = ${linearr[8]}
; Type of "North" that defines the wind direction, "due" (true north) or "grid"
; (grid columns along north-south direction)
north_for_wind_direction = due

[output_files]
; Output GeoTiff image file of footprint
footprint_grid_file = ${fp_file} 
footprint_grid_format = GTiff

[user_runtime_parameters]
; For name of data types to use, see https://gdal.org/user/raster_data_model.html#raster-band
data_type = Float32
scale_factor = ${scale_factor}
add_offset = 0
EOF

    tfile_ini=$(mktemp -p ${tmp_dir} --suffix .ini)
    echo "${INI_STR}" > ${tfile_ini}
    echo ${PYCMD} ${tfile_ini}
    PYTHONPATH=${exe_dir}/fluxfm ${PYCMD} ${tfile_ini}
    
    echo ${fp_file} >> ${fp_file_list}
    rm -f ${tfile_ini}

    cur_fname=${fp_file/%.tif/_lowres.tif}
    gdalwarp -overwrite -tr ${out_res} ${out_res} \
        -r sum \
        ${fp_file} ${cur_fname}
    echo "${cur_fname}" >> ${out_fp_file_list}
    rm -f ${fp_file}
}
export -f doit_ffm_func

# echo ${tfile_csv}
# doit_ffm_func "$(tail -n +5 ${tfile_csv} | head -n 1)"
# exit

parallel --jobs ${njobs} \
    --env _ --env my_importer \
    'my_importer; doit_ffm_func {1}' \
    :::: ${tfile_csv}

rm -f ${tfile_csv}

if [[ $(cat ${fp_file_list} | wc -l) -lt 1 ]]; then
    echo "No footprint generated for the records in the input CSV file."
    rm -f ${fp_file_list}
    exit
fi

# doit_gdalwarp_func () {
#     cur_fp_file=${1}
#     cur_fname=${cur_fp_file/%.tif/_lowres.tif}
#     gdalwarp -overwrite -tr ${out_res} ${out_res} \
#         -r sum \
#         ${cur_fp_file} ${cur_fname}
#     echo "${cur_fname}" >> ${out_fp_file_list}
# }
# export -f doit_gdalwarp_func
# 
# parallel --jobs ${njobs} \
#     --env _ \
#     'doit_gdalwarp_func {1}' \
#     :::: ${fp_file_list}
# 
out_fp_file_arr=($(sort ${out_fp_file_list}))

# New way using xarray and rioxarray to generate a NetCDF file
# Not working very well, the timestamp is not well-displayed in QGIS
python ${exe_dir}/00a_multi_gtiff_to_netcdf.py ${out_fp_file_list} ${out_nc}

# # Old way using GDAL programs to generate a NetCDF file
# out_vrt=${out_nc/%.nc/.vrt}
# gdalbuildvrt -separate ${out_vrt} ${out_fp_file_arr[@]}
# for ((i=0; i<${#out_fp_file_arr[@]}; i++)); do
#     fp_label=$(gdalinfo ${out_fp_file_arr[$i]} | grep Description | cut -d'=' -f2 | tr -d [:blank:])
#     src_str="band=\"$((i+1))\">"
#     read -r -d '' tmp_str <<- EOF
# band="$((i+1))">
#     <Metadata>
#       <MDI key="long_name">footprint_${fp_label}</MDI>
#       <MDI key="NETCDF_VARNAME">${fp_label}</MDI>
#       <MDI key="scale_factor">${scale_factor}</MDI>
#       <MDI key="add_offset">0.0</MDI>
#     </Metadata>
# EOF

#     dst_str=""
#     while IFS= read -r linestr; do
#         dst_str="${dst_str}\n${linestr}"
#     done < <(echo "${tmp_str}")
#     dst_str="${dst_str/#\\n}"

#     read -r -d '' sed_cmd <<- EOF
# s#${src_str}#${dst_str}#
# EOF
#     sed -i "${sed_cmd}" ${out_vrt}
# done

# [[ -f ${out_nc} ]] && rm -f ${out_nc}
# gdal_translate -of "NetCDF" \
#     -co FORMAT=NC4 \
#     -co COMPRESS=DEFLATE \
#     -co ZLEVEL=9 \
#     ${out_vrt} ${out_nc}
# rm -f ${out_vrt}

# (cat ${fp_file_list} | xargs -I{} rm -f {})
rm -f ${fp_file_list}
rm -f ${out_fp_file_arr[@]} && rm -f ${out_fp_file_list}
