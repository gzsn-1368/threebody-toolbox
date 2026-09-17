#!/bin/sh
set -eu

if [ "$#" -ne 4 ]; then
    echo "usage: build_cspice_unix.sh SOURCE_DIR INCLUDE_DIR OUTPUT_LIBRARY OBJECT_DIR" >&2
    exit 2
fi

source_dir=$1
include_dir=$2
output_library=$3
object_dir=$4
compiler=${CC:-cc}

mkdir -p "$object_dir" "$(dirname "$output_library")"

for source_file in "$source_dir"/*.c; do
    stem=$(basename "$source_file" .c)
    object_file="$object_dir/$stem.o"
    if [ ! -f "$object_file" ] || [ "$source_file" -nt "$object_file" ]; then
        echo "Compiling CSPICE: $stem.c"
        "$compiler" -c -O2 -fPIC -std=c89 -DNON_UNIX_STDIO \
            -I"$include_dir" "$source_file" -o "$object_file"
    fi
done

echo "Archiving CSPICE: $output_library"
rm -f "$output_library"
ar rcs "$output_library" "$object_dir"/*.o
