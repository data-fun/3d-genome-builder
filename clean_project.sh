#! /usr/bin/env bash

# Stop at first error or with an unset variable.
set -euo pipefail

confirm() {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case "${response}" in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            echo "Aborting"
            false
            ;;
    esac
}


if [[ $# -eq 0 ]]
then
    echo "Missing argument"
    exit 1
fi

if [[ ! -d "$1" ]]
then
    echo "Directory $1 does not exist"
    exit 1
fi

echo "$1 will be cleaned."
du -csh "$1"/*

confirm "Are you sure? [y/N]"

targets="contact_maps dense_matrix HiC-Pro  logs pastis sequence structure .snakemake"
for name in ${targets}
do
    rm -rf "$1/${name}"
done

echo "$1 is cleaned."
du -csh "$1"/*