set shell := ["bash", "-c"]

# help
default: help

alias h := help
help:
    just -l -f {{justfile()}}

root := justfile_directory()

# meson working directory
work_dir := (root)
# install dir
install_dir := (root) / "install"

alias s := setup
setup:
    cd {{work_dir}}; meson setup --prefix {{install_dir}} build

alias r := resetup
resetup:
    cd {{work_dir}}; meson setup --reconfigure --prefix {{install_dir}} build

alias b := build
build:
    cd {{work_dir}}; meson compile -C build
    cd {{work_dir}}; meson install -C build

alias f := format
format:
    cd {{work_dir}}; clang-format -style=Mozilla -i src/*
