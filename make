#!/usr/bin/python
import argparse
import os
import shutil
import subprocess as s
import sys

parser = argparse.ArgumentParser(
    "a build script using CMake and either Ninja or the OS default generator"
)

parser.add_argument("-t", "--target", default="all", help="the build target")
parser.add_argument(
    "--debug", action="store_true", help="use build type 'Debug' instead of 'Release'"
)
parser.add_argument(
    "--default-generator",
    action="store_true",
    help="use OS default generator instead of ninja",
)

args = parser.parse_args()

target_dir = os.path.dirname(os.path.abspath(__file__))


def generate_build_files():
    cmake_command = ["cmake", "-B", "build", "-S", "source"]

    if not args.default_generator:
        cmake_command.append("-G")
        cmake_command.append("Ninja")

    if args.debug:
        cmake_command.append("-DCMAKE_BUILD_TYPE=Debug")
    else:
        cmake_command.append("-DCMAKE_BUILD_TYPE=Release")

    res = s.run(cmake_command, cwd=target_dir)
    if res.returncode != 0:
        sys.exit(res.returncode)


def build():
    if args.default_generator:
        res = s.run(["make", "-C", "build"], cwd=target_dir)
    else:
        res = s.run(["ninja", "-C", "build"], cwd=target_dir)
    if res.returncode != 0:
        sys.exit(res.returncode)


def remove_build_dir():
    global target_dir
    build_dir = ""
    if target_dir is None:
        build_dir = "build"
    else:
        build_dir = os.path.join(target_dir, "build")
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)


if args.target == "all":
    generate_build_files()
    build()
elif args.target == "cmake":
    generate_build_files()
elif args.target == "rebuild" or args.target == "build":
    build()
elif args.target == "clear":
    remove_build_dir()
    print("Sucessfully removed build directory")
else:
    print('The target "' + args.target + '" is unknown.')
    parser.print_help()
    sys.exit(1)
