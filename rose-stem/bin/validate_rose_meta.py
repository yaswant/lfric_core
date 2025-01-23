#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Validate the rose metadata in lfric_apps
"""

import os
import sys
import re
import subprocess
import argparse

DEFAULT_APPS = {
   "skeleton": ["skeleton"],
   "mesh_tools": ["mesh_tools"]
}

COMPONENTS = [
    "driver",
]

APPLICATIONS = [
    "skeleton",
    "mesh_tools"
]

OTHERS = [
    "infrastructure",
    "mesh_tools"
]


def run_subprocess(command, source_dir):
    """
    Launch a subprocess command and return output
    If source_dir provided, changes directory to source directory top level
    as rose commands need to be launched from here for meta imports to work.
    Then returns to original loc.
    """
    initial_dir = os.getcwd()
    if source_dir:
        os.chdir(source_dir)
    process = subprocess.run(
        command, shell=True, capture_output=True, text=True
    )
    os.chdir(initial_dir)
    return process


def check_metadata(application, source_dir):
    """
    Run rose metadata-check on application metadata
    """
    if application in OTHERS:
        meta_path = os.path.join(".", application, "rose-meta")
    elif application in APPLICATIONS:
        meta_path = os.path.join("applications", application, "rose-meta")
    elif application in COMPONENTS:
        meta_path = os.path.join("components", application, "rose-meta")
    else:
        sys.exit(f"Application {application} not defined in either "
                 "APPLICATIONS or COMPONENTS lists.")
    retcode = 0
    paths = os.listdir(os.path.join(source_dir, meta_path))
    for path in paths:
        meta_dir = os.path.join(meta_path, path, "HEAD")
        validate_command = f"rose metadata-check --verbose -C {meta_dir}"
        process = run_subprocess(validate_command, source_dir)
        print(process.stdout)
        if process.returncode:
            print(f"Failures checking {application}", file=sys.stderr)
            print(process.stderr, file=sys.stderr)
        retcode += process.returncode
    return retcode


def validate_app(app, source_dir, ignore_sc):
    """
    Run rose macro --validate on rose-stem apps
    """
    app_dir = os.path.join(source_dir, "rose-stem", "app", app)
    validate_command = (
        f"rose macro --validate -M {source_dir} -C {app_dir} --no-warn version"
    )
    process = run_subprocess(validate_command, False)
    print(process.stdout)
    if process.returncode:
        err = []
        split_err = process.stderr.split("\n")
        i = 0
        while i < len(split_err):
            line = split_err[i]
            if "opts=suite_controlled" in line and ignore_sc:
                i += 2
                continue
            if line.strip():
                err.append(line)
            i += 1
        if len(err) > 1:
            print(f"Failures validating {app}", file=sys.stderr)
            print("\n".join(err), file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    """
    Main function - Loop over DEFAULT_APPS and call validate_rose_meta script
    on each, looping over custom apps as required.
    """

    failures = False

    parser = argparse.ArgumentParser(
        description="Validate the rose metadata in lfric_apps"
    )
    parser.add_argument(
        "-s",
        "--source",
        help="The top level of lfric directory.",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--ignore_sc",
        help="Ignore failures from rose-app-suite_controlled.conf files.",
        default=True
    )
    args = parser.parse_args()

    print("[INFO] Checking Metadata")
    print("=================")
    for application in DEFAULT_APPS:
        print("\n---------------------------------------------------")
        print(f"[INFO]Running rose metadata-check on {application}")
        retcode = check_metadata(application, args.source)
        if retcode:
            print(f"[FAIL] {application} failed rose metadata-check")
            failures = True
        else:
            print(f"[PASS] {application} passed rose metadata-check")


    print("[INFO] Validating App Configurations")
    if args.ignore_sc:
        print("[WARNING]: Ignoring rose-app-suite_controlled.conf failures")
    print("===========================================================")
    for application in DEFAULT_APPS:
        apps = DEFAULT_APPS[application]
        for app in apps:
            print("\n--------------------------------------------------------")
            print(f"[INFO] Validating {application} with app {app}")
            retcode = validate_app(app, args.source, args.ignore_sc)
            if retcode:
                print(f"[FAIL] {application} with app {app} failed to validate")
                failures = True
            else:
                print(f"[PASS] {application} with app{app} "
                      "validated successfully")

    if failures:
        sys.exit("There were metadata validation failures")
    print("\n[PASS] Metadata Validated Successfully")
