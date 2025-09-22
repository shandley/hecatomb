import os


"""
Database location
"""
if not config["hecatomb"]["args"]["databases"]:
    try:
        assert(os.environ["HECATOMB_DB"]) is not None
        config["hecatomb"]["args"]["databases"] = os.environ["HECATOMB_DB"]
    except (KeyError, AssertionError):
        config["hecatomb"]["args"]["databases"] = os.path.join(workflow.basedir,"..","resources","databases")


"""
Expand output, temp, database dir and file paths to include base directories
"""
for output_file_path in config["hecatomb"]["args"]["output_paths"]:
    config["hecatomb"]["args"]["output_paths"][output_file_path] =  os.path.join(
        config["hecatomb"]["args"]["output"], config["hecatomb"]["args"]["output_paths"][output_file_path]
    )
for temp_file_path in config["hecatomb"]["args"]["temp_paths"]:
    config["hecatomb"]["args"]["temp_paths"][temp_file_path] =  os.path.join(
        config["hecatomb"]["args"]["output_paths"]["temp"], config["hecatomb"]["args"]["temp_paths"][temp_file_path]
    )
for db_file_path in config["hecatomb"]["args"]["database_paths"]:
    config["hecatomb"]["args"]["database_paths"][db_file_path] =  os.path.join(
        config["hecatomb"]["args"]["databases"], config["hecatomb"]["args"]["database_paths"][db_file_path]
    )


"""
Testing databases
"""
if config["hecatomb"]["args"]["testing"]:
    config["hecatomb"]["args"]["database_paths"]["primaryAA"] += "_testing"
    config["hecatomb"]["args"]["database_paths"]["secondaryAA"] += "_testing"
    config["hecatomb"]["args"]["database_paths"]["primaryNT"] += "_testing"
    config["hecatomb"]["args"]["database_paths"]["secondaryNT"] += "_testing"

