#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: plotting facility for L95 and QG
"""

import os
import argparse
import sys

# -----------------------------------------------------------------------------
# Main parser -----------------------------------------------------------------
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser()
subparser_model = parser.add_subparsers(help="model help", required=True, dest="model")

# -----------------------------------------------------------------------------
# L95 model sub-parsers --------------------------------------------------------
# -----------------------------------------------------------------------------
parser_l95 = subparser_model.add_parser("l95", help="L95 parser")

# L95 diagnostics sub-parsers
subparser_diagnostic_l95 = parser_l95.add_subparsers(help="L95 diagnostic help", required=True, dest="diagnostic")
parser_l95_cost = subparser_diagnostic_l95.add_parser("cost", help="L95 cost parser")
parser_l95_fields = subparser_diagnostic_l95.add_parser("fields", help="L95 fields parser")
parser_l95_timeseries = subparser_diagnostic_l95.add_parser("timeseries", help="L95 timeseries parser")
parser_l95_errors = subparser_diagnostic_l95.add_parser("errors", help="L95 errors parser")
parser_l95_increments = subparser_diagnostic_l95.add_parser("increments", help="L95 increments parser")

# L95 cost arguments
parser_l95_cost.add_argument("filepath", type=str, help="File path")
parser_l95_cost.add_argument("--output", help="Output file path")

# L95 fields arguments
parser_l95_fields.add_argument("filepath", help="Analysis file path")
parser_l95_fields.add_argument("-bg", "--bgfilepath", help="Background file path", default=None)
parser_l95_fields.add_argument("-t",  "--truthfilepath", help="Truth file path", default=None)
parser_l95_fields.add_argument("-o",  "--obsfilepath",  help="Obs file path", default=None)
parser_l95_fields.add_argument("--output", help="Output file path")

# L95 RMSE time series arguments
parser_l95_timeseries.add_argument("filepath", help="Files path with %id% template that will be replaced with values in times argument")
parser_l95_timeseries.add_argument("-fKey", "--fileKey", help="Specify a Key in the legend for file 1",default=None)
parser_l95_timeseries.add_argument("-f2", "--filepath2", help="a second File path with %id% template that will be replaced with values in times argument",default=None)
parser_l95_timeseries.add_argument("-f2Key", "--file2Key", help="Specify a Key in the legend for file 2",default=None)
parser_l95_timeseries.add_argument("-f3", "--filepath3", help="a third File path with %id% template that will be replaced with values in times argument",default=None)
parser_l95_timeseries.add_argument("-f3Key", "--file3Key", help="Specify a Key in the legend for file 3",default=None)
parser_l95_timeseries.add_argument("-f4", "--filepath4", help="a fourth File path with %id% template that will be replaced with values in times argument", default=None)
parser_l95_timeseries.add_argument("-f4Key", "--file4Key", help="Specify a Key in the legend for file 4",default=None)
parser_l95_timeseries.add_argument("truthfilepath", help="Truth files path with %id% template that will be replaced with values in times argument")
parser_l95_timeseries.add_argument("times", help="Time series pattern values (will be used to replace %id% in the file names")
parser_l95_timeseries.add_argument("--output", help="Output file path")

# L95 errors arguments
parser_l95_errors.add_argument("filepath", help="Analysis file path")
parser_l95_errors.add_argument("truthfilepath", help="Truth file path")
parser_l95_errors.add_argument("-bg", "--bgfilepath", help="Background file path", default=None)
parser_l95_errors.add_argument("--output", help="Output file path")
parser_l95_errors.add_argument("--recenter", action="store_true", help="Center plot at 1st and last variables?")
parser_l95_errors.add_argument("--title", help="Specify a title for the plot with --title 'my title' ")
#L95 increments arguments
parser_l95_increments.add_argument("filepath", help="Analysis file path")
parser_l95_increments.add_argument("bgfilepath", help="Background file path")
parser_l95_increments.add_argument("-t", "--truthfilepath", help="Truth file path", default=None)
parser_l95_increments.add_argument("--output", help="Output file path")
parser_l95_increments.add_argument("--recenter", action="store_true", help="Center plot at 1st and last variables?")
parser_l95_increments.add_argument("--title", help="Specify a title for the plot with --title 'my title' ")
# -----------------------------------------------------------------------------
# QG model sub-parsers --------------------------------------------------------
# -----------------------------------------------------------------------------
parser_qg = subparser_model.add_parser("qg", help="QG parser")

# QG diagnostics sub-parsers
subparser_diagnostic_qg = parser_qg.add_subparsers(help="QG diagnostic help", required=True, dest="diagnostic")
parser_qg_cost = subparser_diagnostic_qg.add_parser("cost", help="QG cost parser")
parser_qg_fields = subparser_diagnostic_qg.add_parser("fields", help="QG fields parser")
parser_qg_obs = subparser_diagnostic_qg.add_parser("obs", help="QG obs parser")

# QG cost arguments
parser_qg_cost.add_argument("filepath", type=str, help="File path")
parser_qg_cost.add_argument("--output", help="Output file path")

# QG fields arguments
parser_qg_fields.add_argument("filepath", help="File path")
parser_qg_fields.add_argument("basefilepath", nargs="?", help="Base file path", default=None)
parser_qg_fields.add_argument("--plotObsLocations", help="Specify an observation file to plot the obs positions")
parser_qg_fields.add_argument("--plotwind", dest="plotwind", action="store_true", help="Plot wind")
parser_qg_fields.add_argument("--gif", help="Gif pattern values, separated with commas, replacing %id% in the file path")
parser_qg_fields.add_argument("--output", help="Output file path")
parser_qg_fields.add_argument("--title", help="Specify a title for the plot with --title 'my title' ")

# QG obs arguments
parser_qg_obs.add_argument("filepath", type=str, help="File path")
parser_qg_obs.add_argument("--output", help="Output file path")

# -----------------------------------------------------------------------------
# Parse and print arguments ---------------------------------------------------
# -----------------------------------------------------------------------------
# Parse argument
args = parser.parse_args()

# Print arguments
print("Parameters:")
for arg in vars(args):
    if not arg is None:
        print(" - " + arg + ": " + str(getattr(args, arg)))

# -----------------------------------------------------------------------------
# Load appropriate module and run function ------------------------------------
# -----------------------------------------------------------------------------
sys.path.insert(1, "plot")
module = __import__(args.model + "_" + args.diagnostic)
func = getattr(module, "func")
print("Run script")
func(args)
