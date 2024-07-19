#!/bin/bash

# Prompt the user for the installation directory
read -p "Enter the directory where you want to install SRIM (write the full path, don't use wildcards): " SRIMDIR

# Create the directory if it doesn't exist
mkdir -p "$SRIMDIR"

# Download the SRIM installer to the specified directory
wget --output-document="$SRIMDIR/SRIM_INSTALL.exe" http://www.srim.org/SRIM/SRIM-2013-Std.e

# Run the installer using Wine
wine "$SRIMDIR/SRIM_INSTALL.exe"
