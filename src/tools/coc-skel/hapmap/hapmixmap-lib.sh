#!/bin/bash

# A library of common functions for bash scripts.

function print_variable {
	VAR_NAME="$1"
	eval "echo \"${VAR_NAME}\"=\\'\"\$${VAR_NAME}\"\\'"
}

function print_settings {
	print_variable POPULATIONS
}


