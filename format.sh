#!/bin/bash
echo "Formatting Python code"
autopep8 --in-place --recursive . --max-line-length 120
