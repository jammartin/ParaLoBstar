#!/bin/bash

# clear log output
echo "Removing old logfiles ..."
ls | grep -E ".*\.[eo][0-9]+(_feedback)?" | xargs -d"\n" rm

