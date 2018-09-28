#!/bin/bash

## One-liner to remove non ATCG REF or ALT, https://www.biostars.org/p/340440/:

mawk '$1 ~ /^#/ {print $0;next} {if ($4 ~ /A|C|T|G/ && $5 ~ /A|C|T|G/) print $0}' $1
