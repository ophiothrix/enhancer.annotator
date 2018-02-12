#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

R --vanilla --args $1 $2 < annotate.enhancers.with.genes.R
