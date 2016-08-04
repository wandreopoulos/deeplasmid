#!/bin/bash -l
set -ex
rsync -a --exclude=".*" io/ $PREFIX/
