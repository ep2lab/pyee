#!/bin/bash
addgroup --gid $2 user
adduser --disabled-password --gecos '' --uid $1 --gid $2 user
exit