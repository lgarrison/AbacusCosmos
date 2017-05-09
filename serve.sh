#!/bin/bash

# Locally serve the Jekyll site for development
# just uses the first non-loopback IP

set -v
bundle exec jekyll serve -H $(hostname -I |cut -f1 -d' ') -w --config _config.yml,_config_dev.yml
