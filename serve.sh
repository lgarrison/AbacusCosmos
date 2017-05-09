#!/bin/bash

# Locally serve the Jekyll site for development
# serves on all IPs to facilitate SSH port forwarding

set -v
bundle exec jekyll serve -H '*' -P 4001 -w --config _config.yml,_config_dev.yml
