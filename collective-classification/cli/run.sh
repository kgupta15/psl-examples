#!/bin/bash

readonly JAR_URL='https://linqs-data.soe.ucsc.edu/maven/repositories/psl-releases/org/linqs/psl-cli/CANARY/psl-cli-CANARY.jar'
readonly JAR_FILE='psl-cli-CANARY.jar'
readonly BASE_NAME='simple_cc'

FETCH_COMMAND=''

function err() {
   echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: $@" >&2
}

# Check for:
#  - wget or curl (final choice to be set in FETCH_COMMAND)
#  - java
function check_requirements() {
   local hasWget

   type wget > /dev/null 2> /dev/null
   hasWget=$?

   type curl > /dev/null 2> /dev/null
   if [[ "$?" -eq 0 ]]; then
      FETCH_COMMAND="curl -o"
   elif [[ "${hasWget}" -eq 0 ]]; then
      FETCH_COMMAND="wget -O"
   else
      err 'wget or curl required to download jar'
      exit 10
   fi

   type java > /dev/null 2> /dev/null
   if [[ "$?" -ne 0 ]]; then
      err 'java required to run project'
      exit 13
   fi
}

function fetch_jar() {
   if [[ -e "${JAR_FILE}" ]]; then
      echo "Jar found cached, skipping download."
      return
   fi

   echo "Downloading the jar with command: $FETCH_COMMAND"
   $FETCH_COMMAND "${JAR_FILE}" "${JAR_URL}"
   if [[ "$?" -ne 0 ]]; then
      err 'Failed to download jar'
      exit 20
   fi
}

function run() {
   java -jar "${JAR_FILE}" -infer -model "${BASE_NAME}.psl" -data "${BASE_NAME}.data"
   if [[ "$?" -ne 0 ]]; then
      err 'Failed to run'
      exit 60
   fi
}

function main() {
   check_requirements
   fetch_jar
   run
}

main "$@"
