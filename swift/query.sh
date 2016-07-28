#!/bin/bash

sqlite3 -noheader -separator , -batch swift_provenance.db "select script_run_id, duration from script_run;"
