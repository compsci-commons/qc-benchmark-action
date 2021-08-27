const core = require('@actions/core');
const exec = require('@actions/exec');
const github = require('@actions/github');
const fs = require('fs');
const common = require('./common.js');

try {
    const benchmarkName = common.getBenchmarkName();
    console.log(`Benchmark ${benchmarkName}!`);
    
    let prefix = common.initEnv(getBenchmarkFile('download-env.yaml'), 'download-env');

    common.exec(`${prefix}; source download.sh`);
} catch (error) {
    core.setFailed(error.message);
}