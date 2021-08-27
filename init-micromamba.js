const core = require('@actions/core');
const exec = require('@actions/exec');
const github = require('@actions/github');
const fs = require('fs');
const common = require('./common.js');

try {
    let outdir = common.getBenchmarkOutdir();
    fs.mkdirSync(outdir, { recursive: true });
    common.initMicromamba();
} catch (error) {
    core.setFailed(error.message);
}