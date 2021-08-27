const core = require('@actions/core');
const common = require('./common.js');
const fs = require('fs');

try {
    const benchmarkName = common.getBenchmarkName();
    console.log(`Download data for ${benchmarkName}.`);
    
    let prefix = common.initEnv(common.getBenchmarkFile('download-env.yaml'), 'download-env');

    fs.mkdirSync(common.getBenchmarkOutdir(), {recursive: true});

    common.initMicromamba();

    common.exec(`${prefix}; ${common.getOutpathEnvvars()} bash ${common.getBenchmarkFile('download.sh')}`);
} catch (error) {
    core.setFailed(error.message);
}