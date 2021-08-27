const core = require('@actions/core');
const common = require('./common.js');

try {
    const benchmarkName = common.getBenchmarkName();
    console.log(`Evaluate results with ${benchmarkName}.`);

    let prefix = common.initEnv(common.getBenchmarkFile('eval-env.yaml'), 'eval-env');

    let result = core.getInput('result');

    common.exec(`${prefix}; ${common.getOutpathEnvvars()} result="${result}" bash ${common.getBenchmarkFile('eval.sh')}`);

    core.setOutput("report", common.getOutpath('report'));
} catch (error) {
    core.setFailed(error.message);
}