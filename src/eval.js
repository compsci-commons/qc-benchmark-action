const core = require('@actions/core')
const common = require('./common.js')

async function doEval () {
  const benchmarkName = common.getBenchmarkName()
  console.log(`Evaluate results with ${benchmarkName}.`)

  const envpath = common.getBenchmarkFile('eval-env.yaml')
  const envname = 'eval-env'
  const prefix = common.getEnvActivate(envpath, envname)
  await common.initEnv(envpath, envname)

  const result = core.getInput('result')

  await common.exec(`${prefix}; ${common.getOutpathEnvvars()} result="${result}" bash ${common.getBenchmarkFile('eval.sh')}`)

  core.setOutput('report', common.getOutpath('report'))
}

module.exports = doEval
