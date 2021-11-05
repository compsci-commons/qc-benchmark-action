const common = require('./common.js')
const core = require('@actions/core')

async function doInit () {
  const benchmarkName = common.getBenchmarkName()
  console.log(`Download data for ${benchmarkName}.`)

  await common.initMamba()

  const envpath = common.getBenchmarkFile('download-env.yaml')
  const envname = 'download-env'
  const prefix = common.getEnvActivate(envpath, envname)
  await common.initEnv(envpath, envname)

  await common.exec(`${prefix}; ${common.getOutpathEnvvars()} bash ${common.getBenchmarkFile('download.sh')}`)

  core.setOutput('data', common.getBenchmarkOutdir())
  console.log(`Storing output files in ${common.getBenchmarkOutdir()}`)
}

module.exports = doInit
