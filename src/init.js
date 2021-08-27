const common = require('./common.js')

async function doInit () {
  const benchmarkName = common.getBenchmarkName()
  console.log(`Download data for ${benchmarkName}.`)

  await common.initMicromamba()

  const envpath = common.getBenchmarkFile('download-env.yaml')
  const envname = 'download-env'
  const prefix = common.getEnvActivate(envpath, envname)
  await common.initEnv(envpath, envname)

  await common.exec(`${prefix}; ${common.getOutpathEnvvars()} bash ${common.getBenchmarkFile('download.sh')}`)
}

module.exports = doInit
