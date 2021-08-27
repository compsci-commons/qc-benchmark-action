const common = require('./common.js')

function doInit () {
  const benchmarkName = common.getBenchmarkName()
  console.log(`Download data for ${benchmarkName}.`)

  const prefix = common.initEnv(common.getBenchmarkFile('download-env.yaml'), 'download-env')

  common.initMicromamba()

  common.exec(`${prefix}; ${common.getOutpathEnvvars()} bash ${common.getBenchmarkFile('download.sh')}`)
}

module.exports = doInit
