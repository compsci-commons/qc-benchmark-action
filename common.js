const core = require('@actions/core')
const exec = require('@actions/exec')
const yaml = require('js-yaml')
const fs = require('fs')

const benchmarkName = core.getInput('benchmark-name')

const benchmarkOutdir = `benchmark-data/${benchmarkName}`
const bashrc = `${benchmarkOutdir}/.bashrc`
const micromamba = `${benchmarkOutdir}/bin/micromamba`

fs.mkdirSync(benchmarkOutdir, { recursive: true })
fs.closeSync(fs.openSync(bashrc, 'w'))

function _exec (cmd) {
  return exec('bash', ['-l', '-c', `source ${bashrc}; ${cmd}`])
}

let meta = `benchmarks/${benchmarkName}/meta.yaml`
meta = yaml.load(fs.readFileSync(meta, 'utf-8'))

const common = {
  exec: function (cmd) {
    return _exec(cmd)
  },
  initMicromamba: function () {
    _exec(`curl -L https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj ${micromamba}`)
    _exec(`${micromamba} shell init -s bash -p ${benchmarkOutdir}/tmp/micromamba --rc-file ${bashrc}`)
  },
  initEnv: function (envpath, name) {
    if (fs.existsSync(envpath)) {
      _exec(`${micromamba} create -n ${name} -f ${envpath}`)
      return `micromamba activate ${name}`
    }
    return ''
  },
  getBenchmarkName: function () {
    return benchmarkName
  },
  getBenchmarkFile: function (filename) {
    return `benchmarks/${benchmarkName}/${filename}`
  },
  getBenchmarkOutdir: function () {
    return benchmarkOutdir
  },
  getOutpathEnvvars: function () {
    return Object.entries(meta.variables).map(function ([key, path]) { return `${key}="${benchmarkOutdir}/${path}"` }).join(' ')
  },
  getOutpath: function (name) {
    return meta.variables[name]
  }
}

module.exports = common
