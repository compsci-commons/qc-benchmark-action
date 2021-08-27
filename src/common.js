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

async function _exec (cmd) {
  await exec.exec('bash', ['-l', '-c', `source ${bashrc}; ${cmd}`])
}

let meta = `benchmarks/${benchmarkName}/meta.yaml`
meta = yaml.load(fs.readFileSync(meta, 'utf-8'))

const common = {
  exec: async function (cmd) {
    await _exec(cmd)
  },
  initMicromamba: async function () {
    await _exec(`curl -L https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C ${benchmarkOutdir} bin/micromamba`)
    await _exec(`${micromamba} shell hook -s bash -p ${benchmarkOutdir}/tmp/micromamba > ${bashrc}`)
  },
  getEnvActivate: function (envpath, name) {
    if (fs.existsSync(envpath)) {
      return `micromamba activate ${name}`
    }
    return ''
  },
  initEnv: async function (envpath, name) {
    if (fs.existsSync(envpath)) {
      await _exec(`${micromamba} create --root-prefix ${benchmarkOutdir}/tmp/micromamba -n ${name} -f ${envpath}`)
    }
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
