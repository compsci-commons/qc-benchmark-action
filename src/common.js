const core = require('@actions/core')
const exec = require('@actions/exec')
const yaml = require('js-yaml')
const fs = require('fs')
const path = require('path')

const benchmarkName = core.getInput('benchmark_name')


const benchmarkOutdir = `benchmark-data/${benchmarkName}`
const mambaPrefix = `${benchmarkOutdir}/mamba`
const condaInit = `${mambaPrefix}/etc/profile.d/conda.sh`;

async function _exec (cmd) {
  await exec.exec('bash', ['-l', '-c', `test -f ${condaInit} && source ${condaInit}; ${cmd}`])
}

let meta = path.join(__dirname, `../benchmarks/${benchmarkName}/meta.yaml`)
meta = yaml.load(fs.readFileSync(meta, 'utf-8'))

const common = {
  exec: async function (cmd) {
    await _exec(cmd)
  },
  initMamba: async function () {
    await _exec(`curl -L --insecure https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh > mambaforge.sh; bash mambaforge.sh -b -p ${mambaPrefix}`)
  },
  getEnvActivate: function (envpath, name) {
    if (fs.existsSync(envpath)) {
      return `source activate ${name}`
    }
    return ''
  },
  initEnv: async function (envpath, name) {
    if (fs.existsSync(envpath)) {
      await _exec(`mamba env create -n ${name} -f ${envpath}`)
    }
  },
  getBenchmarkName: function () {
    return benchmarkName
  },
  getBenchmarkFile: function (filename) {
    return path.join(__dirname, `../benchmarks/${benchmarkName}/${filename}`)
  },
  getBenchmarkOutdir: function () {
    return benchmarkOutdir
  },
  getOutpathEnvvars: function () {
    return Object.entries(meta.paths).map(function ([key, path]) { return `${key}="${benchmarkOutdir}/${path}"` }).join(' ')
  },
  getOutpath: function (name) {
    return `${benchmarkOutdir}/${meta.paths[name]}`
  }
}

module.exports = common
